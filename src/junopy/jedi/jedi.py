# Copyright (c) 2018 Yash Sarkango
# Licensed under the MIT License

import os
import numpy as np
import pandas as pd
import datetime as dt
import glob

import matplotlib.colors as colors

from pathlib import Path
from scipy.constants import m_e, m_p, e

class JediData:
    """
    Class to read and do things with a single JEDI data file.

    Args:
        filename (Path | None): Path to the JEDI data file.
        start (datetime | None): Start time for the data.
        stop (datetime | None): Stop time for the data.
        species (str | None): Particle species to filter data for.
        units (str): Units of the data, either 'FLUX' or 'PSD'.
        num_pitchangle_bins (int): Number of pitch angle bins.
        do_save_only (bool): If True, only save the data without further processing.
        min_energy (float | None): Minimum energy.
        max_energy (float | None): Maximum energy.
        min_pitchangle (float): Minimum pitch angle.
        max_pitchangle (float): Maximum pitch angle.

    Attributes:
        num_ssd (int): Number of JEDI SSDs.
        num_channels (int): Number of energy channels per SSD.
        has_valid_data (bool): Whether the data contains valid entries.
        has_energy_flux (bool): Whether energy flux data is available.
        df_raw (pd.DataFrame): Raw data DataFrame.
        metadata (pd.DataFrame): Metadata DataFrame containing channel information.
        ssd_data (list[pd.DataFrame]): List containing data split by SSD.
        ssd_sigma (list[pd.DataFrame]): List containing measurement uncertainty for each SSD.
        ssd_pitchangles (list[pd.Series]): List containing pitch angles for each SSD.

    """

    num_ssd = 6
    num_channels = -1
    has_valid_data = False
    has_energy_flux = False

    def __init__(
        self,
        filename: Path = None,
        start: dt.datetime = None,
        stop: dt.datetime = None,
        species: str | None = None,
        units: str = "FLUX",
        num_pitchangle_bins: int = 18,
        do_save_only: bool = False,
        min_energy: float | None = None,
        max_energy: float | None = None,
        min_pitchangle: float = 0,
        max_pitchangle: float = 180,
    ):
        self.filename = filename
        self.start = start
        self.stop = stop
        self.species = species.upper() if species else None
        self.units = units.upper()
        self.num_pitchangle_bins = num_pitchangle_bins
        self.min_energy = min_energy
        self.max_energy = max_energy
        self.min_pitchangle = min_pitchangle
        self.max_pitchangle = max_pitchangle

        if not self.filename:  # Creates empty JediData object
            return None

       # Read the data file and convert it to HDF5 if necessary.
        if ".hdf" in self.filename.suffix:
            self.read_hdf_file()
        else:
            self.read_raw_csv_file()
            self.convert_to_hdf_file()

        # If we only want to save the data without further processing, exit early.
        if do_save_only:
            return

        # This method processes the metadata to filter out unwanted channels
        # and species based on the arguments provided.
        self._process_metadata()

        self._assign_species_mass()

        # Set all negative values to zero in the raw data.
        self._negative_values_to_zero()

        if 'PSD' in self.units:
            self.calc_phase_space_density()

        # Infer geometric factor using the ratio between count-rate and differential
        # intensity.  Geff = CR / J / Energy
        self.infer_geometric_factor()

        # Infer the measurement uncertainty (sigma) for each channel based on Poisson statistics.
        self.infer_measurement_uncertainty()

        # Split the data into separate DataFrames for each SSD.
        self.split_data_by_ssd()


    def read_raw_csv_file(self):
        """
        Read the JEDI CSV file and extract metadata and raw data.
        """
        with open(self.filename, "r") as f:
            channel = [x.strip() for x in f.readline().split(",")]
            energy_low = [x.strip() for x in f.readline().split(",")]
            energy_high = [x.strip() for x in f.readline().split(",")]
            particle = [x.strip() for x in f.readline().split(",")]

        df = pd.DataFrame(
            {
                "Channel": channel[: len(energy_low)],
                "Elow": energy_low,
                "Ehigh": energy_high,
                "Particle": particle,
            }
        )

        self.metadata = df

        self.data_raw = pd.read_csv(self.filename, names=channel, skiprows=5)

        self.data_raw["time"] = self.data_raw.apply(
            lambda x: dt.datetime.strptime(x["UTC"], "%Y-%m-%dT%H:%M:%S"), axis=1
        )

        self.data_raw.drop(labels="UTC", axis=1, inplace=True)
        self.data_raw.set_index(["time"], inplace=True)
        self.data_raw.sort_index(inplace=True)

    def convert_to_hdf_file(self):
        """
        Convert the current data file to HDF5 format for faster access and storage.
        This will create a new file with the same name but with a .hdf5 extension
        and store the raw data and metadata in it.
        If the file already exists, it will be overwritten.
        """
        hdf5_filename: Path = self.filename.with_suffix(".hdf5")
        self.data_raw.to_hdf(hdf5_filename, "df_raw", mode="a")
        self.metadata.to_hdf(hdf5_filename, "metadata", mode="a")
        self.filename.unlink()
        self.filename = hdf5_filename

    def read_hdf_file(self):
        """
        Read JEDI data from a pre-processed HDF5 file.
        """
        self.data_raw = pd.read_hdf(self.filename, key="df_raw")
        self.metadata = pd.read_hdf(self.filename, key="metadata")

    def _process_metadata(self):
        """
        Process the metadata DataFrame to filter out unwanted channels and
        species.
        """
        df = self.metadata.copy()

        # Filter out unwanted units
        df = df[df["Channel"].str.contains(self.units)]

        # Filter out unwanted species
        df["Particle"] = df["Particle"].str.upper()
        df = df[~df["Particle"].str.contains("NONPARTICLE")]
        df = df[df["Particle"].str.contains(self.species)]

        df["Elow"] = df["Elow"].astype("float")
        df["Ehigh"] = df["Ehigh"].astype("float")


        # Filter unwanted energy channels (should we leave this to the user?)
        # if self.species == "ELECTRON":
        #     for i in range(4):
        #         filter = df["Channel"][i][-7:]  # what is the significance of -7?
        #         df = df[~df["Channel"].str.contains(filter)]

        # if self.species == "PROTON":
        #     for i in range(1):
        #         filter = df["Channel"][i][-7:]
        #         df = df[~df["Channel"].str.contains(filter)]

        # df.reset_index(
        #     inplace=True, drop=True
        # )  # Have to reset index with every df change ??

        # Filter out invalid energy channels based on min and max energy
        if not self.min_energy:
            self.min_energy = df["Elow"].min()
        if not self.max_energy:
            self.max_energy = df["Ehigh"].max()

        invalid_channels = \
            (df["Elow"] < self.min_energy) | (df["Ehigh"] > self.max_energy)
        df = df[~invalid_channels]

        df.reset_index(inplace=True, drop=True)
        self.metadata = df

    def _assign_species_mass(self):
        """
        Assign ion mass and rename some strings for plotting.
        """
        # Assign ion mass and rename some strings for plotting.
        if self.species == r"O\+S" or self.species == "HEAVY":
            self.species = "O & S"
            self.mass = 16 * m_p

        elif self.species == "OXYGEN":
            self.species = "O"
            self.mass = 16 * m_p

        elif self.species == "SULFUR":
            self.species = "S"
            self.mass = 32 * m_p

        elif self.species == "PROTON":
            self.species = "H+"
            self.mass = m_p * 1.0

        else:
            self.species = "e-"
            self.mass = m_e * 1.0

    def _calculate_intensity_to_psd_conversion(self):
        """
        Calculate the conversion factor from intensity (J) to phase space
        density (PSD).
        This is based on the formula:
        f = m / v^2 j
        where:
            - f is the phase space density (PSD) in [s^3 cm^-6],
            - m is the mass of the particle,
            - v is the velocity of the particle,
            - j is the intensity in [s-1 cm-2 eV^-1 sr^-1].
        """
        self.metadata["j2psd"] = (self.mass**2 / e**2) \
            / np.sqrt(self.metadata["Ehigh"] * self.metadata["Elow"]) / (2 * 1e14)

    def _negative_values_to_zero(self):
        """
        Specify negative values in the raw data as 0 or no-count.
        This is necessary because negative values in the data do not make sense
        and can cause issues in further processing.
        """
        # Specify negative values as 0 or no-count
        vals = self.data_raw.loc[:, self.metadata["Channel"]].values
        vals[vals < 0] = 0.0
        self.data_raw.loc[:, self.metadata["Channel"]] = vals

        if len(self.metadata.index) % self.num_ssd == 0:
            self.num_channels = int(len(self.metadata.index) / self.num_ssd)
        else:
            raise Exception("Unequal energy channels within each SSD.")

    def infer_geometric_factor(self):
        """
        Infer the geometric factor for each energy channel using the ratio
        between count-rate and differential intensity.
        Geff = CR / J / Energy
        where:
            - Geff is the geometric factor in [cm^2 sr keV/keV],
            - CR is the count-rate in [s^-1],
            - J is the differential intensity in [s^-1 cm^-2 keV^-1 sr^-1],
            - Energy is the energy in [eV] for the channel.
        """
        self.metadata["Geff"] = np.nan

        for _, row in self.metadata.iterrows():
            channel = [c.replace(self.units, "FLUX") for c in row["Channel"]]
            energy = np.sqrt(row["Elow"] * row["Ehigh"])

            intensities = self.data_raw[channel].values
            count_rates = self.data_raw[channel.replace("FLUX", "CPS")].values

            with np.errstate(divide='ignore', invalid='ignore'):
                geff_values = np.where(
                    intensities > 0,
                    count_rates / (intensities * energy),
                    np.nan)

            # Take the median of the inferred geometric factor, ignoring NaNs
            geff_median = np.nanmedian(geff_values)

            self.metadata.loc[self.metadata["Channel"] == channel, "Geff"] = geff_median

    def infer_measurement_uncertainty(self):
        """
        Infer the measurement uncertainty (sigma) for each channel based on
        Poisson statistics.
        Sigma = sqrt(CR) / (Geff * Energy)
        where:
            - Sigma is the measurement uncertainty in [s^-1 cm^-2 keV^-1 sr^-1],
            - CR is the count-rate in [s^-1],
            - Geff is the geometric factor in [cm^2 sr keV/keV],
            - Energy is the energy in [eV] for the channel.
        """
        sigma_columns = [x.replace(self.units, "SIGMA") for x in self.metadata["Channel"]]

        df = pd.DataFrame(index=self.data_raw.index, columns=sigma_columns)
        for _, row in self.metadata.iterrows():
            channel = row["Channel"].replace(self.units, "FLUX")
            cps_channel = row["Channel"].replace(self.units, "CPS")
            energy = np.sqrt(row["Elow"] * row["Ehigh"])
            geff = row["Geff"]

            count_rates = self.data_raw[cps_channel].values

            with np.errstate(divide='ignore', invalid='ignore'):
                sigma_values = np.where(
                    count_rates > 0,
                    np.sqrt(count_rates) / (geff * energy),
                    np.nan)

            df[row["Channel"].replace(self.units, "SIGMA")] = sigma_values

        self.data_raw = pd.concat([self.data_raw, df], axis=1)

    def calc_phase_space_density(self):
        """
        Add columns of PSD to the data_raw DataFrame. Also change the
        metadata DataFrame to reflect the new channel names.
        """
        self._calculate_intensity_to_psd_conversion()

        self.metadata["ChannelPSD"] = self.metadata.apply(
            lambda x: x["Channel"].replace("FLUX", "PSD"), axis=1
        )

        df = pd.DataFrame(index=self.data_raw.index, columns=self.metadata["ChannelPSD"])
        for _, row in self.metadata.iterrows():
            df[row["ChannelPSD"]] = self.data_raw[row["Channel"]] * row["j2psd"]

        self.data_raw = pd.concat([self.data_raw, df], axis=1)

        self.metadata.drop(labels="Channel", axis=1, inplace=True)
        self.metadata.rename(columns={"ChannelPSD": "Channel"}, inplace=True)

    def split_data_by_ssd(self):
        """
        Split the data into separate DataFrames for each SSD.
        """
        self.ssd_data = []
        self.ssd_sigma = []
        self.ssd_pitchangles = []

        # Create a new DataFrame for each SSD
        for i in range(self.num_ssd):

            ssd_columns = self.metadata.loc[
                self.metadata['channel'].str.startswith(f"T{i}"), 'channel'
            ].tolist()

            if not ssd_columns:
                continue

            ssd_df = self.data_raw[ssd_columns].copy()

            # Add the SSD DataFrame to the list
            try:
                self.ssd_data.append(ssd_df)
                self.ssd_pitchangles.append(self.data_raw['Pitch_Angle_T{:d}'.format(i)])
            except Exception as e:
                print(f"Error adding SSD DataFrame for SSD {i}: {e}")

    def create_energytime_df(self):
        """
        Create the energy-time spectra in the form of a dataframe. Each energy
        channel has Elow and Ehigh values, which can be used to create non-contiguous
        energy bins.

        |El1---Eh1|, |El2---Eh2|, |El3---Eh3|, ...

        This method will create separate energy-time spectra for each SSD.
        The first column of the dataframe will be the instantaneous pitch
        """
<<<<<<< Updated upstream
        return None
        if self.efilter:
            self.energybins = np.exp(np.linspace(
                np.log(self.efilter[0]), np.log(self.efilter[1]),
                self.nEbins + 1))
        else:
            self.energybins = np.exp(np.linspace(
                np.log(self.metadata['Elow'].min()),
                np.log(self.metadata['Ehigh'].max()), self.nEbins + 1))

        self.energytime = pd.DataFrame(
            index=self.df_raw.index,
            columns=['Ebin_{:d}'.format(x) for x in range(self.nEbins)])

        for j in range(self.nEbins):

            meta = self.metadata.copy()

            meta['deltaE'] = meta['Ehigh'] - meta['Elow']

            # Prevalence of each channel in this energy bins.
            meta['dE'] = np.minimum(meta['Ehigh'], self.energybins[j + 1]) - \
                         np.maximum(meta['Elow'],  self.energybins[j])


            meta = meta[~(meta['dE'] < 0)]

            # Now we know the energy channels that are present within the j-th
            # energy bin. Knowing this, we can identify the pitch
            # angle associated with each energy channel.

            denom = np.zeros(len(self.df_raw.index))
            val = np.zeros(len(self.df_raw.index))

            for _, row in meta.iterrows():

                scope = row['channel'][:2]    # Get telescope ID for E channel
                pa = self.df_raw.loc[:, 'Pitch_Angle_{:s}'.format(scope)]\
                    .astype(float).values

                # This is the solid angle of the measurement at a particular
                # pitch angle, given that JEDI-I has an angular resolution of
                # approximately 26.7 degrees.
                deltaSigma = 4 * np.pi * np.sin(pa * np.pi / 180) \
                             * np.sin(13.35 * np.pi / 180)

                # Need to ignore NaN values while summing.
                deltaSigma[np.isnan(deltaSigma)] = 0
                intensities = self.df_raw.loc[:, row['channel']].fillna(0)

                # deltaSigma[intensities == 0] = 0

                val += intensities * row['dE'] # * deltaSigma
                denom += row['dE'] # * deltaSigma

                # How to handle nans?

            condition =  (denom > 0)

            self.energytime.loc[condition, 'Ebin_{:d}'.format(j)] = \
                val[condition] / denom[condition]

            self.energytime.loc[~condition, 'Ebin_{:d}'.format(j)] = np.nan


    # To-do:
    # [ ] Change 'Channel' to contain the telescope ID and remove the inherent units. Various data
    #     processing functions can then add the units as necessary.
    # [ ] Calculate the geometric factor using the count-rate and differential intensity.
    # [ ] Calculate the measurement uncertainty (sigma) using Poisson statistics.
    # [ ] Create energy-time spectrograms for each SSD, and a combined one as well.
    # [ ] Create pitch-angle spectrograms for particular energies for individual SSDs and combined.
    # [ ] Implement a method to visualize the spectrograms.
    # [ ] Create a vPerp-vParallel plot for selected time intervals.
=======
        pass
    
>>>>>>> Stashed changes
