# Copyright (c) 2018 Yash Sarkango
# Licensed under the MIT License

import datetime as dt
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.constants import e, m_e, m_p

from junopy import spice, utils

JEDI_ANGULAR_RESOLUTION = 26.7  # degrees


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
        filename: Path | None = None,
        start: dt.datetime | None = None,
        stop: dt.datetime | None = None,
        species: str | None = None,
        units: str = "FLUX",
        # num_pitchangle_bins: int = 18,
        do_save_only: bool = False,
        # min_energy: float | None = None,
        # max_energy: float | None = None,
        # min_pitchangle: float = 0,
        # max_pitchangle: float = 180,
    ):
        self.filename: Path = filename
        self.start: dt.datetime = start
        self.stop: dt.datetime = stop
        self.species: str | None = species.upper() if species else None
        self.units: str = units.upper()
        # self.num_pitchangle_bins: int = num_pitchangle_bins
        # self.min_energy: float | None = min_energy
        # self.max_energy: float | None = max_energy
        # self.min_pitchangle: float = min_pitchangle
        # self.max_pitchangle: float = max_pitchangle

        if not self.filename:  # Creates empty JediData object
            return None

        # Read the data file and convert it to HDF5 if necessary.
        if ".hdf" in self.filename.suffix:
            self.read_hdf_file()
        else:
            self.read_raw_csv_file()
            self.convert_to_hdf_file()

        if do_save_only:
            return

        if start and stop:
            self.data_raw = self.data_raw[start:stop].copy()

        # This method processes the metadata to filter out unwanted channels
        # and species based on the arguments provided.
        self._process_metadata()

        self._assign_species_mass()

        # Set all negative values to zero in the raw data.
        self._negative_values_to_zero()

        if "PSD" in self.units:
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
        Reads the JEDI CSV file and extracts metadata and raw data.

        Notes:
            - Metadata is extracted from the first four lines of the CSV.
            - Raw data is loaded into `self.data_raw` as a DataFrame.
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

        # Remove the unit from the channel name and remove rows that do not
        # correspond to an energy channel.
        df[["Channel", "Unit"]] = df["Channel"].str.split(" ", expand=True)
        df = df[~df["Unit"].isnull()]
        df = df[df["Unit"].str.contains("COUNTS")]
        df.drop(labels="Unit", axis=1, inplace=True)

        self.metadata: pd.DataFrame = df

        self.data_raw: pd.DataFrame = pd.read_csv(
            self.filename, names=channel, skiprows=5
        )

        times = self.data_raw["ETLookDirection"].values
        times = spice.convert_to_utc_v(times)

        self.data_raw["time"] = [
            dt.datetime.strptime(t, "%Y-%jT%H:%M:%S.%f") for t in times
        ]
        self.data_raw.set_index(["time"], inplace=True)
        self.data_raw.sort_index(inplace=True)

    def convert_to_hdf_file(self):
        """
        Converts the current data file to HDF5 format for faster access and storage.

        - Creates a new file with the same name but with a `.hdf5` extension.
        - Stores both raw data and metadata in the HDF5 file.
        - Overwrites the file if it already exists.
        """
        hdf5_filename: Path = self.filename.with_suffix(".hdf5")
        self.data_raw.to_hdf(hdf5_filename, key="df_raw", mode="a")
        self.metadata.to_hdf(hdf5_filename, key="metadata", mode="a")
        self.filename.unlink()
        self.filename = hdf5_filename

    def read_hdf_file(self):
        """
        Reads JEDI data from a pre-processed HDF5 file.

        Loads both raw data and metadata into their respective attributes.
        """
        self.data_raw: pd.DataFrame = pd.read_hdf(self.filename, key="df_raw")
        self.metadata: pd.DataFrame = pd.read_hdf(self.filename, key="metadata")

    def _process_metadata(self):
        """
        Processes the metadata DataFrame to filter out unwanted channels and species.

        - Filters out non-particle entries and unwanted species.
        - Converts energy columns to float.
        - Applies energy range filters if specified.
        """
        df = self.metadata.copy()

        # Filter out unwanted species
        df["Particle"] = df["Particle"].str.upper()
        df = df[~df["Particle"].str.contains("NONPARTICLE")]
        df = df[df["Particle"].str.contains(self.species)]

        df["Elow"] = df["Elow"].astype("float")
        df["Ehigh"] = df["Ehigh"].astype("float")

        # Filter out invalid energy channels based on min and max energy
        # if not self.min_energy:
        #     self.min_energy = df["Elow"].min()
        # if not self.max_energy:
        #     self.max_energy = df["Ehigh"].max()

        # invalid_channels = (df["Elow"] < self.min_energy) | (
        #     df["Ehigh"] > self.max_energy
        # )

        # df = df[~invalid_channels]

        df.reset_index(inplace=True, drop=True)
        self.metadata = df

    def _assign_species_mass(self):
        """
        Assigns ion mass and renames species strings for plotting and calculations.

        Sets the `self.mass` attribute based on the selected species.
        """
        if self.species == r"O\+S" or self.species == "HEAVY":
            self.species = "O & S"
            self.mass = 16.0 * m_p

        elif self.species == "OXYGEN":
            self.species = "O"
            self.mass = 16.0 * m_p

        elif self.species == "SULFUR":
            self.species = "S"
            self.mass = 32.0 * m_p

        elif self.species == "PROTON":
            self.species = "H+"
            self.mass = m_p * 1.0

        else:
            self.species = "e-"
            self.mass = m_e * 1.0

    def _calculate_intensity_to_psd_conversion(self):
        """
        Calculates the conversion factor from intensity (J) to phase space density (PSD).

        Formula:
            f = m / v^2 * j

        Where:
            - f: phase space density (PSD) [s^3 cm^-6]
            - m: mass of the particle
            - v: velocity of the particle
            - j: intensity [s^-1 cm^-2 keV^-1 sr^-1]
        """
        self.metadata["j2psd"] = (
            (self.mass**2 / e**2)
            / np.sqrt(self.metadata["Ehigh"] * self.metadata["Elow"])
            / (2 * 1e14)
        )

    def _negative_values_to_zero(self):
        """
        Sets negative values in the raw data to zero (no-count).

        This ensures that all data values are physically meaningful and prevents issues in further processing.
        """
        channel_names = [
            " ".join([channel, self.units]) for channel in self.metadata["Channel"]
        ]

        # Specify negative values as 0 or no-count
        vals = self.data_raw.loc[:, channel_names].values
        vals[vals < 0] = 0.0
        self.data_raw.loc[:, channel_names] = vals

        if len(self.metadata.index) % self.num_ssd == 0:
            self.num_channels = len(self.metadata.index) // self.num_ssd
        else:
            raise Exception("Unequal energy channels between SSDs.")

    def infer_geometric_factor(self):
        """
        Infers the geometric factor for each energy channel using the ratio between count-rate and differential intensity.

        Formula:
            Geff = CR / (J * Energy)

        Where:
            - Geff: geometric factor [cm^2 sr keV/keV]
            - CR: count-rate [s^-1]
            - J: differential intensity [s^-1 cm^-2 keV^-1 sr^-1]
            - Energy: channel energy [keV]
        """
        self.metadata["Geff"] = (
            np.nan
        )  # There is only one value of Geff per energy channel

        for _, row in self.metadata.iterrows():
            channel_name = row["Channel"] + " FLUX"
            energy = np.sqrt(row["Elow"] * row["Ehigh"])

            intensities = self.data_raw[channel_name].values
            count_rates = self.data_raw[row["Channel"] + " CPS"].values

            with np.errstate(divide="ignore", invalid="ignore"):
                geff_values = np.where(
                    intensities > 0, count_rates / (intensities * energy), np.nan
                )

            # Take the median of the inferred geometric factor, ignoring NaNs
            row["Geff"] = np.nanmedian(geff_values)

    def infer_measurement_uncertainty(self):
        """
        Infers the measurement uncertainty (sigma) for each channel based on Poisson statistics.

        Formula:
            Sigma = sqrt(CR) / (Geff * Energy)

        Where:
            - Sigma: measurement uncertainty [s^-1 cm^-2 keV^-1 sr^-1]
            - CR: count-rate [s^-1]
            - Geff: geometric factor [cm^2 sr keV/keV]
            - Energy: channel energy [eV]
        """
        sigma_columns = [x + " SIGMA" for x in self.metadata["Channel"]]

        df = pd.DataFrame(index=self.data_raw.index, columns=sigma_columns)
        for _, row in self.metadata.iterrows():
            cps_channel = row["Channel"] + " CPS"
            energies = np.sqrt(row["Elow"] * row["Ehigh"])
            geff = row["Geff"]

            count_rates = self.data_raw[cps_channel].values

            with np.errstate(divide="ignore", invalid="ignore"):
                sigma_values = np.where(
                    count_rates > 0, np.sqrt(count_rates) / (geff * energies), np.nan
                )

            df[row["Channel"] + " SIGMA"] = sigma_values

        self.data_raw = pd.concat([self.data_raw, df], axis=1)

    def calc_phase_space_density(self):
        """
        Adds columns of phase space density (PSD) to the `data_raw` DataFrame.

        Also updates the metadata DataFrame to reflect the new channel names.
        """
        self._calculate_intensity_to_psd_conversion()

        channels_psd = [x + " PSD" for x in self.metadata["Channel"]]

        df = pd.DataFrame(index=self.data_raw.index, columns=channels_psd)
        for _, row in self.metadata.iterrows():
            df[row["Channel"] + " PSD"] = (
                self.data_raw[row["Channel"] + " FLUX"] * row["j2psd"]
            )

        self.data_raw = pd.concat([self.data_raw, df], axis=1)

    def split_data_by_ssd(self):
        """
        Splits the data into separate DataFrames for each SSD (Solid State Detector).

        Populates:
            - `self.ssd_data`: list of DataFrames for each SSD
            - `self.ssd_pitchangles`: list of pitch angle Series for each SSD
        """
        self.ssd_data = []

        # Create a new DataFrame for each SSD
        for i in range(self.num_ssd):
            ssd_metadata = self.metadata.loc[
                self.metadata["Channel"].str.startswith(f"T{i}"), :
            ]

            ssd_channels = ssd_metadata["Channel"].to_list()
            if not ssd_channels:
                continue

            ssd_data_columns = [" ".join([x, self.units]) for x in ssd_channels]
            ssd_sigma_columns = [x + " SIGMA" for x in ssd_channels]

            ssd_df = self.data_raw[ssd_data_columns].copy()
            ssd_sigma_df = self.data_raw[ssd_sigma_columns].copy()
            ssd_pitchangles = self.data_raw[f"Pitch_Angle_T{i}"]

            # Filter pitch-angles
            # condition = (ssd_pitchangles >= self.min_pitchangle) & (
            #     ssd_pitchangles <= self.max_pitchangle
            # )
            # ssd_df.loc[~condition, :] = np.nan
            # ssd_sigma_df.loc[~condition, :] = np.nan
            # ssd_pitchangles.loc[~condition] = np.nan

            self.ssd_data.append(
                JediSsdData(
                    ssd_metadata,
                    ssd_df,
                    ssd_sigma_df,
                    ssd_pitchangles,
                    units=self.units,
                )
            )


class JediSsdData:
    def __init__(
        self,
        metadata,
        data,
        sigma,
        pitchangles,
        units,
    ):
        """
        Store data for each SSD separately.

        Attributes:
            metadata (pd.DataFrame): Contains information for the channel name, min. and max.
                energies, and other metadata related to each SSD channel.
            data (pd.DataFrame): All data corresponding to this SSD in the user-specified units.
            sigma (pd.DataFrame): Measurement uncertainty at each time and energy channel.
            pitchangles_low (pd.Series): Lower limit of the pitchangles measured by that SSD.
            pitchangles_high (pd.Series): Upper limit for the pitchangles measured by that SSD.
        """
        self.metadata = metadata
        self.data = data
        self.sigma = sigma

        pitchangles_low = pitchangles - JEDI_ANGULAR_RESOLUTION * 0.5
        pitchangles_high = pitchangles + JEDI_ANGULAR_RESOLUTION * 0.5
        pitchangles_low[pitchangles_low < 0.0] = 0.0
        pitchangles_high[pitchangles_high > 180.0] = 180.0

        self.pitchangles_low = pitchangles_low
        self.pitchangles_high = pitchangles_high
        self.units = units

    def rebin_energy_pitchangle(
        self,
        min_energy: float = np.nan,
        max_energy: float = np.nan,
        min_pitchangle: float = 0.0,
        max_pitchangle: float = 180.0,
        num_energy_bins: int = np.nan,
        num_pitchangle_bins: int = 18,
        weigh_solid_angle: bool = False,
        mask_type: str = "min-to-max",
    ):
        """ """
        if np.isnan(min_energy):
            min_energy = self.metadata["Elow"].min()
        if np.isnan(max_energy):
            max_energy = self.metadata["Ehigh"].max()
        if np.isnan(num_energy_bins):
            num_energy_bins = len(self.metadata["Channel"])

        self.min_energy = min_energy
        self.max_energy = max_energy
        self.num_energy_bins = num_energy_bins
        self.energy_bins = np.geomspace(
            self.min_energy,
            self.max_energy,
            self.num_energy_bins + 1,
        )

        if mask_type == "min-to-max":
            self.num_pitchangle_bins = num_pitchangle_bins

            self.pitchangle_bins = np.linspace(
                min_pitchangle,
                max_pitchangle,
                self.num_pitchangle_bins + 1,
            )

            min_pa = pd.Series(index=self.pitchangles_low.index, data=min_pitchangle)
            max_pa = pd.Series(index=self.pitchangles_low.index, data=max_pitchangle)

        else:
            juno_trace = utils.field.JunoFieldTrace(self.data.index)

            # Lower loss cone is from PA=0 to PA=LC
            if mask_type == "lower-losscone":
                min_pa = pd.Series(index=self.data.index, data=0.0)
                max_pa = pd.Series(
                    index=self.data.index, data=juno_trace.df["losscone_local"]
                )

            # Upper loss cone is from PA=(180-LC) to PA=180
            if mask_type == "upper-losscone":
                min_pa = pd.Series(
                    index=self.data.index, data=180.0 - juno_trace.df["losscone_local"]
                )
                max_pa = pd.Series(index=self.data.index, data=180.0)

        self.mask_pitch_angles(min_pa, max_pa)
        self.calc_energy_spectra(weigh_solid_angle)
        self.calc_pitch_angle_spectra()

    def mask_pitch_angles(self, min_pa: pd.Series, max_pa: pd.Series):
        """
        Set data corresponding to unnecessary pitch-angles to np.nan
        """
        # Edge case when the lower PA is outside bounds, but not the upper
        condition = (self.pitchangles_low < min_pa) & (self.pitchangles_high > min_pa)
        self.pitchangles_low[condition] = min_pa[condition]

        # Edge case when the lower PA is outside bounds, but not the upper
        condition = (self.pitchangles_low < max_pa) & (self.pitchangles_high > max_pa)
        self.pitchangles_high[condition] = max_pa[condition]

        # Actual mask for the majority of the data
        condition = (self.pitchangles_low >= min_pa) & (self.pitchangles_high <= max_pa)

        self.data.loc[~condition, :] = np.nan

    def calc_energy_spectra(self, weigh_solid_angle):
        times = self.data.index
        data = self.data.fillna(0)

        self.energy_spectra = pd.DataFrame(
            index=times,
            columns=[f"energy_bin_{i}" for i in range(self.num_energy_bins)],
        )

        for i in range(self.num_energy_bins):  # For each user-specified bin
            meta = self.metadata.copy()
            meta["bandwidth"] = meta["Ehigh"] - meta["Elow"]

            meta["support"] = np.minimum(
                meta["Ehigh"], self.energy_bins[i + 1]
            ) - np.maximum(meta["Elow"], self.energy_bins[i])

            meta = meta[~(meta["support"] < 0)]

            # Now we know the energy channels that are present within the j-th
            # energy bin. Knowing this, we can identify the pitch
            # angle associated with each energy channel.

            denom = np.zeros(len(times))
            val = np.zeros(len(times))

            for _, row in meta.iterrows():  # For each energy channel
                if weigh_solid_angle:
                    # This is the solid angle of the measurement at a particular
                    # pitch angle, given that JEDI-I has an angular resolution of
                    # approximately 26.7 degrees.
                    delta_sigma = (
                        2
                        * np.pi
                        * (
                            np.cos(self.pitchangles_low * np.pi / 180)
                            - np.cos(self.pitchangles_high * np.pi / 180)
                        )
                    )

                    # Need to ignore NaN values while summing.
                    delta_sigma[np.isnan(delta_sigma)] = 0.0
                else:
                    delta_sigma = 1.0

                col_name = " ".join([row["Channel"], self.units])
                val += row["support"] * data.loc[:, col_name].values.T * delta_sigma
                denom += row["support"] * delta_sigma

            condition = denom > 0

            self.energy_spectra.loc[condition, f"energy_bin_{i}"] = (
                val[condition] / denom[condition]
            )
            self.energy_spectra.loc[~condition, f"energy_bin_{i}"] = np.nan

    def calc_pitch_angle_spectra(self):
        pass

    # To-do:
    # [x] Change 'Channel' to contain the telescope ID and remove the inherent units. Various data
    #     processing functions can then add the units as necessary.
    # [x] Calculate the geometric factor using the count-rate and differential intensity.
    # [x] Calculate the measurement uncertainty (sigma) using Poisson statistics.


def merge_spectra(df_list: list[pd.DataFrame]):
    """
    Merge multiple time-series spectra together, assuming that they share
    index and bins.  The index does not have to be exactly the same but items
    in the index must be present in both time-series. (i.e. its a join).

    Args:
        df_list (list[pd.DataFrame]): List of spectra to be merged.

    Returns:
        df_merged (pd.DataFrame): Returns a merged pandas DataFrame that
            corresponds to the time-series data.
    """
    df_merged = pd.concat(df_list)
    return df_merged.groupby(df_merged.index).mean()
