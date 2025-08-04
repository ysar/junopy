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
        ssd_data (dict[pd.DataFrame]): Dictionary containing data split by SSD.

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
                "channel": channel[: len(energy_low)],
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
        df = df[df["channel"].str.contains(self.units)]

        # Filter out unwanted species
        df["Particle"] = df["Particle"].str.upper()
        df = df[~df["Particle"].str.contains("NONPARTICLE")]
        df = df[df["Particle"].str.contains(self.species)]

        df["Elow"] = df["Elow"].astype("float")
        df["Ehigh"] = df["Ehigh"].astype("float")


        # Filter unwanted energy channels (should we leave this to the user?)
        # if self.species == "ELECTRON":
        #     for i in range(4):
        #         filter = df["channel"][i][-7:]  # what is the significance of -7?
        #         df = df[~df["channel"].str.contains(filter)]

        # if self.species == "PROTON":
        #     for i in range(1):
        #         filter = df["channel"][i][-7:]
        #         df = df[~df["channel"].str.contains(filter)]

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
        vals = self.data_raw.loc[:, self.metadata["channel"]].values
        vals[vals < 0] = 0.0
        self.data_raw.loc[:, self.metadata["channel"]] = vals

        if len(self.metadata.index) % self.num_ssd == 0:
            self.num_channels = int(len(self.metadata.index) / self.num_ssd)
        else:
            raise Exception("Unequal energy channels within each SSD.")
    
    def calc_phase_space_density(self):
        """
        Add columns of PSD to the data_raw DataFrame. Also change the
        metadata DataFrame to reflect the new channel names.
        """
        self._calculate_intensity_to_psd_conversion()

        self.metadata["channelpsd"] = self.metadata.apply(
            lambda x: x["channel"].replace("FLUX", "PSD"), axis=1
        )

        df = pd.DataFrame(index=self.data_raw.index, columns=self.metadata["channelpsd"])
        for _, row in self.metadata.iterrows():
            df[row["channelpsd"]] = self.data_raw[row["channel"]] * row["j2psd"]

        self.data_raw = pd.concat([self.data_raw, df], axis=1)

        self.metadata.drop(labels="channel", axis=1, inplace=True)
        self.metadata.rename(columns={"channelpsd": "channel"}, inplace=True)

    def split_data_by_ssd(self):
        """
        Split the data into separate DataFrames for each SSD.
        """
        self.ssd_data = {}

        # Create a new DataFrame for each SSD
        for i in range(self.num_ssd):

            ssd_columns = self.metadata.loc[
                self.metadata['channel'].str.startswith(f"T{i}"), 'channel'
            ].tolist()

            if not ssd_columns:
                continue

            ssd_columns.append('Pitch_Angle_T{:d}'.format(i))

            ssd_df = self.df_raw[ssd_columns].copy()

            # Add the SSD DataFrame to the dictionary
            self.ssd_data[f"SSD_{i}"] = ssd_df

    def create_energytime_df(self):
        """
        Create the energy-time spectra in the form of a dataframe. Each energy 
        channel has Elow and Ehigh values, which can be used to create non-contiguous
        energy bins.

        |El1---Eh1|, |El2---Eh2|, |El3---Eh3|, ...
        
        This method will create separate energy-time spectra for each SSD.
        The first column of the dataframe will be the instantaneous pitch 
        """
        pass
    