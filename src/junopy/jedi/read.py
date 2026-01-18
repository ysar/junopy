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
        if not filename:  # Creates empty JediData object
            return None

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

    def merge_rebinned_ssd_data(
        self,
        **kwargs,
    ):
        """
        Merge data from all SSDs for this JEDI file.

        Args:
            min_energy (float): Minimum energy to be used across SSDs.
            max_energy (float): Maximum energy to be used across SSDs.
            **kwargs: Arguments passed to `JediSsdData.rebin_data` for each SSD.
        """
        min_energy = kwargs.pop("min_energy", None)
        max_energy = kwargs.pop("max_energy", None)
        if not min_energy or not max_energy:
            # Identify new energy bins based on minimum and maximum energies for each SSD. This is
            # necessary otherwise `rebin_data` would provide different energy bins for each SSD.
            min_energy = 1.0e8
            max_energy = 0.01

            for i_ssd in range(self.num_ssd):
                min_energy = np.minimum(
                    min_energy, self.ssd_data[i_ssd].metadata["Elow"].min()
                )
                max_energy = np.maximum(
                    max_energy, self.ssd_data[i_ssd].metadata["Ehigh"].max()
                )

        for ssd in self.ssd_data:
            ssd.rebin_data(min_energy=min_energy, max_energy=max_energy, **kwargs)

        self.energy_bins = self.ssd_data[0].energy_bins
        self.pitchangle_bins = self.ssd_data[0].pitchangle_bins
        self.num_energy_bins = self.ssd_data[0].num_energy_bins
        self.num_pitchangle_bins = self.ssd_data[0].num_pitchangle_bins

        self.data_rebinned = _merge_dataframe(
            [ssd.data_rebinned for ssd in self.ssd_data]
        )

        self.energy_spectra = _calc_energy_spectra(
            self.data_rebinned,
            self.pitchangle_bins,
            self.num_energy_bins,
            self.num_pitchangle_bins,
            self.ssd_data[0].weigh_solid_angle,
        )

        self.pitchangle_spectra = _calc_pitchangle_spectra(
            self.data_rebinned,
            self.energy_bins,
            self.num_energy_bins,
            self.num_pitchangle_bins,
            self.ssd_data[0].weigh_energy,
        )


def merge_multiple_jedi_data(
    jedi_data_list,
    **kwargs,
):
    """
    Merge data from several JEDI data files.

    Args:
        jedi_data_list (list[JediData]): List of different JEDI data files corresponding to
            different JEDI sensors (090, 180, 270)
        **kwargs: Arguments passed to `JediSsdData.rebin_data` for each SSD.
    """
    min_energy = kwargs.pop("min_energy", None)
    max_energy = kwargs.pop("max_energy", None)
    if not min_energy or not max_energy:
        # Identify new energy bins based on minimum and maximum energies for each SSD. This is
        # necessary otherwise `rebin_data` would provide different energy bins for each SSD.
        min_energy = 1.0e8
        max_energy = 0.01

        for jedi_data in jedi_data_list:
            for i_ssd in range(jedi_data.num_ssd):
                min_energy = np.minimum(
                    min_energy, jedi_data.ssd_data[i_ssd].metadata["Elow"].min()
                )
                max_energy = np.maximum(
                    max_energy, jedi_data.ssd_data[i_ssd].metadata["Ehigh"].max()
                )

    rebinned_data_list = []
    for jedi_data in jedi_data_list:
        for ssd in jedi_data.ssd_data:
            ssd.rebin_data(min_energy=min_energy, max_energy=max_energy, **kwargs)
        rebinned_data_list += [ssd.data_rebinned for ssd in jedi_data.ssd_data]

    merged = JediData()
    merged.start = jedi_data_list[0].start
    merged.stop = jedi_data_list[0].stop
    merged.species = jedi_data_list[0].species
    merged.units = jedi_data_list[0].units
    merged.energy_bins = jedi_data_list[0].ssd_data[0].energy_bins
    merged.pitchangle_bins = jedi_data_list[0].ssd_data[0].pitchangle_bins
    merged.num_energy_bins = jedi_data_list[0].ssd_data[0].num_energy_bins
    merged.num_pitchangle_bins = jedi_data_list[0].ssd_data[0].num_pitchangle_bins
    weigh_energy = jedi_data_list[0].ssd_data[0].weigh_energy
    weigh_solid_angle = jedi_data_list[0].ssd_data[0].weigh_solid_angle

    merged.data_rebinned = _merge_dataframe(rebinned_data_list)
    merged.energy_spectra = _calc_energy_spectra(
        merged.data_rebinned,
        merged.pitchangle_bins,
        merged.num_energy_bins,
        merged.num_pitchangle_bins,
        weigh_solid_angle,
    )

    merged.pitchangle_spectra = _calc_pitchangle_spectra(
        merged.data_rebinned,
        merged.energy_bins,
        merged.num_energy_bins,
        merged.num_pitchangle_bins,
        weigh_energy,
    )
    return merged


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
        Store data for each solid-state detector (SSD) separately.

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

    def rebin_data(
        self,
        min_energy: float = np.nan,
        max_energy: float = np.nan,
        num_energy_bins: int = np.nan,
        min_pitchangle: float = 0.0,
        max_pitchangle: float = 180.0,
        num_pitchangle_bins: int = 18,
        weigh_solid_angle: bool = False,
        weigh_energy: bool = False,
        pitchangle_mask_type: str = "min-to-max",
        cadence: dt.timedelta = dt.timedelta(seconds=1),
    ):
        """
        Re-bin the data into [time, energy, pitchangle] bins.  Store in a pandas DataFrame with
        multiple columns.

        For each energy bin (i) and pitch-angle bin (j),
            - Identify the support of the SSD measurement in that bin (extension from 1d to 2d).
                In 1D the support is defined as `s = min(E_h, E_i+1) - max(E_l, E_i)`.  We extend
                this to 2D by consider different supports for energy and pitch-angle and using the
                product as the support for 2D averaging.
            -

        Args:
            min_energy (float): Minimum energy to consider.
            max_energy (float): Maximum energy to consider.
            num_energy_bins (int): Number of energy bins.
            min_pitchangle (float): Minimum pitch-angle to consider.
            max_pitchangle (float): Maximum pitch-angle to consider.
            num_pitchangle_bins (int): Number of pitch-angle bins.
            weigh_solid_angle (bool): Whether to use solid-angle weights during averaging
                (True) or use unweighted averages (False).
            weigh_energy (bool): Whether to use energy weights during averaging (True) or use
                unweighted averages (False).
            pitchangle_mask_type (str): Type of pitch-angle mask (Options - 'min-to-max', 'lower-losscone', 'upper-losscone')
        """

        self.weigh_solid_angle = weigh_solid_angle
        self.weigh_energy = weigh_energy

        if np.isnan(min_energy):
            min_energy = self.metadata["Elow"].min()
        if np.isnan(max_energy):
            max_energy = self.metadata["Ehigh"].max()
        if np.isnan(num_energy_bins):
            num_energy_bins = len(self.metadata["Channel"])

        # Create energy bins
        self.num_energy_bins = num_energy_bins
        self.energy_bins = np.geomspace(
            min_energy,
            max_energy,
            self.num_energy_bins + 1,
        )

        self.num_pitchangle_bins = num_pitchangle_bins
        self.pitchangle_bins = np.linspace(
            min_pitchangle,
            max_pitchangle,
            self.num_pitchangle_bins + 1,
        )
        mid_pitchangles = 0.5 * (self.pitchangles_high + self.pitchangles_low)
        times = self.data.index
        data = self.data.copy().fillna(0.0)

        # Mask data correponding to unnecessary pitch-angles
        pitchangle_mask, pitchangles_low, pitchangles_high = self._mask_pitch_angles(
            mask_type=pitchangle_mask_type,
            min_pitchangle=min_pitchangle,
            max_pitchangle=max_pitchangle,
        )
        data.loc[~pitchangle_mask, :] = np.nan

        col_names: list[str] = []
        for i in range(self.num_energy_bins):
            col_names.extend(
                [f"energy_{i}_pitchangle_{j}_" for j in range(self.num_pitchangle_bins)]
            )

        df = pd.DataFrame(index=times, columns=col_names, data=np.nan)

        # We calculate the pitch-angle support for this SSD.  This is time-variable, so we
        # store it as a pd.DataFrame with index the same as self.pitchangle_low and self.pitchangle_high
        # and columns as different user-defined pitch-angle bins.  We can also do this beforehand in
        # a vectorized manner.
        support_pitchangle = pd.DataFrame(
            index=times,
            columns=[f"PA_{j}" for j in range(self.num_pitchangle_bins)],
            data=0.0,
        )

        for j in range(self.num_pitchangle_bins):
            support_pitchangle.loc[:, f"PA_{j}"] = np.minimum(
                pitchangles_high, self.pitchangle_bins[j + 1]
            ) - np.maximum(pitchangles_low, self.pitchangle_bins[j])

        # Averaging loop - outer over energies
        for i in range(self.num_energy_bins):  # For each user-specified bin
            # First we have to identify the native energy bins that are valid for this user-defined
            # energy bin.
            meta = self.metadata.copy()
            meta["bandwidth"] = meta["Ehigh"] - meta["Elow"]

            meta["support"] = np.minimum(
                meta["Ehigh"], self.energy_bins[i + 1]
            ) - np.maximum(meta["Elow"], self.energy_bins[i])

            meta = meta[~(meta["support"] < 0)]

            # Averaging loop - inner over pitch angles
            for j in range(self.num_pitchangle_bins):
                # The final loop over times is vectorized.
                denom = np.zeros(len(times))
                val = np.zeros(len(times))

                for _, row in meta.iterrows():  # For each energy channel
                    if weigh_solid_angle:
                        # This is the solid angle of the measurement at a particular
                        # pitch-angle.
                        weight_pitchangle = (
                            2
                            * np.pi
                            * np.sin(mid_pitchangles)
                            * np.sin(0.5 * support_pitchangle.loc[:, f"PA_{j}"])
                        )
                    else:
                        weight_pitchangle = support_pitchangle.loc[:, f"PA_{j}"]

                    if weigh_energy:
                        weight_energy = row["support"] * row["bandwidth"]
                    else:
                        weight_energy = row["support"]

                    # Need to ignore NaN values while summing.
                    weight_pitchangle[np.isnan(weight_pitchangle)] = 0.0
                    if np.isnan(weight_energy):
                        weight_energy = 0.0

                    col_name = " ".join([row["Channel"], self.units])
                    val += (
                        data.loc[:, col_name].values.T
                        * weight_pitchangle
                        * weight_energy
                    )
                    denom += weight_pitchangle * weight_energy

                # End of energy channel loop

                condition = denom > 0

                df.loc[condition, f"energy_{i}_pitchangle_{j}_"] = (
                    val[condition] / denom[condition]
                )
                # df.loc[~condition, f"energy_{i}_pitchangle_{j}_"] = np.nan
            # End of user-defined pitch-angle bins loop.
        # End of user-defined energy bins loop.

        self.data_rebinned = df.resample(cadence).mean()

    def calc_energy_spectra(self):
        self.energy_spectra = _calc_energy_spectra(
            self.data_rebinned,
            self.pitchangle_bins,
            self.num_energy_bins,
            self.num_pitchangle_bins,
            self.weigh_solid_angle,
        )

    def calc_pitchangle_spectra(self):
        self.pitchangle_spectra = _calc_pitchangle_spectra(
            self.data_rebinned,
            self.energy_bins,
            self.num_energy_bins,
            self.num_pitchangle_bins,
            self.weigh_energy,
        )

    def _mask_pitch_angles(
        self,
        mask_type: str = "min-to-max",
        min_pitchangle: float = 0.0,
        max_pitchangle: float = 180.0,
    ):
        """
        Identify pitch-angles to be masked during the energy spectra calculations.

        Args:
            mask_type (str): Type of pitch-angle mask (Options - 'min-to-max', 'lower-losscone', 'upper-losscone')
            min_pitchangle (float): Minimum pitch-angle to consider.
            max_pitchangle (float): Maximum pitch-angle to consider.
        """
        pitchangles_low = self.pitchangles_low.copy()
        pitchangles_high = self.pitchangles_high.copy()

        times = self.pitchangles_low.index

        if mask_type == "min-to-max":
            min_pa = pd.Series(index=times, data=min_pitchangle)
            max_pa = pd.Series(index=times, data=max_pitchangle)

        # Lower loss cone is from PA=0 to PA=LC
        elif mask_type == "lower-losscone":
            juno_trace = utils.field.JunoFieldTrace(self.data.index)
            min_pa = pd.Series(index=times, data=0.0)
            max_pa = pd.Series(index=times, data=juno_trace.df["losscone_local"])

        # Upper loss cone is from PA=(180-LC) to PA=180
        elif mask_type == "upper-losscone":
            juno_trace = utils.field.JunoFieldTrace(self.data.index)
            min_pa = pd.Series(
                index=times, data=180.0 - juno_trace.df["losscone_local"]
            )
            max_pa = pd.Series(index=times, data=180.0)

        else:
            raise Exception("Unknown mask type for energy spectra calculation.")

        # Edge case when the lower PA is outside bounds, but not the upper
        condition = (pitchangles_low < min_pa) & (pitchangles_high > min_pa)
        pitchangles_low[condition] = min_pa[condition]

        # Edge case when the upper PA is outside bounds, but not the lower
        condition = (pitchangles_low < max_pa) & (pitchangles_high > max_pa)
        pitchangles_high[condition] = max_pa[condition]

        # Actual mask for the majority of the data
        condition = (pitchangles_low >= min_pa) & (pitchangles_high <= max_pa)

        return condition, pitchangles_low, pitchangles_high


# When aggregating in pandas, we want to avoid counting nans as zeros.
def _weighted_nanmean(x, weights):
    _valid_indices = ~(np.isnan(x) | np.isnan(weights))
    return np.average(x[_valid_indices], weights=weights[_valid_indices])


def _calc_energy_spectra(
    data: pd.DataFrame,
    pitchangle_bins: np.ndarray,
    num_energy_bins: int,
    num_pitchangle_bins: int,
    weigh_solid_angle: bool = False,
):
    """
    Calculate the energy spectra based on the rebinned data.
    """
    energy_cols = [f"energy_{i}" for i in range(num_energy_bins)]
    df = pd.DataFrame(index=data.index, columns=energy_cols, data=np.nan)

    if weigh_solid_angle:
        pitchangle_delta = pitchangle_bins[1:] - pitchangle_bins[:-1]
        pitchangle_mid = 0.5 * (pitchangle_bins[1:] + pitchangle_bins[:-1])
        weights = 2 * np.pi * np.sin(pitchangle_mid) * np.sin(0.5 * pitchangle_delta)
    else:
        weights = np.full(num_pitchangle_bins, 1.0)

    for i in range(num_energy_bins):
        df.iloc[:, i] = (
            data.filter(like=f"energy_{i}_", axis=1)
            .fillna(0.0)
            .agg(lambda x: _weighted_nanmean(x, weights), axis=1)
        )
    return df


def _calc_pitchangle_spectra(
    data: pd.DataFrame,
    energy_bins: np.ndarray,
    num_energy_bins: int,
    num_pitchangle_bins: int,
    weigh_energy: bool = False,
):
    """
    Calculate the pitch-angle spectra based on the rebinned data.
    """
    pitchangle_cols = [f"pitchangle_{j}_" for j in range(num_pitchangle_bins)]
    df = pd.DataFrame(index=data.index, columns=pitchangle_cols, data=np.nan)

    if weigh_energy:
        weights = energy_bins[1:] - energy_bins[:-1]
    else:
        weights = np.full(num_energy_bins, 1.0)

    for j in range(num_pitchangle_bins):
        df.iloc[:, j] = (
            data.filter(like=f"pitchangle_{j}_", axis=1)
            .fillna(0.0)
            .agg(lambda x: _weighted_nanmean(x, weights), axis=1)
        )

    return df


def _merge_dataframe(df_list: list[pd.DataFrame]):
    """
    Merge multiple time-series DataFrame together, assuming that they share
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
