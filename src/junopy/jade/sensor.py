# Copyright (c) 2018 Yash Sarkango
# Licensed under the MIT License

import datetime as dt
from copy import deepcopy

import numpy as np
import pandas as pd
from scipy.constants import e, m_e

from .. import mag
from . import calib
from .read import JadeBinaryData


def _create_energy_bins(energies: np.ndarray) -> np.ndarray:
    """
    Given an array of energies, create the energy bins, assuming that the
    energies are the geometric mean of the bins.  By definition, the bins are
    contiguous.

    Args:
        energies (np.ndarray): The mean energies (assumed to be geometric mean)
            of the bins. Shape (num_energies,).

    Returns:
        np.ndarray: The energy bins spaced logarithmically with shape
            (num_energies + 1,).
    """
    energy_bins = np.zeros(len(energies) + 1)
    energy_bins[1:-1] = np.sqrt(energies[:-1] * energies[1:])
    energy_bins[0] = energy_bins[1] ** 2 / energy_bins[2]
    energy_bins[-1] = energy_bins[-2] ** 2 / energy_bins[-3]
    return energy_bins


class SensorData(object):
    def __init__(self, df: pd.DataFrame, id: int):
        """
        A class to hold the data from a single JADE sensor or packet-ID.

        Args:
            df (pd.DataFrame): The data frame containing the sensor data.
            id (int): The packet ID of the sensor.

        Attributes:
            id (int): The packet ID of the sensor.
            name (str): The name of the sensor.
            times (np.ndarray): The times of the measurements.
            data_raw (dict): The raw data from the sensor as a dictionary.
            geom_factor_eff (np.ndarray): The effective geometric factor of the sensor.
            accumulation_time (np.ndarray): The accumulation time of the measurements.
            data_rebinned (np.ndarray): The rebinned data.
            uncertainty_rebinned (np.ndarray): The rebinned uncertainty.
            energy_bins (np.ndarray): The energy bins used for rebinning.
            pitchangle_bins (np.ndarray): The pitch angle bins used for rebinning.
            energy_spectra (pd.DataFrame): The energy spectra of the sensor.
            pitchangle_spectra (pd.DataFrame): The pitch angle spectra of the sensor.
        """

        self.id = id

        if self.id == 104:
            self.name = "E060"
        elif self.id == 106:
            self.name = "E180"
        elif self.id == 131:
            self.name = "Species-3"
        elif self.id == 132:
            self.name = "Species-4"
        elif self.id == 133:
            self.name = "Species-5"

        self.times = df.index.to_pydatetime()

        self.data_raw = JadeBinaryData.create_dict_from_df(df)

    def calc_geometric_factor_electrons(self, is_high_res=True):
        """
        Calculate the effective geometric factor for the electron sensors.
        This is done using the calibration functions in calib.py using the
        equations described in the docstring of calib.py.

        Args:
            is_high_res (bool): Whether the data is high resolution or not. If True,
                the anode and deflection corrections are applied. If False, they
                are set to 1.
        """

        delta = self.data_raw["DIM2_ELEVATION"]
        energies = self.data_raw["DIM1_E"]
        times = self.times

        sensor_const = calib._get_sensor_constants(self.id)

        if is_high_res:
            _F = calib._calc_anode_correction_F(sensor_const)[
                np.newaxis, np.newaxis, :
            ]  # shape (nTime, nE, 16)
            _D = calib._calc_deflection_correction_D(
                delta, sensor_const
            )  # shape (nTime, nE, 16)
        else:
            _F = _D = 1.0

        tau = np.zeros(len(times))
        for i in range(len(times)):
            tau[i] = calib._calc_time_correction_tau(times[i], sensor_const)

        tau = tau[:, np.newaxis, np.newaxis]  # shape (nTime, nE, 16)

        _geometric_factor_sim = calib._calc_simulation_geom_factor(energies)
        _x = calib._calc_sensor_dependence_exp(sensor_const, energies)

        # Energy Geometric Factor in units of cm2 sr eV/eV
        self.geom_factor_eff = _geometric_factor_sim * _x * tau * _D * _F

    def calc_geometric_factor_ion(self):
        """
        Calculate the effective geometric factor for the ion sensors.
        This is done using the empirical formula provided in the JADE
        documentation.
        """
        energies = self.data_raw["DIM1_E"]
        self.geom_factor_eff = (
            6.84e-5 + (2.52e-6 + np.log10(energies)) + (6.84e-5 * energies ** (-0.28))
        )

        if self.id == 105 or self.id == 121:
            self.geom_factor_eff *= 12

        self.geom_factor_eff *= 1e-4  # Convert from m2 to cm2

    def convert_units(
        self,
        units_in: str = "countrate",
        units_out: str = "dnf",
        m: float = m_e,
        q: int = 1,
    ):
        """
        Convert the units of the JADE `raw_data` dictionary from "units_in" to
        "units_out". The conversion is done in-place.

        Args:
            units_in (str): The input units of the data. Can be 'countrate' [s-1],
                'def' (differential energy flux) [cm-2 sr-1 s-1],
                'dnf' (differential number flux) [cm-2 sr-1 s-1 keV-1],
                or 'psd' (phase space density) [s3 cm-6].
            units_out (str): The output units of the data. Can be 'def', 'dnf',
                or 'psd'.
            m (float): The mass of the particle in kg. Default is electron mass.
            q (int): The charge of the particle in elementary charges. Default is 1.
        """

        to_convert = ["DATA", "DATA_SIGMA", "BACKGROUND", "BACKGROUND_SIGMA"]

        if self.geom_factor_eff is None:
            raise ValueError(
                "Geometric factor not calculated."
                "Call calc_geometric_factor_electrons() or calc_geometric_factor_ion() first."
            )

        if "countrate" in units_in:
            assert (self.data_raw["DATA_UNITS"] == 2).all()
            if "def" in units_out:
                for k in to_convert:
                    self.data_raw[k] /= self.geom_factor_eff  # cm-2 s-1 sr-1
                self.data_raw["DATA_UNITS"][:] = 6

            elif "dnf" in units_out:
                for k in to_convert:
                    # Differential Intensity [cm-2 s-1 sr-1 keV-1]
                    self.data_raw[k] /= (
                        self.geom_factor_eff * q * self.data_raw["DIM1_E"] * 1e-3
                    )
                self.data_raw["DATA_UNITS"][:] = 7

            elif "psd" in units_out:
                for k in to_convert:
                    # Phase space density [s3 cm-6]
                    speed_squared = 2 * self.data_raw["DIM1_E"] * q * e / m  # m2 s-2
                    self.data_raw[k] /= self.geom_factor_eff * 0.5e8 * speed_squared**2
        else:
            raise NotImplementedError(
                "Only conversion from `countrate` is implemented."
            )

    def calc_pitchangles(self, use_internal_mag: bool = False):
        """
        Calculate the pitch angles of the measurements using either the internal
        magnetometer data (from JADE) or the Juno magnetometer investigation.
        The pitch angle is defined as the angle between the negative of the
        look direction and the magnetic field vector.

        Args:
            use_internal_mag (bool): Whether to use the internal magnetometer
                data from JADE. If False, the Juno MAG data is used.
        """
        elevation = self.data_raw["DIM2_ELEVATION"] * np.pi / 180
        azimuth = self.data_raw["DIM2_AZIMUTH_DESPUN"] * np.pi / 180

        if use_internal_mag:
            magvector = self.data_raw["MAG_VECTOR"]
            Bmag = np.sqrt(np.sum(magvector**2, axis=0))
            Bx = magvector[0, :] / Bmag
            By = magvector[1, :] / Bmag
            Bz = magvector[2, :] / Bmag

        else:
            magdata = mag.load_mag_data(
                self.times[0], self.times[-1], coordsystem="Juno_JADE_despun"
            ).df.loc[:, ["Bx", "By", "Bz", "Bmag"]]

            magdata = magdata.reindex(self.times, method="nearest")

            for c in ["Bx", "By", "Bz"]:
                magdata[c] /= magdata["Bmag"]

            Bx = magdata["Bx"].values
            By = magdata["By"].values
            Bz = magdata["Bz"].values

        # fmt: off
        self.data_raw["PITCH_ANGLE"] = (
            np.arccos(
                -(
                    np.cos(elevation) * np.cos(azimuth) * Bx[:, np.newaxis, np.newaxis]
                    + np.cos(elevation) * np.sin(azimuth) * By[:, np.newaxis, np.newaxis]
                    + np.sin(elevation) * Bz[:, np.newaxis, np.newaxis]
                )
            )
            * 180
            / np.pi
        )
        # fmt: on

    def rebin_sensor_data(
        self,
        min_energy: float | None = None,
        max_energy: float | None = None,
        min_pitchangle: float = 0.0,
        max_pitchangle: float = 180.0,
        num_pitchangle_bins: int = 12,
        add_background: bool = False,
    ):
        """
        Rebin the sensor data from (nTime, nE, nLookDir) -> (nTime, nE, nPA).
        This is useful for averaging later on. This method creates a new data
        array `dataR` and `uncertaintyR` in the object to hold the rebinned
        data. This method works with both LRS data (whose first axis is look-direction)
        and HRS data (whose first axis is anode).

        Args:
            min_energy (float | None): Minimum energy to consider.
            max_energy (float | None): Maximum energy to consider.
            min_pitchangle (float): Minimum pitch angle to consider (default=0.).
            max_pitchangle (float): Maximum pitch angle to consider (default=180.).
            num_pitchangle_bins (int): Number of pitch angle bins.
            add_background (bool): Whether to add background data.
        """
        d = self.data_raw

        energies = d["DIM1_E"]
        data = d["DATA"]
        pitchangles = d["PITCH_ANGLE"]
        uncertainty = d["DATA_SIGMA"]

        if add_background:
            data += d["BACKGROUND"]

        # Set invalid values to zero or NaN
        data[data < -999990] = np.nan
        data[(data < 0) & (data > -999990)] = 0.0

        uncertainty[uncertainty < -999990] = np.nan
        uncertainty[(uncertainty < 0) & (uncertainty > -999990)] = 0.0

        # Filter out the appropriate energies
        _energies = energies[0, :, 0]
        condition = np.full(_energies.shape, True)
        if min_energy:
            condition *= _energies >= min_energy

        if max_energy:
            condition *= _energies <= max_energy

        energies = energies[:, condition, :]
        _energies = _energies[condition]
        data = data[:, condition, :]
        uncertainty = uncertainty[:, condition, :]
        pitchangles = pitchangles[:, condition, :]

        energy_bins = _create_energy_bins(_energies)
        num_energy_bins = len(energy_bins) - 1

        pitchangle_bins = np.linspace(
            min_pitchangle, max_pitchangle, num_pitchangle_bins + 1
        )

        # Create masks -
        # 1. Mask invalid pitch angle values (fortunately this is only on one axis)
        # 2. Mask any NaN values prior to averaging

        num_times = data.shape[0]
        new_shape = (num_times, num_energy_bins, num_pitchangle_bins)
        data_rebinned = np.zeros(new_shape)
        uncertainty_rebinned = np.zeros(new_shape)

        nan_mask = ~np.isnan(data)
        sigma_mask = ~np.isnan(uncertainty)

        for i in range(num_pitchangle_bins):
            pitchangle_mask = (pitchangle_bins[i] <= pitchangles) & (
                pitchangle_bins[i + 1] > pitchangles
            )
            mask = pitchangle_mask & nan_mask
            with np.errstate(divide="ignore", invalid="ignore"):
                data_rebinned[:, :, i] = np.sum(data * mask, axis=-1) / np.sum(
                    mask, axis=-1
                )

                uncertainty_rebinned[:, :, i] = np.sum(
                    uncertainty**2 * mask * sigma_mask, axis=-1
                ) / np.sum(mask * sigma_mask, axis=-1)

        data_rebinned[data_rebinned == 0.0] = np.nan

        uncertainty_rebinned = np.sqrt(uncertainty_rebinned)
        uncertainty_rebinned[uncertainty_rebinned == 0.0] = np.nan

        self.data_rebinned = data_rebinned
        self.uncertainty_rebinned = uncertainty_rebinned

        self.energy_bins = energy_bins
        self.pitchangle_bins = pitchangle_bins

    def calc_energy_spectra(self, weigh_solid_angle: bool = False):
        """
        Given some data that has already been rebinned into (nPA, nE, nTime),
        create the energy spectra (nE, nTime) by averaging over the different
        pitch angles.  Creates a new DataFrame `energy_spectra` in the object to hold
        the energy spectra that has index as time and columns as different energy bins.

        Args:
            weigh_solid_angle (bool): Whether to weigh the pitch angle bins by
                their solid angle. If False, all pitch angle bins are weighed equally.
        """
        pitchangle_bins = self.pitchangle_bins

        if weigh_solid_angle:
            delta_pitchangle = (pitchangle_bins[1] - pitchangle_bins[0]) * np.pi / 180
            mid_pitchangle = (
                0.5 * (pitchangle_bins[1:] + pitchangle_bins[:-1]) * np.pi / 180
            )
            solid_angles = (
                np.sin(mid_pitchangle) * np.sin(delta_pitchangle * 0.5) * 4 * np.pi
            )

            weights = solid_angles
        else:
            weights = np.zeros(len(pitchangle_bins) - 1) + 1.0

        # We do this for numpy broadcasting magic
        weights = weights[np.newaxis, np.newaxis, :]

        nanmask = ~np.isnan(self.data_rebinned)
        self.data_rebinned[~nanmask] = 0.0

        # We use this context as there is no guarantee that nanmask has
        # atleast one True value, which would make the denominator 0.
        # 0 in the denominator would lead to a NaN, which is perfectly
        # fine in our case. Physically, it means that no measurements
        # were made for that particular bin.

        with np.errstate(divide="ignore", invalid="ignore"):
            energy_spectra = np.sum(
                self.data_rebinned * weights * nanmask, axis=-1
            ) / np.sum(weights * nanmask, axis=-1)

        energy_spectra[energy_spectra == 0.0] = np.nan

        self.energy_spectra = pd.DataFrame(
            index=self.times,
            columns=[f"Ebin_{x}" for x in range(len(self.energy_bins) - 1)],
            data=energy_spectra,
        )

    def calc_pitch_angle_spectra(self, weigh_energy: bool = True):
        """
        Given some data that has already been rebinned into (nPA, nE, nTime),
        create the pitch angle spectra. Default weighing along the energy axis is
        the energy bandwidth.  Creates a new DataFrame `pitchangle_spectra`
        in the object to hold the pitch angle spectra that has index as time
        and columns as different pitch angle bins.

        Args:
            weigh_energy (bool): Whether to weigh the energy bins by their
                bandwidth. If False, all energy bins are weighed equally.
        """
        energy_bins = self.energy_bins

        if weigh_energy:
            weights = np.sqrt(energy_bins[1:] * energy_bins[:-1])
        else:
            weights = np.zeros(len(energy_bins) - 1) + 1.0

        # We do this for numpy broadcasting magic
        weights = weights[np.newaxis, :, np.newaxis]

        nanmask = ~np.isnan(self.data_rebinned)
        self.data_rebinned[~nanmask] = 0.0

        with np.errstate(divide="ignore", invalid="ignore"):
            pitchangle_spectra = np.sum(
                self.data_rebinned * weights * nanmask, axis=1
            ) / np.sum(weights * nanmask, axis=1)

        pitchangle_spectra[pitchangle_spectra == 0.0] = np.nan

        self.pitchangle_spectra = pd.DataFrame(
            index=self.times,
            columns=[
                "PAbin_{:d}".format(x) for x in range(len(self.pitchangle_bins) - 1)
            ],
            data=pitchangle_spectra,
        )

    def calc_accumulation_time_seconds(self):
        """
        Calculate the accumulation time in seconds from the DIM0_UTC_UPPER and
        DIM0_UTC_LOWER columns in the raw data. This is done by converting the
        strings to datetime objects and taking the difference. Creates a new
        array `accumulation_time` in the object to hold the accumulation times
        in seconds.
        """
        fmt = "%Y-%jT%H:%M:%S.%f"
        zipped = zip(self.data_raw["DIM0_UTC_UPPER"], self.data_raw["DIM0_UTC_LOWER"])

        self.accumulation_time = np.array(
            [
                (
                    dt.datetime.strptime(x, fmt) - dt.datetime.strptime(y, fmt)
                ).total_seconds()
                for x, y in zipped
            ]
        )

    def valid_times_copy(self, times):
        """
        Create a deep copy of the object with only the data corresponding to the
        given times. This is useful for creating a new object with only the
        data corresponding to a specific times. The times should be a list
        of datetime objects.

        Args:
            times (list of datetime): The times to keep in the new object.
        """
        x = deepcopy(self)

        indices = self.times != self.times
        for t in times:
            indices += t == self.times

        x.times = x.times[indices]

        for key in ["DATA", "DATA_SIGMA", "BACKGROUND", "BACKGROUND_SIGMA"]:
            x.data_raw[key] = x.data_raw[key][indices, :, :]

        x.accumulation_time = x.accumulation_time[indices]

        x.geom_factor_eff = x.geom_factor_eff[indices, :, :]

        try:
            x.energy_spectra = x.energy_spectra.iloc[indices].copy()
            x.pitchangle_spectra = x.pitchangle_spectra.iloc[indices, :].copy()
        except AttributeError:
            pass

        try:
            x.data_rebinned = x.data_rebinned[indices, :, :]
            x.uncertainty_rebinned = x.uncertainty_rebinned[indices, :, :]
        except AttributeError:
            pass  # ??? what is this. Why do we expect AttributeErrors?

        return x
