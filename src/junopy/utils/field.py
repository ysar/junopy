"""
Methods to work with Jupiter's magnetic field models to calculate the
equatorial pitch-angle, M-shell, etc. using iupitermag.
"""

import datetime as dt

import iupitermag as im
import numpy as np
import pandas as pd

from junopy import spice

DEFAULT_INTERNAL_FIELD = "JRM33"
DEFAULT_EXTERNAL_FIELD = "CON2020"


class JunoFieldTrace:
    def __init__(
        self,
        times: list[dt.datetime],
        internal_field: str = DEFAULT_INTERNAL_FIELD,
        external_field: str = DEFAULT_EXTERNAL_FIELD,
    ):
        """
        Trace field lines from Juno's position and calculate M-shell,
        equatorial and polar magnetic field strengths.

        Attributes:
            times (list[dt.datetime]): List of times traversed by Juno in UTC.
            internal_field (str): Internal field identifier passed to iupiterMag.
            external_field (str): External field identifier passed to iupiterMag.
            df (pd.DataFrame): DataFrame containing the M-shell, Bmin, Bmax
                obtained from tracing magnetic field lines.
        """
        self.times = times

        self.internal_field = im.InternalField(internal_field)
        self.external_field = im.CurrentSheetField(external_field)

        self.df = pd.DataFrame(
            index=times, columns=["m_shell", "b_min", "b_max"], data=np.nan
        )

        self._do_trace()
        self._calc_m_shell()
        self._calc_b()
        self._calc_bmin_bmax()
        self._calc_loss_cone()

    def _do_trace(self):
        """
        Trace magnetic field lines from Juno's position at times defined in the
        `times` attribute for this class.
        """
        self.positions, _ = spice.get_juno_position(self.times, "IAU_JUPITER")

        self.traces = im.trace_field_to_planet(
            self.positions, self.internal_field, self.external_field
        )

    def _calc_m_shell(self):
        """
        Calculate Juno's M-shell and populate the `df` attribute for this class.

        Note: This libraries defines M-shell as the maximum radial extent of the field line. This
        is a reasonable assumption given the axial asymmetry of the field.
        """
        m_shell = np.full(len(self.times), np.nan)
        for i in range(len(self.times)):
            rad = np.sqrt(np.sum(self.traces[i] ** 2, axis=1))
            m_shell[i] = rad[np.argmax(rad)]

        self.df["m_shell"] = m_shell

    def _calc_b(self):
        """
        Calculate the magnetic field vector for each point in each trace corresponding to times in
        self.times.
        """
        self.traces_b = []
        for i in range(len(self.times)):
            self.traces_b.append(
                self.internal_field.parmap_calc_field_xyz(self.traces[i])
                + self.external_field.parmap_calc_field_xyz(self.traces[i])
            )

    def _calc_bmin_bmax(self):
        """
        Calculate the minimum and maximum magnetic field along the field line
        that Juno is on at times `self.times`.
        """
        b_min = np.full(len(self.times), np.nan)
        b_max = np.full(len(self.times), np.nan)
        for i in range(len(self.times)):
            b_mag = np.sqrt(np.sum(self.traces_b[i] ** 2, axis=1))
            b_min[i], b_max[i] = np.min(b_mag), np.max(b_mag)

        self.df["b_min"], self.df["b_max"] = b_min, b_max

    def _calc_loss_cone(self):
        """
        Calculate the local and equatorial loss cones using the local magnetic field
        strength and the equatorial and polar magnetic field strength.
        """
        b_local = np.full(len(self.times), np.nan)

        b_local = self.internal_field.parmap_calc_field_xyz(
            self.positions
        ) + self.external_field.parmap_calc_field_xyz(self.positions)

        self.df["b_local"] = np.sqrt(np.sum(b_local**2, axis=1))
        self.df["losscone_local"] = (
            np.arcsin(np.sqrt(self.df["b_local"] / self.df["b_max"])) * 180.0 / np.pi
        )
        self.df["losscone_eq"] = (
            np.arcsin(np.sqrt(self.df["b_min"] / self.df["b_max"])) * 180.0 / np.pi
        )


def get_model_field_at_times(
    times: list[dt.datetime],
    coordsystem="JUNO_JSS",
    internal_field=DEFAULT_INTERNAL_FIELD,
    external_field=DEFAULT_EXTERNAL_FIELD,
) -> np.ndarray:
    """
    Returns the Jovian internal + current sheet field calculated at Juno's position at a given
    list of times.

    Args:
        times (list[dt.datetime]): List of times at which to calculate B.
        coordsystem (str): Coordinate system recognized by SPICE.
        internal_field (str): Internal field (default="jrm33")
        external_field (str): External field (default="con2020")

    Returns:
        b_model (np.ndarray): Magnetic field vector in `coordsystem` frame for each of the
            requested times.
    """
    internal_field = im.InternalField(internal_field)
    external_field = im.CurrentSheetField(external_field)

    juno_position, _ = spice.get_juno_position(times, coordsystem="IAU_JUPITER")

    num_times = len(times)
    if num_times > 20:
        b_model_internal = internal_field.parmap_calc_field_xyz(juno_position)
        b_model_external = external_field.parmap_calc_field_xyz(juno_position)
    else:
        b_model_internal = internal_field.map_calc_field_xyz(juno_position)
        b_model_external = external_field.map_calc_field_xyz(juno_position)

    b_model = b_model_internal + b_model_external
    rotation_matrix = spice.calc_rotation_matrix(times, "IAU_JUPITER", coordsystem)
    return np.matvec(rotation_matrix, b_model)
