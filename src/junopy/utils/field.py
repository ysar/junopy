"""
Methods to work with Jupiter's magnetic field models to calculate the
equatorial pitch-angle, M-shell, etc. using JupiterMag.

> Note: To calculate Jupiter's internal and current sheet field, use the appropriate
methods from [`JupiterMag`](https://github.com/mattkjames7/JupiterMag) directly.
This module only provides a simpler interface to calculate higher-level quantities
like M-shell, loss-cone, etc., using `JupiterMag`.
"""

import datetime as dt

import JupiterMag as jm
import numpy as np
import pandas as pd

from junopy import spice

DEFAULT_INTERNAL_FIELD = "jrm33"
DEFAULT_EXTERNAL_FIELD = "Con2020"


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
            internal_field (str): Internal field identifier passed to JupiterMag.
            external_field (str): External field identifier passed to JupiterMag.
            df (pd.DataFrame): DataFrame containing the M-shell, Bmin, Bmax
                obtained from tracing magnetic field lines.
        """
        self.times = times
        self.internal_field = internal_field
        self.external_field = external_field

        jm.Con2020.Config(equation_type="analytic")
        jm.Internal.Config(Model=internal_field, CartesianIn=True, CartesianOut=True)
        self.df = pd.DataFrame(
            index=times, columns=["m_shell", "b_min", "b_max"], data=np.nan
        )

        self._do_trace()
        self._calc_m_shell()
        self._calc_bmin_bmax()
        self._calc_loss_cone()

    def _do_trace(self):
        """
        Trace magnetic field lines from Juno's position at times defined in the
        `times` attribute for this class.
        """
        self.positions, _ = spice.get_juno_position(self.times, "IAU_JUPITER")
        trace_analytic = jm.TraceField(
            self.positions[:, 0],
            self.positions[:, 1],
            self.positions[:, 2],
            IntModel=self.internal_field,
            ExtModel=self.external_field,
        )
        self._trace = trace_analytic.TraceDict()

    def _calc_m_shell(self):
        """
        Calculate Juno's M-shell and populate the `df` attribute for this class.
        """
        m_shell = np.full(len(self.times), np.nan)
        for i in range(len(self.times)):
            b_mag = np.sqrt(
                self._trace["Bx"][i] ** 2
                + self._trace["By"][i] ** 2
                + self._trace["Bz"][i] ** 2
            )
            rad = np.sqrt(
                self._trace["x"][i] ** 2
                + self._trace["y"][i] ** 2
                + self._trace["z"][i] ** 2
            )
            m_shell[i] = rad[np.argmin(b_mag)]

        self.df["m_shell"] = m_shell

    def _calc_bmin_bmax(self):
        """
        Calculate the minimum and maximum magnetic field along the field line
        that Juno is on at times `self.times`.
        """
        b_min = np.full(len(self.times), np.nan)
        b_max = np.full(len(self.times), np.nan)
        for i in range(len(self.times)):
            b_mag = np.sqrt(
                self._trace["Bx"][i] ** 2
                + self._trace["By"][i] ** 2
                + self._trace["Bz"][i] ** 2
            )
            b_min[i], b_max[i] = np.min(b_mag), np.max(b_mag)

        self.df["b_min"], self.df["b_max"] = b_min, b_max

    def _calc_loss_cone(self):
        """
        Calculate the local and equatorial loss cones using the local magnetic field
        strength and the equatorial and polar magnetic field strength.
        """
        b_local = np.full(len(self.times), np.nan)

        for i in range(len(self.times)):
            b_local_internal = jm.Internal.Field(
                self.positions[i, 0],
                self.positions[i, 1],
                self.positions[i, 2],
            )

            b_local_external = jm.Con2020.Field(
                self.positions[i, 0],
                self.positions[i, 1],
                self.positions[i, 2],
            )

            b_local[i] = np.sqrt(
                (b_local_internal[0] + b_local_external[0]) ** 2
                + (b_local_internal[1] + b_local_external[1]) ** 2
                + (b_local_internal[2] + b_local_external[2]) ** 2
            )

        self.df["b_local"] = b_local
        self.df["losscone_local"] = (
            np.arcsin(np.sqrt(self.df["b_local"] / self.df["b_max"])) * 180.0 / np.pi
        )
        self.df["losscone_eq"] = (
            np.arcsin(np.sqrt(self.df["b_min"] / self.df["b_max"])) * 180.0 / np.pi
        )


# def get_model_field_at_times(
#     times: list[dt.datetime],
#     coordsystem="JUNO_JSS",
#     internal_field=DEFAULT_INTERNAL_FIELD,
#     external_field=DEFAULT_EXTERNAL_FIELD,
# ) -> np.ndarray:
#     """
#     Returns the Jovian internal + current sheet field calculated at Juno's position at a given
#     list of times.

#     Args:
#         times (list[dt.datetime]): List of times at which to calculate B.
#         coordsystem (str): Coordinate system recognized by SPICE.
#         internal_field (str): Internal field (default="jrm33")
#         external_field (str): External field (default="con2020")
#     """
#     juno_position, _ = spice.get_juno_position(times, coordsystem="IAU_JUPITER")
#     jm.Con2020.Config(equation_type="analytic")
#     jm.Internal.Config(Model=internal_field, CartesianIn=True, CartesianOut=True)

#     b_model_internal = np.array(
#         jm.Internal.Field(juno_position[:, 0], juno_position[:, 1], juno_position[:, 2])
#     )
#     b_model_external = np.array(
#         jm.Con2020.Field(juno_position[:, 0], juno_position[:, 1], juno_position[:, 2])
#     )

#     b_model = b_model_internal + b_model_external

#     rotation_matrix = spice.calc_rotation_matrix(times, "IAU_JUPITER", coordsystem)
