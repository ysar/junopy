# junopy - A Python API for the Juno spacecraft's particle and fields datasets.
# Copyright (C) 2018  Yash Sarkango

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

import numpy as np
import datetime as dt
import pandas as pd

from junopy import spice
from . import convert
# from . import field


def find_ticks_between(
    start: dt.datetime, end: dt.datetime, nticks: int = 6
) -> tuple[pd.DatetimeIndex, list[str]]:
    """
    Given an axes object make the ticks appropriately for subplots.
    """
    time_array = (
        pd.date_range(start=start, end=end, periods=nticks).to_pydatetime().tolist()
    )

    position, _ = spice.get_juno_position(time_array)

    tickstring = []
    for i, pos in enumerate(position):
        R, LT, Latitude = convert.convert_pos_spherical(pos)

        if i == 0:
            tickstring.append(
                "Time - {:s}".format(time_array[i].strftime("%H:%M:%S"))
                + "\n"
                + "R  - {:7.2f}".format(R)
                + "\n"
                + "LT - {:7.2f}".format(LT)
                + "\n"
                + "Latitude - {:7.2f}".format(Latitude)
            )
        else:
            tickstring.append(
                time_array[i].strftime("%H:%M:%S")
                + "\n"
                + "{:7.2f}".format(R)
                + "\n"
                + "{:7.2f}".format(LT)
                + "\n"
                + "{:7.2f}".format(Latitude)
            )

    return time_array, tickstring


# def find_ticks_between_Mshell(
#     start, end, nticks=6, timeformat="%H:%M:%S", appendlabels=True
# ):
#     """
#     Given an axes object make the ticks appropriately for subplots.
#     We would like to have :
#     M shell, Latitude, SIII-Longitude
#     """
#     time_array = (
#         pd.date_range(start=start, end=end, periods=nticks).to_pydatetime().tolist()
#     )

#     position, _ = spice.get_juno_position(time_array, coordsystem="IAU_JUPITER")
#     position_mag, _ = spice.get_juno_position(time_array, coordsystem="JUNO_MAG_VIP4")

#     tickstring = []
#     for i, pos in enumerate(position):
#         M, _, lon = field.get_M_shell(pos[0], pos[1], pos[2])

#         mlat = (
#             np.arctan2(
#                 position_mag[i, 2],
#                 np.sqrt(position_mag[i, 0] ** 2 + position_mag[i, 1] ** 2),
#             )
#             * 180
#             / np.pi
#         )

#         r = np.sqrt(pos[0] ** 2 + pos[1] ** 2 + pos[2] ** 2)

#         if i == 0 and appendlabels:
#             tickstring.append(
#                 "Time   {:s}".format(time_array[i].strftime(timeformat))
#                 + "\n"
#                 + "R     {:7.2f}".format(r)
#                 + "\n"
#                 + "M     {:7.2f}".format(M)
#                 + "\n"
#                 + "MLat   {:7.2f}".format(mlat)
#                 + "\n"
#                 + "Lon   {:7.2f}".format(lon)
#             )
#         else:
#             tickstring.append(
#                 time_array[i].strftime(timeformat)
#                 + "\n"
#                 + "{:7.2f}".format(r)
#                 + "\n"
#                 + "{:7.2f}".format(M)
#                 + "\n"
#                 + "{:7.2f}".format(mlat)
#                 + "\n"
#                 + "{:7.2f}".format(lon)
#             )

#     return time_array, tickstring


# def find_ticks_between_deltaAlfven(
#     start, end, satellite="IO", nticks=6, timeformat="%H:%M:%S"
# ):
#     """
#     Given an axes object make the ticks appropriately for subplots.
#     We would like to have :
#     M shell, Latitude, SIII-Longitude
#     """
#     time_array = (
#         pd.date_range(start=start, end=end, periods=nticks).to_pydatetime().tolist()
#     )

#     position, _ = spice.get_juno_position(time_array, coordsystem="IAU_JUPITER")
#     position_mag, _ = spice.get_juno_position(time_array, coordsystem="JUNO_MAG_VIP4")

#     dAlfven = field.calc_delta_alfven_juno(time_array, satellite=satellite)

#     tickstring = []
#     for i, pos in enumerate(position):
#         M, _, _ = field.get_M_shell(pos[0], pos[1], pos[2])

#         mlat = (
#             np.arctan2(
#                 position_mag[i, 2],
#                 np.sqrt(position_mag[i, 0] ** 2 + position_mag[i, 1] ** 2),
#             )
#             * 180
#             / np.pi
#         )

#         if i == 0:
#             tickstring.append(
#                 "Time   {:s}".format(time_array[i].strftime(timeformat))
#                 + "\n"
#                 + "M     {:7.2f}".format(M)
#                 + "\n"
#                 + "MLat   {:7.2f}".format(mlat)
#                 + "\n"
#                 + r"$\Delta\lambda$   {:7.2f}".format(dAlfven[i])
#             )
#         else:
#             tickstring.append(
#                 time_array[i].strftime(timeformat)
#                 + "\n"
#                 + "{:7.2f}".format(M)
#                 + "\n"
#                 + "{:7.2f}".format(mlat)
#                 + "\n"
#                 + "{:7.2f}".format(dAlfven[i])
#             )

#     return time_array, tickstring
