import datetime as dt

import matplotlib as mpl
import numpy as np
import pandas as pd
from matplotlib import axes, colors
from matplotlib.collections import PatchCollection
from matplotlib.patches import Ellipse, Wedge

from junopy import spice

from . import convert, misc

# from . import field


def find_ticks_between(
    start: dt.datetime, end: dt.datetime, nticks: int = 6
) -> tuple[pd.DatetimeIndex, list[str]]:
    """
    Given an axes object make the ticks appropriately for subplots.
    """
    times = pd.date_range(start=start, end=end, periods=nticks).to_pydatetime().tolist()

    position, _ = spice.get_juno_position(times)

    tickstring = []
    for i, pos in enumerate(position):
        radius, localtime, latitude = convert.convert_pos_spherical(pos)

        if i == 0:
            tickstring.append(
                "Time - "
                + times[i].strftime("%H:%M:%S")
                + f"\nR - {radius:7.2f}"
                + f"\nLT - {localtime:7.2f}"
                + f"\nLat - {latitude:7.2f}"
            )
        else:
            tickstring.append(
                times[i].strftime("%H:%M:%S")
                + f"\n{radius:7.2f}"
                + f"\n{localtime:7.2f}"
                + f"\n{latitude:7.2f}"
            )

    return times, tickstring


def plot_jupiter(ax: axes.Axes, **kwargs):
    """
    Plot Jupiter this axes.

    Args:
        ax (axes.Axes): Matplotlib axes on which to make this plot.
    """
    jupiter = Ellipse((0.0, 0.0), 2.0, 2.0 * (1.0 - 1.0 / 15.4), **kwargs)
    ax.add_patch(jupiter)


def plot_distribution_perp_parallel(
    ax: mpl.axes.Axes,
    energy_bins: np.ndarray,
    pitchangle_bins: np.ndarray,
    data: np.ndarray,
    vmin: None | float = None,
    vmax: None | float = None,
    cmap: str = "turbo",
) -> PatchCollection:
    """
    Given energy bins and pitch-angle bins, and 2D data, visualize the 2D
    distribution function in perp and parallel coordinates.
    """

    if not vmin:
        vmin = misc.get_percentile(data, 1)

    if not vmax:
        vmax = misc.get_percentile(data, 99)

    def _get_patch(energy1, energy2, pa1, pa2):
        return Wedge(
            (0, 0),
            np.log10(energy2),
            pa1 - 90,
            pa2 - 90,
            width=(np.log10(energy2) - np.log10(energy1)),
        )

    patches = []
    for i in range(len(energy_bins) - 1):
        for j in range(len(pitchangle_bins) - 1):
            patches.append(
                _get_patch(
                    energy_bins[i],
                    energy_bins[i + 1],
                    pitchangle_bins[j],
                    pitchangle_bins[j + 1],
                )
            )

    patch_collection = PatchCollection(patches)
    patch_collection.set_array(data.T.ravel())
    patch_collection.set_cmap(cmap)
    patch_collection.set_norm(colors.LogNorm(vmin=vmin, vmax=vmax))
    ax.add_collection(patch_collection)
    return patch_collection


# def find_ticks_between_Mshell(
#     start, end, nticks=6, timeformat="%H:%M:%S", appendlabels=True
# ):
#     """
#     Given an axes object make the ticks appropriately for subplots.
#     We would like to have :
#     M shell, Latitude, SIII-Longitude
#     """
#     times = pd.date_range(start=start, end=end, periods=nticks).to_pydatetime().tolist()

#     position, _ = spice.get_juno_position(times, coordsystem="IAU_JUPITER")
#     position_mag, _ = spice.get_juno_position(times, coordsystem="JUNO_MAG_VIP4")

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
#                 "Time   {:s}".format(times[i].strftime(timeformat))
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
#                 times[i].strftime(timeformat)
#                 + "\n"
#                 + "{:7.2f}".format(r)
#                 + "\n"
#                 + "{:7.2f}".format(M)
#                 + "\n"
#                 + "{:7.2f}".format(mlat)
#                 + "\n"
#                 + "{:7.2f}".format(lon)
#             )

#     return times, tickstring


# def find_ticks_between_deltaAlfven(
#     start, end, satellite="IO", nticks=6, timeformat="%H:%M:%S"
# ):
#     """
#     Given an axes object make the ticks appropriately for subplots.
#     We would like to have :
#     M shell, Latitude, SIII-Longitude
#     """
#     times = (
#         pd.date_range(start=start, end=end, periods=nticks).to_pydatetime().tolist()
#     )

#     position, _ = spice.get_juno_position(times, coordsystem="IAU_JUPITER")
#     position_mag, _ = spice.get_juno_position(times, coordsystem="JUNO_MAG_VIP4")

#     dAlfven = field.calc_delta_alfven_juno(times, satellite=satellite)

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
#                 "Time   {:s}".format(times[i].strftime(timeformat))
#                 + "\n"
#                 + "M     {:7.2f}".format(M)
#                 + "\n"
#                 + "MLat   {:7.2f}".format(mlat)
#                 + "\n"
#                 + r"$\Delta\lambda$   {:7.2f}".format(dAlfven[i])
#             )
#         else:
#             tickstring.append(
#                 times[i].strftime(timeformat)
#                 + "\n"
#                 + "{:7.2f}".format(M)
#                 + "\n"
#                 + "{:7.2f}".format(mlat)
#                 + "\n"
#                 + "{:7.2f}".format(dAlfven[i])
#             )

#     return times, tickstring
