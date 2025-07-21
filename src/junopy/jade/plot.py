# Copyright (c) 2018 Yash Sarkango
# Licensed under the MIT License

import numpy as np
import datetime as dt
import pandas as pd
import matplotlib as mpl
from matplotlib import cm
from matplotlib import colors

from .sensor import SensorData


def plot_timeseries_spectra(
    df_spectra: pd.DataFrame,
    bins: np.ndarray,
    ax: mpl.axes.Axes,
    vmin: float | None = None,
    vmax: float | None = None,
    cmap: str = "turbo",
    zorder: int = 1,
    show_colorbar: bool = True,
    cbar_label: str | None = None,
) -> tuple[mpl.axes.Axes, mpl.collections.QuadMesh]:
    """
    Plot a time series of spectra from a dataframe.
    """
    # x data is time bins
    times = df_spectra.index.to_pydatetime()
    delta_t = times[1] - times[0]
    time_bins = times - delta_t
    time_bins = np.append(time_bins, time_bins[-1] + (time_bins[1] - time_bins[0]))

    # y data is just passed bins

    # z data are spectral values
    z = df_spectra.values.T

    # Correct for non-contigous data.
    for i in range(len(times)):
        if (time_bins[i + 1] - time_bins[i]).total_seconds() > 301:
            z[:, i] = np.nan

    if not vmin and not vmax:

        def get_percentile(arr, p):
            return np.exp(np.percentile(np.log(arr[arr > 0 & ~np.isnan(arr)]), p))

        vmin = get_percentile(z, 1)
        vmax = get_percentile(z, 99)

    norm = colors.LogNorm(vmin=vmin, vmax=vmax)

    cmap = mpl.colormaps.get_cmap(cmap)
    cmap.set_bad(color="0.5")

    mp = ax.pcolormesh(time_bins, bins, z, norm=norm, cmap=cmap, zorder=zorder)

    ax.set_ylim([bins[0], bins[-1]])

    if show_colorbar:
        bbox = ax.get_position()
        cax = ax.figure.add_axes(
            [bbox.x1 + 0.01, bbox.y0, 0.02, 0.8 * (bbox.y1 - bbox.y0)]
        )
        cbar = ax.figure.colorbar(mp, cax=cax)
        cbar.set_label(cbar_label)

    return ax, mp


def plot_energy_pa_spectra(
    ax: mpl.axes.Axes,
    sensor_data: SensorData,
    start: dt.datetime | None = None,
    stop: dt.datetime | None = None,
    normalize_energy_bins: bool = False,
    cmap: str = "turbo",
    show_colorbar: bool = True,
    vmin: float | None = None,
    vmax: float | None = None,
) -> tuple[mpl.axes.Axes, mpl.collections.QuadMesh]:
    """
    Plot a pitch angle vs energy spectrogram from JADE sensor data.
    """
    if start and stop:
        indices = (sensor_data.times >= start) & (sensor_data.times <= stop)
    else:
        indices = np.full((sensor_data.data_rebinned.shape[-1],), True)

    # z data are spectral values
    z = np.nanmean(sensor_data.data_rebinned[:, :, indices], axis=2).T

    if normalize_energy_bins:
        z = np.apply_along_axis(lambda x: (x - x.min()) / (x.max() - x.min()), 1, z)
        norm = colors.Normalize
    else:
        norm = colors.LogNorm

    if not vmin and not vmax:

        def get_percentile(arr, p):
            return np.exp(np.percentile(np.log(arr[arr > 0 & ~np.isnan(arr)]), p))

        vmin = get_percentile(z, 1)
        vmax = get_percentile(z, 99)

    mp = ax.pcolormesh(
        sensor_data.pitchangle_bins,
        sensor_data.energy_bins,
        z,
        norm=norm(vmin=vmin, vmax=vmax),
        cmap=cmap,
    )

    ax.set_ylim([sensor_data.energy_bins[0], sensor_data.energy_bins[-1]])
    ax.set_xlim([0, 180])
    ax.set_yscale("log")
    ax.set_xticks([0, 45, 90, 135, 180])

    if show_colorbar:
        bbox = ax.get_position()
        cax = ax.figure.add_axes(
            [bbox.x1 + 0.01, bbox.y0, 0.02, 0.8 * (bbox.y1 - bbox.y0)]
        )
        _ = ax.figure.colorbar(mp, cax=cax)

    return ax, mp


def plot_anode_energy_time_spectra(
    axs: np.array[mpl.axes.Axes],
    sensor_data: SensorData,
    data: str = "data",
    vmin: float | None = None,
    vmax: float | None = None,
    cmap: str = "turbo",
    zorder: int = 1,
    show_colorbar: bool = True,
    locator: mpl.ticker.Locator | None = None,
):
    """
    Plot JADE anode energy-time spectrograms.
    """
    # x data is time bins
    times = sensor_data.times
    delta_t = times[1] - times[0]
    time_bins = times - delta_t
    time_bins = np.append(time_bins, time_bins[-1] + (time_bins[1] - time_bins[0]))

    # y data are the energy bins
    y = sensor_data.energy_bins

    # z data are spectral values
    if data == "data":
        d = sensor_data.data_raw["DATA"]
    elif data == "data+bg":
        sensor_data.data_raw["DATA"][np.isnan(sensor_data.data_raw["DATA"])] = 0.0
        d = sensor_data.data_raw["DATA"] + sensor_data.data_raw["BACKGROUND"]
    elif data == "pitchangle":
        d = sensor_data.data_raw["PITCH_ANGLE"]
    elif data == "s/n":
        sensor_data.data_raw["DATA"][np.isnan(sensor_data.data_raw["DATA"])] = 0.0
        d = sensor_data.data_raw["DATA"] / sensor_data.data_raw["BACKGROUND"]

    if not data == "pitchangle":
        if not vmin and not vmax:

            def get_percentile(arr, p):
                return np.exp(np.percentile(np.log(arr[arr > 0 & ~np.isnan(arr)]), p))

            vmin = get_percentile(d, 1)
            vmax = get_percentile(d, 99)

        norm = colors.LogNorm(vmin=vmin, vmax=vmax)
    else:
        vmin = 0
        vmax = 180
        norm = colors.Normalize(vmin=vmin, vmax=vmax)

    for i in range(4):
        for j in range(4):
            iAnode = 4 * i + j
            axs[i, j].set_facecolor("0.7")
            axs[i, j].pcolormesh(time_bins, y, d[iAnode, :, :], norm=norm, cmap=cmap)
            axs[i, j].set_yscale("log")
            axs[i, j].set_xlabel("Time")
            axs[i, j].set_ylabel("Energy [eV]")
            axs[i, j].text(
                0.1,
                0.8,
                f"Anode={iAnode}",
                backgroundcolor="w",
                transform=axs[i, j].transAxes,
                fontsize=7,
            )
            if locator:
                axs[i, j].xaxis.set_major_locator(locator)
            axs[i, j].label_outer()

    mp = cm.ScalarMappable(norm=norm, cmap=cmap)

    ax = axs[3, 3]
    bbox = ax.get_position()
    cax = ax.figure.add_axes([bbox.x1 + 0.01, bbox.y0, 0.02, 2.5 * (bbox.y1 - bbox.y0)])
    _ = ax.figure.colorbar(mp, cax=cax)
    # cax.text(0.9, 1.1, "Counts s$^{-1}$", transform=cax.transAxes)

    return axs, mp
