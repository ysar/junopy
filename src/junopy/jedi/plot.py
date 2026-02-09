import matplotlib as mpl
import numpy as np
import pandas as pd
from matplotlib import axes, collections, colors

from junopy.utils import misc


def plot_energy_spectra_ssd_native(
    df: pd.DataFrame,
    metadata: pd.DataFrame,
    ax: axes.Axes,
    vmin: float | None = None,
    vmax: float | None = None,
    cmap: str | colors.Colormap = "turbo",
    zorder: int = 2,
    show_colorbar: bool = True,
    cbar_label: str = "",
):
    times = df.index.to_pydatetime()
    delta_t = times[1] - times[0]
    time_bins = times - delta_t * 0.5
    time_bins = np.append(time_bins, time_bins[-1] + delta_t * 0.5)

    bins = np.zeros(df.values.shape[1] * 2)
    bins[::2] = metadata["Elow"]
    bins[1::2] = metadata["Ehigh"]

    z_in = df.values.T
    z = np.zeros((z_in.shape[0] * 2 - 1, z_in.shape[1])) * np.nan
    z[::2, :] = z_in

    if not vmin:
        vmin = misc.get_percentile(z, 1)

    if not vmax:
        vmax = misc.get_percentile(z, 99)

    norm = colors.LogNorm(vmin=vmin, vmax=vmax)

    cmap = mpl.colormaps.get_cmap(cmap)
    cmap.set_bad(color="0.5")

    mp = ax.pcolormesh(time_bins, bins, z, norm=norm, cmap=cmap, zorder=zorder)

    if show_colorbar:
        bbox = ax.get_position()
        cax = ax.figure.add_axes(
            (bbox.x1 + 0.01, bbox.y0, 0.02, 0.8 * (bbox.y1 - bbox.y0)),
        )
        cbar = ax.figure.colorbar(mp, cax=cax)
        cbar.set_label(cbar_label)

    ax.set_yscale("log")


def plot_timeseries_spectra(
    df_spectra: pd.DataFrame,
    bins: np.ndarray,
    ax: axes.Axes,
    vmin: float | None = None,
    vmax: float | None = None,
    cmap: str = "turbo",
    zorder: int = 1,
    show_colorbar: bool = True,
    cbar_label: str = "",
) -> tuple[axes.Axes, collections.QuadMesh]:
    """
    Plot a time series of spectra from a dataframe.

    Args:
        df_spectra (pd.DataFrame): DataFrame containing spectral data with time as index.
        bins (np.ndarray): Array of spectral bins.
        ax (axes.Axes): Matplotlib axes to plot on.
        vmin (float | None): Minimum value for color normalization.
        vmax (float | None): Maximum value for color normalization.
        cmap (str): Colormap to use for the plot.
        zorder (int): Z-order for the plot.
        show_colorbar (bool): Whether to show the colorbar.
        cbar_label (str): Label for the colorbar.

    Returns:
        tuple: Tuple containing the axes and the QuadMesh object.
    """
    # x data is time bins
    times = df_spectra.index.to_pydatetime()
    delta_t = times[1] - times[0]
    time_bins = times - delta_t
    time_bins = np.append(time_bins, time_bins[-1] + (time_bins[1] - time_bins[0]))

    # y data is just passed bins

    # z data are spectral values
    z = df_spectra.values.astype(float).T

    # # Correct for non-contigous data.
    # for i in range(len(times)):
    #     if (time_bins[i + 1] - time_bins[i]).total_seconds() > 301:
    #         z[:, i] = np.nan

    if not vmin:
        vmin = misc.get_percentile(z, 1)

    if not vmax:
        vmax = misc.get_percentile(z, 99)

    norm = colors.LogNorm(vmin=vmin, vmax=vmax)

    _cmap = mpl.colormaps.get_cmap(cmap)
    _cmap.set_bad(color="0.5")

    mp = ax.pcolormesh(time_bins, bins, z, norm=norm, cmap=_cmap, zorder=zorder)

    ax.set_ylim(ymin=bins[0], ymax=bins[-1])

    if show_colorbar:
        bbox = ax.get_position()
        cax = ax.figure.add_axes(
            (bbox.x1 + 0.01, bbox.y0, 0.02, 0.8 * (bbox.y1 - bbox.y0))
        )
        cbar = ax.figure.colorbar(mp, cax=cax)
        cbar.set_label(cbar_label)

    return ax, mp
