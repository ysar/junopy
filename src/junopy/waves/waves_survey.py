import numpy as np
import copy
import pandas as pd
import scipy
import datetime as dt
import glob

import matplotlib.colors as colors

from typing import Optional

from junopy import config
from junopy import index


class WavesData(object):
    def __init__(
        self,
        start: Optional[dt.datetime] = None,
        stop: Optional[dt.datetime] = None,
        filename: Optional[str] = None,
        freq_limits: tuple[Optional(float), Optional(float)] = (None, None),
    ):
        """
        Class to parse data from the Waves instrument.

        Args:
            start (dt.datetime | None): Start time for observations.
            stop (dt.datetime | None): Stop time for observations.
            filename (str | None): Name of the data file to be read.
            freq_limits (tuple[Optional(float), Optional(float)]): Limit analysis to these
                frequencies in Hz. Default: (None, None)
        """
        if filename:
            self.filename = filename
            self.parse_waves_file()

        if start and stop:
            self.df = self.df.loc[
                start.strftime("%Y-%m-%d %H:%M:%S") : stop.strftime("%Y-%m-%d %H:%M:%S")
            ].copy()

        if freq_limits[0]:
            self.channels = self.channels[self.channels["freq"] >= freq_limits[0]]

        if freq_limits[1]:
            self.channels = self.channels[self.channels["freq"] <= freq_limits[1]]

    def parse_waves_file(self):
        """
        Parse a Waves-survey file.
        2 lines contain nonsense.
        3rd line contains channel names.
        4th line contains channel frequency centers.
        5th line contains channel bandwidths.
        Everything else is data.
        """

        with open(self.filename, "r") as f:
            if "Magnetic Field" in f.readline():
                self.yvar = "B\nFreq"
                self.units = " (Hz)"
                self.cbarlabel = "nT$^2$/Hz"
            else:
                self.yvar = "E\nFreq"
                self.units = " (Hz)"
                self.cbarlabel = "(V/m)$^2$/Hz"

            f.readline()

            channel_names = [x.replace('"', "") for x in f.readline()[:-2].split(",")]
            channel_freq = [x.replace('"', "") for x in f.readline()[:-2].split(",")]
            mixed_list = [x.replace('"', "") for x in f.readline()[:-2].split(",")]

        assert len(channel_names) == len(channel_freq) == len(mixed_list)

        ichannels = [i for i, x in enumerate(channel_names) if x]

        self.channels = pd.DataFrame()

        self.channels["names"] = [x for x in channel_names if x]
        self.channels["freq"] = [float(channel_freq[i]) for i in ichannels]
        self.channels["bandwidth"] = [float(mixed_list[i]) for i in ichannels]

        # ctemp = mixed_list[:ichannels[0]] + self.channels['names']
        columns = []
        for cname in mixed_list[: ichannels[0]] + self.channels["names"].to_list():
            if cname not in columns:
                columns.append(cname)
            else:
                columns.append(
                    cname
                    + "_{:d}".format(
                        len(list(filter(lambda x: cname in x, columns))) + 1
                    )
                )

        self.df = pd.read_csv(self.filename, names=columns, skiprows=5)

        self.df["time"] = self.df.apply(
            lambda x: dt.datetime.strptime(x["SCET"], "%Y-%jT%H:%M:%S.%f"), axis=1
        )

        self.df.set_index(["time"], inplace=True)
        self.df.drop(inplace=True, columns=["SCET", "SCLK (sec)"])

    def limit_freq(self, lower, upper):
        self.channels = self.channels[
            (self.channels["freq"] >= lower) & (self.channels["freq"] <= upper)
        ]

    @staticmethod
    def calc_cutoff(x, freq):
        intensity = x.values

        maxima = scipy.signal.argrelextrema(intensity, np.greater, order=4)[0]
        minima = scipy.signal.argrelextrema(intensity, np.less, order=4)[0]
        maxima = maxima[freq[maxima] > 500][0]
        minima = minima[freq[minima] > 500][0]

        return 0.5 * (freq[maxima] + freq[minima])

    def find_cutoff(self):
        """Find the low-frequency cutoff for the continuum radiation"""
        self.df["cutoff"] = self.df[self.channels["names"]].apply(
            lambda x: WavesData.calc_cutoff(x, self.channels["freq"].values), axis=1
        )

    @staticmethod
    def calc_upperhybrid(x, freq):
        intensity = x.values

        maxima = scipy.signal.argrelextrema(np.log(intensity), np.greater, order=10)[0]
        if len(maxima) > 0:
            return freq[maxima[-1]]
        else:
            return np.nan

    def find_upperhybrid(self):
        """Find the low-frequency cutoff for the continuum radiation"""
        self.df["UH"] = self.df[self.channels["names"]].apply(
            lambda x: WavesData.calc_upperhybrid(x, self.channels["freq"].values),
            axis=1,
        )

    def plot(
        self,
        ax,
        mag_data=None,
        freq_low=5e1,
        freq_high=20e3,
        vmin=None,
        vmax=None,
        use_log_scale=True,
        cbar_label="",
        show_cutoff=False,
        **kwargs,
    ):
        """
        Plot Waves-Survey Spectra  (Freq vs time).

        Args:
            ax (mpl.Axes): The matplotlib axes on which to make this plot.
            mag_data (Optional[junopy.MagData]): The magnetic field data from the
                Juno magnetometer (optional).
            freq_low (float): The lowest frequency to plot (Hz).
            freq_high (float): The higher frequency to plot (Hz).
            vmin (Optional[float]): The lowest PSD on the colorbar.
            vmax (Optional[float]): The highest PSD on the colorbar.
            use_log_scale (bool): Whether to use log scale on the Y axis
            cbar_label (str): Label for the color bar.
            show_cutoff (bool): Whether to show the cutoff (if calculated).
            **kwargs: Other arguments to pass to [plt.pcolormesh][]
        """

        channels = self.channels[
            (self.channels["freq"] <= freq_high) & (self.channels["freq"] > freq_low)
        ]

        timearray = self.df.index.to_pydatetime()

        dtfirst = timearray[1] - timearray[0]

        timearray = np.append(timearray - dtfirst / 2, timearray[-1] + dtfirst / 2)

        # Create low and high freq bins.
        ybins = []

        for freq, bw in zip(channels["freq"], channels["bandwidth"]):
            ybins.extend([freq - bw, freq + bw])

        ybins = np.asarray(ybins)

        # Create cdata

        cdata = self.df[channels["names"]].values.T

        c = np.zeros((len(ybins) - 1, len(timearray) - 1))
        c[:, :] = np.nan
        c[::2, :] = cdata[:, :]

        mp = ax.pcolormesh(
            timearray, ybins, c, norm=colors.LogNorm(vmin=vmin, vmax=vmax), **kwargs
        )  # , vmin=1e-15, vmax=1e-12)
        # 14, 10

        if show_cutoff:
            ax.plot(
                self.df.index,
                self.df["cutoff"].values,
                ".",
                markersize=1.0,
                color="m",
                lw=0.5,
            )

        if mag_data:
            ax = mag_data.plot_electron_cyc_freq(ax)

        if use_log_scale:
            ax.set_ylim([freq_low, freq_high])
            ax.set_yscale("log")
            ax.set_yticks([1e2, 1e3, 1e4])
        else:
            ax.set_ylim([100, 1500])
            ax.set_yticks([500, 1000, 1500])

        ax.set_ylabel(self.yvar + self.units)

        bbox = ax.get_position()
        cax = ax.figure.add_axes(
            [bbox.x1 + 0.01, bbox.y0, 0.02, 0.8 * (bbox.y1 - bbox.y0)]
        )
        cbar = ax.figure.colorbar(mp, cax=cax)
        cbar.set_label(cbar_label, rotation=270, labelpad=10.0)

        label = kwargs.pop("label", None)
        if label:
            ax.text(0.1, 1.05, label, transform=ax.transAxes)

        # bar.set_ticks(ticker.LogLocator(base=10, numticks=3, subs='all'))

        return ax


def get_waves_data(
    start: dt.datetime, end: dt.datetime
) -> tuple[Optional(WavesData), Optional(WavesData)]:
    """
    Get the electric and magnetic field data from the Waves survey dataset.

    Args:
        start (dt.datetime): Start time for the observations.
        stop (dt.datetime): Stop time for the observations.

    Returns:
        electric_spectra, magnetic_spectra (tuple[Optional(WavesData),
        Optional(WavesData)]):
            Tuple where the first item is the electric field spectra (from all
            antennae), and the second item is the magnetic field spectra.

    """

    datafolder = config.ConfigVars.get("WAVES_SURVEY_DATA_DIR")
    wavesfiles = glob.glob(datafolder + "/*{:s}*.CSV".format(start.strftime("%Y%j")))

    if len(wavesfiles) != 2:
        wavesindex = index.IndexFile(dataset="waves-survey")
        _ = wavesindex.get_files_between_dates(
            start, end, DoDownload=True, sid="SURVEY"
        )

        wavesfiles = glob.glob(
            datafolder + "/*{:s}*.CSV".format(start.strftime("%Y%j"))
        )

        # assert(len(wavesfiles) == 2)

    for wfile in wavesfiles:
        print("Reading Waves file - {:s}".format(wfile))
        if "_E_" in wfile:
            wavesdataE = WavesData(filename=wfile, start=start, end=end)
        else:
            wavesdataB = WavesData(filename=wfile, start=start, end=end)

    return wavesdataE, wavesdataB


def merge_wavesdata(w1: WavesData, w2: WavesData) -> WavesData:
    """
    Merge two instances of WavesData:

    Args:
        w1 (WavesData): First instance
        w2 (WavesData): Second instance

    Returns:
        WavesData: The merged dataset.
    """
    w3 = copy.copy(w1)
    w3.df = pd.concat([w1.df, w2.df])
    return w3
