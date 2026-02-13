import datetime as dt
import logging
from pathlib import Path
from typing import Optional

import matplotlib.colors as colors
import numpy as np
import pandas as pd
from matplotlib import axes

from junopy import index


class WavesSurveyData(object):
    def __init__(
        self,
        start: Optional[dt.datetime] = None,
        stop: Optional[dt.datetime] = None,
        filename: Optional[Path] = None,
        freq_low: Optional[float] = None,
        freq_high: Optional[float] = None,
        is_hdf: bool = False,
    ):
        """
        Class to parse data from the Waves instrument.

        Args:
            start (dt.datetime | None): Start time for observations.
            stop (dt.datetime | None): Stop time for observations.
            filename (Path | None): Name of the data file to be read.
            freq_low (float | None): Lower frequency limit for analysis in Hz.
            freq_high (float | None): Upper frequency limit for analysis in Hz.
            is_hdf (bool): True if the data file is an HDF5 file.
        """
        if filename:
            logging.info(f"Reading Waves file - {str(filename)}")
            if is_hdf:
                self.read_hdf5_file(filename)
            else:
                self.filename = filename
                self._parse_waves_file()
                self._convert_to_hdf()

        if start and stop:
            self.df = self.df.loc[
                start.strftime("%Y-%m-%d %H:%M:%S") : stop.strftime("%Y-%m-%d %H:%M:%S")
            ].copy()

        if freq_low:
            self.channels = self.channels[self.channels["freq"] >= freq_low]

        if freq_high:
            self.channels = self.channels[self.channels["freq"] <= freq_high]

    def _parse_waves_file(self):
        """
        Parse a Waves-survey file.
        2 lines can be ignored.
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

    def _convert_to_hdf(self):
        """
        Store all attributes of the WavesSurveyData class into an HDF5 file.
        DataFrames are stored as datasets, scalars/strings as metadata.
        """
        hdf_file = self.filename.with_suffix(".hdf5")
        with pd.HDFStore(hdf_file, mode="w") as store:
            # Store the main dataframe
            store.put("df", self.df)
            # Store the channels dataframe if it exists
            if hasattr(self, "channels") and isinstance(self.channels, pd.DataFrame):
                store.put("channels", self.channels)
            # Store scalar/string attributes as metadata on the main dataset
            attrs = ["yvar", "units", "cbarlabel", "filename"]
            for attr in attrs:
                if hasattr(self, attr):
                    store.get_storer("df").attrs[attr] = getattr(self, attr)
        # Remove the original file and update filename
        self.filename.unlink()
        self.filename = hdf_file

    def read_hdf5_file(self, filename: Path) -> None:
        """
        Load all attributes from an HDF5 file into this WavesSurveyData instance.
        """
        with pd.HDFStore(filename, mode="r") as store:
            self.df = store["df"]
            self.channels = store["channels"] if "channels" in store else None
            attrs = store.get_storer("df").attrs
            self.filename = filename
            if hasattr(attrs, "yvar"):
                self.yvar = attrs.yvar
            if hasattr(attrs, "units"):
                self.units = attrs.units
            if hasattr(attrs, "cbarlabel"):
                self.cbarlabel = attrs.cbarlabel

    def limit_freq(self, lower: float, upper: float) -> None:
        """
        Limit the frequency range of the data to the specified lower and upper bounds.

        Args:
            lower (float): Lower frequency limit in Hz.
            upper (float): Upper frequency limit in Hz.
        """
        self.channels = self.channels[
            (self.channels["freq"] >= lower) & (self.channels["freq"] <= upper)
        ]

    def plot(
        self,
        ax: axes.Axes,
        freq_low: float | None = None,
        freq_high: float | None = None,
        vmin: float | None = None,
        vmax: float | None = None,
        use_log_scale: bool = True,
        cbar_label: str = "",
        **kwargs,
    ) -> None:
        """
        Plot Waves-Survey Spectra  (Freq vs time).

        Args:
            ax (mpl.Axes): The matplotlib axes on which to make this plot.
            freq_low (float): The lowest frequency to plot (Hz).
            freq_high (float): The higher frequency to plot (Hz).
            vmin (Optional[float]): The lowest PSD on the colorbar.
            vmax (Optional[float]): The highest PSD on the colorbar.
            use_log_scale (bool): Whether to use log scale on the Y axis
            cbar_label (str): Label for the color bar.
            **kwargs: Other arguments to pass to [plt.pcolormesh][]
        """

        valid_channels = np.ones(len(self.channels), dtype=bool)

        if freq_low:
            valid_channels &= self.channels["freq"] >= freq_low

        if freq_high:
            valid_channels &= self.channels["freq"] <= freq_high

        channels = self.channels[valid_channels]

        timearray = self.df.index.to_pydatetime()

        dtfirst = timearray[1] - timearray[0]

        timearray = np.append(timearray - dtfirst / 2, timearray[-1] + dtfirst / 2)

        # Create low and high freq bins.
        ybins = []

        for freq, bw in zip(channels["freq"], channels["bandwidth"]):
            ybins.extend([freq - bw / 2, freq + bw / 2])

        ybins = np.asarray(ybins)

        # Create cdata

        cdata = self.df[channels["names"]].values.T

        c = np.zeros((len(ybins) - 1, len(timearray) - 1))
        c[:, :] = np.nan
        c[::2, :] = cdata[:, :]

        mp = ax.pcolormesh(
            timearray, ybins, c, norm=colors.LogNorm(vmin=vmin, vmax=vmax), **kwargs
        )

        if freq_low and freq_high:
            ax.set_ylim([freq_low, freq_high])

        if use_log_scale:
            ax.set_yscale("log")

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


def get_waves_data_for_date(
    date: dt.datetime,
    start: dt.datetime | None = None,
    stop: dt.datetime | None = None,
) -> tuple[WavesSurveyData | None, WavesSurveyData | None]:
    """
    Get the electric and magnetic field data from the Waves survey dataset.

    Args:
        start (dt.datetime): Start time for the observations.
        stop (dt.datetime): Stop time for the observations.

    Returns:
        electric_spectra (WavesSurveyData | None): The electric field spectra.
        magnetic_spectra (WavesSurveyData | None): The magnetic field spectra.
    """

    # data_folder = config.ConfigVars.get("WAVES_SURVEY_DATA_DIR")

    if not start:
        start = dt.datetime(date.year, date.month, date.day, 0, 0, 0, 1)
    if not stop:
        stop = dt.datetime(date.year, date.month, date.day, 23, 59, 59, 99999)

    waves_index = index.IndexFile(dataset="waves-survey", do_download=False)
    waves_files = waves_index.get_files_between_dates(
        start, stop, do_download=True, sid="SURVEY"
    )

    waves_files = [
        f.with_suffix(".hdf5") if f.with_suffix(".hdf5").exists() else f
        for f in waves_files
    ]

    for wfile in waves_files:
        is_hdf = wfile.suffix.lower() == ".hdf5"
        if "_E_" in str(wfile):
            waves_psd_electric = WavesSurveyData(
                filename=wfile, start=start, stop=stop, is_hdf=is_hdf
            )
        else:
            waves_psd_magnetic = WavesSurveyData(
                filename=wfile, start=start, stop=stop, is_hdf=is_hdf
            )

    return waves_psd_electric, waves_psd_magnetic


# def merge_waves_survey_data(w1: WavesSurveyData, w2: WavesSurveyData) -> WavesSurveyData:
#     """
#     Merge two instances of WavesSurveyData:

#     Args:
#         w1 (WavesSurveyData): First instance
#         w2 (WavesSurveyData): Second instance

#     Returns:
#         WavesSurveyData: The merged dataset.
#     """
#     w3 = copy.copy(w1)
#     w3.df = pd.concat([w1.df, w2.df])
#     return w3
