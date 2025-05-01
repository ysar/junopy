
import numpy as np
import copy
import pandas as pd
from pandas.io.formats.style_render import Optional
import scipy
import datetime as dt
import glob

import matplotlib.colors as colors
from matplotlib import ticker

from junopy import config
from junopy import index

class WavesData(object):

    def __init__(
        self,
        start: Optional[dt.datetime] = None,
        stop: Optional[dt.datetime] = None,
        filename: Optional[str] = None,
        freq_limits: tuple[float, float] | None = None
    ):
        """
        Class to parse data from the Waves instrument.

        Args:
            start (dt.datetime | None): Start time for observations.
            stop (dt.datetime | None): Stop time for observations.
            filename (str | None): Name of the data file to be read.
            freq_limits (tuple[float, float]): Limit analysis to these
                frequencies in Hz.
        """
        if filename:
            self.filename = filename
            self.parse_waves_file()

        if start and stop:
            self.df = self.df.loc[start.strftime('%Y-%m-%d %H:%M:%S'):
                                  stop.strftime('%Y-%m-%d %H:%M:%S')].copy()

        if freq_limits:
            self.channels = self.channels[(self.channels['freq'] >= freq_limits[0]) &
                                          (self.channels['freq'] <= freq_limits[1])]

    def parse_waves_file(self):
        """
        Parse a Waves-survey file.
        2 lines contain nonsense.
        3rd line contains channel names.
        4th line contains channel frequency centers.
        5th line contains channel bandwidths.
        Everything else is data.
        """

        with open(self.filename, 'r') as f:

            if 'Magnetic Field' in f.readline():
                self.yvar = 'B\nFreq'
                self.units  = ' (Hz)'
                self.cbarlabel = 'nT$^2$/Hz'
            else:
                self.yvar = 'E\nFreq'
                self.units = ' (Hz)'
                self.cbarlabel = '(V/m)$^2$/Hz'

            f.readline()

            channel_names = [
                x.replace('"', '') for x in  f.readline()[:-2].split(',')]
            channel_freq = [
                x.replace('"', '') for x in  f.readline()[:-2].split(',')]
            mixed_list = [
                x.replace('"', '') for x in  f.readline()[:-2].split(',')]

        assert (len(channel_names) == len(channel_freq) == len(mixed_list))

        ichannels = [i for i, x in enumerate(channel_names) if x]


        self.channels = pd.DataFrame()

        self.channels['names'] = [x for x in channel_names if x]
        self.channels['freq'] = [float(channel_freq[i]) for i in ichannels]
        self.channels['bandwidth'] = [float(mixed_list[i]) for i in ichannels]

        # ctemp = mixed_list[:ichannels[0]] + self.channels['names']
        columns=[]
        for cname in (mixed_list[:ichannels[0]] + self.channels['names'].to_list()):
            if cname not in columns:
                columns.append(cname)
            else:
                columns.append(
                    cname + '_{:d}'.format(
                        len(list(filter(lambda x: cname in x, columns))) + 1))

        self.df = pd.read_csv(self.filename, names=columns, skiprows=5)

        self.df['time'] = self.df.apply(
            lambda x: dt.datetime.strptime(x['SCET'], '%Y-%jT%H:%M:%S.%f'),
            axis=1)

        self.df.set_index(['time'], inplace=True)
        self.df.drop(
            inplace=True,
            columns=['SCET', 'SCLK (sec)']
            )

    def limit_freq(self, lower, upper):
        self.channels = self.channels[(self.channels['freq'] >= lower) &
                                      (self.channels['freq'] <= upper)]

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
        self.df['cutoff'] = self.df[self.channels['names']].apply(
            lambda x: WavesData.calc_cutoff(x, self.channels['freq'].values), axis=1)

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
        self.df['UH'] = self.df[self.channels['names']].apply(
            lambda x: WavesData.calc_upperhybrid(x, self.channels['freq'].values), axis=1)

    def plot(self, ax=None, magdata=None, extendmag=False,
             flow=50, fhigh=20e3, vmin=None, vmax=None,
             logscale=True, cbarlabel=None, showcutoff=False, **kwargs):
        """
        Plot Waves-Survey Spectra  (Freq vs time).
        """

        channels = self.channels[(self.channels['freq'] <= fhigh) & \
                                      (self.channels['freq'] > flow)]


        timearray = self.df.index.to_pydatetime()

        dtfirst = (timearray[1] - timearray[0])

        timearray = np.append(timearray - dtfirst / 2,
                              timearray[-1] + dtfirst / 2)

        # Create low and high freq bins.
        ybins = []

        for freq, bw in zip(channels['freq'], channels['bandwidth']):
            ybins.extend([freq - bw, freq + bw])

        ybins = np.asarray(ybins)

        # Create cdata

        cdata = self.df[channels['names']].values.T

        c = np.zeros((len(ybins) - 1, len(timearray) - 1))
        c[:, :] = np.nan
        c[::2, :] = cdata[:, :]

        mp = ax.pcolormesh(timearray, ybins, c,
                           norm=colors.LogNorm(vmin=vmin, vmax=vmax), **kwargs) # , vmin=1e-15, vmax=1e-12)
        # 14, 10

        if showcutoff:
            ax.plot(self.df.index, self.df['cutoff'].values, '.', markersize=1., color='m', lw=0.5)

        if magdata:
            ax = magdata.plot_electron_cyc_freq(ax)

            if extendmag:    # Perform rolling-average FFT for the magnetometer data.
                pass

        if logscale:
            ax.set_ylim([flow, fhigh])
            ax.set_yscale('log')
            ax.set_yticks([1e2, 1e3, 1e4])
        else:
            ax.set_ylim([100, 1500])
            ax.set_yticks([500, 1000, 1500])

        ax.set_ylabel(self.yvar + self.units)

        size = 0.02
        bbox = ax.get_position()
        cax = ax.figure.add_axes([bbox.x1 + 0.01, bbox.y0, 0.02, 0.8 * (bbox.y1- bbox.y0)])
        cbar = ax.figure.colorbar(mp, cax=cax)

        if cbarlabel:
            cbar.set_label(cbarlabel, rotation=270, labelpad=10.0)

        label = kwargs.pop('label', None)
        if label:
            ax.text(0.1, 1.05, label, transform=ax.transAxes)

        # bar.set_ticks(ticker.LogLocator(base=10, numticks=3, subs='all'))

        return ax



def get_waves_data(start, end):
    """Get the electric and magnetic field data from the Waves survey."""

    datafolder = config.ConfigVars.get("WAVES_SURVEY_DATA_DIR")
    wavesfiles = glob.glob(datafolder + '/*{:s}*.CSV'.format(start.strftime('%Y%j')))

    if len(wavesfiles) != 2:
        wavesindex = index.IndexFile(dataset='waves-survey')
        _ = wavesindex.get_files_between_dates(
            start, end, DoDownload=True, sid='SURVEY')

        wavesfiles = glob.glob(datafolder + '/*{:s}*.CSV'.format(start.strftime('%Y%j')))

        # assert(len(wavesfiles) == 2)

    for wfile in wavesfiles:
        print('Reading Waves file - {:s}'.format(wfile))
        if '_E_' in wfile:
            wavesdataE = WavesData(filename=wfile, start=start, end=end)
        else:
            wavesdataB = WavesData(filename=wfile, start=start, end=end)

    return wavesdataE, wavesdataB


def merge_wavesdata(w1, w2):

    w3 = copy.copy(w1)

    w3.df = pd.concat([w1.df, w2.df])

#     m3.filename = os.path.join(
#         os.path.split(m1.filename)[0],
#         os.path.split(m1.filename)[-1][:-4] + '_' + os.path.split(m2.filename)[-1])


    return w3
