
import numpy as np
import pandas as pd
import datetime as dt


class WavesData(object):

    def __init__(self, filename=None):
        """
        Class to parse data from the Waves instrument.
        """
        if filename:
            self.filename = filename
            self.parse_waves_file()


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
                self.yvar = 'B-Field Spectra\nFreq'
                self.units  = ' (Hz)'
                self.cbarlabel = 'nT$^2$/Hz'
            else:
                self.yvar = 'E-Field Spectra\nFreq'
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
        self.channels = dict()
        self.channels['names'] = [x for x in channel_names if x]
        self.channels['freq'] = np.array(
            [float(channel_freq[i]) for i in ichannels])
        self.channels['bandwidth'] = np.array(
            [float(mixed_list[i]) for i in ichannels])

        # ctemp = mixed_list[:ichannels[0]] + self.channels['names']
        columns=[]
        for cname in (mixed_list[:ichannels[0]] + self.channels['names']):
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


    def plot(self, ax=None, magdata=None):
        """
        Plot Waves-Survey Spectra  (Freq vs time).
        """
        import matplotlib.colors as colors

        x = self.df.index
        y = self.channels['freq'][self.channels['freq'] <= 20e3]
        X, Y = np.meshgrid(x, y)
        cnames = [x for i, x in enumerate(self.channels['names']) \
                  if self.channels['freq'][i] <= 20e3]
        C = self.df[cnames].values.T

        mp = ax.pcolormesh(X, Y, C,
                           cmap='cividis',
                           norm=colors.LogNorm())

        if magdata:
            ax = magdata.plot_electron_cyc_freq(ax)
            ax.set_ylim([50, 41e6])

        ax.set_yscale('log')
        ax.set_ylabel(self.yvar + self.units)

        size = 0.02
        bbox = ax.get_position()
        cax = ax.figure.add_axes([bbox.x1, bbox.y0, 0.02, bbox.y1- bbox.y0])
        cbar = ax.figure.colorbar(mp, cax=cax)
        cbar.set_label(self.cbarlabel, rotation=270, labelpad=10.0)

        return ax
