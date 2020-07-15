
import numpy as np
import pandas as pd
import datetime as dt

import matplotlib.colors as colors

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

    def plot(self, ax=None, magdata=None, freqlimit=False, **kwargs):
        """
        Plot Waves-Survey Spectra  (Freq vs time).
        """



        cnames = [x for i, x in enumerate(self.channels['names'])]
                  # if (self.channels['freq'][i] <= 20e3 and self.channels['freq'][i] >= 200)]
        cfreq = [f for f in self.channels['freq']] # if (f <= 20e3 and f >= 200)]

        df = self.df[cnames]

        x = df.index
        y = cfreq
        X, Y = np.meshgrid(x, y)
        C = df.values.T

        mp = ax.pcolormesh(X, Y, C,
                           cmap='jet',
                           norm=colors.LogNorm(), **kwargs) # , vmin=1e-15, vmax=1e-12)
        # 14, 10

        if magdata:
            ax = magdata.plot_electron_cyc_freq(ax)
       
        if freqlimit:
            ax.set_ylim([50, 20e3])
        # ax.set_ylim([200, 20e3]) # self.channels['freq'][-1]]) # 20e3])
        # ax.set_ylim([self.channels['freq'][1], self.channels['freq'][2]])
    
        ax.set_yscale('log')
        ax.set_ylabel(self.yvar + self.units)

        size = 0.02
        bbox = ax.get_position()
        cax = ax.figure.add_axes([bbox.x1, bbox.y0, 0.02, bbox.y1- bbox.y0])
        cbar = ax.figure.colorbar(mp, cax=cax)
        cbar.set_label(self.cbarlabel, rotation=270, labelpad=10.0)

        return ax


class WavesDataDas(object):
    
    DAS_SERVER = 'site:/uiowa/juno/wav/survey/das2'
    
    def __init__(self, start, end):
        """Given a start and end time in datetime, fetch data from Iowa DAS server"""
        import das2
        
        src = das2.get_source(self.DAS_SERVER)
        datasets = src.get({'time' : (start.strftime('%Y-%m-%dT%H:%M:%S'), 
                                     end.strftime('%Y-%m-%dT%H:%M:%S'))})
        
        for d in datasets:
            if 'LFRL' in d.__str__():
                
                self.time = np.zeros(d['Time']['center'].array.shape, dtype=object)
                for i in range(len(self.time[:, 0])):
                    for j in range(len(self.time[0, :])):
                        self.time[i,j] = dt.datetime.utcfromtimestamp(
                            d['Time']['center'].array[i, j].astype(int) * 1e-9)
            
                self.frequency = d['frequency']['center'].array
                self.spec_dens = d['spec_dens']['center'].array
        
        self.cbarlabel = '(V/m)$^2$/Hz'
        
    
    def plot(self, ax=None, magdata=None, freqlimit=False, **kwargs):
        """
        Plot Waves-Survey Spectra  (Freq vs time).
        """
        scalenorm = colors.LogNorm(**kwargs)
        
        mp = ax.pcolormesh(self.time, self.frequency, self.spec_dens,
                           cmap='jet',
                           norm=scalenorm, shading='auto') # , vmin=1e-15, vmax=1e-12)
        # 14, 10

        if magdata:
            ax = magdata.plot_electron_cyc_freq(ax)
       
        if freqlimit:
            ax.set_ylim([50, 20e3])
        # ax.set_ylim([200, 20e3]) # self.channels['freq'][-1]]) # 20e3])
        # ax.set_ylim([self.channels['freq'][1], self.channels['freq'][2]])
    
        ax.set_yscale('log')
        ax.set_ylabel('Freq (Hz)')

        size = 0.02
        bbox = ax.get_position()
        cax = ax.figure.add_axes([bbox.x1, bbox.y0, 0.02, bbox.y1- bbox.y0])
        cbar = ax.figure.colorbar(mp, cax=cax)
        cbar.set_label(self.cbarlabel, rotation=270, labelpad=10.0)

        return ax