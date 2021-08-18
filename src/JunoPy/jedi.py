
import os
import numpy as np
import pandas as pd
import datetime as dt
import glob

from matplotlib import ticker
import matplotlib.colors as colors

from JunoPy.index import IndexFile


class JediData:
    
    nSSD = 6
    pitchangles = np.linspace(0, 170, 18)
    
    def __init__(self, filename=None, start=None, end=None, species=None, dtbin='T'):
        """
        Class to read and do things with a SINGLE JEDI data file.
        """
        
        if not filename:   # Useful for creating empty JediData object.
            return None
        
        self.filename = filename
        
        self.species = species
        
        self.dtbin = dtbin
        
        self.read_csv_file()
        
        if start and end:
            
            
            self.df_raw = self.df_raw.loc[start.strftime('%Y-%m-%d %H:%M:%S'):
                                          end.strftime('%Y-%m-%d %H:%M:%S')]
            
        self.process_data()
            
    
    def read_csv_file(self, dtype='FLUX'):
        
         with open(self.filename, 'r') as f:
        
            channel_names = \
                [x.strip() for x in f.readline().split(',')]
            ebin_lows = [x.strip() for x in f.readline().split(',')]
            ebin_highs = [x.strip() for x in f.readline().split(',')]
            particle_type = [x.strip() for x in f.readline().split(',')]
            
            f.readline()
            
            self.create_jedi_metadata(
                channel_names, ebin_lows, ebin_highs, particle_type, dtype=dtype)
            
            self.df_raw = pd.read_csv(self.filename, names=channel_names, skiprows=5)
            
            self.df_raw['time'] = self.df_raw.apply(
                lambda x: dt.datetime.strptime(x['UTC'], '%Y-%m-%dT%H:%M:%S'), axis=1)
        
            self.df_raw.drop('UTC', 1, inplace=True)
            self.df_raw.set_index(['time'], inplace=True)
            self.df_raw.sort_index(inplace=True)
            
            
    def create_jedi_metadata(self, columns, lowE, highE, particletype, dtype='FLUX'):
        """Return metadata for the datafile."""
        
        df = pd.DataFrame()
        df['channels'] = columns
        df['indices'] = df.index.values
        
        df = df[df['channels'].str.contains(dtype)]
            
        indices = [int(x) for x in df['indices'].values]
        
        df['particletype'] = [particletype[x] for x in indices]
        df['lowEbins'] = [float(lowE[x]) for x in indices]
        df['highEbins'] = [float(highE[x]) for x in indices]
            
        df = df[~df['particletype'].str.contains('NONPARTICLE')]
        
        df.reset_index(inplace=True)
        
        
        if self.species == 'ELECTRON':
            filter_string = df['channels'][0][-7:]    # Filter lowest E channel for e-.
            df = df[~df['channels'].str.contains(filter_string)]
            filter_string = df['channels'][1][-7:]    # Filter lowest2 E channel for e-.
            df = df[~df['channels'].str.contains(filter_string)]
            filter_string = df['channels'][2][-7:]    # Filter lowest3 E channel for e-.
            df = df[~df['channels'].str.contains(filter_string)]
            filter_string = df['channels'][3][-7:]    # Filter lowest3 E channel for e-.
            df = df[~df['channels'].str.contains(filter_string)]
            
        if self.species == 'PROTON':
            filter_string = df['channels'][0][-7:]    # Filter lowest E channel for H+.
            df = df[~df['channels'].str.contains(filter_string)]
            
        if self.species in ['O\+S', 'OXYGEN', 'SULFUR']:
            df = df[df['particletype'].str.contains(self.species)]
            
        if self.species == 'O\+S':
            self.species = 'O & S'
        elif self.species == 'OXYGEN':
            self.species = 'O+/S++'
        elif self.species == 'SULFUR':
            self.species = 'S+'
        elif self.species == 'PROTON':
            self.species = 'H+'
        else:
            self.species = 'e-'

            
        self.metadata = df
        
    
    def process_data(self, start=None, end=None):
        """From the raw data, process the data into a cleaner format."""
        
        # Create pandas dataframe with Energy / time 
        collect_durations = self.df_raw['Duration'].values
        
#         if self.dtbin == 'T':
#             self.energytime = self.df_raw[self.metadata['channels']]\
#                 .apply(lambda x: x * collect_durations, axis=0)\
#                 .resample('T', offset='30S')\
#                 .sum()\
#                 .apply(lambda x: x / 60 / (4 * np.pi), axis=0)
            
#         elif self.dtbin == '30S':
#             self.energytime = self.df_raw[self.metadata['channels']]\
#                 .apply(lambda x: x * collect_durations, axis=0)\
#                 .resample('30S', offset='15S')\
#                 .sum()\
#                 .apply(lambda x: x / 30 / (4 * np.pi), axis=0)
            
        if self.dtbin == 'T':
            self.energytime = self.df_raw[self.metadata['channels']]\
                .resample('T', offset='30S')\
                .mean()
            
        elif self.dtbin == '30S':
            self.energytime = self.df_raw[self.metadata['channels']]\
                .resample('30S', offset='15S')\
                .mean()

        
        # Create energy bins manually based on the largest allowable bins based on min and max.
        
        if len(self.metadata) % self.nSSD != 0:
            raise Exception('Unequal channels among SSDs')
            
        n_energybins = len(self.metadata) // self.nSSD
        
        lowE  = self.metadata['lowEbins'].to_numpy() 
        highE = self.metadata['highEbins'].to_numpy()
        
        self.energybins_low = lowE[:n_energybins]
        self.energybins_high = highE[:n_energybins]
            
        for j in range(1, self.nSSD):

            self.energybins_low  = np.minimum(
                self.energybins_low,   lowE[j * n_energybins: (j + 1) * n_energybins])
            
            self.energybins_high = np.maximum(
                self.energybins_high, highE[j * n_energybins: (j + 1) * n_energybins])
        
        self.n_energybins = n_energybins
        
        diffE = self.energybins_high - self.energybins_low
        
        # Combine flux from different SSDs.
        
        array = self.energytime.values
        
        for i in range(self.n_energybins):
            self.energytime['Channel{:d}_Flux'.format(i)] = array[:, i::self.n_energybins].mean(axis=1)
            
        self.energytime.drop(labels=self.metadata['channels'], axis=1, inplace=True)
        
        
        # Create pitch-angle data (flux units)
        self.pitchangletime = self.df_raw\
            .loc[:, self.metadata['channels']].copy()
        
        pitchanglecols = self.df_raw\
            .loc[:, [x for x in self.df_raw.columns if 'Pitch_Angle' in x]].copy()
        
        for j, pa in enumerate(self.pitchangles):

            self.pitchangletime['Intensity_PA_bin_{:d}'.format(j)] = 0
            self.pitchangletime['nSSDs_PA_bin_{:d}'.format(j)] = 0

            for k in range(self.nSSD):

                condition = (pitchanglecols['Pitch_Angle_T{:d}'.format(k)] >= pa) & \
                            (pitchanglecols['Pitch_Angle_T{:d}'.format(k)] < pa + 10)

#                 self.pitchangletime.loc[condition, 'Intensity_PA_bin_{:d}'.format(j)] += \
#                     self.pitchangletime.loc[condition, 'SSD{:d} Counts'.format(k)].values

                self.pitchangletime.loc[condition, 'Intensity_PA_bin_{:d}'.format(j)] += \
                    self.pitchangletime.loc[
                        condition, 
                        self.pitchangletime.columns.str.contains('T{:d}'.format(k))].mean(axis=1)
                    # .apply(lambda x: (x.values).mean(), axis=1)  # *diffE).sum()
    
                self.pitchangletime.loc[condition, 'nSSDs_PA_bin_{:d}'.format(j)] += 1
            
            self.pitchangletime['Intensity_PA_bin_{:d}'.format(j)] /= \
                self.pitchangletime['nSSDs_PA_bin_{:d}'.format(j)]
        
            
#         if self.dtbin == 'T':
#             self.pitchangletime = self.pitchangletime\
#                 .loc[:, self.pitchangletime.columns.str.contains('Intensity_PA_bin')]\
#                 .apply(lambda x: x * collect_durations, axis=0)\
#                 .resample('T', offset='30S')\
#                 .sum()\
#                 .apply(lambda x: x / 60, axis=0)
            
#         elif self.dtbin == '30S':
#             self.pitchangletime = self.pitchangletime\
#                 .loc[:, self.pitchangletime.columns.str.contains('Intensity_PA_bin')]\
#                 .apply(lambda x: x * collect_durations, axis=0)\
#                 .resample('30S', offset='15S')\
#                 .sum()\
#                 .apply(lambda x: x / 30, axis=0)
            
        if self.dtbin == 'T':
            self.pitchangletime = self.pitchangletime\
                .loc[:, self.pitchangletime.columns.str.contains('Intensity_PA_bin')]\
                .resample('T', offset='30S')\
                .mean()

        elif self.dtbin == '30S':
            self.pitchangletime = self.pitchangletime\
                .loc[:, self.pitchangletime.columns.str.contains('Intensity_PA_bin')]\
                .resample('30S', offset='15S')\
                .mean()

        
        for key in self.pitchangletime.keys():
            if not 'Intensity_PA_bin' in key:
                self.pitchangletime.drop(label=key, axis=1)
            
    
        
    def plot_energy_spectra(self, ax, vmin=None, vmax=None, cmap=None):
        """Energy / time dynamic spectra using pcolormesh."""
        
        # Create x data containing start and end
    
        if self.dtbin == 'T':
            timearray = (self.energytime.index - dt.timedelta(seconds=30)).to_pydatetime()
        elif self.dtbin == '30S':
            timearray = (self.energytime.index - dt.timedelta(seconds=15)).to_pydatetime() 
            
        timearray = np.append(timearray, timearray[-1] + (timearray[1] - timearray[0]))

        # Create y data (pitch-angle or energy bins)
        bins = []
        for lbin, hbin in zip(self.energybins_low, self.energybins_high):
                bins.extend([lbin, hbin])    
        bins = np.asarray(bins)

        # Create c data. Mask appropriately. 

        c0 = self.energytime.values.T
        c0 = np.ma.masked_where(((c0 < 0) | (c0 > 1e15)), c0)

        c = np.zeros((len(bins) - 1, len(timearray) - 1))
        c[:, :] = np.NaN
        c[::2, :] = c0[:, :]

        mp = ax.pcolormesh(timearray, bins, c, 
                           norm=colors.LogNorm(vmin=vmin, vmax=vmax), cmap=cmap)
        
        ax.set_yscale('log')
        
        if self.species in ['oxygen', 'sulfur']:
            ax.set_ylim([100, 2e4])
            ax.set_yticks([100, 1000, 10000])
        else:
            ax.set_ylim([30, 1e3])
            ax.set_yticks([100, 1000])
            
        
        ax.set_ylabel('{:s}\nEnergy'.format(self.species) + '\nkeV')

        size = 0.02
        bbox = ax.get_position()
        cax = ax.figure.add_axes([bbox.x1 + 0.01, bbox.y0, 0.02, 0.8 * (bbox.y1- bbox.y0)])
        cbar = ax.figure.colorbar(mp, cax=cax)
        cbar.set_ticks(ticker.LogLocator(numticks=7, subs='all'))
        
        # cbar.set_label('counts/s/cm$^2$/str/keV', rotation=270, labelpad=25.0, fontsize=9)
    
        
        
        return ax
    
    
    def plot_pitchangle_spectra(self, ax, vmin=None, vmax=None, cmap=None, combine=None):
        
        import matplotlib.colors as colors 
        
        # Create x data containing start and end
        
        if self.dtbin == 'T':
            timearray = (self.energytime.index - dt.timedelta(seconds=30)).to_pydatetime()
        elif self.dtbin == '30S':
            timearray = (self.energytime.index - dt.timedelta(seconds=15)).to_pydatetime() 
            
        timearray = np.append(timearray, timearray[-1] + (timearray[1] - timearray[0]))

        # Create y data (pitch-angle or energy bins)  
        bins = np.append(self.pitchangles, 180)

        # Create c data. Mask appropriately. 
        
        if combine:
            c0 = self.pitchangletime.values.T \
                 + combine[0].pitchangletime.values.T + combine[1].pitchangletime.values.T
            
        c0 = self.pitchangletime.values.T
        c0 = np.ma.masked_where(((c0 < 0) | (c0 > 1e15)), c0)
        
        mp = ax.pcolormesh(timearray, bins, c0, 
                           norm=colors.LogNorm(vmin=vmin, vmax=vmax), cmap=cmap)
        
        ax.set_ylim([0, 180])
        
        if combine:
            ax.set_ylabel('Non-proton\n' + 'PA $^\circ$')
        else:
            ax.set_ylabel('{:s}\n'.format(self.species) + 'PAD $^\circ$')
        
        ax.set_yticks([0, 90, 180])


        size = 0.02
        bbox = ax.get_position()
        cax = ax.figure.add_axes([bbox.x1 + 0.01, bbox.y0, 0.02, 0.8 * (bbox.y1- bbox.y0)])
        
        cbar = ax.figure.colorbar(mp, cax=cax)
        # cbar.set_label('counts/s/cm$^2$/str/keV', rotation=270, labelpad=25.0, fontsize=9)
        # cbar.set_label('intensity', rotation=270, labelpad=20.0)
        cbar.set_ticks(ticker.LogLocator(numticks=7, subs='all'))
    
    
    def plot_flux_timeseries(self, ax, disable_legend=False):
        
        from matplotlib import cm
        
        start = 0.1
        stop = 0.9
        cm_subsection = np.linspace(start, stop, self.n_energybins)
        colors = [ cm.turbo(x) for x in cm_subsection ]
            
        energytime_masked = np.ma.masked_where(((self.energytime.values <= 0) | 
                                                (self.energytime.values > 1e15)), self.energytime.values)
        
        for i in range(self.n_energybins):
            ax.plot(self.energytime.index, energytime_masked[:, i], # '.',
                    lw=2.0, color=colors[i], markersize=4.,
                    label='{:.0f}-{:.0f} keV'.format(self.energybins_low[i], self.energybins_high[i]))
        
        if not disable_legend:
            ax.legend(loc=(1.01, 0.01), fontsize=6, framealpha=1.0, frameon=False)
            
        ax.set_yscale('log')
        # ax.grid('on', which='minor', lw='0.5', ls='dotted')
        # ax.set_xlim([self.df.index[0], self.df.index[-1]])
        
        ax.set_ylabel('{:s}\nDiff.\nIntensity'.format(self.species))
        ax.yaxis.set_major_locator(ticker.LogLocator(base=10, numticks=15))
    
    
def merge_two_jedi_data(j1, j2):
    """
    Merges two jeDI data files of same type & species. 
    Returns a JediData object.
    """
        
    j3 = JediData()    # Instantiate empty for merge. 

    
    if not (j1.n_energybins == j2.n_energybins): 
        raise Exception('Unequal Energy bins.')
    else:
        j3.n_energybins = j1.n_energybins
        
#     if not (j1.energytime.index.equals(j2.energytime.index)):
#         raise Exception('Unequal time tags.')
        
    if not (j1.species == j2.species):
        raise Exception('Merging files of different species.')
    else:
        j3.species = j1.species
        
    j3.dtbin = j1.dtbin

    j3.pitchangles = j1.pitchangles.copy()
    
    df_concat = pd.concat([j1.energytime, j2.energytime])
    j3.energytime = df_concat.groupby(df_concat.index).mean()
    
    df_concat = pd.concat([j1.pitchangletime, j2.pitchangletime])
    j3.pitchangletime = df_concat.groupby(df_concat.index).mean()
        
    j3.energybins_low = np.minimum(j1.energybins_low, j2.energybins_low)
    j3.energybins_high = np.maximum(j1.energybins_high, j2.energybins_high)
    
    return j3

        

def get_jedi_data(start, end, dtbin='T'):
    
    electron_data = check_or_fetch('ERSESP', start, end, 'ELECTRON', dtbin=dtbin)
    
    proton_data = check_or_fetch('ERSTOFXER', start, end, 'PROTON', dtbin=dtbin)
    
    nonproton_data = []
    nonproton_data.append(check_or_fetch('NONPTOFXER', start, end, 'O\+S', dtbin=dtbin))
    nonproton_data.append(check_or_fetch('NONPTOFXER', start, end, 'OXYGEN', dtbin=dtbin))
    nonproton_data.append(check_or_fetch('NONPTOFXER', start, end, 'SULFUR', dtbin=dtbin))
    
    return electron_data, proton_data, nonproton_data
        

def check_or_fetch(sid, start, end, species, dtbin='T'):
    
    

    jediindex = IndexFile(dataset='jedi')
    
    files = glob.glob('data/jedi/*{:s}_CDR_{:s}_*.TAB'.format(sid, start.strftime('%Y%j')))
    
    if not files:
        
        jediindex.get_files_between_dates(start, end, DoDownload=True, sid=sid)
        
        files = glob.glob('data/jedi/*{:s}_CDR_{:s}_*.TAB'.format(sid, start.strftime('%Y%j')))
    
    
    if len(files) == 0:
        return None
    else:
        print(files)
        
    merged = JediData(filename=files[0], start=start, end=end, species=species, dtbin=dtbin)
    
    if len(files) == 1:
        return merged
    
    
    for i in range(1, len(files)):
        
        merged = merge_two_jedi_data(
            merged, JediData(filename=files[i], start=start, end=end, species=species, dtbin=dtbin)) 
        
    return merged