import os
import numpy as np
import pandas as pd
import datetime as dt
import glob

import matplotlib.colors as colors

from scipy.constants import c as cspeed
from scipy.constants import m_e, m_p, e

from JunoPy.spice import get_juno_position


class JediData:
    """
    Class to read and do things with a single JEDI data file.
    
    :param str filename: Path to file to be read.
    :param datetime.datetime start: Starting time to slice data.
    :param datetime.datetime end: Ending time to slice data.
    :param str species: Plasma species (Options - 'O\+S', 'OXYGEN', 'SULFUR',
                        'ELECTRON')
    :param str dtbin: String describing the binning frequency (datetime code).
    :param list(float) efilter: List with upper and lower energy limits in keV.
    """
    nSSD = 6
    hasEQPAD = False
    hasValidData = False

    def __init__(self, **kwargs):

        self.filename = kwargs.pop('filename', None)
        if not self.filename:    # Creates empty JediData object 
            return None

        self.species  = kwargs.pop('species', None)
        self.dtbin    = kwargs.pop('dtbin', dt.timedelta(seconds=60))
        self.efilter  = kwargs.pop('efilter', None)
        self.units    = kwargs.pop('units', 'FLUX')
        self.dbfile   = kwargs.pop('dbfile', None) 
        self.nEbins   = kwargs.pop('nEbins', 20)
        self.nPAbins  = kwargs.pop('nPAbins', 18)
        self.start    = kwargs.pop('start', None)
        self.end      = kwargs.pop('end', None)
        self.nEQPAbins = kwargs.pop('nEQPAbins', 9)
        self.dosaveonly = kwargs.pop('dosaveonly', False)

        if '.hdf' in self.filename:
            self.read_hdf_file()
        else:
            self.read_raw_csv_file()
            self.convert_to_hdf_file()

        if self.dosaveonly: return
        else:
            self._filter_metadata()

        if self.start and self.end:
            self.df_raw = self.df_raw.loc[
                self.start.strftime('%Y-%m-%d %H:%M:%S'):
                self.end.strftime('%Y-%m-%d %H:%M:%S')]

        if self.df_raw.empty: return
        else:   
            self.hasValidData = True

        if self.units == 'PSD':
            self.calc_phase_space_density()

        self.create_energytime_df()
        self.create_pitchangle_df()

        if self.dbfile:
            self.calc_equatorial_pad()
            self.create_pitchangle_df(eqpad=True)
            self.hasEQPAD = True

        self.resample()


    def read_raw_csv_file(self):
        """
        Internal function to read the JEDI CSV file.
        :param str dtype: Type of JEDI file to be read 
                          (Options='FLUX'(Default), 'COUNTS', 'CPS')
        """
        with open(self.filename, 'r') as f:
            channel = [x.strip() for x in f.readline().split(',')]
            Elow = [x.strip() for x in f.readline().split(',')]
            Ehigh = [x.strip() for x in f.readline().split(',')]
            Particle = [x.strip() for x in f.readline().split(',')]

        df = pd.DataFrame({
            'channel': channel[:len(Elow)],
            'Elow': Elow,
            'Ehigh': Ehigh,
            'Particle': Particle})

        self.metadata = df

        self.df_raw = pd.read_csv(
            self.filename, names=channel, skiprows=5)

        self.df_raw['time'] = self.df_raw.apply(
                lambda x: dt.datetime.strptime(
                    x['UTC'], '%Y-%m-%dT%H:%M:%S'), axis=1)
        
        self.df_raw.drop(labels='UTC', axis=1, inplace=True)
        self.df_raw.set_index(['time'], inplace=True)
        self.df_raw.sort_index(inplace=True)


    def _filter_metadata(self):
        """
        Impose additional filters on the metadata for later analysis. 
        """
        df = self.metadata.copy()
        df['Particle'] = df['Particle'].str.upper()

        df = df[df['channel'].str.contains('FLUX')]
        df = df[~df['Particle'].str.contains('NONPARTICLE')]
        df = df[df['Particle'].str.contains(self.species)]

        df['Elow'] = df['Elow'].astype('float')
        df['Ehigh'] = df['Ehigh'].astype('float')
        
        df.reset_index(inplace=True)
        

        # Filter unwanted energy channels
        if self.species == 'ELECTRON':
            for i in range(4):
                filter = df['channel'][i][-7:]
                df = df[~df['channel'].str.contains(filter)]      

        if self.species == 'PROTON':
            filter = df['channel'][0][-7:]
            df = df[~df['channel'].str.contains(filter)]

        if self.efilter:
            df = df[(df['Elow'] >= self.efilter[0]) &
                    (df['Ehigh'] <= self.efilter[1])]
        
        # Assign ion mass and rename some strings for plotting.
        if self.species == 'O\+S' or self.species == 'HEAVY':
            self.species = 'O & S'
            m = 16 * m_p

        elif self.species == 'OXYGEN':
            self.species = 'O'
            m = 16 * m_p

        elif self.species == 'SULFUR':
            self.species = 'S'
            m = 32 * m_p

        elif self.species == 'PROTON':
            self.species = 'H+'
            m = m_p * 1.0

        else:
            self.species = 'e-'
            m = m_e * 1.0

        # Have to convert energy units to MeV so p2c2 is in Mev^2
        # df['p2c2'] = np.sqrt( \
        #     df['Elow'] * 1e-3  * (
        #         df['Elow'] * 1e-3 + (2 * m * cspeed**2) / e / 1e6) * \
        #      df['Ehigh'] * 1e-3 * (
        #         df['Ehigh']* 1e-3 + (2 * m * cspeed**2) / e / 1e6))

        # df['p2c2'] = 1 / 3 * 1e-6 * \
        #     (df['Ehigh']**3 - df['Elow']**3) / (df['Ehigh'] - df['Elow']) \
        #         - m**2 * cspeed**4 / e**2 / 1e12
        df['j2psd'] = \
             (m**2 / e**2) * 10 / 4 / \
                np.exp(
                (np.log(df['Ehigh']) + np.log(df['Elow'])) / 2) * 1e-3 * 1e-12
        # (m / m_p)**2 / np.exp(
        #     (np.log(df['Ehigh']) + np.log(df['Elow'])) / 2) * 0.53707 * 1e-18

             # Units of MeV^2

        self.metadata = df.copy()

        vals = self.df_raw.loc[:, self.metadata['channel']].values
        vals[vals < 0] = np.NaN
        self.df_raw.loc[:, self.metadata['channel']] = vals

    def convert_to_hdf_file(self):
        """
        Save the dataframe information into hdf5 files.
        """
        hdf5filename = self.filename[:-4] + '.hdf5'
        self.df_raw.to_hdf(hdf5filename, 'df_raw', mode='a')
        self.metadata.to_hdf(hdf5filename, 'metadata', mode='a')
        os.remove(self.filename)
        self.filename = hdf5filename


    def read_hdf_file(self):
        """
        Read the dataframe information from hdf5 files.
        """
        self.df_raw = pd.read_hdf(self.filename, key='df_raw')
        self.metadata = pd.read_hdf(self.filename, key='metadata')


    def calc_phase_space_density(self):
        """We add columns of PSD to the df_raw DataFrame. Also change the 
        metadata DataFrame to reflect the new channel names."""

        self.metadata['channelpsd'] = self.metadata.apply(
            lambda x: x['channel'].replace('FLUX', 'PSD'), axis=1)

        df = pd.DataFrame(
            index=self.df_raw.index, columns=self.metadata['channelpsd'])
        for i, row in self.metadata.iterrows():
            # df[row['channelpsd']] = \
            #     self.df_raw[row['channel']] / row['p2c2'] * 1.66e-10 * 1e-18
            df[row['channelpsd']] = \
                self.df_raw[row['channel']] * row['j2psd']

        self.df_raw = pd.concat([self.df_raw, df], axis=1)

        self.metadata.drop(labels='channel', axis=1, inplace=True)
        self.metadata.rename(columns={'channelpsd': 'channel'}, inplace=True)


    def create_energytime_df(self):
        """
        Create a new DataFrame with time as the index and different columns
        corresponding to user-defined energy bins.
        """
        if self.efilter:
            self.energybins = np.exp(np.linspace(
                np.log(self.efilter[0]), np.log(self.efilter[1]), 
                self.nEbins + 1))
        else:
            self.energybins = np.exp(np.linspace(
                np.log(self.metadata['Elow'].min()), 
                np.log(self.metadata['Ehigh'].max()), self.nEbins + 1))

        self.energytime = pd.DataFrame(
            index=self.df_raw.index, 
            columns=['Ebin_{:d}'.format(x) for x in range(self.nEbins)])

        for i in range(self.nEbins):

            df = self.metadata.copy()

            df['deltaE'] = df['Ehigh'] - df['Elow']

            df['dE'] = np.minimum(df['Ehigh'], self.energybins[i + 1]) - \
                       np.maximum(df['Elow'] , self.energybins[i])

            df = df[~((df['Ehigh'] < self.energybins[i]) & 
                      (df['Elow']  < self.energybins[i]))]

            df = df[~((df['Ehigh'] > self.energybins[i + 1]) & 
                      (df['Elow']  > self.energybins[i + 1]))]
            
            self.energytime.loc[:, 'Ebin_{:d}'.format(i)] = \
                np.matmul(self.df_raw.loc[:, df['channel']].values, 
                          df['dE'].values) / df['dE'].sum()


    def create_pitchangle_df(self, eqpad=False):
        """
        Create a new DataFrame with time as the index and different columns
        corresponding to user-defined energy bins.
        """
        if eqpad:
            s = 'EQ'
            nbins = self.nEQPAbins
            bins = np.linspace(0, 90, nbins + 1)
            assert(self.units == 'PSD')
        else:
            s = ''
            nbins = self.nPAbins
            bins = np.linspace(0, 180, nbins + 1)
        
        df = pd.DataFrame(
            index=self.df_raw.index, 
            columns=\
                ['PAbin_{:d}'.format(x) for x in range(nbins)] + 
                ['nSSDbin_{:d}'.format(x) for x in range(nbins)], 
            data=0.)

        datacols = []    # Identify data columns for each SSD beforehand
        for j in range(self.nSSD):
            datacols.append(
                self.metadata.loc[
                    self.metadata['channel'].str.contains('T{:d}'.format(j)),
                    'channel'])

        for i in range(nbins):
            for j in range(self.nSSD):

                condition = (
                (self.df_raw['{:s}Pitch_Angle_T{:d}'.format(s, j)] > bins[i]) & 
                (self.df_raw['{:s}Pitch_Angle_T{:d}'.format(s, j)] <=  
                    bins[i + 1]))

                df.loc[     # Average all channels for one SSD
                    condition, 'PAbin_{:d}'.format(i)] += \
                    np.nanmean(
                        self.df_raw.loc[condition, datacols[j]].values, axis=1)

                df.loc[condition, 'nSSDbin_{:d}'.format(i)] += 1

            df['PAbin_{:d}'.format(i)] /= df['nSSDbin_{:d}'.format(i)]

        df.drop(
            labels=['nSSDbin_{:d}'.format(x) for x in range(nbins)],
            axis=1, inplace=True)

        if eqpad:
            self.eqpad = df
            self.eq_pitchanglebins = bins
            self.hasEQPAD = True
        else:
            self.pitchangletime = df
            self.pitchanglebins = bins


    def calc_equatorial_pad(self):
        """
        Create new columns called EQPitch_Angles corresponding to the 
        equatorial pitch angle for each Pitch_Angle column.
        """
        self.pa_database = pd.read_hdf(self.dbfile)

        df = pd.DataFrame(
            index=self.df_raw.index, 
            columns=['EQPitch_Angle_T{:d}'.format(x) for x in range(self.nSSD)])

        self.df_raw = pd.concat([self.df_raw, df], axis=1)

        positionJuno, _ = get_juno_position(
            self.df_raw.index, coordsystem='JUNO_MAG_VIP4')

        pacols = ['Pitch_Angle_T{:d}'.format(x) for x in range(self.nSSD)]
        eqpacols = ['EQPitch_Angle_T{:d}'.format(x) for x in range(self.nSSD)]
        
        pavalues = self.df_raw.loc[:, pacols].values.astype(float)

        pavalues[((positionJuno[:, 2] >= 0) & (pavalues <= 90).T).T] = np.NaN
        pavalues[((positionJuno[:, 2] <  0) & (pavalues >= 90).T).T] = np.NaN

        count = 0
        for t in self.df_raw.index:
            
            Beq = np.nanmean(
                self.pa_database.loc[
                    t - dt.timedelta(minutes=2.5): 
                    t + dt.timedelta(minutes=2.5), 'Beq'].values)

            BJuno = np.nanmean(
                self.pa_database.loc[
                    t - dt.timedelta(minutes=2.5): 
                    t + dt.timedelta(minutes=2.5), 'BJuno'].values)
            
            if BJuno != BJuno:
                continue
            elif Beq > BJuno:
                self.df_raw.loc[t, eqpacols] = np.arcsin(
                    np.sin(pavalues[count, :] * np.pi / 180)) * 180 / np.pi
            else:
                self.df_raw.loc[t, eqpacols] = np.arcsin(
                    np.sin(pavalues[count, :] * np.pi / 180) * \
                    np.sqrt(Beq / BJuno)) * 180 / np.pi
            
            count += 1
        

    def resample(self):
        """
        Resample the time-series DataFrames according to the dtbin specified 
        by the user. 
        """
        self.energytime = self.energytime\
            .resample(self.dtbin, offset=self.dtbin/2)\
            .agg(lambda x: np.nanmean(x))

        self.pitchangletime = self.pitchangletime\
            .resample(self.dtbin, offset=self.dtbin/2)\
            .agg(lambda x: np.nanmean(x))
        
        if self.dbfile:
            self.eqpad = self.eqpad\
                .resample(self.dtbin, offset=self.dtbin/2)\
                .agg(lambda x: np.nanmean(x))


    def plot_spectra(
        self, ax, data='Energy', vmin=None, vmax=None, cmap=None, logplot=True):
        """
        Generic method to plot a timeseries spectra.
        """
        if data == 'Energy':
            ax.set_yscale('log')
            ax.set_yticks([100, 1000, 10000])
            df = self.energytime
            bins = self.energybins
        elif data == 'PAD':
            ax.set_yticks([0, 90, 180])
            df = self.pitchangletime
            bins = self.pitchanglebins 
        elif data == 'EQPAD':
            ax.set_yticks([0, 45, 90])
            df = self.eqpad
            bins = self.eq_pitchanglebins 
        
        # x-data , shifted left by half of dtbin 
        timearray = (df.index - self.dtbin / 2).to_pydatetime()
        timearray = np.append(
            timearray, timearray[-1] + (timearray[1] - timearray[0]))

        # y-data is just bins 

        # z-data is the values of the timeseries DataFrame (transposed)
        z = df.values.T

        if logplot:
            norm = colors.LogNorm(vmin=vmin, vmax=vmax)
        else:
            norm = colors.Normalize(vmin=vmin, vmax=vmax)

        mp = ax.pcolormesh(timearray, bins, z, norm=norm, cmap=cmap)

        ax.set_ylim([bins[0], bins[-1]])

        bbox = ax.get_position()
        cax = ax.figure.add_axes(
            [bbox.x1 + 0.01, bbox.y0, 0.02, 0.8 * (bbox.y1- bbox.y0)])
        cbar = ax.figure.colorbar(mp, cax=cax)

        return ax


def merge_two_timeseries(df1, df2):
    df = pd.concat([df1, df2])
    return df.groupby(df.index).agg(lambda x: np.nanmean(x))


def jedi_merge(j1, j2):
    """
    Merge two Jedi datafiles that have been resampled and created using the 
    same parameters.
    """
    j3 = JediData()

    if j1.hasValidData and not j2.hasValidData: return j1
    if j2.hasValidData and not j1.hasValidData: return j2
    if not j1.hasValidData and not j2.hasValidData:
        return j3

    assert (j1.species == j2.species)     ; j3.species = j1.species
    assert (j1.dtbin == j1.dtbin)         ; j3.dtbin = j1.dtbin
    assert (j1.efilter == j2.efilter)     ; j3.efilter = j1.efilter
    assert (j1.units == j2.units)         ; j3.units = j1.units
    assert (j1.dbfile == j2.dbfile)       ; j3.dbfile = j1.dbfile
    assert (j1.nEbins == j2.nEbins)       ; j3.nEbins = j1.nEbins
    assert (j1.nPAbins == j2.nPAbins)     ; j3.nPAbins = j1.nPAbins
    assert (j1.start == j2.start and j1.end == j2.end)
    assert (j1.nEQPAbins == j2.nEQPAbins) ; j3.nEQPAbins = j1.nEQPAbins

    j3.start = j1.start
    j3.end = j1.end

    j3.energybins = j1.energybins
    j3.pitchanglebins = j1.pitchanglebins
    j3.energytime = merge_two_timeseries(j1.energytime, j2.energytime)
    j3.pitchangletime = merge_two_timeseries(
        j1.pitchangletime, j2.pitchangletime)
    j3.hasValidData = True 
    
    if j1.hasEQPAD and j2.hasEQPAD:
        j3.eq_pitchanglebins = j1.eq_pitchanglebins
        j3.eqpad = merge_two_timeseries(j1.eqpad, j2.eqpad)
        j3.hasEQPAD = True

    return j3

    
def get_jedi_data(**kwargs):
    """
    This script will return JEDI data from all three detectors for all times 
    between start and end. The timeseries data from different days will be 
    combined. 
    """
    if 'ELECTRON' in kwargs['species']:
        idstring = 'ERSESP'
    elif 'PROTON' in kwargs['species']:
        idstring = 'ERSTOFXER'
    elif 'HEAVY' in kwargs['species']:
        idstring = 'NONPTOFXPHR'
    else:
        idstring = 'NONPTOFXER'

    AllData = []

    for day in pd.date_range(kwargs['start'], kwargs['end']):

        dataFiles = glob.glob(
            'data/jedi/*{:s}_CDR_{:s}_*.hdf5'.format(
                idstring, day.strftime('%Y%j')))

        hrsfiles = [x for x in dataFiles if 'HI' in x or 'NONP' in x]
        lrsfiles = [x for x in dataFiles if 'LO' in x]

        if len(hrsfiles) > 0:
            data_hrs = JediData(filename=hrsfiles[0], **kwargs)
            for i in range(1, len(hrsfiles)):
                data_hrs = jedi_merge(
                    data_hrs, JediData(filename=hrsfiles[i], **kwargs))
        else:
            data_hrs = JediData()

        if len(lrsfiles) > 0:
            data_lrs = JediData(filename=lrsfiles[0], **kwargs)
            for i in range(1, len(lrsfiles)):
                data_lrs = jedi_merge(
                    data_lrs, JediData(filename=lrsfiles[i], **kwargs))
        else:
            data_lrs = JediData()

        AllData.append(jedi_merge(data_hrs, data_lrs))

    # Merge different days together. 
    data = AllData[0]
    for i in range(1, len(AllData)):
        data = jedi_merge(data, AllData[i])

    return data

    
  
                
                        



            



            


        
        

        


    



    


