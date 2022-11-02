from JunoPy.jedi import merge_two_timeseries
import numpy as np
import datetime as dt
import pandas as pd
import scipy
import numpy.ma as ma
import pidly
import contextlib

import matplotlib.colors as colors

from scipy.io import readsav
from JunoPy.spice import get_juno_position
from JunoPy.library import find_ticks_between_Mshell
from JunoPy.mag import get_mag_data


class JadeData(object):
    """
    This class reads in the .sav files written out by the IDL tools for JADE.
    The .sav files are read into numpy recarray (by default). The purpose of 
    this class is to convert from recarray into pandas DataFrames.
    """
    hasEQPAD = False
    hasValidData = False

    def __init__(self, **kwargs):
        
        self.filename = kwargs.pop('filename', None)
        if not self.filename:
            return None

        self.magfolder = kwargs.pop('magfolder', 'data/mag')
        self.efilter   = kwargs.pop('efilter', None)
        self.dtbin     = kwargs.pop('dtbin', dt.timedelta(seconds=120))
        self.dbfile    = kwargs.pop('dbfile', None)
        self.nEbins    = kwargs.pop('nEbins', 20)
        self.nPAbins   = kwargs.pop('nPAbins', 18)
        self.nEQPAbins = kwargs.pop('nEQPAbins', 9)

        # units, start and end are defined within the temp .sav file, not here

        check = self.read_sav_file()
        if check == -1:
            return
        
        self.hasValidData = True
        self.create_energytime_df()
        self.create_pitchangle_df()

        self.resample()


    def read_sav_file(self):
        """Read in the .sav file stored using IDL."""
        self.data_raw = readsav(self.filename)['data']

        if isinstance(self.data_raw, np.int16) and self.data_raw == -1:
            return -1

        df = pd.DataFrame()
        df['UTC'] = self.data_raw['DIM0_UTC'][0]

        df['time'] = df.apply(
            lambda x: dt.datetime.strptime(
                x['UTC'].decode('UTF-8'), '%Y-%jT%H:%M:%S.%f'), axis=1)

        self.times = df['time'].squeeze()
        return 0

        
    def create_energytime_df(self):
        """
        Create a dataframe with time in index and energy bins as columns.
        """
        data = self.data_raw['DATA'][0][:, :, :]
        energies = self.data_raw['DIM1_E'][0][:, :, :]

        # Filter the invalid data
        self.NaNMask = (data < 0)

        # Create mask based on Efilter and specify new energy bins
        if self.efilter:
            self.Emask = (energies > self.efilter[0]) & \
                         (energies <= self.efilter[1])
            self.energybins = np.exp(np.linspace(
                np.log(self.efilter[0]), np.log(self.efilter[1]), 
                self.nEbins + 1))
        else:
            self.Emask = np.zeros(energies.shape, dtype=bool)
            self.energybins = np.exp(np.linspace(
                np.log(energies.min()), np.log(energies.max()), 
                self.nEbins + 1))

        self.energytime = pd.DataFrame(
            index=self.times,
            columns=['Ebin_{:d}'.format(x) for x in range(self.nEbins)])

        for i in range(self.nEbins):

            tmpEmask = (energies > self.energybins[i]) & \
                       (energies <= self.energybins[i + 1])

            self.energytime.loc[:, 'Ebin_{:d}'.format(i)] = \
                ma.masked_where(
                    ~(self.Emask & tmpEmask) | self.NaNMask, 
                    data).mean(axis=(0, 1))


    def create_pitchangle_df(self):
        """
        Create a dataframe with time in index and pitch-angle bins as columns.
        """
        data = self.data_raw['DATA'][0][:, :, :]
        elevation = self.data_raw['DIM2_ELEVATION'][0][:, :, :] * np.pi / 180
        azimuth = self.data_raw['DIM2_AZIMUTH_DESPUN'][0][:, :, :] * np.pi/180

        magdata = get_mag_data(
            self.times.iloc[0], self.times.iloc[-1], 
            self.magfolder, coordsystem='Juno_JADE_despun')

        # Interpolate MAG data onto times for JADE
        df = pd.DataFrame(index=self.times, columns=['Bx', 'By', 'Bz'])

        for t, row in df.iterrows():
            df.loc[t, :] = magdata.df.loc[
                t - dt.timedelta(seconds=0.6):
                t + dt.timedelta(seconds=0.6), 
                ['Bx', 'By', 'Bz']].mean(axis=0).values

        Bx = df['Bx'].values.astype('float')
        By = df['By'].values.astype('float')
        Bz = df['Bz'].values.astype('float')

        Bmag = np.sqrt(Bx**2 + By**2 + Bz**2)
        Bx /= Bmag; By /= Bmag; Bz /= Bmag

        self.pitchangles = np.arccos(-(
            np.cos(elevation) * np.cos(azimuth) * Bx + \
            np.cos(elevation) * np.sin(azimuth) * By + \
            np.sin(elevation) * Bz)) * 180 / np.pi

        self.magdata = df
        self.pitchanglebins = np.linspace(0, 180, self.nPAbins + 1)

        self.pitchangletime = pd.DataFrame(
            index=self.times,
            columns=['PAbin_{:d}'.format(x) for x in range(self.nPAbins)])

        for i in range(self.nPAbins):

            tmpPAmask = (self.pitchangles > self.pitchanglebins[i]) & \
                        (self.pitchangles <= self.pitchanglebins[i + 1])

            self.pitchangletime.loc[:, 'PAbin_{:d}'.format(i)] = \
                ma.masked_where(
                    ~(self.Emask & tmpPAmask) | self.NaNMask, 
                    data).mean(axis=(0, 1))

        # Calculate equatorial pitch-angles if database file is provided.
        if not self.dbfile:
            return None

        self.pa_database = pd.read_hdf(self.dbfile)
        self.eq_pitchanglebins = np.linspace(0, 90, self.nEQPAbins + 1)
        self.eqpad = pd.DataFrame(
            index=self.times,
            columns=['PAbin_{:d}'.format(x) for x in range(self.nEQPAbins)])
        
        self.eqpa = self.pitchangles * 1.0

        count = 0
        for t in self.times:

            Beq = np.nanmean(
                self.pa_database.loc[
                    t - dt.timedelta(minutes=2.5): 
                    t + dt.timedelta(minutes=2.5), 'Beq'].values)

            BJuno = Bmag[count]

            if BJuno != BJuno:
                self.eqpa[:, :, count] = -1.
            elif Beq > BJuno:
                self.eqpa[:, :, count] = np.arcsin(
                    np.sin(self.eqpa[:, :, count] * np.pi / 180)) * 180 / np.pi
            else:
                self.eqpa[:, :, count] = np.arcsin(
                    np.sin(self.eqpa[:, :, count] * np.pi / 180) * \
                    np.sqrt(Beq / BJuno)) * 180 / np.pi

            count += 1

        invalidEQPA = (self.eqpa < 0)

        for i in range(self.nEQPAbins):

            tmpPAmask = (self.eqpa > self.eq_pitchanglebins[i]) & \
                        (self.eqpa <= self.eq_pitchanglebins[i + 1])

            self.eqpad.loc[:, 'PAbin_{:d}'.format(i)] = \
                ma.masked_where(
                    ~(self.Emask & tmpPAmask) | self.NaNMask | invalidEQPA, 
                    data).mean(axis=(0, 1))

        self.hasEQPAD = True


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


def jade_merge(j1, j2):
    j3 = JadeData()

    assert (j1.hasValidData and j2.hasValidData) ; j3.hasValidData = True
    assert (j1.efilter == j2.efilter)       ; j3.efilter = j1.efilter 
    assert (j1.dtbin == j2.dtbin)           ; j3.dtbin = j1.dtbin
    assert (j1.dbfile == j2.dbfile)         ; j3.dbfile = j1.dbfile 
    assert (j1.nEbins == j2.nEbins)         ; j3.nEbins = j1.nEbins
    assert (j1.nPAbins == j2.nPAbins)       ; j3.nPAbins = j1.nPAbins 
    assert (j1.nEQPAbins == j2.nEQPAbins)   ; j3.nEQPAbins = j1.nEQPAbins

    j3.energybins = j1.energybins
    j3.pitchanglebins = j1.pitchanglebins
    j3.energytime = merge_two_timeseries(j1.energytime, j2.energytime)
    j3.pitchangletime = merge_two_timeseries(
        j1.pitchangletime, j2.pitchangletime)

    if j1.hasEQPAD and j2.hasEQPAD:
        j3.eq_pitchanglebins = j1.eq_pitchanglebins
        j3.eqpad = merge_two_timeseries(j1.eqpad, j2.eqpad)
        j3.hasEQPAD = True

    return j3


def get_jade_data(start, end, species=3, **kwargs):
    """
    This method will return JADE data combined for high/low res between the 
    start and end times. 
    """
    if isinstance(start, dt.datetime):
        ts = start.strftime('%Y-%jT%H:%M')
        te = end.strftime('%Y-%jT%H:%M')
    else:
        ts = start
        te = end

    kwargs['species'] = species

    with contextlib.redirect_stdout(None):
        idl = pidly.IDL('/Users/ysar/Applications/harris/idl88/bin/idl')
        idl(".run set_jade_vars")
        idl("set_jade_vars")
        idl("CD, '/Users/ysar/Documents/code/JADE/'")

        filename = '/Users/ysar/Documents/code/JADE/jade-yash/tmp/temp.sav'

        idl("data = jade_js_get_sci(['{:s}', '{:s}']"
            ",type='ion',hrs=1,species={:d})".format(ts, te, species))
        idl("data = jade_convert_data_units_from_cnt(data, UNITS='PSD', m=1, q=1)")
        idl("save, data, FILENAME = '{:s}'".format(filename))
        
        data_highres = JadeData(filename=filename, **kwargs)
        
        idl("data = jade_js_get_sci(['{:s}', '{:s}']"
            ",type='ion',species={:d})".format(ts, te, species))
        idl("data = jade_convert_data_units_from_cnt(data, UNITS='PSD', m=1, q=1)")
        idl("save, data, FILENAME = '{:s}'".format(filename))

        data_lowres = JadeData(filename=filename, **kwargs)

    if data_highres.hasValidData and data_lowres.hasValidData:
        return jade_merge(data_highres, data_lowres)
    elif data_highres.hasValidData:
        return data_highres
    elif data_lowres.hasValidData:
        return data_lowres
    else:
        return JadeData()

        
        

        
