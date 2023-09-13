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
from scipy.constants import e
from JunoPy.spice import get_juno_position
from JunoPy.library import find_ticks_between_Mshell
from JunoPy.library import cxcs_Hplus_O, cxcs_Hplus_H2
from JunoPy.mag import get_mag_data


class JadeData(object):
    """
    This class reads in the .sav files written out by the IDL tools for JADE.
    The .sav files are read into numpy recarray (by default). The purpose of 
    this class is to convert from recarray into pandas DataFrames.
    """
    hasEQPAD = False
    hasValidData = False
    hasCXrate = False
    hasEflux = False

    def __init__(self, **kwargs):
        
        self.filename = kwargs.pop('filename', None)
        if not self.filename:
            return None

        self.magfolder = kwargs.pop('magfolder', 'data/mag')
        self.efilter   = kwargs.pop('efilter', None)
        self.dtbin     = kwargs.pop('dtbin', dt.timedelta(seconds=120))
        self.dbfile    = kwargs.pop('dbfile', None)
        # self.nEbins    = kwargs.pop('nEbins', 20)
        self.isHighRate= kwargs.pop('isHighRate', False)
        self.nPAbins   = kwargs.pop('nPAbins', 12)
        self.nEQPAbins = kwargs.pop('nEQPAbins', 9)
        self.units     = kwargs.pop('units', 'PSD')
        self.calcEflux = kwargs.pop('calcEflux', False)
        self.paMin     = kwargs.pop('paMin', 0)
        self.paMax     = kwargs.pop('paMax', 180)
        self.addBackground = kwargs.pop('addBackground', False)
        self.useOnlyBackground = kwargs.pop('useOnlyBackground', False)

        docalcCX = kwargs.pop('docalcCX', False)
        self.cxtype = kwargs.pop('cxtype', 'Hplus-O')

        # units, start and end are defined within the temp .sav file, not here

        check = self.read_sav_file()
        if check == -1:
            return
        
        self.hasValidData = True
        self.reorganize_data()

        if docalcCX:
            self.calc_cx_rate()

        if self.calcEflux:
            self.calc_energy_flux(
                pa1=self.paMin, pa2=self.paMax)
            
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

        
    def reorganize_data(self):
        """
        We first rebin the data from look direction to pitch angle and then 
        average over the relevant axes to get the energy or pitch angle spectra.
        Similarly, we will rebin into equatorial pitch angles. 
        """
        energies  = self.data_raw['DIM1_E'][0] * 1.0
        flux      = self.data_raw['DATA'][0] * 1.0
        elevation = self.data_raw['DIM2_ELEVATION'][0] * np.pi / 180
        azimuth   = self.data_raw['DIM2_AZIMUTH_DESPUN'][0] * np.pi / 180
        
        self.data_sigma = self.data_raw['DATA_SIGMA'][0] * 1.0

        # Add in the background (CAUTION: DO THIS ONLY WHEN NEEDED)
        if self.addBackground:
            background = self.data_raw['BACKGROUND'][0] * 1.0
            self.data_sigma = np.sqrt(
                self.data_raw['DATA_SIGMA'][0]**2 + \
                self.data_raw['BACKGROUND_SIGMA'][0]**2)
            flux += background

        if self.useOnlyBackground:
            flux = self.data_raw['BACKGROUND'][0] * 1.0
            self.data_sigma = self.data_raw['BACKGROUND_SIGMA'][0] * 1.0
        

        # Change all negative values to zero
        flux[flux < -999990] = np.NaN
        flux[(flux < 0) & (flux > -999990)] = 0 

        # Create energybins based on log-space intervals 
        # (and)
        # Filter out those energies requested for by user

        energylin = energies[0, :, 0]

        if self.efilter:

            indices = np.asarray(
                (energylin >= self.efilter[0]) & 
                (energylin <= self.efilter[1])).nonzero()[0]

            energies    = energies[:, indices, :]
            elevation   = elevation[:, indices, :]
            azimuth     = azimuth[:, indices, :]
            flux        = flux[:, indices, :]
            self.data_sigma = self.data_sigma[:, indices, :]
            energylin   = energylin[indices]


        self.energybins = np.zeros(len(energylin) + 1)
        self.energybins[1:-1] = np.sqrt(energylin[:-1] * energylin[1:])
        self.energybins[0] = self.energybins[1]**2 / self.energybins[2]
        self.energybins[-1] = self.energybins[-2]**2 / self.energybins[-3]
        self.nEbins = len(self.energybins) - 1
        energybinsize = self.energybins[1:] - self.energybins[:-1]

        # If high rate, remove data from the first 4 anodes
        # if self.isHighRate:
        #     assert (len(elevation[:, 0, 0] == 12))
        #     energies = energies[4:, :, :]
        #     elevation = elevation[4:, :, :]
        #     azimuth = azimuth[4:, :, :]
        #     flux = flux[4:, :, :]


        # Get magnetometer data for pitch angle calculations
        magdata = get_mag_data(
            self.times.iloc[0], self.times.iloc[-1], 
            datafolder=self.magfolder, coordsystem='Juno_JADE_despun')

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

        pitchangles = np.arccos(-(
            np.cos(elevation) * np.cos(azimuth) * Bx + \
            np.cos(elevation) * np.sin(azimuth) * By + \
            np.sin(elevation) * Bz)) * 180 / np.pi

        self.pitchanglebins = np.linspace(0, 180, self.nPAbins + 1)
        pitchanglemid = 0.5 * (
            self.pitchanglebins[:-1] + self.pitchanglebins[1:])

        nEorig = flux.shape[1]
        nTimes = flux.shape[2]

        # Mask invalid pitch angles
        mask = ((pitchangles < self.paMin) | (pitchangles > self.paMax))
        flux = ma.masked_where(mask, flux)

        # Rebin from LookDirection to PitchAngle
        flux_rebin = np.zeros((self.nPAbins, nEorig, nTimes))
        
        nanmask = np.isnan(flux)

        for ii in range(self.nPAbins):
            
            mask = (self.pitchanglebins[ii] <= pitchangles) & \
                   (self.pitchanglebins[ii+1]> pitchangles)
            
            flux_rebin[ii, :, :] = ma.masked_where(
                ~mask | nanmask, flux).mean(axis=0).filled(np.NaN)

        dPA = self.pitchanglebins[1] - self.pitchanglebins[0]
        solidAngle = np.sin(pitchanglemid * np.pi / 180) *\
                     np.sin(dPA * np.pi / 180 / 2) * 4 * np.pi

        # We need to mask the nan values for averaging
        mask = np.isnan(flux_rebin)
        masked_arr = ma.masked_where(mask, flux_rebin)
        self.data_pro = masked_arr.copy()

        # print(self.data_pro.shape, len(self.times), solidAngle.shape)
        self.energytime = pd.DataFrame(
            index=self.times,
            columns=['Ebin_{:d}'.format(x) for x in range(self.nEbins)],
            data=np.average(self.data_pro, axis=0).T)
        
        self.energytime_sigma = pd.DataFrame(
            index=self.times,
            columns=['Ebin_{:d}'.format(x) for x in range(self.nEbins)],
            data=np.sqrt(np.nansum(self.data_sigma**2, axis=0)).T / flux.shape[0])
            # Combine uncertainty over all look-directions

        self.pitchangletime = pd.DataFrame(
            index=self.times,
            columns=['PAbin_{:d}'.format(x) for x in range(self.nPAbins)],
            data=np.average(self.data_pro, axis=1, weights=energybinsize).T)


        # Calculate equatorial pitch angles if database file is provided 
        if not self.dbfile:
            return None

        self.pa_database = pd.read_hdf(self.dbfile)
        self.eq_pitchanglebins = np.linspace(0, 90, self.nEQPAbins + 1)
        
        self.eqpa = pitchangles * 1.0

        count = 0
        for t in self.times:

            Beq = np.nanmean(
                self.pa_database.loc[
                    t - dt.timedelta(minutes=2.5): 
                    t + dt.timedelta(minutes=2.5), 'Beq'].values)

            # Blocal = Bmag[count]

            Blocal = np.nanmean(
                self.pa_database.loc[
                    t - dt.timedelta(minutes=2.5): 
                    t + dt.timedelta(minutes=2.5), 'Blocal_model'].values)

            if Blocal != Blocal:
                self.eqpa[:, :, count] = -1.
            elif Beq > Blocal:
                self.eqpa[:, :, count] = np.arcsin(
                    np.sin(self.eqpa[:, :, count] * np.pi / 180)) * 180 / np.pi
            else:
                self.eqpa[:, :, count] = np.arcsin(
                    np.sin(self.eqpa[:, :, count] * np.pi / 180) * \
                    np.sqrt(Beq / Blocal)) * 180 / np.pi

            count += 1

        flux_rebin = np.zeros((self.nEQPAbins, nEorig, nTimes))

        positionJuno, _ = get_juno_position(
            self.times, coordsystem='JUNO_MAG_VIP4')
        higherPAmask = \
            (((positionJuno[:, 2] >=0) & (pitchangles <= 90)) | \
            ((positionJuno[:, 2] < 0) & (pitchangles >= 90))) 
            # This expands over the time axis

        for ii in range(self.nEQPAbins):
            
            mask = (self.eq_pitchanglebins[ii] <= self.eqpa) & \
                   (self.eq_pitchanglebins[ii+1]> self.eqpa)
            
            flux_rebin[ii, :, :] = ma.masked_where(
                ~mask | higherPAmask, flux).mean(axis=0).filled(np.NaN)


        # We need to mask the nan values for averaging
        mask = np.isnan(flux_rebin)
        masked_eqpa = ma.masked_where(mask, flux_rebin)
        self.data_eqpa = masked_eqpa.copy()
    
        self.eqpad = pd.DataFrame(
            index=self.times,
            columns=['EQPAbin_{:d}'.format(x) for x in range(self.nEQPAbins)],
            data=np.average(self.data_eqpa, axis=1, weights=energybinsize).T)

        self.hasEQPAD = True


    def resample(self):
        """
        Resample the time-series DataFrames according to the dtbin specified 
        by the user. 
        """
        self.energytime = self.energytime\
            .resample(self.dtbin, offset=self.dtbin/2)\
            .agg(np.nanmean)

        self.pitchangletime = self.pitchangletime\
            .resample(self.dtbin, offset=self.dtbin/2)\
            .agg(np.nanmean)
        
        if self.dbfile:
            self.eqpad = self.eqpad\
                .resample(self.dtbin, offset=self.dtbin/2)\
                .agg(np.nanmean)

        if self.hasCXrate:
            self.cxrate = self.cxrate\
                .resample(self.dtbin, offset=self.dtbin/2)\
                .agg(np.nanmean)
            # if self.isHighRate:
            #     self.cxrate = self.cxrate\
            #         .resample(dt.timedelta(seconds=30), offset=dt.timedelta(seconds=15))\
            #         .agg(np.nansum)\
            #         .resample(self.dtbin, offset=self.dtbin/2)\
            #         .agg(np.nanmean)
            # else:
            #     self.cxrate = self.cxrate\
            #         .resample(self.dtbin, offset=self.dtbin/2)\
            #         .agg(np.nanmean)

        # if self.hasEflux: 
        #     self.energyflux = self.energyflux\
        #         .resample(self.dtbin, offset=self.dtbin/2)\
        #         .agg(np.nanmean)


    def calc_cx_rate(self):
        """
        Integrate the intensity over solid angles
        i.e. convert from 
        [s-1 cm-2 str-1 keV-1] to [s-1 cm-2 keV-1].
        """
        assert (self.units == 'FLUX')
        assert (self.hasValidData and self.hasEQPAD)
        
        if self.cxtype == 'Hplus-O':
            cxcs = cxcs_Hplus_O
        elif self.cxtype == 'Hplus-H2':
            cxcs = cxcs_Hplus_H2

        self.cxrate = pd.DataFrame(
            index=self.times,
            columns=['Ebin_{:d}'.format(x) for x in range(self.nEbins)])

        solidAngle = 4 * np.pi * \
            np.sin(
                (self.eq_pitchanglebins[:-1] + 
                 self.eq_pitchanglebins[1:]) * 0.5 * np.pi / 180) * \
            np.sin(
                (self.eq_pitchanglebins[1] - 
                 self.eq_pitchanglebins[0]) * 0.5 * np.pi / 180)

        
        for i in range(self.nEbins):

            Eeff = np.sqrt(self.energybins[i] * self.energybins[i + 1])
            
            dE = self.energybins[i + 1] - self.energybins[i]

            intensity2cxrate = \
                solidAngle * cxcs(Eeff * 1e-3) * 1e-16 * 1e-3 * dE * 2
                # This has shape (nEQPAbins)

                # Multiply by 2 to get full sky

            self.cxrate.loc[:, 'Ebin_{:d}'.format(i)] = \
                np.nansum(
                    self.data_eqpa[:, i, :].T * intensity2cxrate, axis=1)

                # np.dot(intensity2cxrate, 
                #            self.data_eqpa[:, i, :].filled(0))
                    # filling nans in this step is not accurate because it 
                    # undercounts the equatorial pitch angles that are not 
                    # observed. 

        self.cxrate = self.cxrate.astype(float)
        
        self.hasCXrate = True


    def plot_spectra(
        self, ax, data='Energy', donorm=False,
        vmin=None, vmax=None, cmap=None, logplot=True,
        label=None, zorder=1, plotcolorbar=True):
        """
        Generic method to plot a timeseries spectra.
        """
        if data == 'Energy':
            ax.set_yscale('log')
            ax.set_yticks([0.1, 1, 10])
            df = self.energytime
            bins = self.energybins / 1e3
        elif data == 'PAD':
            ax.set_yticks([0, 90, 180])
            df = self.pitchangletime
            bins = self.pitchanglebins 
        elif data == 'EQPAD':
            ax.set_yticks([0, 45, 90])
            df = self.eqpad
            bins = self.eq_pitchanglebins 
        elif data == 'cxrate':
            ax.set_yscale('log')
            ax.set_yticks([0.1, 1, 10])
            df = self.cxrate
            bins = self.energybins / 1e3
        
        # x-data , shifted left by half of dtbin 
        timearray = (df.index - self.dtbin / 2).to_pydatetime()
        timearray = np.append(
            timearray, timearray[-1] + (timearray[1] - timearray[0]))

        # y-data is just bins 

        # z-data is the values of the timeseries DataFrame (transposed)
        z = df.values.T

        if donorm:
            z /= np.nanmean(z, axis=0)
            
        if logplot:
            norm = colors.LogNorm(vmin=vmin, vmax=vmax)
        else:
            norm = colors.Normalize(vmin=vmin, vmax=vmax)

        mp = ax.pcolormesh(
            timearray, bins, z, norm=norm, cmap=cmap, zorder=zorder)

        ax.set_ylim([bins[0], bins[-1]])

        if plotcolorbar:
            bbox = ax.get_position()
            cax = ax.figure.add_axes(
                [bbox.x1 + 0.01, bbox.y0, 0.02, 0.8 * (bbox.y1- bbox.y0)])
            cbar = ax.figure.colorbar(mp, cax=cax)

        if label:
            ax.text(0.1, 1.05, label, transform=ax.transAxes)

        return ax, mp
    
    
    def find_valid_data_time(self):
        """Find times with valid data."""
        val = self.eqpad.values
        count = 0
        for i in range(val.shape[0]):
            if np.count_nonzero(np.isnan(val[i, :])) < 2:
                count += 1

        return (count * (self.eqpad.index[1] - self.eqpad.index[0])).total_seconds()


    def calc_energy_flux(self, pa1=0, pa2=180):
        """
        Calculate the upward/downward energy flux based on the formula given
        in Mauk et al. (2014).
        """
        midpa = 0.5 * (self.pitchanglebins[1:] + self.pitchanglebins[:-1])
        condition = (midpa >= pa1) & (midpa <= pa2)
        
        self.energyflux = pd.Series(index=self.times, data=0.)
        for i in range(self.nEbins):

            Eeff = np.sqrt(self.energybins[i] * self.energybins[i + 1])
            dE = self.energybins[i + 1] - self.energybins[i]

            conversionfactor = Eeff * dE * e * 1e4 * np.pi

            self.energyflux += np.sum(
                self.data_pro[condition, i, :].T * conversionfactor, axis=1)
            
        self.energyflux[self.energyflux == 0] = np.NaN 
        self.hasEflux = True

    

def merge_two_timeseries(df1, df2):
    df = pd.concat([df1, df2])
    return df.groupby(df.index).agg(lambda x: np.nanmean(x))

def merge_two_timeseries_quadrature(df1, df2):
    df = pd.concat([df1, df2])
    return df.groupby(df.index).agg(lambda x: np.sqrt(np.nansum(x**2)) / len(x))

def jade_merge(j1, j2):
    j3 = JadeData()

    assert (j1.hasValidData and j2.hasValidData) ; j3.hasValidData = True
    assert (j1.efilter == j2.efilter)       ; j3.efilter = j1.efilter 
    assert (j1.dtbin == j2.dtbin)           ; j3.dtbin = j1.dtbin
    assert (j1.dbfile == j2.dbfile)         ; j3.dbfile = j1.dbfile 
    assert (j1.nEbins == j2.nEbins)         ; j3.nEbins = j1.nEbins
    assert (j1.nPAbins == j2.nPAbins)       ; j3.nPAbins = j1.nPAbins 
    assert (j1.nEQPAbins == j2.nEQPAbins)   ; j3.nEQPAbins = j1.nEQPAbins
    assert (j1.units == j2.units)           ; j3.units = j1.units

    j3.energybins = j1.energybins
    j3.pitchanglebins = j1.pitchanglebins
    j3.energytime = merge_two_timeseries(j1.energytime, j2.energytime)
    j3.pitchangletime = merge_two_timeseries(
        j1.pitchangletime, j2.pitchangletime)
    
    # Merge data_sigma in quadrature
    j3.energytime_sigma = merge_two_timeseries_quadrature(
        j1.energytime_sigma, j2.energytime_sigma)

    if j1.hasEQPAD and j2.hasEQPAD:
        j3.eq_pitchanglebins = j1.eq_pitchanglebins
        j3.eqpad = merge_two_timeseries(j1.eqpad, j2.eqpad)
        j3.hasEQPAD = True

    if j1.hasCXrate and j2.hasCXrate:
        j3.cxrate = merge_two_timeseries(j1.cxrate, j2.cxrate)
        j3.cxrate.replace(0, np.NaN, inplace=True)
        j3.hasCXrate = True

    if j1.hasEflux and j2.hasEflux:
        j3.energyflux = merge_two_timeseries(j1.energyflux, j2.energyflux)
        j3.hasEflux = True

    return j3


def get_jade_data(start, end, species=3, **kwargs):
    """
    This method will return JADE data combined for high/low res between the 
    start and end times. 
    """
    if isinstance(start, dt.datetime):
        ts = start.strftime('%Y-%jT%H:%M:%S')
        te = end.strftime('%Y-%jT%H:%M:%S')
    else:
        ts = start
        te = end

    kwargs['species'] = species
    getHRS = kwargs.pop('getHRS', True)
    getLRS = kwargs.pop('getLRS', True)

    # with contextlib.redirect_stdout(None):
    idl = pidly.IDL('/Users/ysar/Applications/harris/idl88/bin/idl')
    idl(".run set_jade_vars")
    idl("set_jade_vars")
    idl("CD, '/Users/ysar/Documents/code/JADE/'")

    filename = '/Users/ysar/Documents/code/JADE/jade-yash/tmp/temp.sav'

    if species == -1:
        datatype = 'elc'
    else:
        datatype = 'ion'

    if getHRS:
        idl("data = jade_js_get_sci(['{:s}', '{:s}']"
            ",type='{:s}',hrs=1,species={:d})".format(ts, te, datatype, species))
        
        if kwargs['units'] == 'PSD':
            if species == -1:
                idl("data = jade_convert_data_units_from_cnt(data, UNITS='PSD')")
            else:
                idl("data = jade_convert_data_units_from_cnt(data, UNITS='PSD', m=1, q=1)")
        elif kwargs['units'] == 'FLUX':
            if species == -1:
                idl("data = jade_convert_data_units_from_cnt(data, UNITS='DNF')")
            else:
                idl("data = jade_convert_data_units_from_cnt(data, UNITS='DNF', q=1)")
        elif kwargs['units'] == 'DEF':
            idl("data = jade_convert_data_units_from_cnt(data, UNITS='DEF')")

        idl("save, data, FILENAME = '{:s}'".format(filename))
        
        data_highres = JadeData(filename=filename, isHighRate=True, **kwargs)
    
    if getLRS:
        idl("data = jade_js_get_sci(['{:s}', '{:s}']"
            ",type='{:s}',species={:d})".format(ts, te, datatype, species))
        
        if kwargs['units'] == 'PSD':
            if species == -1:
                idl("data = jade_convert_data_units_from_cnt(data, UNITS='PSD')")
            else:
                idl("data = jade_convert_data_units_from_cnt(data, UNITS='PSD', m=1, q=1)")
        elif kwargs['units'] == 'FLUX':
            if species == -1:
                idl("data = jade_convert_data_units_from_cnt(data, UNITS='DNF')")
            else:
                idl("data = jade_convert_data_units_from_cnt(data, UNITS='DNF', q=1)")
        elif kwargs['units'] == 'DEF':
            idl("data = jade_convert_data_units_from_cnt(data, UNITS='DEF')")


        idl("save, data, FILENAME = '{:s}'".format(filename))

        data_lowres = JadeData(filename=filename, **kwargs)
    
    if getHRS and getLRS:
        if data_highres.hasValidData and data_lowres.hasValidData:
            return jade_merge(data_highres, data_lowres)
        elif data_highres.hasValidData:
            return data_highres
        elif data_lowres.hasValidData:
            return data_lowres
        else:
            return JadeData()
    elif getHRS and data_highres.hasValidData:
        return data_highres
    elif getLRS and data_lowres.hasValidData:
        return data_lowres
    else:
        return JadeData()
    
# -----------------------------------------------------------------------------
class JadeMomentsData(object):
    """
    This class reads in the .sav files written out by the IDL tools for JADE
    moments. The .sav files are read into numpy recarray (by default). The 
    purpose of this class is to convert from recarray into pandas DataFrames.
    """
    def __init__(self, **kwargs):
        
        self.filename = kwargs.pop('filename', None)
        if not self.filename:
            return None
        
        check = self.read_moments_file()
        if check == -1:
            return
        
        self.hasValidData = True
        self.reorganize_data()


    def read_moments_file(self):
        """Read in the IDL .sav file."""
        self.data_raw = readsav(self.filename)['momentsdata']

        if isinstance(self.data_raw, np.int16) and self.data_raw == -1:
            return -1

        df = pd.DataFrame()
        df['UTC'] = self.data_raw['UTC'][0]

        df['time'] = df.apply(
            lambda x: dt.datetime.strptime(
                x['UTC'].decode('UTF-8'), '%Y-%jT%H:%M:%S.%f'), axis=1)

        self.times = df['time'].squeeze()
        return 0
    

    def reorganize_data(self):
        """Reorganize data from recarrays to DataFrames."""

        self.data = pd.DataFrame(
            index=self.times,
            columns=['Density', 'Pressure', 'Temp', 'Vmag', 'Vx_JSS', 'Vy_JSS', 'Vz_JSS']
            )
    
        self.data['Density']  = self.data_raw['N_CC'][0] * 1.0
        self.data['Pressure'] = self.data_raw['PRESSURE_PA'][0] * 1.0
        self.data['Temp']     = self.data_raw['TEMP_EV'][0] * 1.0
        self.data['Vmag']     = self.data_raw['VMAG_JSS_KMPS'][0] * 1.0
        self.data['Vx_JSS']   = self.data_raw['V_JSSXYZ_KMPS'][0][0, :] * 1.0
        self.data['Vy_JSS']   = self.data_raw['V_JSSXYZ_KMPS'][0][1, :] * 1.0
        self.data['Vz_JSS']   = self.data_raw['V_JSSXYZ_KMPS'][0][2, :] * 1.0

        # Uncertainties
        self.data['SigmaDensity']  = self.data_raw['N_SIGMA_CC'][0] * 1.0
        self.data['SigmaPressure'] = self.data_raw['PRESSURE_SIGMA_PA'][0] * 1.0
        self.data['SigmaTemp']     = self.data_raw['TEMP_SIGMA_EV'][0] * 1.0
        self.data['SigmaVmag']     = \
            self.data_raw['VMAG_JSS_SIGMA_KMPS'][0] * 1.0
        self.data['SigmaVx_JSS']   = \
            self.data_raw['V_JSSXYZ_SIGMA_KMPS'][0][0, :] * 1.0
        self.data['SigmaVy_JSS']   = \
            self.data_raw['V_JSSXYZ_SIGMA_KMPS'][0][1, :] * 1.0
        self.data['SigmaVz_JSS']   = \
            self.data_raw['V_JSSXYZ_SIGMA_KMPS'][0][2, :] * 1.0

 
def merge_moments_data(j1, j2):
    """Merge two instances of JadeMomentsData"""
    assert (j1.hasValidData and j2.hasValidData) ; j3.hasValidData = True
    
    j3 = JadeMomentsData()
    j3.data = merge_two_timeseries(j1.data, j2.data)


def get_jade_moments(start, end, species=3, **kwargs):
    """
    Use the IDL JADE tools to calculate moments for specific species.
    """
    if isinstance(start, dt.datetime):
        ts = start.strftime('%Y-%jT%H:%M')
        te = end.strftime('%Y-%jT%H:%M')
    else:
        ts = start
        te = end

    kwargs['species'] = species
    getHRS = kwargs.pop('getHRS', True)
    getLRS = kwargs.pop('getLRS', True)

    # with contextlib.redirect_stdout(None):
    idl = pidly.IDL('/Users/ysar/Applications/harris/idl88/bin/idl')
    idl(".run set_jade_vars")
    idl("set_jade_vars")
    idl("CD, '/Users/ysar/Documents/code/JADE/'")

    filename = '/Users/ysar/Documents/code/JADE/jade-yash/tmp/temp.sav'

    if species == -1:
        datatype = 'elc'
    else:
        datatype = 'ion'

    if getHRS:
        idl("data = jade_js_get_sci(['{:s}', '{:s}']"
            ",type='{:s}',hrs=1,species={:d})".format(ts, te, datatype, species))
        
        if species == -1:
            idl("momentsdata = jade_numerical_moments_elc_pa(data)")
        else:
            idl("momentsdata = jade_numerical_moments(data)")

        idl("save, momentsdata, FILENAME = '{:s}'".format(filename))
        
        moments_highres = JadeMomentsData(filename=filename, isHighRate=True, **kwargs)
    
    if getLRS:
        idl("data = jade_js_get_sci(['{:s}', '{:s}']"
            ",type='{:s}',species={:d})".format(ts, te, datatype, species))
        
        if species == -1:
            idl("momentsdata = jade_numerical_moments_elc_pa(data)")
        else:
            idl("momentsdata = jade_numerical_moments(data)")

        idl("save, momentsdata, FILENAME = '{:s}'".format(filename))
        
        moments_lowres = JadeMomentsData(filename=filename, isHighRate=True, **kwargs)
    
        
    
    if getHRS and getLRS:
        if moments_highres.hasValidData and moments_lowres.hasValidData:
            return merge_moments_data(moments_highres, moments_lowres)
        elif moments_highres.hasValidData:
            return moments_highres
        elif moments_lowres.hasValidData:
            return moments_lowres
        else:
            return JadeMomentsData()
    elif getHRS and moments_highres.hasValidData:
        return moments_highres
    elif getLRS and moments_lowres.hasValidData:
        return moments_lowres
    else:
        return JadeMomentsData()
        
# -----------------------------------------------------------------------------