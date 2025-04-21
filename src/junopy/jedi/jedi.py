"""
junopy - A Python API for the Juno spacecraft's particle and fields datasets.
Copyright (C) 2018  Yash Sarkango

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
"""
import os
import numpy as np
import pandas as pd
import datetime as dt
import glob

import numpy.ma as ma

import matplotlib.colors as colors

from scipy.constants import m_e, m_p, e

from JunoPy.spice import get_juno_position
from JunoPy.library import cxcs_Hplus_O, cxcs_Hplus_H2
from JunoPy.library import cxcs_Hplus_Oplus, cxcs_Hplus_Oplusplus
from JunoPy.library import cxcs_Hplus_Splus, cxcs_Hplus_Splusplus


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
    hasCXrate = False
    nChannels = -1
    hasEflux = False

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
        self.stop     = kwargs.pop('stop', None)
        self.nEQPAbins = kwargs.pop('nEQPAbins', 9)
        self.dosaveonly = kwargs.pop('dosaveonly', False)
        self.docalcCX   = kwargs.pop('docalcCX', False)
        self.cxtype    = kwargs.pop('cxtype', 'Hplus-O')
        self.calcEflux = kwargs.pop('calcEflux', False)
        self.paMin     = kwargs.pop('paMin', 0)
        self.paMax     = kwargs.pop('paMax', 180)

        if self.species == 'PROTON_TOFxPH':
            self.species = 'PROTON'

        if '.hdf' in self.filename:
            self.read_hdf_file()
        else:
            self.read_raw_csv_file()
            self.convert_to_hdf_file()

        if self.dosaveonly: return
        else:
            self._filter_metadata()

        if self.start and self.stop:
            self.df_raw = self.df_raw.loc[
                self.start.strftime('%Y-%m-%d %H:%M:%S'):
                self.stop.strftime('%Y-%m-%d %H:%M:%S')]

        self.mask_invalid_pa()

        if self.df_raw.empty:
            return
        else:
            self.hasValidData = True

        if self.units == 'PSD':
            self.calc_phase_space_density()

        self.create_energytime_df()
        self.create_pitchangle_df()

        if self.dbfile:
            self.calc_equatorial_pad()
            self.create_pitchangle_df(eqpad=True)

        if self.docalcCX:
            self.calc_cx_rate()

        if self.calcEflux:
            self.calc_energy_flux(
                pa1=self.paMin, pa2=self.paMax)

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

        df.reset_index(inplace=True, drop=True)

        # Filter unwanted energy channels
        if self.species == 'ELECTRON':
            for i in range(4):
                filter = df['channel'][i][-7:]
                df = df[~df['channel'].str.contains(filter)]

        if self.species == 'PROTON':
            for i in range(1):
                filter = df['channel'][i][-7:]
                df = df[~df['channel'].str.contains(filter)]

        df.reset_index(inplace=True, drop=True) # Have to reset index with every df change

        # A better way to filter energy channels so it doesnt
        # lead to unequal distribution across SSDs
        if self.efilter:

            nChannels = int(len(df.index) / self.nSSD)

            # Calculate mean energy across all SSDs

            Emean = np.zeros((6, nChannels))

            for i in range(self.nSSD):

                dftmp = df[df['channel'].str[:2] == 'T{:d}'.format(i)]
                Emean[i, :] = np.sqrt(dftmp['Elow'] * dftmp['Ehigh']).values

            Emean = Emean.mean(axis=0)

            invalidChannels = np.asarray(
                (Emean < self.efilter[0]) | (Emean >= self.efilter[1]))\
                .nonzero()[0]

            for i in invalidChannels:
                filter = df['channel'][i][-7:]
                df = df[~df['channel'].str.contains(filter)]

        df.reset_index(inplace=True, drop=True)

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

        # Conversion faction from J to PSD
        # f = m / v^2 j
        #   = [kg] / []

        df['j2psd'] = \
            (m**2 / e**2) / np.sqrt(df['Ehigh'] * df['Elow']) \
                / (2 * 1e14)


        self.metadata = df.copy()

        # Specify negative values as 0 or no-count
        vals = self.df_raw.loc[:, self.metadata['channel']].values
        vals[vals < 0] = 0.
        self.df_raw.loc[:, self.metadata['channel']] = vals

        if len(df.index) % self.nSSD == 0:
            self.nChannels = int(len(df.index) / self.nSSD)
        else:
            raise Exception('Unequal energy channels within each SSD.')


    def mask_invalid_pa(self):
        """
        Mask all values outside the valid pitch angle range.
        """
        for i in range(self.nSSD):
            channels = self.metadata.loc[
                self.metadata['channel'].str.contains('T{:d}'.format(i)),
                'channel']

            condition = (
                (self.df_raw['Pitch_Angle_T{:d}'.format(i)] < self.paMin) |
                (self.df_raw['Pitch_Angle_T{:d}'.format(i)] > self.paMax))

            self.df_raw.loc[condition, 'Pitch_Angle_T{:d}'.format(i)] = np.nan



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
        for _, row in self.metadata.iterrows():
            df[row['channelpsd']] = \
                self.df_raw[row['channel']] * row['j2psd']

        self.df_raw = pd.concat([self.df_raw, df], axis=1)

        self.metadata.drop(labels='channel', axis=1, inplace=True)
        self.metadata.rename(columns={'channelpsd': 'channel'}, inplace=True)


    def create_energytime_df(self):
        """
        Create the energy-time spectra in the form of a dataframe.
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

        for j in range(self.nEbins):

            meta = self.metadata.copy()

            meta['deltaE'] = meta['Ehigh'] - meta['Elow']

            # Prevalence of each channel in this energy bins.
            meta['dE'] = np.minimum(meta['Ehigh'], self.energybins[j + 1]) - \
                         np.maximum(meta['Elow'],  self.energybins[j])


            meta = meta[~(meta['dE'] < 0)]

            # Now we know the energy channels that are present within the j-th
            # energy bin. Knowing this, we can identify the pitch
            # angle associated with each energy channel.

            denom = np.zeros(len(self.df_raw.index))
            val = np.zeros(len(self.df_raw.index))

            for _, row in meta.iterrows():

                scope = row['channel'][:2]    # Get telescope ID for E channel
                pa = self.df_raw.loc[:, 'Pitch_Angle_{:s}'.format(scope)]\
                    .astype(float).values

                # This is the solid angle of the measurement at a particular
                # pitch angle, given that JEDI-I has an angular resolution of
                # approximately 26.7 degrees.
                deltaSigma = 4 * np.pi * np.sin(pa * np.pi / 180) \
                             * np.sin(13.35 * np.pi / 180)

                # Need to ignore NaN values while summing.
                deltaSigma[np.isnan(deltaSigma)] = 0
                intensities = self.df_raw.loc[:, row['channel']].fillna(0)

                # deltaSigma[intensities == 0] = 0

                val += intensities * row['dE'] # * deltaSigma
                denom += row['dE'] # * deltaSigma

                # How to handle nans?

            condition =  (denom > 0)

            self.energytime.loc[condition, 'Ebin_{:d}'.format(j)] = \
                val[condition] / denom[condition]

            self.energytime.loc[~condition, 'Ebin_{:d}'.format(j)] = np.nan


    def create_pitchangle_df(self, eqpad=False):
        """
        Create pitchangle spectra by averaging over energy channels and
        accounting for their width.
        """
        if eqpad:
            s = 'EQ'
            nbins = self.nEQPAbins
            bins = np.linspace(0, 90, nbins + 1)
        else:
            s = ''
            nbins = self.nPAbins
            bins = np.linspace(0, 180, nbins + 1)

        df = pd.DataFrame(
            data=0.,
            index=self.df_raw.index,
            columns=\
                ['PAbin_{:d}'.format(x) for x in range(nbins)])

        df_nSSD = pd.DataFrame(
            index=self.df_raw.index,
            columns=['nSSDbin_{:d}'.format(x) for x in range(nbins)],
            data=0.)

        datacols = []    # Identify data columns for each SSD beforehand
        energybinsize = [] # np.zeros((len(self.index), self.nSSD))
        for j in range(self.nSSD):

            condition = self.metadata['channel'].str.contains('T{:d}'.format(j))
            datacols.append(
                self.metadata.loc[condition, 'channel'])

            energybinsize.append(
                self.metadata.loc[condition, 'Ehigh'] -
                self.metadata.loc[condition, 'Elow'])

        for i in range(nbins):
            for j in range(self.nSSD):

                condition = (
                (self.df_raw['{:s}Pitch_Angle_T{:d}'.format(s, j)] > bins[i]) &
                (self.df_raw['{:s}Pitch_Angle_T{:d}'.format(s, j)] <=
                    bins[i + 1]))

                # Average all channels for one SSD and consider width of
                # energy bin as weight.
                val = self.df_raw.loc[condition, datacols[j]].values

                df.loc[condition, 'PAbin_{:d}'.format(i)] += \
                    np.average(val, axis=1, weights=energybinsize[j])

                # Count number of SSDs in this PA bin
                df_nSSD.loc[condition, 'nSSDbin_{:d}'.format(i)] += 1

            # If multiple SSDs in the bin, divide by this number
            condition = (df_nSSD['nSSDbin_{:d}'.format(i)] > 0)

            df.loc[condition, 'PAbin_{:d}'.format(i)] /= \
                df_nSSD.loc[condition, 'nSSDbin_{:d}'.format(i)]

            df.loc[~condition, 'PAbin_{:d}'.format(i)] = np.nan

        if eqpad:
            self.eqpad = df
            self.eq_pitchanglebins = bins
            self.hasEQPAD = True
        else:
            self.pitchangletime = df
            self.pitchanglebins = bins


    def create_energy_pitchangle_array(self):
        """
        Create an array analogous to the JADE data structure organized by
        (pitchangle, energy, time)
        """
        data = np.zeros((
            len(self.pitchanglebins) - 1,
            len(self.nEbins) - 1,
            len(self.df_raw.index)))

        for i, row in self.df_raw.iterrows():
            pass




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

        pavalues[((positionJuno[:, 2] >= 0) & (pavalues <= 90).T).T] = np.nan
        pavalues[((positionJuno[:, 2] <  0) & (pavalues >= 90).T).T] = np.nan

        count = 0
        for t in self.df_raw.index:

            db = self.pa_database.loc[
                t - dt.timedelta(minutes=2.5):
                t + dt.timedelta(minutes=2.5)].copy()

            if len(db.index) == 0:
                self.df_raw.loc[t, eqpacols] = np.nan
                continue

            Beq = np.nanmean(db['Beq'].values)

            Blocal = np.nanmean(db['Blocal_model'].values)

            if Beq > Blocal:
                self.df_raw.loc[t, eqpacols] = np.arcsin(
                    np.sin(pavalues[count, :] * np.pi / 180)) * 180 / np.pi
            else:
                self.df_raw.loc[t, eqpacols] = np.arcsin(
                    np.sin(pavalues[count, :] * np.pi / 180) * \
                    np.sqrt(Beq / Blocal)) * 180 / np.pi

            count += 1


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

                # .resample(dt.timedelta(seconds=30), offset=dt.timedelta(seconds=15))\
                # .agg(np.nansum)\

                # We assume that every JEDI measurement measures the full-sky
                # distribution. Unlike JADE, JEDI-90 and 270 sample the full-sky
                # provided the magnetic field is in the spin plane.

        # if self.hasEflux:
        #     self.energyflux = self.energyflux\
        #         .resample(self.dtbin, offset=self.dtbin/2)\
        #         .agg(np.nanmean)


    def calc_cx_rate(self):
        """
        Integrate the intensity over solid angles
        i.e. convert from
        [s-1 cm-2 str-1 keV-1] to [s-1 cm-2 keV-1].
        We will need to use df_raw as the pitchangle dataframe is reduced in
        energy space already.
        """
        assert (self.units == 'FLUX')

        cxrate = pd.DataFrame(
            0,
            index=self.df_raw.index,
            columns=['Ebin_{:d}'.format(x) for x in range(self.nEbins)])

        if self.cxtype == 'Hplus-O':
            cxcs = cxcs_Hplus_O
        elif self.cxtype == 'Hplus-H2':
            cxcs = cxcs_Hplus_H2
        elif self.cxtype == 'Hplus-Oplus':
            cxcs = cxcs_Hplus_Oplus
        elif self.cxtype == 'Hplus-Oplusplus':
            cxcs = cxcs_Hplus_Oplusplus
        elif self.cxtype == 'Hplus-Splus':
            cxcs = cxcs_Hplus_Splus
        elif self.cxtype == 'Hplus-Splusplus':
            cxcs = cxcs_Hplus_Splusplus

        for j in range(self.nEbins):
            # For each energy bin, calculate Flux

            meta = self.metadata.copy()

            meta['deltaE'] = meta['Ehigh'] - meta['Elow']
            meta['dE'] = np.minimum(meta['Ehigh'], self.energybins[j + 1]) - \
                         np.maximum(meta['Elow'],  self.energybins[j])

            meta = meta[~(meta['dE'] < 0)]

            # Now we know the energy channels that are present within the j-th
            # energy bin. Knowing this, we can identify the equatorial pitch
            # angle associated with each energy channel.

            dEsum = np.zeros(len(self.df_raw.index))

            for _, row in meta.iterrows():

                scope = row['channel'][:2]    # Get telescope ID for E channel
                eqpa = self.df_raw.loc[:, 'EQPitch_Angle_{:s}'.format(scope)]\
                    .astype(float).values

                # This is the solid angle of the measurement at a particular
                # pitch angle, given that JEDI-I has an angular resolution of
                # approximately 26.7 degrees.
                deltaSigma = 4 * np.pi * np.sin(eqpa * np.pi / 180) \
                             * np.sin(13.35 * np.pi / 180)

                # Need to ignore NaN values while summing
                # (that is done automatically)
                dEsum[~np.isnan(deltaSigma)] += row['dE']
                deltaSigma[np.isnan(deltaSigma)] = 0

                intensities = self.df_raw.loc[:, row['channel']] # .fillna(0)

                # deltaFlux = Intensity * deltaSolidAngle
                Eeff = np.sqrt(row['Ehigh'] * row['Elow'])
                cxrate['Ebin_{:d}'.format(j)] += intensities * deltaSigma * \
                    cxcs(Eeff) * 1e-16 * row['deltaE'] * row['dE'] * 2
                    # Multiply by 2 to account for full-sky

            cxrate['Ebin_{:d}'.format(j)] /= dEsum # meta['dE'].sum()

        cxrate.replace(0, np.nan, inplace=True)
        self.cxrate = cxrate
        self.hasCXrate = True


    def plot_spectra(
        self, ax, data='Energy', donorm=False,
        vmin=None, vmax=None, cmap=None, logplot=True, label=None):
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
        elif data == 'cxrate':
            ax.set_yscale('log')
            ax.set_yticks([100, 1000, 10000])
            df = self.cxrate
            bins = self.energybins


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

        mp = ax.pcolormesh(timearray, bins, z, norm=norm, cmap=cmap)

        ax.set_ylim([bins[0], bins[-1]])

        bbox = ax.get_position()
        cax = ax.figure.add_axes(
            [bbox.x1 + 0.01, bbox.y0, 0.02, 0.8 * (bbox.y1- bbox.y0)])
        cbar = ax.figure.colorbar(mp, cax=cax)

        if label:
            ax.text(0.1, 1.05, label, transform=ax.transAxes)

        return ax


    def find_valid_data_time(self):
        """Find times with valid data."""
        val = self.eqpad.values
        count = 0
        for i in range(val.shape[0]):
            if np.count_nonzero(np.isnan(val[i, :])) < 2:
                count += 1

        return (count * (self.eqpad.index[1] - self.eqpad.index[0])).total_seconds()


    def calc_energy_flux(self, pa1=0, pa2=15):
        """
        Calculate the upward/downward energy flux based on the formula given
        in Mauk et al. (2014).
        pa1, pa2 are in degrees.
        """

        self.energyflux = pd.Series(index=self.df_raw.index, data=0.)

        for j in range(self.nSSD):

            condition = (
                (self.df_raw['Pitch_Angle_T{:d}'.format(j)] >= pa1) &
                (self.df_raw['Pitch_Angle_T{:d}'.format(j)] <= pa2)
            )

            channels = self.metadata['channel'].str.contains('T{:d}'.format(j))
            datacols = self.metadata.loc[channels, 'channel']

            deltaE = self.metadata.loc[channels, 'Ehigh'] - \
                     self.metadata.loc[channels, 'Elow']

            midE = np.sqrt(
                self.metadata.loc[channels, 'Ehigh'] * \
                self.metadata.loc[channels, 'Elow']
            )

            for i, col in enumerate(datacols):

                val = self.df_raw.loc[:, col].values
                val[~condition] = 0.

                self.energyflux += val \
                    * deltaE.iloc[i] * midE.iloc[i] \
                    * np.pi * e * 1e10     # Unit conversion and
                                           # "Area-projected weighted size" ?

        self.energyflux[self.energyflux == 0] = np.nan
        self.hasEflux = True



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
    assert (j1.start == j2.start and j1.stop == j2.stop)
    assert (j1.nEQPAbins == j2.nEQPAbins) ; j3.nEQPAbins = j1.nEQPAbins

    j3.start = j1.start
    j3.stop = j1.stop

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

    if j1.hasCXrate and j2.hasCXrate:
        j3.cxrate = merge_two_timeseries(j1.cxrate, j2.cxrate)
        j3.hasCXrate = True

    if j1.hasEflux and j2.hasEflux:
        j3.energyflux = merge_two_timeseries(j1.energyflux, j2.energyflux)
        j3.hasEflux = True

    return j3


def get_jedi_data(**kwargs):
    """
    This script will return JEDI data from all three detectors for all times
    between start and end. The timeseries data from different days will be
    combined.
    """
    if 'ELECTRON' in kwargs['species']:
        idstring = 'ERSESP'
    elif kwargs['species'] == 'PROTON':
        idstring = 'ERSTOFXER'
    elif 'HEAVY' in kwargs['species']:
        idstring = 'NONPTOFXPHR'
    elif kwargs['species'] == 'PROTON_TOFxPH':
        idstring = 'ERSTOFXPHR'
    else:
        idstring = 'NONPTOFXER'

    datafolder = kwargs.pop('datafolder', 'data/jedi/')
    sensor = kwargs.pop('sensor', ['JED_090', 'JED_270', 'JED_180'])

    AllData = []

    for day in pd.date_range(kwargs['start'], kwargs['stop']):
        print(datafolder + day.strftime("/%Y/%j/") + '/*{:s}_CDR_{:s}_*.hdf5'.format(idstring, day.strftime('%Y%j')))
        dataFiles = glob.glob(
            datafolder + \
            day.strftime("/%Y/%j/") + \
            '/*{:s}_CDR_{:s}_*.hdf5'.format(idstring, day.strftime('%Y%j')))

        if 'JEDI-90' in sensor:
            dataFiles = [x for x in dataFiles if 'JED_090' in x]
        if 'JEDI-270' in sensor:
            dataFiles = [x for x in dataFiles if 'JED_270' in x]
        if 'JEDI-180' in sensor:
            dataFiles = [x for x in dataFiles if 'JED_180' in x]

        print(dataFiles)
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
