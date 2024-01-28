"""
This module can be used to read level-3 data files from the JADE instrument.

Usage:
    from .read import read_jade_binary_file
    filename = 'JAD_L30_HRS_ELC_TWO_CNT_2021201_V04.DAT.gz'
    df = read_jade_binary_file(filename)

The `read_jade_binary_file` method returns a `pandas` DataFrame that contains 
all information as in the .DAT file. Each row of the DataFrame is a different 
record (i.e. measurement time), and each column is a different value, e.g. 
'DATA', 'BACKGROUND', 'DIM1_E', etc. Note that the values within a cell of the 
DataFrame can be numpy arrays when applicable. 

The above example reads a file with the .DAT.gz expansion (compressed), but 
reading the .DAT (uncompressed) file is also possible.
 
"""
import datetime as dt
import numpy as np
import pandas as pd
import struct
import re
import gzip
import logging

from matplotlib import pyplot as plt
from matplotlib import colors
from numpy import ma

from JunoPy import mag
from .calib import get_Geff_electrons

def open_file(filename, *args):
    if '.gz' == filename[-3:]:
        return gzip.open(filename, *args)
    else:
        return open(filename, *args)


def read_jade_lbl_file(lblfile):
    """
    Read a JADE .LBL file to infer the column names, types, and byte sizes.
    
    Args:
        lblfile (str): Path to the .LBL file to be read.

    Returns:
        dfLabel (pandas.DataFrame): DataFrame with the contents of the .LBL 
            file containing the name, type, and dimensions of each column.
    """
    dfLabel = pd.DataFrame(columns=['name', 'dataID-np', 'nDim', 'shape'])

    with open_file(lblfile, 'rt') as f:

        i = 0

        for line in f:
            if 'FILE_RECORDS' in line:
                nRecords = int(re.search(r'(\d+)', line)[0])
                break
        
        for line in f:
            if '/* RJW' in line:
                
                if 'BYTES_PER_RECORD' in line:
                    continue
                elif 'OBJECTS_PER_RECORD' in line:
                    continue

                a = line.strip()[3:-3].split(',')[1:]
                dfLabel.at[i, 'name'] = a[0].strip()
                dfLabel.at[i, 'dataID-np'] = a[1].strip()
                dfLabel.at[i, 'nDim'] = int(a[2].strip())
                dfLabel.at[i, 'shape'] = np.array([int(x) for x in a[3:]])
                i += 1


    dfLabel['dataID-np'].replace('c', 's', inplace=True)

    dfLabel['nbytes'] = dfLabel.apply(
        lambda x: np.prod(x['shape']) * struct.calcsize(x["dataID-np"]), 
        axis=1)
        
    return dfLabel, nRecords


def _read_jade_bytes(filename, dfLabel, nRecords):
    """
    Reads in various bytes from the .DAT file and inserts them into a DataFrame.
    Does not decode or unpack the bytes, which is faster when done as an apply 
    to the entire column.

    Args:
        filename (str): The JADE .DAT file to be read
        dfLabel (`pandas.DataFrame`): 
            DataFrame with the contents of the .LBL file.
        nRecords (int): Total number of records ("rows") in the .DAT file.

    Returns:
        A `pandas.DataFrame` containing the contents of the .DAT file. Each row
        corresponds to a measurement interval and each column has the data
        associated with that interval. Contents of each cell of the DataFrame
        may be 1D arrays. For e.g., each value in the 'BACKGROUND' column is an
        array with shape `(m,n)` with m and n defined in the JADE SIS. 
        Currently, this is stored as a 1D array with size `(m*n)` (it is 
        reshaped later).
    """
    data = pd.DataFrame(columns=dfLabel['name'].to_list())

    # We will create numpy arrays to store data and then copy over the array to
    # the DataFrame. This is much faster than using .at() or .loc().
    datadict = {}
    for i, col in dfLabel.iterrows():
        datadict[col['name']] = np.empty(shape=(nRecords,), dtype='object')

    # Pre-making these lists is much faster than using .iterrows(). 
    n = [np.prod(x) for x in dfLabel['shape']]
    fmt = dfLabel['dataID-np'].to_list()
    nbytes = dfLabel['nbytes'].to_list()
    names = dfLabel['name'].to_list()


    # Loop over every record (i) and every column (j)
    # Unfortunately, so far I have not discovered a way to vectorize this.

    f = open_file(filename, 'rb')

    i = 0

    while True:
        for j in range(len(dfLabel)):
            readdata = f.read(nbytes[j])
            
            if readdata == b'':    # Check for end of file.
                for name in names: 
                    data[name] = datadict[name]
                f.close()
                return data

            datadict[names[j]][i] = struct.unpack(f'{n[j]}{fmt[j]}', readdata)
            
        i += 1
    

def read_jade_binary_file(filename):
    """
    Reads the JADE .DAT file into a single pandas DataFrame. This method assumes
    that the corresponding .LBL file is stored in the same directory as the .DAT
    file.

    Args:
        filename (str): Path to the .DAT file being read.

    Returns:
        df (pandas.DataFrame): DataFrame containing the information stored in
            the JADE .DAT file. All strings have been decoded using UTF-8 and 
            arrays have been reshaped according to the format specified in the 
            .LBL files.
    """
    if filename[-3:] == '.gz':
        lblfile = filename[:-7] + '.LBL'
    else:
        lblfile = filename[:-4] + '.LBL'

    print(lblfile)

    dfLabel, nRecords = read_jade_lbl_file(lblfile)
    
    df = _read_jade_bytes(filename, dfLabel, nRecords)

    for _, column in dfLabel.iterrows():

        if column['dataID-np'] == 's':
            df[column['name']] = df[column['name']].apply(
                lambda x: x[0].decode('UTF-8').strip())

        elif column['nDim'] > 1:
            df[column['name']] = df[column['name']].apply(
                lambda x: np.array(x).reshape(column['shape']))

        else:
            df[column['name']] = df[column['name']].apply(lambda x: x[0])

    return df



    # Return numpy array for compatibility with IDL
    datadict = {}
    for key in df.keys():
        datadict[key] = np.stack(df[key])

    return datadict


class JadeData(object):
    """
    This class reads in the JADE data from the .gz compressed files. It works 
    roughly like this - 
    1. Extract bytes from .gz file
    2. Read in bytes to a dictionary with each key corresponding to the column
       label
    3. Interprets this dictionary to make pitch angle and energy time spectra
    """
    hasEQPAD = False
    hasValidData = False
    hasCXrate = False
    hasEflux = False

    def __init__(self, **kwargs):
        
        self.filename = kwargs.pop('filename', None)
        if not self.filename:
            return None

        self.magfolder = kwargs.pop('magfolder', None)
        self.efilter   = kwargs.pop('efilter', None)
        self.dtbin     = kwargs.pop('dtbin', dt.timedelta(seconds=120))
        # self.nEbins    = kwargs.pop('nEbins', 20)    # Decrepated for now
        self.nPAbins   = kwargs.pop('nPAbins', 12)
        self.units     = kwargs.pop('units', 'COUNTS')
        self.calcEflux = kwargs.pop('calcEflux', False)
        self.paMin     = kwargs.pop('paMin', 0)
        self.paMax     = kwargs.pop('paMax', 180)

        if not self.magfolder:
            logging.error(
            "%name: JadeData needs magfolder= keyword to calculate pitch angles.")

        self.read_binary_file()
        self.calc_pitchangles()
        self.detect_sensor()
        self.Geff = get_Geff_electrons(
            self.data_raw['DIM1_E'], self.times, self.data_raw['DIM2_ELEVATION'])
        
        # self.convert_to_units()
        
        self.hasValidData = True
        self.calc_energy_spectra()
        self.calc_pa_spectra()

        if self.calcEflux:
            self.calc_energy_flux(
                pa1=self.paMin, pa2=self.paMax)
            
        # self.resample()

        
    def read_binary_file(self):
        # Read the binary files
        df = read_jade_binary_file(self.filename)
        self.data_raw = {}
        for key in df.keys():
            self.data_raw[key] = np.stack(df[key]).T
            # if self.data_raw[key].ndim > 1:
            #     self.data_raw[key] = self.data_raw[key].T
    
        self.times = pd.Series(
            map(lambda x: dt.datetime.strptime(x, '%Y-%jT%H:%M:%S.%f'), 
                self.data_raw['DIM0_UTC']))
        
        self.data_raw['DATA'][self.data_raw['DATA'] < -999990] = np.NaN
        self.data_raw['DATA'][
            (self.data_raw['DATA'] > -999990) & 
            (self.data_raw['DATA'] < 0)] = np.NaN
        
        self.data_sigma = self.data_raw['DATA_SIGMA']


    def calc_pitchangles(self):
        """
        Calculate pitch angles from the azimuth and elevation values.
        """
        logging.info('Calculating pitch angles ...')
        elevation = self.data_raw['DIM2_ELEVATION'] * np.pi / 180
        azimuth   = self.data_raw['DIM2_AZIMUTH_DESPUN'] * np.pi / 180

        # Get magnetometer data for pitch angle calculations
        magdata = mag.get_mag_data(
            self.times.iloc[0], self.times.iloc[-1], 
            datafolder=self.magfolder, coordsystem='Juno_JADE_despun')

        # Interpolate MAG data onto times for JADE
        df = pd.DataFrame(index=self.times, columns=['Bx', 'By', 'Bz'])
        print("interpolating mag data ...")
        for t, _ in df.iterrows():
            df.loc[t, :] = magdata.df.loc[
                t - dt.timedelta(seconds=0.6):
                t + dt.timedelta(seconds=0.6), 
                ['Bx', 'By', 'Bz']].mean(axis=0).values

        Bx = df['Bx'].values.astype('float')
        By = df['By'].values.astype('float')
        Bz = df['Bz'].values.astype('float')

        Bmag = np.sqrt(Bx**2 + By**2 + Bz**2)
        Bx /= Bmag; By /= Bmag; Bz /= Bmag

        self.data_raw['PITCH_ANGLE'] = np.arccos(-(
            np.cos(elevation) * np.cos(azimuth) * Bx + \
            np.cos(elevation) * np.sin(azimuth) * By + \
            np.sin(elevation) * Bz)) * 180 / np.pi

        self.pitchanglebins = np.linspace(0, 180, self.nPAbins + 1)


    def detect_sensor(self):
        """
        Returns the sensor in use.
        """
        packetid = self.data_raw['PACKETID']
        if packetid == 105:
            return 'NONE'
        elif packetid == 104:
            return '060'
        elif packetid == 106:
            return '180'
        elif packetid == 107:
            return '300'
        elif packetid == 255:
            return 'UNKNOWN'
        

    
    def calc_energy_spectra(self):
        """
        Calculate the energy spectra.
        """
        print('Calculating energy spectra ...')
        energies = self.data_raw['DIM1_E']
        flux = self.data_raw['DATA']
        pitchangles = self.data_raw['PITCH_ANGLE']


        # Create energybins based on log-space intervals 
        # (and)
        # Filter out those energies requested for by user
        energylin = energies[0, :, 0]

        if self.efilter:
            indices = np.asarray(
                (energylin >= self.efilter[0]) & 
                (energylin <= self.efilter[1])).nonzero()[0]

            energies    = energies[:, indices, :]
            pitchangles = pitchangles[:, indices, :]
            flux        = flux[:, indices, :]
            self.data_sigma = self.data_sigma[:, indices, :]
            energylin   = energylin[indices]

        self.energybins = np.zeros(len(energylin) + 1)
        self.energybins[1:-1] = np.sqrt(energylin[:-1] * energylin[1:])
        self.energybins[0] = self.energybins[1]**2 / self.energybins[2]
        self.energybins[-1] = self.energybins[-2]**2 / self.energybins[-3]
        self.nEbins = len(self.energybins) - 1
        
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
        pitchanglemid = 0.5 * (
            self.pitchanglebins[1:] + self.pitchanglebins[:-1])
        
        solidAngle = np.sin(pitchanglemid * np.pi / 180) *\
                     np.sin(dPA * np.pi / 180 / 2) * 4 * np.pi

        # We need to mask the nan values before averaging
        mask = np.isnan(flux_rebin)
        masked_arr = ma.masked_where(mask, flux_rebin)
        self.data_pro = masked_arr.copy()

        self.energytime = pd.DataFrame(
            index=self.times,
            columns=['Ebin_{:d}'.format(x) for x in range(self.nEbins)],
            data=np.average(self.data_pro, axis=0).T)
        
        self.energytime_sigma = pd.DataFrame(
            index=self.times,
            columns=['Ebin_{:d}'.format(x) for x in range(self.nEbins)],
            data=np.sqrt(np.nansum(self.data_sigma**2, axis=0)).T / flux.shape[0])
            # Combine uncertainty over all look-directions

        self.energytime = self.energytime\
            .resample(self.dtbin, offset=self.dtbin/2)\
            .agg(np.nanmean)
        

    def calc_pa_spectra(self):

        print('Calculating pitch angle spectra ...')
        
        energybinsize = self.energybins[1:] - self.energybins[:-1]

        self.pitchangletime = pd.DataFrame(
            index=self.times,
            columns=['PAbin_{:d}'.format(x) for x in range(self.nPAbins)],
            data=np.average(self.data_pro, axis=1, weights=energybinsize).T)
        
        # When averagin over energy channels, weigh the data by the energy bin
        # width.
        # We need to swap axes to enable numpy broadcasting
        # pasigma = np.swapaxes(self.data_sigma, 1, 2) * energybinsize
        # pasigma = np.sqrt(np.nansum(pasigma**2, axis=2)) / np.sum(energybinsize)

        # self.pitchangletime_sigma = pd.DataFrame(
        #     index=self.times,
        #     columns=['Ebin_{:d}'.format(x) for x in range(self.nEbins)],
        #     data=pasigma.T
        # )


    def resample(self):
        """
        Resample the time-series DataFrames according to the dtbin specified 
        by the user. 
        """
        

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
        
        # x-data , shifted left by half of dtbin 
        timearray = (df.index - self.dtbin / 2).to_pydatetime()
        timearray = np.append(
            timearray, timearray[-1] + (timearray[1] - timearray[0]))

        # y-data is just bins 

        # z-data is the values of the timeseries DataFrame (transposed)
        z = df.values.T

        def get_percentile(df, p):
            return np.exp(np.percentile(np.log(df[df > 0].dropna()), p))
    
        if not vmin:
            vmin = get_percentile(df, 1)
        if not vmax:
            vmax = get_percentile(df, 99)

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
    


    def plot_energy_vs_pitchangle_spectra(self, ax):
        """
        Plot energy vs pitch angle spectra for some time interval.
        """
        pass

    

    def calc_energy_flux(self, pa1=0, pa2=180):
        """
        Calculate the upward/downward energy flux based on the formula given
        in Mauk et al. (2014).
        """
        midpa = 0.5 * (self.pitchanglebins[1:] + self.pitchanglebins[:-1])
        condition = (midpa >= pa1) & (midpa <= pa2)
        
        solidAngle = 4 * np.pi * \
            np.sin(
                (self.pitchanglebins[:-1] + 
                 self.pitchanglebins[1:]) * 0.5 * np.pi / 180) * \
            np.sin(
                (self.pitchanglebins[1] - 
                 self.pitchanglebins[0]) * 0.5 * np.pi / 180)
        
        self.energyflux = pd.Series(index=self.times, data=0.)
        for i in range(self.nEbins):
            Eeff = np.sqrt(self.energybins[i] * self.energybins[i + 1])
            dE = self.energybins[i + 1] - self.energybins[i]

            # We need to convert (cm-2 . s-1 . keV-1 . sr-1) to (mW . m-2)
            # This needs a 1e4 factor  and multiplication with e 
            # = (1e4 to get from cm-2 to m-2) 
            # * (1e-3 to get from keV-1 to eV-1)
            # * (e to convert eV to J)
            # * (1e3 to convert W to mW).

            conversionfactor = \
                solidAngle[condition] * Eeff * dE * e * 1e4
        
            self.energyflux += \
                np.nansum(
                    self.data_pro[condition, i, :].T * conversionfactor, axis=1)
            
            
        self.energyflux[self.energyflux == 0] = np.NaN 
        self.hasEflux = True
