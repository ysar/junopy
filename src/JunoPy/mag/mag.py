import os
import datetime as dt
import numpy as np
import pandas as pd
import copy
import spiceypy
import logging

from pathlib import Path
from JunoPy import spice
from JunoPy import index

class MagData(object):

    columns = [
        "year", "DOY", "hour", "minute", "second", "msec", "DecD",
        "Bx", "By", "Bz", "Brange", "X", "Y", "Z"
        ]
    
    coordsystem = 'JUNO_JSS'

    def __init__(
            self, filename=None, ishdf=True, start=None, stop=None, 
            coordsystem=None):
        """
        Class to parse magnetometer data files. Stores data into a pandas
        dataframe.
        """
        if filename:
            logging.debug(f"Reading mag file - {filename}")
            if ishdf:
                self.filename = filename
                self.df = pd.read_hdf(filename)
            else:
                self.filename = filename
                self.parse_mag_file()
                self.convert_to_hdf()
        
        if start and stop:
            self.df = self.df.loc[start: stop].copy()

        if coordsystem:
            if not self.df.empty:
                self.convert_to_coordsystem(outCoord=coordsystem)
            else:
                self.coordsystem = coordsystem
            

    def __repr__(self):
        """String representation for class. Writes the dataframe head"""
        return self.df.head().__repr__()
    
    
    def _count_header_lines(self):
        """Count header lines in the ASCII magnetometer file."""
        
        i = -1
        with open(self.filename, 'r') as f:

            n_start = 0
            n_end = 0
            for i, line in enumerate(f):
                if "OBJECT" in line:
                    if "END_OBJECT" in line:
                        n_end += 1
                    else:
                        n_start += 1

                if n_end == n_start:
                    break
        
        return i + 1


    def parse_mag_file(self):
        """
        Parse magnetometer file into a pandas dataframe.
        """
        logging.debug(f"Parsing ASCII mag file - {self.filename}")
        nheader = self._count_header_lines()

        self.df = pd.read_csv(
            self.filename,
            names=self.columns,
            delim_whitespace=True,
            skiprows=nheader
            )

        self.df['time'] = self.df.apply(
            lambda x: dt.datetime.strptime(
                '{:04.0f}-{:03.0f}T{:02.0f}:{:02.0f}:{:02.0f}.{:03.0f}'.format(
                    x['year'],
                    x['DOY'],
                    x['hour'],
                    x['minute'],
                    x['second'],
                    x['msec']
                ), '%Y-%jT%H:%M:%S.%f'
            ), axis=1)

        self.df.set_index(['time'], inplace=True)
        self.df.drop(
            inplace=True,
            columns=['year', 'DOY', 'hour', 'minute', 'second', 'msec']
            )

        self.df['Bmag'] = np.sqrt(
            self.df['Bx']**2 + self.df['By']**2 + self.df['Bz']**2)
        

    def read_hdf(self):
        return pd.read_hdf(self.filename)
    
    def convert_to_hdf(self):
        logging.debug(f"Converting file to HDF - {self.filename}")
        assert('.sts' in self.filename)
        self.df.to_hdf(self.filename[:-4] + '.hdf5', 'data')
        os.remove(self.filename)
        self.filename = self.filename[:-4] + '.hdf5'


    def convert_to_spherical(self, usemagcoords=True):
        """
        Calculate magnetic field values in spherical coordinates based
        on the XYZ location.
        """
        logging.debug("Converting to spherical coordinates ...")
        # Calculate R and LT for axes
        if usemagcoords:
            x = self.df['X']
            y = self.df['Y']
            z = self.df['Z']
        else:
            pos, _ = spice.get_juno_position(
                self.df.index, coordsystem=self.coordsystem)
            
            x = pos[:, 0] * 71492.
            y = pos[:, 1] * 71492.
            z = pos[:, 2] * 71492.

        
        self.df['R'] = np.sqrt(x**2 + y**2 + z**2)

        self.df['LT'] = 6 / 90 * np.arctan2(y, x) * 180 / np.pi + 12

        self.df['Latitude'] = np.arctan2(
            z, np.sqrt(x**2 + y**2)) * 180 / np.pi

        phi = np.arctan2(y, x)
        theta = np.pi / 2 - self.df['Latitude'] * np.pi / 180

        self.df['Br'] = np.sin(theta) * np.cos(phi) * self.df['Bx'] + \
            np.sin(theta) * np.sin(phi) * self.df['By'] + \
            np.cos(theta) * self.df['Bz']

        self.df['Btheta'] = np.cos(theta) * np.cos(phi) * self.df['Bx'] + \
            np.cos(theta) * np.sin(phi) * self.df['By'] - \
            np.sin(theta) * self.df['Bz']

        self.df['Bphi'] = - np.sin(phi) * self.df['Bx'] + \
            np.cos(phi) * self.df['By']


    @spice.withspiceloaded()
    def convert_to_coordsystem(self, outCoord='JUNO_JSS', usemagcoords=True):
        logging.debug("Converting from/to coordinate system - {:s}/{:s} ..."
                      .format(self.coordsystem, outCoord))
        
        # spice_load_kernels()

        def tmp_func_convert(row, outCoord=outCoord):

            time = spiceypy.str2et(row.name.strftime('%Y-%m-%dT%H:%M:%S.%f'))
            R = spiceypy.pxform(self.coordsystem, outCoord, time)
            B = np.array([row['Bx'], row['By'], row['Bz']])
            X = np.array([row['X'], row['Y'], row['Z']])

            X_new = np.matmul(R, X)
            B_new = np.matmul(R, B)

            return B_new[0], B_new[1], B_new[2],\
                   X_new[0], X_new[1], X_new[2]
        
        t = self.df.apply(tmp_func_convert, axis=1, result_type='expand')

        try:
            self.df['Bx'], self.df['By'], self.df['Bz'], \
                self.df['X'], self.df['Y'], self.df['Z'] \
                = t[0], t[1], t[2], t[3], t[4], t[5]
            
        except KeyError:
            print(t)
            print(self.filename)
            raise Exception()
        
        self.coordsystem = outCoord
        self.convert_to_spherical(usemagcoords=usemagcoords)

        # spice_unload_kernels()
        
    
    def calculate_bendback(self):
        """
        Calculate bendback angle alpha.
        Note : tan(alpha) / R = (B_phi / (R * B_R)) is often plotted.
        alpha is the angle between the magnetic field in the (R, phi) plane and the radial direction.
        """
        self.df['alpha'] = np.arctan2(
            self.df['Bphi'], self.df['Br']) * 180 / np.pi
    
    
    def plot(self, vars=None, axs=None, color='k'):
        """
        Plot magnetometer data. Returns an axes object.
        """
        
        if 'Btheta' in vars:
            var_string = {
                'Btheta': r'$B_\theta$',
                'Bphi': r'$B_\phi$',
                'Br': r'$B_r$',
                'Bmag': r'|B|',
            }
        
        if 'alpha' in vars:
            var_string.update({'alpha': r'$\alpha'})

        for i, ax in enumerate(axs):
            
            ax.plot(self.df.index, self.df[vars[i]], color=color, lw=0.75)
            ax.grid(which='both', lw='0.5', ls='dotted')

                
            _, _, ymin, ymax = ax.axis()
            if ymin < 0 < ymax:
                ax.axhline(y=0, color='k', ls='dashed', lw=1)
            
            if 'Btheta' in vars:
                ax.set_ylabel(var_string[vars[i]] + ' (nT)')
            elif 'alpha' in vars:
                ax.set_ylabel(var_string[vars[i]] + r'$^\circ$')
            else:
                ax.set_ylabel(vars[i] + ' (nT)')

        return axs
    

def merge_magdata(m1, m2):
    """Merge two magdata files. Usually used for concatenation."""
    m3 = copy.copy(m1)
    m3.df = pd.concat([m1.df, m2.df])    
    m3.df = m3.df.groupby(m3.df.index).agg(lambda x: np.nanmean(x))

    return m3


def get_mag_data(
        start, stop, doDownload=True,
        datafolder='data/mag/', coordsystem=None, sid='SS 1 SECOND'):
    """
    Get magnetometer data between certain start and stop times.
    Downloads files if they do not exist.
    """
    magindex = index.IndexFile(
        dataset='mag', outfolder=datafolder + '/../../index/')

    logging.debug(
        f"Getting magnetometer files. doDownload = {doDownload}")
    
    files = magindex.get_files_between_dates( 
        start, stop, 
        doDownload=doDownload,
        sid=sid, outfolder=datafolder)
    
    if not files:
        logging.error("No relevant magnetometer file found on the PDS.")

    for f in files:
        if not f.with_suffix(".hdf5").exists():
            _ = MagData(filename=str(f), ishdf=False)
        elif not f.exists():
            logging.error("Magnetometer file missing. Use doDownload=True")

    files = [str(x.with_suffix('.hdf5')) for x in files]

    magdatalist = []

    for f in files:
        logging.debug(f"Appending MAG file - {f}")
        magdatalist.append(
            MagData(
                filename=f,
                start=start,
                stop=stop,
                coordsystem=coordsystem
            )
        )

    merged = magdatalist[0]
    for magdata in magdatalist[1:]:
        merged = merge_magdata(merged, magdata)

    return merged