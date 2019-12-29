import os
import numpy as np
import pandas as pd
import datetime as dt
import h5py

from matplotlib.ticker import AutoMinorLocator


class MagData(object):

    columns = [
        "year", "DOY", "hour", "minute", "second", "msec", "DecD",
        "Bx", "By", "Bz", "Brange", "X", "Y", "Z"
    ]

    def __init__(self, filename=None, ishdf=True):
        """
        Class to parse magnetometer data files. Stores data into a pandas
        dataframe.
        """
        if filename:
            if ishdf:
                self.filename = filename
                self.df = pd.read_hdf(filename)
            else:
                self.filename = filename
                self.parse_mag_file()


    def count_header_lines(self):
        """Count header lines in the magnetometer file."""
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
        nheader = self.count_header_lines()

        self.df = pd.read_csv(self.filename,
                              names=self.columns,
                              delim_whitespace=True,
                              skiprows=nheader)

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


    def convert_to_spherical(self):
        """
        Calculate magnetic field values in spherical coordinates based
        on the XYZ location.
        """
        # Calculate R and LT for axes
        self.df['R'] = np.sqrt(
            self.df['X']**2 + self.df['Y']**2 + self.df['Z']**2)

        self.df['LT'] = 6 / 90 * \
            np.arctan2(self.df['Y'], self.df['X']) * 180 / np.pi + 12

        self.df['Latitude'] = np.arctan2(
            self.df['Z'],
            np.sqrt(self.df['X']**2 + self.df['Y']**2)
            ) * 180 / np.pi

        phi = np.arctan2(self.df['Y'], self.df['X'])
        theta = np.pi / 2 - self.df['Latitude'] * np.pi / 180

        self.df['Br'] = np.sin(theta) * np.cos(phi) * self.df['Bx'] + \
            np.sin(theta) * np.sin(phi) * self.df['By'] + \
            np.cos(theta) * self.df['Bz']

        self.df['Btheta'] = np.cos(theta) * np.cos(phi) * self.df['Bx'] + \
            np.cos(theta) * np.sin(phi) * self.df['By'] - \
            np.sin(theta) * self.df['Bz']

        self.df['Bphi'] = - np.sin(phi) * self.df['Bx'] + \
            np.cos(phi) * self.df['By']


    def plot(self, vars=None, axs=None):
        """
        Plot magnetometer data. Returns an axes object.
        """
        for i, ax in enumerate(axs):

            self.df.plot(
                y=vars[i],
                ax=ax,
                grid=True,
                color='k',
                legend=False,
                lw=0.75
            )

            ax.xaxis.set_minor_locator(AutoMinorLocator())
            ax.grid('on', which='minor', lw='0.5', ls='dotted')

            xmin, xmax, ymin, ymax = ax.axis()
            if ymin < 0 < ymax:
                ax.axhline(y=0, color='k', ls='dashed', lw=1)
            ax.set_ylabel(vars[i] + ' (nT)')
            ax.set_xlim([self.df.index[0], self.df.index[-1]])

        return axs


    def plot_electron_cyc_freq(self, ax=None):
        """
        Calculate and plot the electron cyclotron frequency.
        """
        self.df['omega_ce'] = 1.60217e-19 * self.df['Bmag'] * 1e-9 \
                              / 9.1e-31 / 2 / np.pi
        self.df.plot(
            y='omega_ce',
            ax=ax,
            color='k',
            legend=False,
            lw=0.4
        )
        return ax
    
    
    def save_to_hdf5(self):
        """
        Save the current instance to an hdf5 file.
        """
        hdf5_filename = self.filename[:-4] + '.hdf5'
        self.df.to_hdf(hdf5_filename, 
                       os.path.split(self.filename)[-1][:-4], 
                       mode='w')
        os.remove(self.filename)
    
    
def mva(magdata):
    """
    Perform minimum variance analysis (MVA).
    Input - MagData object containing pandas dataframe in self.df
    Ouput - MagData object with additional variables
            B1, B2, B3, eigenval, x1, x2, x3
    """
    magdata.x1 = np.zeros(3)
    magdata.x2 = np.zeros(3)
    magdata.x3 = np.zeros(3)
    magdata.eigvals = np.zeros(3)
    magdata.df['B1'] = 0
    magdata.df['B2'] = 0
    magdata.df['B3'] = 0

    M = np.zeros((3, 3))
    for m, mcom in enumerate(['Bx', 'By', 'Bz']):
        for n, ncom in enumerate(['Bx', 'By', 'Bz']):
            M[m, n] = (magdata.df[mcom] * magdata.df[ncom]).mean() - \
                      magdata.df[mcom].mean() * magdata.df[ncom].mean()

    eigval, eigvec = np.linalg.eig(M)

    imax = np.argmax(eigval)
    imin = np.argmin(eigval)
    imid = list(set([0, 1, 2]) - set([imin, imax]))[0]

    magdata.eigval[0] = eigval[imin]
    magdata.eigval[1] = eigval[imid]
    magdata.eigval[2] = eigval[imax]
    magdata.x1 = eigvec[:, imin]
    magdata.x2 = eigvec[:, imid]
    magdata.x3 = eigvec[:, imax]

    magdata.df['B1'] = magdata.df.apply(
        lambda x: x['Bx'] * magdata.x1[0]
                  + x['By'] * magdata.x1[1]
                  + x['Bz'] * magdata.x1[2], axis=1)
    magdata.df['B2'] = magdata.df.apply(
        lambda x: x['Bx'] * magdata.x2[0]
                  + x['By'] * magdata.x2[1]
                  + x['Bz'] * magdata.x2[2], axis=1)
    magdata.df['B3'] = magdata.df.apply(
        lambda x: x['Bx'] * magdata.x3[0]
                  + x['By'] * magdata.x3[1]
                  + x['Bz'] * magdata.x3[2], axis=1)

    return magdata
