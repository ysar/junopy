
import numpy as np
import pandas as pd
import datetime as dt

from matplotlib.ticker import AutoMinorLocator


class MagData(object):

    columns = [
        "year", "DOY", "hour", "minute", "second", "msec", "DecD",
        "Bx", "By", "Bz", "Brange", "X", "Y", "Z"
    ]


    def __init__(self, filename=None):
        """
        Class to parse magnetometer data files. Stores data into a pandas
        dataframe.
        """
        if filename:
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
                lw=0.4
            )

            ax.xaxis.set_minor_locator(AutoMinorLocator())
            ax.grid('on', which='minor', lw='0.5', ls='dotted')

            xmin, xmax, ymin, ymax = ax.axis()
            if ymin < 0 < ymax:
                ax.axhline(y=0, color='k', ls='dashed', lw=1)
            ax.set_ylabel(vars[i] + ' (nT)')

        return axs


    def plot_electron_cyc_freq(self, ax=None):
        """
        Calculate and plot the electron cyclotron frequency.
        """
        self.df['omega_ce'] = 1.60217e-19 * self.df['Bmag'] * \
                              1e-9 / 9.1e-31 / 2 / np.pi
        self.df.plot(
            y='omega_ce',
            ax=ax,
            color='k',
            legend=False,
            lw=0.4
        )

        return ax
