# junopy - A Python API for the Juno spacecraft's particle and fields datasets.
# Copyright (C) 2018  Yash Sarkango

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

import os
import datetime as dt
import numpy as np
import pandas as pd
import copy
import spiceypy
import logging
import matplotlib as mpl
from typing import Optional, Any

from junopy import spice
from junopy import index
from junopy import utils

class MagData(object):

    columns = [
        "year", "DOY", "hour", "minute", "second", "msec", "DecD",
        "Bx", "By", "Bz", "Brange", "X", "Y", "Z"
        ]

    coordsystem = 'JUNO_JSS'

    def __init__(
            self,
            filename: Optional[str] = None,
            is_hdf: bool = True,
            start: Optional[dt.datetime] = None,
            stop: Optional[dt.datetime] = None,
            coordsystem: Optional[str] = None):
        """
        Class to parse magnetometer data files. Stores data into a pandas
        dataframe.

        Args:
            filename (str | None): Name of the data file to be read.
            is_hdf (bool): True if the data file is an HDF5 file.
            start: (dt.datetime | None): Filter result from this starting time.
            stop: (dt.datetime | None): Filter result until this stopping time.
            coordsystem (str): Coordinate system recognized by the SPICE kernels.

        Returns:
            MagData: An instance of the `MagData` class that contains the parsed
                data stored in the file.  If `filename is None`, then an empty
                `MagData` instance is returned.
        """
        if filename:
            logging.info(f"Reading mag file - {filename}")
            if is_hdf:
                self.filename = filename
                self.df = pd.read_hdf(filename)
            else:
                self.filename = filename
                self._parse_mag_file()
                self._convert_to_hdf()

        assert isinstance(self.df, pd.DataFrame)

        # Filter the result based on the starting and stopping times.
        if start and stop:
            self.df = self.df.loc[start: stop].copy()

        # Convert to a coordinate system recognized by SPICE.
        if coordsystem:
            if not self.df.empty:
                self.convert_to_coordsystem(out_coordsystem=coordsystem)
            else:
                self.coordsystem = coordsystem


    def __repr__(self):
        """String representation for class. Writes out the dataframe head"""
        return self.df.head().__repr__()


    def _count_header_lines(self) -> int:
        """
        Count header lines in the ASCII magnetometer file.

        Returns:
            int: Number of header lines
        """

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


    def _parse_mag_file(self):
        """
        Parse magnetometer file into a pandas dataframe.
        """
        logging.info(f"Parsing ASCII mag file - {self.filename}")
        nheader = self._count_header_lines()

        self.df: pd.DataFrame = pd.read_csv(
            self.filename,
            names=self.columns,
            sep=r'\s+',
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


    def _read_hdf(self):
        return pd.read_hdf(self.filename)

    def _convert_to_hdf(self):
        """
        Convert the ASCII magnetometer data file associated with this
        instance to HDF5.  The original ASCII file will be deleted.
        """
        logging.info(f"Converting file to HDF - {self.filename}")
        assert('.sts' in self.filename)
        self.df.to_hdf(self.filename[:-4] + '.hdf5', key='data')
        os.remove(self.filename)
        self.filename = self.filename[:-4] + '.hdf5'


    def convert_to_spherical(self, use_mag_coords: bool = False):
        """
        Calculate magnetic field values in spherical coordinates based
        on the XYZ location.  The following fields are added to the magnetometer
        dataframe - ['R', 'LT', 'Latitude', 'Br', 'Btheta', 'Bphi']
        """
        logging.info("Converting to spherical coordinates ...")
        # Calculate R and LT for axes
        if use_mag_coords:
            pos = np.zeros((len(self.df), 3)) * np.nan
            pos[:, 0] = self.df['X']
            pos[:, 1] = self.df['Y']
            pos[:, 2] = self.df['Z']
        else:
            pos, _ = spice.get_juno_position(
                self.df.index, coordsystem=self.coordsystem)


        self.df['R'], self.df['LT'], self.df['Latitude'] = \
            utils.convert.convert_pos_spherical(pos.T)

        phi = np.arctan2(pos[:, 1], pos[:, 0])
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
    def convert_to_coordsystem(
        self,
        out_coordsystem: str = 'JUNO_JSS',
        use_mag_coords: bool = True
    ) -> None :
        logging.info("Converting from/to coordinate system - {:s}/{:s} ..."
                      .format(self.coordsystem, out_coordsystem))

        # spice_load_kernels()

        def tmp_func_convert(row, outCoord=out_coordsystem):

            et = spiceypy.str2et(row.name.strftime('%Y-%m-%dT%H:%M:%S.%f'))

            # Origin offset between the two coordinate systems
            origin_offset =



            rotation_matrix = spiceypy.pxform(self.coordsystem, out_coordsystem, et)
            b_prev = np.array([row['Bx'], row['By'], row['Bz']])
            pos = np.array([row['X'], row['Y'], row['Z']])

            pos_rotated = np.matmul(rotation_matrix, pos)
            b_rotated = np.matmul(rotation_matrix, b_prev)

            return b_rotated, pos_rotated

        t = self.df.apply(tmp_func_convert, axis=1, result_type='expand')

        # Apply translation.



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


    def plot(
        self,
        vars: list[str],
        axs: list[mpl.axes.Axes],
        color: Any = 'k'):
        """
        Plot magnetometer data. Modifies and returns an axes object.

        Args:
            vars (list[str]): Variables from the allowable keys of
                the `MagData.df` dataframe to be plotted.
            axs (list[mpl.axes.Axes]: `matplotlib` axes on which the variables
                are to be plotted.
            color (Any): Color used for plotting the data.

        Returns:
            mpl.axes.Axes: The axes objects are modified and returned.
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


def merge_magdata(m1: MagData, m2: MagData) -> MagData:
    """
    Merge two magdata files. Usually used for concatenation.

    Args:
        m1 (MagData): 1st instance to be merged
        m2 (MagData): 2nd instance to be merged

    Returns:
        Magdata: Merged instance with a nanmean operation to aggregate the same
            indices.
    """
    m3 = copy.copy(m1)
    m3.df = pd.concat([m1.df, m2.df])
    m3.df = m3.df.groupby(m3.df.index).agg(np.nanmean)

    return m3


def get_mag_data(
        start: dt.datetime,
        stop: dt.datetime,
        do_download: bool = True,
        coord_system: Optional[str] = None,
        sid: str = 'SS 1 SECOND'
    ) -> MagData:
    """
    Get magnetometer data between certain start and stop times.
    Downloads files if they do not exist.

    Args:
        start (dt.datetime): Get data starting from this date-time.
        stop (dt.datetime): Get data until this date-time.
        do_download (bool): Whether to download missing files.
        coord_system (str | None): Coordinate system for resulting data.
        sid (str): SID corresponding to the data file, see Usage>Magnetometer
            for different options.

    """

    magindex = index.IndexFile(dataset='mag')

    logging.info(
        f"Getting magnetometer files. do_download = {do_download}")

    files = magindex.get_files_between_dates(
        start, stop,
        do_download=do_download,
        sid=sid)

    if not files:
        logging.error("No relevant magnetometer file found on the PDS.")

    for f in files:
        if f.exists():
            logging.info("Converting to HDF5 file ...")
            _ = MagData(filename=str(f), is_hdf=False)
        elif not f.with_suffix(".hdf5").exists():
            logging.error("Magnetometer file missing. Use do_download=True")

    files = [str(x.with_suffix('.hdf5')) for x in files]

    magdatalist = []

    for f in files:
        logging.info(f"Appending MAG file - {f}")
        magdatalist.append(
            MagData(
                filename=f,
                start=start,
                stop=stop,
                coordsystem=coord_system
            )
        )

    merged = magdatalist[0]
    for magdata in magdatalist[1:]:
        merged = merge_magdata(merged, magdata)

    return merged
