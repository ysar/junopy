# Copyright (c) 2018 Yash Sarkango
# Licensed under the MIT License

import datetime as dt
import numpy as np
import pandas as pd
import copy
import spiceypy
import logging
from functools import reduce

from pathlib import Path
from typing import Optional, Any
from matplotlib import axes

from junopy import spice
from junopy import index
from junopy import utils


class MagData(object):
    columns = [
        "year",
        "DOY",
        "hour",
        "minute",
        "second",
        "msec",
        "DecD",
        "Bx",
        "By",
        "Bz",
        "Brange",
        "X",
        "Y",
        "Z",
    ]

    coordsystem = "JUNO_JSS"

    def __init__(
        self,
        filename: Optional[Path] = None,
        is_hdf: bool = True,
        start: Optional[dt.datetime] = None,
        stop: Optional[dt.datetime] = None,
        coordsystem: Optional[str] = None,
    ):
        """
        Class to parse magnetometer data files. Stores data into a pandas
        dataframe.

        Args:
            filename (Path | None): Name of the data file to be read.
            is_hdf (bool): True if the data file is an HDF5 file.
            start: (dt.datetime | None): Filter result from this starting time.
            stop: (dt.datetime | None): Filter result until this stopping time.
            coordsystem (str): Coordinate system recognized by the SPICE kernels.

        Returns:
            MagData: An instance of the `MagData` class that contains the parsed
                data stored in the file.  If `filename is None`, then an empty
                `MagData` instance is returned.

        Attributes:
            df (pd.DataFrame): Dataframe containing the magnetometer data.  The index
                is a datetime index and the columns are:
                ['Bx', 'By', 'Bz', 'Brange', 'X', 'Y', 'Z', 'Bmag'] and optionally
                ['R', 'LT', 'Latitude', 'Br', 'Btheta', 'Bphi'] if spherical
                coordinates are calculated.
            coordsystem (str): Coordinate system of the data.  Default is 'JUNO_JSS'.
            filename (Path | None): Path to the data file associated with this instance
            origin (str): Origin of the coordinate system in SPICE terminology.
                Default is 'JUPITER'.
        """
        if filename:
            logging.info(f"Reading mag file - {filename}")
            if is_hdf:
                self.filename = filename
                _df = pd.read_hdf(filename)
                assert isinstance(_df, pd.DataFrame)
                self.df = _df
            else:
                self.filename = filename
                self._parse_mag_file()
                self._convert_to_hdf()

        assert isinstance(self.df, pd.DataFrame)

        # Filter the result based on the starting and stopping times.
        if start and stop:
            self.df = self.df.loc[start:stop].copy()

        # Default origin is JUPITER
        self.origin = "JUPITER"

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
        with open(self.filename, "r") as f:
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
        logging.info(f"Parsing ASCII mag file - {str(self.filename)}")
        nheader = self._count_header_lines()

        self.df: pd.DataFrame = pd.read_csv(
            self.filename, names=self.columns, sep=r"\s+", skiprows=nheader
        )

        self.df["time"] = self.df.apply(
            lambda x: dt.datetime.strptime(
                "{:04.0f}-{:03.0f}T{:02.0f}:{:02.0f}:{:02.0f}.{:03.0f}".format(
                    x["year"], x["DOY"], x["hour"], x["minute"], x["second"], x["msec"]
                ),
                "%Y-%jT%H:%M:%S.%f",
            ),
            axis=1,
        )

        self.df.set_index(["time"], inplace=True)
        self.df.drop(
            inplace=True, columns=["year", "DOY", "hour", "minute", "second", "msec"]
        )

        self.df["Bmag"] = np.sqrt(
            self.df["Bx"] ** 2 + self.df["By"] ** 2 + self.df["Bz"] ** 2
        )

    def _read_hdf(self):
        return pd.read_hdf(self.filename)

    def _convert_to_hdf(self):
        """
        Convert the ASCII magnetometer data file associated with this
        instance to HDF5.  The original ASCII file will be deleted.
        """
        logging.info(f"Converting file to HDF - {self.filename}")

        assert self.filename.suffix == ".sts"

        hdf5_filename = self.filename.with_suffix(".hdf5")
        self.df.to_hdf(hdf5_filename, key="data")
        self.filename.unlink()
        self.filename = hdf5_filename

    def convert_to_spherical(self, use_mag_coords: bool = False):
        """
        Calculate magnetic field values in spherical coordinates based
        on the XYZ location.  The following fields are added to the magnetometer
        dataframe - ['R', 'LT', 'Latitude', 'Br', 'Btheta', 'Bphi']

        Args:
            use_mag_coords (bool): Whether to use the coordinates stored in the
                magnetometer file.
        """
        logging.info("Converting to spherical coordinates ...")
        # Calculate R and LT for axes
        if use_mag_coords:
            pos = np.zeros((len(self.df), 3)) * np.nan
            pos[:, 0] = self.df["X"]
            pos[:, 1] = self.df["Y"]
            pos[:, 2] = self.df["Z"]
        else:
            pos, _ = spice.get_juno_position(
                [x.to_pydatetime() for x in self.df.index], coordsystem=self.coordsystem
            )

        self.df["R"], self.df["LT"], self.df["Latitude"] = (
            utils.convert.convert_pos_spherical(pos.T)
        )

        phi = np.arctan2(pos[:, 1], pos[:, 0])
        theta = np.pi / 2 - self.df["Latitude"] * np.pi / 180

        self.df["Br"] = (
            np.sin(theta) * np.cos(phi) * self.df["Bx"]
            + np.sin(theta) * np.sin(phi) * self.df["By"]
            + np.cos(theta) * self.df["Bz"]
        )

        self.df["Btheta"] = (
            np.cos(theta) * np.cos(phi) * self.df["Bx"]
            + np.cos(theta) * np.sin(phi) * self.df["By"]
            - np.sin(theta) * self.df["Bz"]
        )

        self.df["Bphi"] = -np.sin(phi) * self.df["Bx"] + np.cos(phi) * self.df["By"]

    @spice.withspiceloaded()
    def convert_to_coordsystem(
        self,
        out_coordsystem: str = "JUNO_JSS",
        use_mag_coords: bool = True,
        new_origin: str = "JUPITER",
    ) -> None:
        """
        Converts the magnetometer data (X, Y, Z, Bx, By, Bz) and (R, Theta,
        Phi, Br, Btheta, Bphi) to another coordinate system.  Data within the
        `.df` attribute will be converted to this new system.

        Args:
            out_coordsystem (str): Coordinate system to convert to.  Must be
                recognizable by SPICE.
            use_mag_coords (bool): Whether to use the coordinates provided in
                the magnetometer data file instead of SPICE.
            new_origin (str): The origin of the new coordinate system.  Must be
                recognizable by SPICE.
        """
        logging.info(
            "Coordinate system conversion from {:s} to {:s}".format(
                self.coordsystem, out_coordsystem
            )
        )

        def _convert(row):
            et = spiceypy.str2et(row.name.strftime("%Y-%m-%dT%H:%M:%S.%f"))
            assert isinstance(et, float)

            # Origin offset between the two coordinate systems
            origin_offset, _ = spiceypy.spkpos(
                new_origin, et, self.coordsystem, "NONE", self.origin
            )

            # Apply the translation operation to the new origin
            pos = np.array([row["X"], row["Y"], row["Z"]])
            pos = pos - origin_offset

            # Rotate the translated system.
            rotation_matrix = spiceypy.pxform(self.coordsystem, out_coordsystem, et)
            b_prev = np.array([row["Bx"], row["By"], row["Bz"]])

            pos_rotated = np.matmul(rotation_matrix, pos)
            b_rotated = np.matmul(rotation_matrix, b_prev)

            return tuple(b_rotated) + tuple(pos_rotated)

        _t = self.df.apply(_convert, axis=1, result_type="expand")

        try:
            self.df.loc[:, ["Bx", "By", "Bz", "X", "Y", "Z"]] = _t.values
        except Exception:
            print(_t)
            print(self.filename)
            raise

        self.coordsystem = out_coordsystem
        self.origin = new_origin
        self.convert_to_spherical(use_mag_coords=use_mag_coords)

    def calculate_bendback(self):
        """
        Calculate bendback angle alpha.
        Note : tan(alpha) / R = (B_phi / (R * B_R)) is often plotted.
        alpha is the angle between the magnetic field in the (R, phi) plane and
        the radial direction.
        """
        self.df["alpha"] = np.arctan2(self.df["Bphi"], self.df["Br"]) * 180 / np.pi

    def plot(
        self, vars: list[str], axs: list[axes.Axes], color: Any = "k", **kwargs
    ) -> list[axes.Axes]:
        """
        Plot magnetometer data. Modifies and returns an axes object.

        Args:
            vars (list[str]): Variables from the allowable keys of
                the `MagData.df` dataframe to be plotted.
            axs (list[mpl.axes.Axes]: `matplotlib` axes on which the variables
                are to be plotted.
            color (Any): Color used for plotting the data.
            **kwargs: Additional keyword arguments passed to `ax.plot`.

        Returns:
            mpl.axes.Axes: The axes objects are modified and returned.
        """

        if "Btheta" in vars:
            var_string = {
                "Btheta": r"$B_\theta$",
                "Bphi": r"$B_\phi$",
                "Br": r"$B_r$",
                "Bmag": r"|B|",
            }
        else:
            var_string = dict()

        if "alpha" in vars:
            var_string.update({"alpha": r"$\alpha"})

        for i, ax in enumerate(axs):
            ax.plot(self.df.index, self.df.loc[:, vars[i]], color=color, **kwargs)
            ax.grid(which="both", lw="0.5", ls="dotted")

            _, _, ymin, ymax = ax.axis()
            if ymin < 0 < ymax:
                ax.axhline(y=0, color="k", ls="dashed", lw=1)

            if "Btheta" in vars:
                ax.set_ylabel(var_string[vars[i]] + " (nT)")
            elif "alpha" in vars:
                ax.set_ylabel(var_string[vars[i]] + r"$^\circ$")
            else:
                ax.set_ylabel(vars[i] + " (nT)")

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
    _df = m3.df.groupby(m3.df.index).agg(np.nanmean)
    assert isinstance(_df, pd.DataFrame)
    m3.df = _df
    return m3


def find_mag_data(
    start: dt.datetime,
    stop: dt.datetime,
    do_download: bool = True,
    do_parallel: bool = False,
    sid: str = "SS 1 SECOND",
) -> list[Path]:
    """
    Finds and downloads magnetometer data between certain start and stop times

    This function will convert the files to
    HDF5 format if they are not already in that format.  This method returns a list
    of the relevant HDF5 files that were downloaded or converted or were already
    present.

    Args:
        start (dt.datetime): Get data starting from this date-time.
        stop (dt.datetime): Get data until this date-time.
        do_download (bool): Whether to download missing files.
        do_parallel (bool): Whether to download files in parallel.
        sid (str): SID corresponding to the data file options are "SS 1 SECOND",
            "SS 1 MINUTE", "SUNSTATE", "PC 1 SECOND", "PC 1 MINUTE", etc.

    Returns:
        list[str]: List of relevant HDF5 files that were downloaded or converted
            or were already present in the junopy root directory.
    """
    magindex = index.IndexFile(dataset="mag")

    logging.info(f"Getting magnetometer files. do_download = {do_download}")

    files = magindex.get_files_between_dates(
        start, stop, do_download=do_download, sid=sid
    )

    if not files:
        logging.error("No relevant magnetometer file found on the PDS.")

    for f in files:
        if f.exists():
            logging.info("Converting to HDF5 file ...")
            _ = MagData(filename=f, is_hdf=False)
        elif not f.with_suffix(".hdf5").exists():
            logging.error("Magnetometer file missing. Try do_download=True")

    files = [x.with_suffix(".hdf5") for x in files]
    return files


def load_mag_data(
    start: dt.datetime,
    stop: dt.datetime,
    do_download: bool = True,
    do_parallel: bool = False,
    coordsystem: Optional[str] = None,
    sid: str = "SS 1 SECOND",
) -> MagData:
    """
    Loads the magnetometer data into memory as a single `MagData` instance.

    Warning:
        This method will merge all the files that are within the start
        and stop times, so it is not recommended to use this method if you need to
        work with a large number of files or if the time range is too large, as
        you may have insufficient memory to hold the data.

    Args:
        start (dt.datetime): Load data starting from this date-time.
        stop (dt.datetime): Load data until this date-time.
        do_download (bool): Whether to download missing files.
        do_parallel (bool): Whether to download files in parallel.
        coord_system (str | None): Coordinate system for resulting data.
        sid (str): SID corresponding to the data file options are "SS 1 SECOND",
            "SS 1 MINUTE", "SUNSTATE", "PC 1 SECOND", "PC 1 MINUTE", etc.

    Returns:
        MagData: An instance of the `MagData` class that contains the merged
            data between the start and stop times.
    """
    files = find_mag_data(
        start,
        stop,
        do_download=do_download,
        do_parallel=do_parallel,
        sid=sid,
    )

    magdatalist = []

    for f in files:
        logging.info(f"Appending MAG file - {f}")
        magdatalist.append(
            MagData(filename=f, start=start, stop=stop, coordsystem=coordsystem)
        )

    merged = reduce(merge_magdata, magdatalist)

    return merged
