# Copyright (c) 2018 Yash Sarkango
# Licensed under the MIT License
"""
This module can be used to read level-3 data files from the JADE instrument.

Usage:
    ```python
    from junopy.jade.read import read_jade_binary_file
    filename = 'JAD_L30_HRS_ELC_TWO_CNT_2021201_V04.DAT.gz'
    df = read_jade_binary_file(filename)
        # This will merely load the data into a pandas DataFrame.  It will NOT
        # create a `JadeBinaryData` object that has more functionality.
        # In practice, consider using higher-level functions to get data for a
        # specific time range, e.g. junopy.jade.get.get_jade_data_files
    ```

The `read_jade_binary_file` method returns a `pandas` DataFrame that contains
all information as in the .DAT file. Each row of the DataFrame is a different
record (i.e. measurement time), and each column is a different value, e.g.
'DATA', 'BACKGROUND', 'DIM1_E', etc. Note that the values within a cell of the
DataFrame can be numpy arrays when applicable.

The above example reads a file with the .DAT.gz expansion (compressed), but
reading the .DAT (uncompressed) file is also possible.

"""

import datetime as dt
import logging
import re
import struct
from pathlib import Path

import numpy as np
import pandas as pd
from isal import igzip


def open_file(filename: Path, *args):
    """
    Opens a file, using gzip if the filename ends with .gz

    Args:
        filename (str): The name of the file to open.
        *args: Additional arguments to pass to the open function.
    """
    if filename.suffix == ".gz":
        return igzip.open(filename, *args)
    else:
        return open(filename, *args)


def read_jade_lbl_file(lblfile: Path) -> tuple[pd.DataFrame, int]:
    """
    Read a JADE .LBL file to infer the column names, types, and byte sizes.

    Args:
        lblfile (Path): Path to the .LBL file to be read.

    Returns:
        df_label (pandas.DataFrame): DataFrame with the contents of the .LBL
            file containing the name, type, and dimensions of each column.
        num_records (int): Total number of records in the .DAT file.
    """
    logging.debug(f"Reading JADE label file - {str(lblfile)}")
    df_label = pd.DataFrame(columns=["name", "dataID-np", "nDim", "shape"])

    with open_file(lblfile, "rt") as f:
        i = 0

        for line in f:
            if "FILE_RECORDS" in line:
                num_records = int(re.search(r"(\d+)", line)[0])
                break

        for line in f:
            if "/* RJW" in line:
                if "BYTES_PER_RECORD" in line:
                    continue
                elif "OBJECTS_PER_RECORD" in line:
                    continue

                a = line.strip()[3:-3].split(",")[1:]
                df_label.at[i, "name"] = a[0].strip()
                df_label.at[i, "dataID-np"] = a[1].strip()
                df_label.at[i, "nDim"] = int(a[2].strip())
                df_label.at[i, "shape"] = np.array([int(x) for x in a[3:]])
                i += 1

    df_label["dataID-np"] = df_label["dataID-np"].replace("c", "s")

    df_label["nbytes"] = df_label.apply(
        lambda x: np.prod(x["shape"]) * struct.calcsize(x["dataID-np"]), axis=1
    )

    return df_label, num_records


class JadeBinaryData(object):
    """
    This class reads in the data from the .gz compressed files and stores the
    data for different packetIDs into a dictionary.

    Args:
        filename (Path, optional): Path to the .DAT or .DAT.gz file to be read.
            If not provided, an empty object is created.
        start (dt.datetime, optional): Start time for slicing the data. If not
            provided, all data is read.
        stop (dt.datetime, optional): Stop time for slicing the data. If not
            provided, all data is read.

    Attributes:
        times (pd.DatetimeIndex): Index of the DataFrame containing the times
            of the measurements.
        data_by_packetid (dict[int, pd.DataFrame]): Dictionary with keys as
            packet IDs and values as DataFrames containing the data for that
            packet ID.
        packetids (list[int]): List of packet IDs present in the data.
        filename (Path): Path to the .DAT or .DAT.gz file being read.
    """

    def __init__(
        self,
        filename: Path | None = None,
        start: dt.datetime | None = None,
        stop: dt.datetime | None = None,
    ):
        if not filename:
            return None

        self.filename = filename
        df = JadeBinaryData.read_jade_binary_file(filename, start, stop)

        self.times = df.index

        self.data_by_packetid = JadeBinaryData.separate_packetids(df)
        self.packetids = [int(x) for x in self.data_by_packetid.keys()]

    @staticmethod
    def read_jade_binary_file(
        file_path: Path, start: dt.datetime | None, stop: dt.datetime | None
    ) -> pd.DataFrame:
        """
        Reads the JADE .DAT file into a single pandas DataFrame. This method assumes
        that the corresponding .LBL file is stored in the same directory as the .DAT
        file.

        Args:
            filename (Path): Path to the .DAT file being read.
            start (dt.datetime | None): Start time to read from.
            stop (dt.datetime | None): Stop time to read until.

        Returns:
            df (pandas.DataFrame): DataFrame containing the information stored in
                the JADE .DAT file. All strings have been decoded using UTF-8 and
                arrays have been reshaped according to the format specified in the
                .LBL files.
        """
        if ".gz" in file_path.suffixes:
            lbl_file = file_path.with_suffix("").with_suffix(".LBL")
        else:
            lbl_file = file_path.with_suffix(".LBL")

        df_label, num_records = read_jade_lbl_file(lbl_file)

        df = JadeBinaryData._read_jade_bytes(
            file_path, df_label, num_records, start, stop
        )

        for i, column in df_label.iloc[1:].iterrows():
            # Skips the first column (DIM0_UTC) since we already processed it while reading in the
            # file.
            if column["dataID-np"] == "s":
                df[column["name"]] = df[column["name"]].apply(
                    lambda x: x[0].decode("UTF-8").strip()
                )

            elif column["nDim"] > 1:
                df[column["name"]] = df[column["name"]].apply(
                    lambda x: np.array(x).reshape(column["shape"])
                )

            elif column["nDim"] == 1 and column["shape"] == 1:
                df[column["name"]] = df[column["name"]].apply(lambda x: x[0])

        # df["time"] = pd.Series(
        #     map(lambda x: dt.datetime.strptime(x, "%Y-%jT%H:%M:%S.%f"), df["DIM0_UTC"])
        # )
        df["time"] = df["DIM0_UTC"]

        df.set_index("time", inplace=True)

        return df

    @staticmethod
    def _read_jade_bytes(
        filename: Path,
        df_label: pd.DataFrame,
        num_records: int,
        start: dt.datetime | None,
        stop: dt.datetime | None,
    ) -> pd.DataFrame:
        """
        Reads in various bytes from the .DAT file and inserts them into a DataFrame.
        Does not decode or unpack the bytes, which is faster when done as an apply
        to the entire column.

        Args:
            filename (str): The JADE .DAT file to be read
            df_label (`pandas.DataFrame`):
                DataFrame with the contents of the .LBL file.
            num_records (int): Total number of records ("rows") in the .DAT file.

        Returns:
            A `pandas.DataFrame` containing the contents of the .DAT file. Each row
            corresponds to a measurement interval and each column has the data
            associated with that interval. Contents of each cell of the DataFrame
            may be arrays. For e.g., each value in the 'BACKGROUND' column is an
            array with shape `(m,n)` with m and n defined in the JADE SIS.
            Currently, this is stored as a 1D array with size `(m*n)` (it is
            reshaped later).
        """
        logging.debug(f"Reading JADE binary file - {filename}")
        # data = pd.DataFrame(columns=df_label["name"].to_list())

        # We will create numpy arrays to store data and then copy over the array to
        # the DataFrame. This is much faster than using .at() or .loc().
        datadict = {}
        for i, col in df_label.iterrows():
            datadict[col["name"]] = np.empty(shape=(num_records,), dtype="object")

        # Pre-making these lists is much faster than using .iterrows().
        n = [np.prod(x) for x in df_label["shape"]]
        fmt = df_label["dataID-np"].to_list()
        nbytes = df_label["nbytes"].to_list()
        names = df_label["name"].to_list()

        # Loop over every record (i) and every column (j)
        # Unfortunately, so far I have not discovered a way to vectorize this.

        f = open_file(filename, "rb")

        def __read_record_item(f, num_bytes, num_shape, fmt):
            # Read one item in the record.
            readdata = f.read(num_bytes)
            return struct.unpack(f"{num_shape}{fmt}", readdata)

        def __read_utc_time(f, num_bytes, num_shape, fmt):
            # Read the UTC time that marks the beginning of each record.
            val = __read_record_item(f, num_bytes, num_shape, fmt)
            return dt.datetime.strptime(
                val[0].decode("UTF-8").strip(), "%Y-%jT%H:%M:%S.%f"
            )

        num_bytes_to_skip = df_label["nbytes"].sum() - nbytes[0]
        valid_indices = []

        for i in range(num_records):
            t = __read_utc_time(f, nbytes[0], n[0], fmt[0])

            if start and t < start:
                # Need to skip the rest of the records but advance the file pointer.
                # Here "1" means skip from the current file position.
                f.seek(num_bytes_to_skip, 1)
                continue

            if stop and t > stop:
                # If we have gone beyond what is required, we can close the file.
                f.close()
                break

            valid_indices.append(i)
            datadict[names[0]][i] = t

            # Read the rest of the data for this record.
            for j in range(1, len(df_label)):
                datadict[names[j]][i] = __read_record_item(f, nbytes[j], n[j], fmt[j])

        # Remove invalid indices that correspond to times outside of the requested range.
        for name in names:
            datadict[name] = datadict[name][valid_indices]

        # One all relevant data is read, convert the dictionary into the dataframe.
        return pd.DataFrame(datadict, copy=False)

    @staticmethod
    def separate_packetids(df: pd.DataFrame) -> dict[int, pd.DataFrame]:
        """
        Here we separate the data for different packet IDs.

        ```md
        |-----------|--------------------------|
        | PacketIds | Data                     |
        |-----------|--------------------------|
        | 142       | Electron high-rate       |
        | 104       | E-060 low-rate           |
        | 106       | E-180 low-rate           |
        | 107       | E-270 low-rate           |
        | 128-135   | Ion Species data (0..7)  |
        | 105       | Ion TOF high-rate        |
        | 121       | Ion TOF low-rate         |
        |-----------|--------------------------|
        ```

        Args:
            df (pandas.DataFrame): DataFrame containing the entire JADE level-3 data.

        Returns:
            df_out (dict[int, pd.DataFrame]): A dictionary with keys as the packet IDs and values as
                DataFrames containing the data for that packet ID.
        """
        df_out = {}
        for pid in df["PACKETID"].unique():
            if pid == 142:  # Handle electron HRS data separately
                d = JadeBinaryData.create_dict_from_df(df[df["PACKETID"] == pid])

                num_sensors = d["DATA"].shape[0] // 16
                if num_sensors == 2:
                    sensors = {104: d.copy(), 106: d.copy()}
                elif num_sensors == 3:
                    sensors = {104: d.copy(), 106: d.copy(), 107: d.copy()}
                else:
                    raise Exception("Invalid num of sensors for packetid=142")

                for i, sensor in enumerate(sensors):
                    dict_sensor = sensors[sensor]
                    for k in dict_sensor.keys():
                        if dict_sensor[k].ndim > 1 and dict_sensor[k].shape[0] > 3:
                            dict_sensor[k] = dict_sensor[k][16 * i : 16 * (i + 1), :, :]

                    df_out[sensor] = JadeBinaryData.create_df_from_dict(dict_sensor)

            else:
                df_out[pid] = df[df["PACKETID"] == pid]

        return df_out

    @staticmethod
    def create_dict_from_df(df: pd.DataFrame) -> dict:
        """
        Stacks the JADE data DataFrame into multi-dimensional numpy arrays to
        facilitate array slicing and faster operations.

        Args:
            df (pandas.DataFrame): DataFrame containing the JADE data.

        Returns:
            data_dict (dict): A dictionary with keys as the column names and values
                as numpy arrays. Each array is stacked along the first axis.
        """
        data_dict = {}
        for key in df.keys():
            data_dict[key] = np.stack(df[key]).T

        return data_dict

    @staticmethod
    def create_df_from_dict(d: dict) -> pd.DataFrame:
        """
        Creates a DataFrame from a dict. This method is the reverse of the
        `create_dict_from_df` method.

        Args:
            d (dict): A dictionary with keys as the column names and values
                as numpy arrays. Each array is stacked along the first axis.
        Returns:
            df (pandas.DataFrame): DataFrame containing the data from the dict.
        """
        times = pd.to_datetime(d["DIM0_UTC"], format="%Y-%jT%H:%M:%S.%f")
        df = pd.DataFrame(index=times, columns=d.keys(), dtype=object)

        for key in d.keys():
            if d[key].ndim > 1:
                df[key] = list(d[key].T)
            else:
                df[key] = d[key]

        return df
