# # junopy - A Python API for the Juno spacecraft's particle and fields datasets.
# # Copyright (C) 2018  Yash Sarkango

# # This program is free software: you can redistribute it and/or modify
# # it under the terms of the GNU General Public License as published by
# # the Free Software Foundation, either version 3 of the License, or
# # (at your option) any later version.

# # This program is distributed in the hope that it will be useful,
# # but WITHOUT ANY WARRANTY; without even the implied warranty of
# # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# # GNU General Public License for more details.

# # You should have received a copy of the GNU General Public License
# # along with this program.  If not, see <https://www.gnu.org/licenses/>.
# """
# This module can be used to read level-3 data files from the JADE instrument.

# Usage:
#     from .read import read_jade_binary_file
#     filename = 'JAD_L30_HRS_ELC_TWO_CNT_2021201_V04.DAT.gz'
#     df = read_jade_binary_file(filename)

# The `read_jade_binary_file` method returns a `pandas` DataFrame that contains
# all information as in the .DAT file. Each row of the DataFrame is a different
# record (i.e. measurement time), and each column is a different value, e.g.
# 'DATA', 'BACKGROUND', 'DIM1_E', etc. Note that the values within a cell of the
# DataFrame can be numpy arrays when applicable.

# The above example reads a file with the .DAT.gz expansion (compressed), but
# reading the .DAT (uncompressed) file is also possible.

# """
# import datetime as dt
# import numpy as np
# import pandas as pd
# import struct
# import re
# import logging

# from isal import igzip


# def open_file(filename, *args):
#     if '.gz' == filename[-3:]:
#         return igzip.open(filename, *args)
#     else:
#         return open(filename, *args)


# def read_jade_lbl_file(lblfile):
#     """
#     Read a JADE .LBL file to infer the column names, types, and byte sizes.

#     Args:
#         lblfile (str): Path to the .LBL file to be read.

#     Returns:
#         dfLabel (pandas.DataFrame): DataFrame with the contents of the .LBL
#             file containing the name, type, and dimensions of each column.
#     """
#     logging.debug(f"Reading JADE label file - {lblfile}")
#     dfLabel = pd.DataFrame(columns=['name', 'dataID-np', 'nDim', 'shape'])

#     with open_file(lblfile, 'rt') as f:

#         i = 0

#         for line in f:
#             if 'FILE_RECORDS' in line:
#                 nRecords = int(re.search(r'(\d+)', line)[0])
#                 break

#         for line in f:
#             if '/* RJW' in line:

#                 if 'BYTES_PER_RECORD' in line:
#                     continue
#                 elif 'OBJECTS_PER_RECORD' in line:
#                     continue

#                 a = line.strip()[3:-3].split(',')[1:]
#                 dfLabel.at[i, 'name'] = a[0].strip()
#                 dfLabel.at[i, 'dataID-np'] = a[1].strip()
#                 dfLabel.at[i, 'nDim'] = int(a[2].strip())
#                 dfLabel.at[i, 'shape'] = np.array([int(x) for x in a[3:]])
#                 i += 1


#     dfLabel['dataID-np'] = dfLabel['dataID-np'].replace('c', 's')

#     dfLabel['nbytes'] = dfLabel.apply(
#         lambda x: np.prod(x['shape']) * struct.calcsize(x["dataID-np"]),
#         axis=1)

#     return dfLabel, nRecords


# class JadeBinaryData(object):
#     """
#     This class reads in the data from the .gz compressed files and stores the
#     data for different packetIDs into a dictionary.
#     """
#     def __init__(self, filename=None, start=None, stop=None):

#         if not filename:
#             return None

#         self.filename = filename
#         df = JadeBinaryData.read_jade_binary_file(filename)

#         if start and stop:
#             df = df.loc[start: stop].copy()

#         self.times = df.index

#         self.dataForPacket = JadeBinaryData.separate_packetids(df)

#     @staticmethod
#     def read_jade_binary_file(filename):
#         """
#         Reads the JADE .DAT file into a single pandas DataFrame. This method assumes
#         that the corresponding .LBL file is stored in the same directory as the .DAT
#         file.

#         Args:
#             filename (str): Path to the .DAT file being read.

#         Returns:
#             df (pandas.DataFrame): DataFrame containing the information stored in
#                 the JADE .DAT file. All strings have been decoded using UTF-8 and
#                 arrays have been reshaped according to the format specified in the
#                 .LBL files.
#         """
#         if filename[-3:] == '.gz':
#             lblfile = filename[:-7] + '.LBL'
#         else:
#             lblfile = filename[:-4] + '.LBL'

#         dfLabel, nRecords = read_jade_lbl_file(lblfile)

#         df = JadeBinaryData._read_jade_bytes(filename, dfLabel, nRecords)

#         for _, column in dfLabel.iterrows():

#             if column['dataID-np'] == 's':
#                 df[column['name']] = df[column['name']].apply(
#                     lambda x: x[0].decode('UTF-8').strip())

#             elif column['nDim'] > 1:
#                 df[column['name']] = df[column['name']].apply(
#                     lambda x: np.array(x).reshape(column['shape']))

#             elif column['nDim'] == 1 and column['shape'] == 1:
#                 df[column['name']] = df[column['name']].apply(lambda x: x[0])


#         df['time'] = pd.Series(
#                 map(lambda x: dt.datetime.strptime(x, '%Y-%jT%H:%M:%S.%f'),
#                     df['DIM0_UTC']))

#         df.set_index('time', inplace=True)

#         return df

#     @staticmethod
#     def _read_jade_bytes(filename, dfLabel, nRecords):
#         """
#         Reads in various bytes from the .DAT file and inserts them into a DataFrame.
#         Does not decode or unpack the bytes, which is faster when done as an apply
#         to the entire column.

#         Args:
#             filename (str): The JADE .DAT file to be read
#             dfLabel (`pandas.DataFrame`):
#                 DataFrame with the contents of the .LBL file.
#             nRecords (int): Total number of records ("rows") in the .DAT file.

#         Returns:
#             A `pandas.DataFrame` containing the contents of the .DAT file. Each row
#             corresponds to a measurement interval and each column has the data
#             associated with that interval. Contents of each cell of the DataFrame
#             may be arrays. For e.g., each value in the 'BACKGROUND' column is an
#             array with shape `(m,n)` with m and n defined in the JADE SIS.
#             Currently, this is stored as a 1D array with size `(m*n)` (it is
#             reshaped later).
#         """
#         logging.debug(f"Reading JADE binary file - {filename}")
#         data = pd.DataFrame(columns=dfLabel['name'].to_list())

#         # We will create numpy arrays to store data and then copy over the array to
#         # the DataFrame. This is much faster than using .at() or .loc().
#         datadict = {}
#         for i, col in dfLabel.iterrows():
#             datadict[col['name']] = np.empty(shape=(nRecords,), dtype='object')

#         # Pre-making these lists is much faster than using .iterrows().
#         n = [np.prod(x) for x in dfLabel['shape']]
#         fmt = dfLabel['dataID-np'].to_list()
#         nbytes = dfLabel['nbytes'].to_list()
#         names = dfLabel['name'].to_list()


#         # Loop over every record (i) and every column (j)
#         # Unfortunately, so far I have not discovered a way to vectorize this.

#         f = open_file(filename, 'rb')

#         i = 0

#         while True:
#             for j in range(len(dfLabel)):
#                 readdata = f.read(nbytes[j])

#                 if readdata == b'':    # Check for end of file.
#                     for name in names:
#                         data[name] = datadict[name]
#                     f.close()
#                     return data

#                 datadict[names[j]][i] = struct.unpack(f'{n[j]}{fmt[j]}', readdata)

#             i += 1

#     @staticmethod
#     def separate_packetids(df):
#         """
#         Here we separate the data for different packet IDs.
#         |-----------|--------------------------|
#         | PacketIds | Data                     |
#         |-----------|--------------------------|
#         | 142       | Electron high-rate       |
#         | 104       | E-060 low-rate           |
#         | 106       | E-180 low-rate           |
#         | 107       | E-270 low-rate           |
#         | 128-135   | Ion Species data (0..7)  |
#         | 105       | Ion TOF high-rate        |
#         | 121       | Ion TOF low-rate         |
#         |-----------|--------------------------|
#         """
#         dfOut = {}
#         for pid in df['PACKETID'].unique():

#             if pid == 142:   # Handle electron HRS data separately
#                 d = JadeBinaryData.create_dict_from_df(df[df['PACKETID'] == pid])

#                 nSensors = d['DATA'].shape[0] // 16
#                 if nSensors == 2:
#                     sensors = {104: d.copy(), 106: d.copy()}
#                 elif nSensors == 3:
#                     sensors = {104: d.copy(), 106: d.copy(), 107: d.copy()}
#                 else:
#                     raise Exception("Invalid num of sensors for packetid=142")

#                 for i, sensor in enumerate(sensors):
#                     dSensor = sensors[sensor]
#                     for k in dSensor.keys():
#                         if dSensor[k].ndim > 1 and dSensor[k].shape[0] > 3:
#                             dSensor[k] = dSensor[k][16 * i: 16 * (i + 1), :, :]

#                     dfOut[sensor] = JadeBinaryData.create_df_from_dict(dSensor)

#             else:
#                 dfOut[pid] = df[df['PACKETID'] == pid]

#         return dfOut


#     @staticmethod
#     def create_dict_from_df(df):
#         """
#         Stacks the JADE data DataFrame into multi-dimensional numpy arrays to
#         facilitate array slicing and faster operations.
#         """
#         dataDict = {}
#         for key in df.keys():
#             dataDict[key] = np.stack(df[key]).T

#         return dataDict


#     @staticmethod
#     def create_df_from_dict(dic):
#         """
#         Creates a DataFrame from a dict. This method is the reverse of the
#         `create_dict_from_df` method.
#         """
#         times = pd.to_datetime(dic['DIM0_UTC'], format = '%Y-%jT%H:%M:%S.%f')
#         df = pd.DataFrame(index=times, columns=dic.keys(), dtype=object)

#         for key in dic.keys():
#             if dic[key].ndim > 1:
#                 df[key] = list(dic[key].T)
#             else:
#                 df[key] = dic[key]

#         return df
