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

import numpy as np
import pandas as pd
import struct
import re
import gzip


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
                nRecords = int(re.search('(\d+)', line)[0])
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
        
        lblfile = filename[:-7] + '.LBL.gz'
    else:
        lblfile = filename[:-4] + '.LBL'

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
