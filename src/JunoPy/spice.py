"""
Module to work with Juno spice kernels.

For now we will only use spice to obtain spacecraft position.
"""
import os
import numpy as np
import datetime as dt
import pandas as pd
import spiceypy as spice


from JunoPy.download import download_file

SPICE_DIR = 'spice'
KERNEL_DIR = os.path.join(SPICE_DIR, 'kernels')
JUNO_SPICE_URL = "https://naif.jpl.nasa.gov/pub/naif/pds/data/jno-j_e_ss-spice-6-v1.0/jnosp_1000/"
METAFILE = os.path.join(SPICE_DIR, 'JunoMetaKernel.txt')
JUNO_SPECIFIC_KERNEL = "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/pck00010.tpc"

class InitSpiceKernels(object):

    isKernelsExists = False
    indexfile = os.path.join(SPICE_DIR, 'index_spice.csv')
    
    def __init__(self, dodownload=False):
        """
        Creating a JunoSpiceKernel object downloads the index file and
        identifies which kernels are useful. If not previously done,
        the relevant kernels are downloaded.
        """
        
        self.isKernelsExists = os.path.exists(KERNEL_DIR)
        if not self.isKernelsExists or dodownload:
            print("Kernels don't exist, downloading ...")
            self.download_index_file()
            self.find_kernels()
            self.download_kernels()
            self.download_juno_specific_kernel()
        else:
            self.find_kernels()
            
        self.create_meta_kernel_file()

    def download_index_file(self):
        """Downloads the spice kernel index file from NAIF."""
        INDEX_URL = JUNO_SPICE_URL + "index/index.tab"
        download_file(INDEX_URL, outfolder=SPICE_DIR)
        os.rename(
            os.path.join(SPICE_DIR, 'index.tab'),
            os.path.join(SPICE_DIR, 'index_spice.csv'))
        

    def find_kernels(self):
        """Finds the relevant kernels by parsing the index file."""
        columns = [
            'START_TIME',
            'STOP_TIME',
            'FILE_SPECIFICATION_NAME',
            'DATASET_ID',
            'CREATION_TIME',
            'IDENTIFIER',
            'RELEASE_TIME',
            'KERNEL_TYPE',
            'KERNEL_FILENAME',
            'VOLUME'
        ]

        df = pd.read_csv(self.indexfile, names=columns)

        for col in ['KERNEL_TYPE', 'KERNEL_FILENAME', 'FILE_SPECIFICATION_NAME']:
            df[col] = df[col].apply(lambda x: x.strip())


        self.kernelsdf = df.loc[df['KERNEL_TYPE'].isin(
            ['SPK', 'FK', 'LSK', 'PCK', 'SCLK'])]


    def download_kernels(self):
        """Download the relevant kernels from NAIF."""
        for i, row in self.kernelsdf.iterrows():
            file_url = os.path.join(
                JUNO_SPICE_URL,
                os.path.split(row['FILE_SPECIFICATION_NAME'])[0],
                row['KERNEL_FILENAME']
            )
            download_file(file_url, outfolder=KERNEL_DIR)
        self.isKernelsExists = True
        
    
    def download_juno_specific_kernel(self):
        download_file(JUNO_SPECIFIC_KERNEL, outfolder=KERNEL_DIR)
        

    def create_meta_kernel_file(self):
        """Create the meta spice kernel file to be used with spiceypy"""
        with open(METAFILE, 'w') as f:
            f.write("\\begindata\n")
            f.write("KERNELS_TO_LOAD=(\n")


            for i, kernel in self.kernelsdf.iterrows():
                filename = os.path.join(KERNEL_DIR, kernel['KERNEL_FILENAME'])
                f.write("'{:s}',\n".format(filename))
            
            f.write("'{:s}',\n".format(
                os.path.join(KERNEL_DIR, os.path.split(JUNO_SPECIFIC_KERNEL)[-1])))
            f.write(")\n")
            f.write("\\begintext\n")


def get_juno_position(times, coordsystem='JUNO_JSS', center='JUPITER'):
    """
    Fetch juno position for a list of spacecraft times.
    Input should be in datetime.
    """
    metafile = os.path.join(SPICE_DIR, 'JunoMetaKernel.txt')
    spice.furnsh(metafile)
    spice.furnsh(os.path.join(SPICE_DIR, 'kernels', 'juno_phi_o_v1.tf'))
    
    times_converted = []
    for time in times:
        times_converted.append(
            spice.str2et(dt.datetime.strftime(time, '%Y-%m-%dT%H:%M:%S.%f')))

    position, lightTimes = spice.spkpos(
        'Juno', times_converted, coordsystem, 'NONE', center)

    spice.kclear()
    
    if 'JUPITER' in center:
        position = np.asarray(position) / 71492
    elif '503' in center or 'GANYMEDE' in center:
        position = np.asarray(position) / 2634

    return position, lightTimes


def get_ganymede_pos(time, coordsystem='IAU_JUPITER'):
    
    metafile = os.path.join(SPICE_DIR, 'JunoMetaKernel.txt')
    spice.furnsh(metafile)
    
    time_converted = spice.str2et(time.strftime('%Y-%m-%dT%H:%M:%S.%f'))
    
    position, lighttime = spice.spkpos(
        'GANYMEDE', time_converted, coordsystem, 'NONE', 'JUPITER')
    
    spice.kclear()
    
    return position / 71492