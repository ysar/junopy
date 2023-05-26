"""
Module to work with Juno spice kernels.

For now we will only use spice to obtain spacecraft position.
"""
import os
import numpy as np
import datetime as dt
import pandas as pd
import spiceypy as spice
import ftplib
import JupiterMag

from JunoPy.download import download_file, download_multiple_files
from contextlib import ContextDecorator

SPICE_DIR = 'spice'
KERNEL_DIR = os.path.join(SPICE_DIR, 'kernels')
JUNO_SPICE_URL = "https://naif.jpl.nasa.gov/pub/naif/pds/data/jno-j_e_ss-spice-6-v1.0/jnosp_1000/"
METAFILE = os.path.join(SPICE_DIR, 'JunoMetaKernel.txt')


class withspiceloaded(ContextDecorator):
    def __enter__(self):
        metafile = os.path.join(SPICE_DIR, 'JunoMetaKernel.txt')
        spice.furnsh(metafile)
        spice.furnsh(os.path.join(SPICE_DIR, 'kernels', 'juno_phi_o_v1.tf'))
    
    def __exit__(self, *exc):
        metafile = os.path.join(SPICE_DIR, 'JunoMetaKernel.txt')
        spice.unload(metafile)
        spice.unload(os.path.join(SPICE_DIR, 'kernels', 'juno_phi_o_v1.tf'))


def download_juno_kernels(doDownload=True):
    """Downloads all relevant Juno related SPICE kernels into a directory."""
    
    kernel_urls = [
        'https://naif.jpl.nasa.gov/pub/naif/JUNO/kernels/fk/juno_v12.tf',
        'https://naif.jpl.nasa.gov/pub/naif/JUNO/kernels/lsk/naif0012.tls',
        'https://naif.jpl.nasa.gov/pub/naif/JUNO/kernels/pck/pck00010.tpc',
        'https://naif.jpl.nasa.gov/pub/naif/JUNO/misc/jszalay/juno_phi_o_v1.tf',
        'https://naif.jpl.nasa.gov/pub/naif/JUNO/kernels/sclk/JNO_SCLKSCET.00134.tsc',    
    ]

    with ftplib.FTP("naif.jpl.nasa.gov") as ftp:

        ftp.login()
        ftp.cwd('/pub/naif/JUNO/kernels/spk')

        for fname in ftp.nlst():
            if not ('.bsp' in fname[-4:]):
                continue

            if (fname[0:3] in ['de4', 'jup']) or \
                ('spk_rec' in fname[0:7]):
                kernel_urls.append(
                    'https://naif.jpl.nasa.gov' + ftp.pwd() + '/' + fname)

        ftp.cwd('/pub/naif/JUNO/kernels/ck')

        for fname in ftp.nlst():
            if 'juno_sc' in fname[:7]:
                kernel_urls.append(
                    'https://naif.jpl.nasa.gov' + ftp.pwd() + '/' + fname)
                
    if doDownload:
        print("Downloading SPICE kernels ...")
        download_multiple_files(
            kernel_urls, os.path.join(SPICE_DIR, 'kernels'))
    
    print("Creating Meta Kernel file ...")
    with open(METAFILE, 'w') as f:
        f.write("\\begindata\n")
        f.write("KERNELS_TO_LOAD=(\n")
        
        for kernel in kernel_urls:
            f.write("'spice/kernels/{:s}',\n".format(kernel.split('/')[-1]))

        f.write(")\n")
        f.write("\\begintext\n")

@withspiceloaded()
def get_juno_position(times, coordsystem='JUNO_JSS', center='JUPITER'):
    """
    Fetch juno position for a list of spacecraft times.
    Input should be in datetime.
    """
    # metafile = os.path.join(SPICE_DIR, 'JunoMetaKernel.txt')
    # spice.furnsh(metafile)
    # spice.furnsh(os.path.join(SPICE_DIR, 'kernels', 'juno_phi_o_v1.tf'))
    
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

    # spice.unload(metafile)
    # spice.unload(os.path.join(SPICE_DIR, 'kernels', 'juno_phi_o_v1.tf'))

    return position, lightTimes



@withspiceloaded()
def get_ganymede_pos(time, coordsystem='IAU_JUPITER'):
    
    metafile = os.path.join(SPICE_DIR, 'JunoMetaKernel.txt')
    spice.furnsh(metafile)
    
    time_converted = spice.str2et(time.strftime('%Y-%m-%dT%H:%M:%S.%f'))
    
    position, _ = spice.spkpos(
        'GANYMEDE', time_converted, coordsystem, 'NONE', 'JUPITER')
    
    spice.kclear()
    
    return position / 71492


# def spice_load_kernels():
#     """Load SPICE kernels into cache (?)"""
#     metafile = os.path.join(SPICE_DIR, 'JunoMetaKernel.txt')
#     spice.furnsh(metafile)
#     spice.furnsh(os.path.join(SPICE_DIR, 'kernels', 'juno_phi_o_v1.tf'))


# def spice_unload_kernels():
#     """Unload SPICE kernels."""
#     metafile = os.path.join(SPICE_DIR, 'JunoMetaKernel.txt')
#     spice.unload(metafile)
#     spice.unload(os.path.join(SPICE_DIR, 'kernels', 'juno_phi_o_v1.tf'))