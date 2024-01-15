"""
Module to work with Juno spice kernels.

For now we will only use spice to obtain spacecraft position.
"""
import os
import re
import numpy as np
import datetime as dt
import spiceypy as spice
import urllib

from JunoPy.download import download_multiple_files, download_file
from contextlib import ContextDecorator

R_JUPITER = 71492    # km
R_GANYMEDE = 2634    # km

SPICE_DIR = 'spice'
KERNEL_DIR = os.path.join(SPICE_DIR, 'kernels')
JUNO_SPICE_URL = "https://naif.jpl.nasa.gov/" + \
                 "pub/naif/pds/data/jno-j_e_ss-spice-6-v1.0/jnosp_1000/"
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

    URL_TOP = "https://naif.jpl.nasa.gov/pub/naif/JUNO/kernels/"
    
    kernel_urls = []

    # Identify latest files from the fk, lsk, pck, and sclk folders
    extensions = {
        'fk': '.tf',
        'pck': '.tpc',
        'lsk': '.tls',
        'sclk': 'tsc'
    }

    for folder in extensions:
        with urllib.request.urlopen(URL_TOP + '{:s}/'.format(folder)) as page:
            pattern = r'<a href="(\S+)' + '{:s}">'.format(extensions[folder])
            fname = re.findall(
                pattern, page.read().decode('UTF-8'))[-1]
            kernel_urls.append(
                URL_TOP + '{:s}/{:s}{:s}'.format(
                    folder, fname, extensions[folder]))
            

    # Identify all reconstructed ck files for download
    folder = 'ck'
    pattern = r'<a href="juno_sc_rec_(1[6-9]|2[0-5])(\S+).bc"'
        # Note: This Regex pattern is complicated because we only need the
        # reconstructed ck kernels for years 2016-2025

    with urllib.request.urlopen(URL_TOP + 'ck/') as page:
        files = re.findall(
            pattern, page.read().decode('UTF-8'))
    
    for f in files:
        kernel_urls.append(
            URL_TOP + 'ck/' + 'juno_sc_rec_' + f[0] + f[1] + '.bc')
        
    # Identify all relevant reconstructed spk files for download
    folder = 'spk'
    pattern = r'<a href="(\S+).bsp">'
    with urllib.request.urlopen(URL_TOP + 'spk/') as page:
        files = re.findall(
            pattern, page.read().decode('UTF-8'))
        
    for f in files:
        if f[:3] in ['de4', 'jup'] or \
            'spk_rec' in f[:7]:
            kernel_urls.append(
                URL_TOP + 'spk/' + f + '.bsp')
            
    # Add the PhiO coordinate system kernel
    kernel_urls.append(
        URL_TOP + '../misc/jszalay/juno_phi_o_v1.tf')
                
    if doDownload:
        print("Downloading SPICE kernels ...")

        for url in kernel_urls:
            print(url)
            download_file(
                url, 
                os.path.join(SPICE_DIR, 'kernels', os.path.split(url)[-1]))
            
        # download_multiple_files(
        #     kernel_urls, os.path.join(SPICE_DIR, 'kernels'))
    
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
    times_converted = []
    for time in times:
        times_converted.append(
            spice.str2et(dt.datetime.strftime(time, '%Y-%m-%dT%H:%M:%S.%f')))

    position, lightTimes = spice.spkpos(
        'Juno', times_converted, coordsystem, 'NONE', center)

    spice.kclear()
    
    if 'JUPITER' in center:
        position = np.asarray(position) / R_JUPITER
    elif '503' in center or 'GANYMEDE' in center:
        position = np.asarray(position) / R_GANYMEDE

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
