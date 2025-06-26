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

"""
Module to work with Juno SPICE kernels.
"""
import re
import numpy as np
import datetime as dt
import spiceypy
import urllib.request
import logging
from typing import Sequence
from contextlib import ContextDecorator

from junopy import config
from junopy import utils
from junopy.constants import R_JUPITER, R_GANYMEDE

URL_JUNO_NAIF = "https://naif.jpl.nasa.gov/pub/naif/JUNO/kernels/"
# This is where the final kernels are stored? Use the other NAIF URL
# JUNO_SPICE_URL = "https://naif.jpl.nasa.gov/" + \
#                  "pub/naif/pds/data/jno-j_e_ss-spice-6-v1.0/jnosp_1000/"


class withspiceloaded(ContextDecorator):
    """
    This is a context decorator to wrap functions that require loading of the
    SPICE kernels.  It provides an ergonomic way to load kernels either by
    using it as a decorator on functions or as a `with` clause.

    Usage:
        ```python
        # Option 1 (as decorator)
        @withspiceloaded()
        def method_that_uses_spice():
            ...

        # Option 2 (inline)
        with withspiceloaded() as context:
            ...
        ```
    """
    def __enter__(self):
        spice_dir = config.ConfigVars.get("SPICE_DIR")
        metafile = spice_dir / 'JunoMetaKernel.txt'
        spiceypy.furnsh(str(metafile))

    def __exit__(self, *exc):
        spice_dir = config.ConfigVars.get("SPICE_DIR")
        metafile = spice_dir / 'JunoMetaKernel.txt'
        spiceypy.unload(str(metafile))


def download_juno_kernels(do_download: bool = True):
    """
    Downloads all relevant Juno related SPICE kernels into a directory.

    Args:
        do_download (bool): Whether to (re)download the Juno kernels. If False,
            no kernels will be downloaded but the meta-kernel file will be
            re-written.  This is useful for changing the working directory
            while reusing kernels stored elsewhere.  In that case, the relative
            links in the meta-kernel file need to be re-defined relative to the
            new working directory. See `::: junopy.config` for more details.
    """

    spice_dir = config.ConfigVars.get("SPICE_DIR")

    kernel_urls = []

    # Identify latest files from the fk, lsk, pck, and sclk folders
    extensions = {
        'fk': '.tf',
        'pck': '.tpc',
        'lsk': '.tls',
        'sclk': 'tsc'
    }

    for folder in extensions:
        with urllib.request.urlopen(URL_JUNO_NAIF + f'{folder:s}/') as page:
            pattern = r'<a href="(\S+)' + '{:s}">'.format(extensions[folder])
            fname = re.findall(pattern, page.read().decode("UTF-8"))[-1]
            kernel_urls.append(
                URL_JUNO_NAIF
                + "{:s}/{:s}{:s}".format(folder, fname, extensions[folder])
            )


    # Identify all reconstructed ck files for download
    folder = 'ck'
    pattern = r'<a href="juno_sc_rec_(1[6-9]|2[0-5])(\S+).bc"'
        # Note: This Regex pattern is a hack because we only need the
        # reconstructed ck kernels for years 2016-2025.  May need to modify
        # this for later years.

    with urllib.request.urlopen(URL_JUNO_NAIF + 'ck/') as page:
        files = re.findall(
            pattern, page.read().decode('UTF-8'))

    for f in files:
        kernel_urls.append(
            URL_JUNO_NAIF + 'ck/' + 'juno_sc_rec_' + f[0] + f[1] + '.bc')

    # Identify all relevant reconstructed spk files for download
    folder = 'spk'
    pattern = r'<a href="(\S+).bsp">'
    with urllib.request.urlopen(URL_JUNO_NAIF + 'spk/') as page:
        files = re.findall(
            pattern, page.read().decode('UTF-8'))

    for f in files:
        if f[:3] in ['de4', 'jup'] or \
            'spk_rec' in f[:7]:
            kernel_urls.append(
                URL_JUNO_NAIF + 'spk/' + f + '.bsp')

    # Add the PhiO coordinate system kernel
    kernel_urls.append(URL_JUNO_NAIF + "../misc/jszalay/juno_phi_o_v1.tf")

    if do_download:
        logging.info("Downloading SPICE kernels ...")

        output_folder = spice_dir / "kernels"
        output_files = [output_folder / url.split("/")[-1] for url in kernel_urls]

        utils.download.download_multiple_files(kernel_urls, output_files=output_files)

    _create_meta_kernel_file(kernel_urls)


def _create_meta_kernel_file(kernel_urls: list[str]) -> None:
    """
    Create the SPICE meta kernel.

    Args:
        kernels_urls (list[str]): List of kernel URLs to be downloaded/used.
    """
    spice_dir = config.ConfigVars.get("SPICE_DIR")
    logging.info("Creating Meta Kernel file ...")
    metafile = spice_dir / 'JunoMetaKernel.txt'

    buffer = "\\begindata\n"
    buffer += "KERNELS_TO_LOAD=(\n"
    for kernel in kernel_urls:
        buffer += "'"
        buffer += str(spice_dir / "kernels" / kernel.split('/')[-1])
        buffer += "',\n"

    buffer += ")\n\\begintext\n"

    with open(metafile, 'w') as f:
        f.write(buffer)


@withspiceloaded()
def get_juno_position(
    times: Sequence[dt.datetime],
    coordsystem: str = "JUNO_JSS",
    center: str ="JUPITER"
) -> tuple[np.ndarray, float | np.ndarray]:
    """
    Fetch juno position for a list of spacecraft times.

    Args:
        times (list[dt.datetime]): List of datetimes
        coordsystem (str): The coordinate system of the output coordinates
        center (str): The center of the coordinate system, as passed to SPICE

    Returns:
        positions (np.ndarray): The coordinates of Juno at the requested times
        light_times (float | np.ndarray): The light time correction
    """
    num_times = len(times)
    times_et = np.zeros(num_times)
    for i, time in enumerate(times):
        times_et[i] = \
            spiceypy.str2et(dt.datetime.strftime(time, '%Y-%m-%dT%H:%M:%S.%f'))

    position, light_times = spiceypy.spkpos(
        'Juno', times_et, coordsystem, 'NONE', center)

    spiceypy.kclear()

    if 'JUPITER' in center:
        position = np.asarray(position) / R_JUPITER
    elif '503' in center or 'GANYMEDE' in center:
        position = np.asarray(position) / R_GANYMEDE

    return position, light_times


@withspiceloaded()
def get_object_position(
    times: list[dt.datetime],
    object: str,
    coordsystem: str = 'JUNO_JSS',
    center: str = 'JUPITER'):
    """
    Fetch object position for a list of spacecraft times.
    Similar to `get_juno_position` but for any object recognized by SPICE.

    Args:
        times (list[dt.datetime]): List of datetimes
        object (str): Name or SPICE identifier for target object
        coordsystem (str): The coordinate system of the output coordinates
        center (str): The center of the coordinate system, as passed to SPICE

    Returns:
        positions (np.ndarray): The coordinates of the object at the requested times
        light_times (float | np.ndarray): The light time correction
    """
    num_times = len(times)
    times_et = np.zeros(num_times)
    for i, time in enumerate(times):
        times_et[i] = \
            spiceypy.str2et(dt.datetime.strftime(time, '%Y-%m-%dT%H:%M:%S.%f'))

    position, light_times = spiceypy.spkpos(
        object, times_et, coordsystem, 'NONE', center)

    spiceypy.kclear()

    if 'JUPITER' in center:
        position = np.asarray(position) / R_JUPITER
    elif '503' in center or 'GANYMEDE' in center:
        position = np.asarray(position) / R_GANYMEDE

    return position, light_times


@withspiceloaded()
def convert_to_utc(et):
    """
    Convert time from ephemeris time to UTC format.

    Args:
        et (float): Ephemeris time

    Returns:
        utc (np.datetime64 | np.ndarray): Time in UTC
    """
    return spiceypy.et2utc(et, "ISOD", 6)
