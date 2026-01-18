# Copyright (c) 2018 Yash Sarkango
# Licensed under the MIT License

"""
Module to work with Juno SPICE kernels.
"""

import datetime as dt
import logging
import re
import urllib.request
from contextlib import ContextDecorator
from pathlib import Path
from typing import Sequence

import numpy as np
import spiceypy

from junopy import config
from junopy.constants import R_GANYMEDE, R_JUPITER
from junopy.utils import download

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
        metafile = spice_dir / "JunoMetaKernel.txt"
        try:
            spiceypy.furnsh(str(metafile))

        except spiceypy.utils.exceptions.SpiceNOSUCHFILE:
            logging.info("A SPICE kernel is missing. Attemping to download ...")
            download_juno_kernels()

            logging.info("SPICE kernel update successful.")
            logging.info("Retrying to load SPICE kernels ...")
            spiceypy.furnsh(str(metafile))

        except Exception as e:
            logging.error("Failed to load SPICE kernels: %s", e)
            raise

    def __exit__(self, *exc):
        spice_dir = config.ConfigVars.get("SPICE_DIR")
        metafile = spice_dir / "JunoMetaKernel.txt"
        spiceypy.unload(str(metafile))


def _split_string_for_spiceypy(input_string: str) -> str:
    """
    Split a string into a format that is compatible with SPICEypy.

    Args:
        input_string (str): The input string to be split.

    Returns:
        str: The input string split into a format compatible with SPICEypy.
    """
    string_len = 70

    if len(input_string) <= string_len:
        return "'" + input_string + "',\n"

    output_string = ""
    i = 0
    segment = input_string[:string_len]
    while len(segment) == string_len:
        output_string += "\n    '" + segment + "+',"
        i += 1
        segment = input_string[i * string_len : (i + 1) * string_len]

    output_string += "\n    '" + segment + "'"

    output_string += "\n,\n"

    return output_string


def identify_juno_kernels() -> list[str]:
    """
    Identify the Juno SPICE kernels that are currently available at NAIF.  This
    method by itself does not download the kernels, only identifies them.

    Returns:
        list[str]: URLs of the latest Juno SPICE kernels that are available for
        download.
    """
    kernel_urls = []

    # Identify latest files from the fk, lsk, pck, and sclk folders
    extensions = {"fk": ".tf", "pck": ".tpc", "lsk": ".tls", "sclk": "tsc"}

    for folder in extensions:
        with urllib.request.urlopen(URL_JUNO_NAIF + f"{folder:s}/") as page:
            pattern = r'<a href="(\S+)' + '{:s}">'.format(extensions[folder])
            fname = re.findall(pattern, page.read().decode("UTF-8"))[-1]
            kernel_urls.append(
                URL_JUNO_NAIF
                + "{:s}/{:s}{:s}".format(folder, fname, extensions[folder])
            )

    # Identify all reconstructed ck files for download
    folder = "ck"
    pattern = r'<a href="juno_sc_rec_(1[6-9]|2[0-5])(\S+).bc"'
    # Note: This Regex pattern is a hack because we only need the
    # reconstructed ck kernels for years 2016-2025.  May need to modify
    # this for later years.

    with urllib.request.urlopen(URL_JUNO_NAIF + "ck/") as page:
        files = re.findall(pattern, page.read().decode("UTF-8"))

    for f in files:
        kernel_urls.append(URL_JUNO_NAIF + "ck/" + "juno_sc_rec_" + f[0] + f[1] + ".bc")

    # Identify all relevant reconstructed spk files for download
    folder = "spk"
    pattern = r'<a href="(\S+).bsp">'
    with urllib.request.urlopen(URL_JUNO_NAIF + "spk/") as page:
        files = re.findall(pattern, page.read().decode("UTF-8"))

    for f in files:
        if f[:3] in ["de4", "jup"] or "spk_rec" in f[:7]:
            kernel_urls.append(URL_JUNO_NAIF + "spk/" + f + ".bsp")

    # Add the PhiO coordinate system kernel
    kernel_urls.append(URL_JUNO_NAIF + "../misc/jszalay/juno_phi_o_v1.tf")

    return kernel_urls


def download_juno_kernels(do_download: bool = True) -> None:
    """
    Downloads all relevant Juno related SPICE kernels into the SPICE_DIR directory.

    Args:
        do_download (bool): If True, forces download of Juno kernels even if
            they already exist. If False, only checks for existing kernels,
            updates the meta kernel file and does not download anything.
    """

    spice_dir = config.ConfigVars.get("SPICE_DIR")

    logging.info("Checking for missing SPICE kernels ...")

    kernel_urls = identify_juno_kernels()
    output_folder = spice_dir / "kernels"
    output_folder.mkdir(parents=True, exist_ok=True)

    # Filter out existing files to avoid re-downloading
    kernels_to_download = [
        url for url in kernel_urls if not (output_folder / url.split("/")[-1]).exists()
    ]

    output_files = [
        output_folder / url.split("/")[-1]
        for url in kernels_to_download
        if not (output_folder / url.split("/")[-1]).exists()
    ]

    if do_download and len(kernels_to_download) > 0:
        download.download_multiple_files(kernels_to_download, output_files=output_files)
    else:
        logging.info("No new kernels to download")

    logging.info("Remaking the meta kernel file ...")
    _remake_meta_kernel_file(kernel_urls)


def _remake_meta_kernel_file(kernel_urls: list[str]) -> None:
    """
    Create the SPICE meta kernel.

    Args:
        kernels_urls (list[str]): List of kernel URLs to be downloaded/used.
    """
    spice_dir = config.ConfigVars.get("SPICE_DIR")

    custom_frame_file = _create_custom_frame_file()

    logging.info("Creating Meta Kernel file ...")
    metafile = spice_dir / "JunoMetaKernel.txt"

    buffer = "\\begindata\n"
    buffer += "KERNELS_TO_LOAD=(\n"
    for kernel in kernel_urls:
        kernel_name = str((spice_dir / "kernels" / kernel.split("/")[-1]).absolute())
        buffer += _split_string_for_spiceypy(kernel_name)

    buffer += _split_string_for_spiceypy(str(custom_frame_file.absolute()))
    buffer += ")\n\\begintext\n"

    with open(metafile, "w") as f:
        f.write(buffer)


def _create_custom_frame_file() -> Path:
    """
    This is a hack to add frame definitions to the SPICE kernel directory that
    are not included in the Juno SPICE kernels.  This is necessary for some of
    the coordinate systems used in the JADE data.
    """
    spice_dir = config.ConfigVars.get("SPICE_DIR")
    custom_frame_file = spice_dir / "kernels" / "_custom_frames.tf"

    with open(custom_frame_file, "w") as f:
        f.write(r"""
        \begindata
            FRAME_JUNO_JADE_DESPUN         = -61951
            FRAME_-61951_NAME              = 'JUNO_JADE_DESPUN'
            FRAME_-61951_CLASS             = 5
            FRAME_-61951_CLASS_ID          = -61951
            FRAME_-61951_CENTER            = -61
            FRAME_-61951_RELATIVE          = 'J2000'
            FRAME_-61951_DEF_STYLE         = 'PARAMETERIZED'
            FRAME_-61951_FAMILY            = 'TWO-VECTOR'
            FRAME_-61951_PRI_AXIS          = 'Z'
            FRAME_-61951_PRI_VECTOR_DEF    = 'CONSTANT'
            FRAME_-61951_PRI_FRAME         = 'JUNO_SPACECRAFT'
            FRAME_-61951_PRI_SPEC          = 'RECTANGULAR'
            FRAME_-61951_PRI_VECTOR        = ( 0, 0, 1 )
            FRAME_-61951_SEC_AXIS          = 'X'
            FRAME_-61951_SEC_VECTOR_DEF    = 'CONSTANT'
            FRAME_-61951_SEC_FRAME         = 'ECLIPJ2000'
            FRAME_-61951_SEC_SPEC          = 'RECTANGULAR'
            FRAME_-61951_SEC_VECTOR        = ( 0, 0, 1 )
        \begintext
        """)

    return custom_frame_file


@withspiceloaded()
def get_juno_position(
    times: Sequence[dt.datetime], coordsystem: str = "JUNO_JSS", center: str = "JUPITER"
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
        times_et[i] = spiceypy.str2et(
            dt.datetime.strftime(time, "%Y-%m-%dT%H:%M:%S.%f")
        )

    position, light_times = spiceypy.spkpos(
        "Juno", times_et, coordsystem, "NONE", center
    )

    spiceypy.kclear()

    if "JUPITER" in center:
        position = np.asarray(position) / R_JUPITER
    elif "503" in center or "GANYMEDE" in center:
        position = np.asarray(position) / R_GANYMEDE

    return position, light_times


@withspiceloaded()
def get_object_position(
    times: list[dt.datetime],
    object: str,
    coordsystem: str = "JUNO_JSS",
    center: str = "JUPITER",
):
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
        times_et[i] = spiceypy.str2et(
            dt.datetime.strftime(time, "%Y-%m-%dT%H:%M:%S.%f")
        )

    position, light_times = spiceypy.spkpos(
        object, times_et, coordsystem, "NONE", center
    )

    spiceypy.kclear()

    if "JUPITER" in center:
        position = np.asarray(position) / R_JUPITER
    elif "503" in center or "GANYMEDE" in center:
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


@withspiceloaded()
def convert_to_utc_v(et):
    """
    Convert time from ephemeris time to UTC format. (Vectorized)

    Args:
        et (np.ndarray): Array of ephemeris time.

    Returns:
        utc (np.ndarray): Array of times in UTC
    """
    return spiceypy.et2utc(et, "ISOD", 6)
