import datetime as dt
import glob
from pathlib import Path

import pandas as pd

PATTERNS = {
    "ELECTRONS": "S_ELC_",
    "ION_SPECIES": "_ION_ANY_",
    "ION_LOGICALS": "_ION_LOG_",
    "ION_TOF": "_ION_TOF_",
}


def find_jade_data_files(
    start: dt.datetime,
    stop: dt.datetime,
    get_hrs: bool = True,
    get_lrs: bool = True,
    datatype: str = "ELECTRONS",
) -> list[Path]:
    """
    Return combined JADE low and high-res data.

    Args:
        start (str or datetime-like): Start date of data to retrieve.
        stop (str or datetime-like): End date of data to retrieve.
        get_hrs (bool, optional): Whether to include high-res data. Default is True
        get_lrs (bool, optional): Whether to include low-res data. Default is True
        datatype (str, optional): Type of data to retrieve. Options are
            "ELECTRONS", "ION_SPECIES", "ION_LOGICALS", "ION_TOF".
            Default is "ELECTRONS".

    Returns:
        list of Path: List of file paths matching criteria.
    """
    raise Exception("This is unimplemented")
    # if datatype not in PATTERNS:
    #     raise ValueError(
    #         f"Invalid datatype for JADE '{datatype}'. Must be one of {list(PATTERNS.keys())}."
    #     )

    # sid = PATTERNS[datatype]

    # # Download files that are missing

    # relevant_files = []
    # for day in pd.date_range(start, stop):
    #     relevant_files += glob.glob(
    #         config.ConfigVars.get("JADE_DATA_DIR")
    #         + day.strftime("/%Y/%j/")
    #         + f"*{PATTERNS[datatype]}*.DAT*"
    #     )

    #     if not get_lrs:
    #         relevant_files = [
    #             f for f in relevant_files if "_LRS" not in f and "_LLS" not in f
    #         ]

    #     if not get_hrs:
    #         relevant_files = [
    #             f for f in relevant_files if "_HRS" not in f and "_HLS" not in f
    #         ]

    # relevant_files = [Path(f) for f in relevant_files]

    # return relevant_files
