import datetime as dt
import pandas as pd
import glob

from JunoPy import config

PATTERNS = {
    "ELECTRONS": "_ELC_",
    "ION_SPECIES": "_ION_ANY_",
    "ION_LOGICALS": "_ION_LOG_",
    "ION_TOF": "_ION_TOF_",
}


def get_jade_data_files(start, stop, **kwargs):
    """
    Return combined JADE low and high-res data.
    """
    getHRS = kwargs.pop("getHRS", True)
    getLRS = kwargs.pop("getLRS", True)
    datatype = kwargs.pop("datatype", "ELECTRONS")

    for day in pd.date_range(start, stop):
        
        relevantFiles = glob.glob(
            config.ConfigVars.get("JADE_DATA_DIR")\
            + day.strftime("/%Y/%j/") \
            + f"*{PATTERNS[datatype]}*.DAT*")
        
        if not getLRS:
            relevantFiles = [
                f for f in relevantFiles if "_LRS" not in f and "_LLS" not in f]
            
        if not getHRS:
            relevantFiles = [
                f for f in relevantFiles if "_HRS" not in f and "_HLS" not in f]
            
        return relevantFiles