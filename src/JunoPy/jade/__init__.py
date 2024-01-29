# from .__download import JADE_SID_TYPES
# from .__download import jade_download_data
from .read import read_jade_binary_file, create_dict_from_jade_df

# from .get import get_jade_data_files
# from .read import read_jade_binary_file, JadeData
# from .calib import get_Geff_electrons

# Let's think a bit about how the JADE processing can be improved. Firstly, 
# after reading in the raw data, we should separate it by sensor. 
# JadeData
#   |--- JadeSensorData
#   |--- JadeSensorData
#   |--- JadeSensorData
# Another question to ask is if it is necessary to separate high-res vs low-res
# data?
import datetime as dt
import numpy as np
import logging


class JadeData(object):
    """
    This class should serve as a common base for the read-in Jade data files,
    but it will have some organization within to separate by each sensor.
    """
    def __init__(self, **kwargs):
        
        self.filename = kwargs.pop('filename', None)
        if not self.filename:
            return None

        start          = kwargs.pop('start', None)
        stop           = kwargs.pop('stop', None)
        doProcess      = kwargs.pop('doProcess', False)
        self.magfolder = kwargs.pop('magfolder', None)
        self.dtbin     = kwargs.pop('dtbin', dt.timedelta(seconds=120))
        self.nPAbins   = kwargs.pop('nPAbins', 12)
        self.units     = kwargs.pop('units', 'COUNTS')
        self.paMin     = kwargs.pop('paMin', 0)
        self.paMax     = kwargs.pop('paMax', 180)
        self.energyMin = kwargs.pop('energyMin', 1)
        self.energyMax = kwargs.pop('energyMax', 100e3)

        self.df_raw = read_jade_binary_file(self.filename)

        # Truncate data based on what is required
        if start and stop:
            self.df_raw = self.df_raw[start: stop]

        if not doProcess:
            return None
        
        self.data_raw = create_dict_from_jade_df(self.df_raw)
        
        self.correct_invalid_values()

    
    def correct_invalid_values(self):
        
        # if "_V04" in self.filename:
        #     missingConstant = -999999
        # else:
        #     missingConstant = -1
        for k in self.data_raw.keys():

        self.data_raw[self.data_raw < 0] = np.NaN



    def 


