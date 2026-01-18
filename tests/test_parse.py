import datetime as dt
import numpy as np
from junopy import mag
from junopy import waves

def test_parse_mag():
    """
    Test to parse a magnetometer file.
    """
    start = dt.datetime(2020, 1, 1, 2)
    stop = dt.datetime(2020, 1, 1, 23)

    mag_data = mag.load_mag_data(start, stop, coordsystem="JUNO_JSS")
    assert np.isclose(mag_data.df['Bx'].iloc[0], 9.46) 
    assert np.isclose(mag_data.df['Btheta'].iloc[0], 2.11828) 

def test_parse_wavessurvey():
    """
    Test to parse a waves-survey data file.
    """
    start = dt.datetime(2020, 1, 1, 2)
    _ = dt.datetime(2020, 1, 1, 23)

    electric_spectra, magnetic_spectra = waves.get_waves_data_for_date(date=start.date()) 
    assert np.isclose(electric_spectra.channels['freq'][1], 97.656)
    assert np.isclose(magnetic_spectra.channels['freq'][1], 97.656)
    assert np.isclose(electric_spectra.df['HFR_HI_30'].iloc[0], 1.0865e-13) 
    assert np.isclose(magnetic_spectra.df['LFR_LO_B_01'].iloc[0], 2.6888e-06) 


# import os
# import test_index

# from JunoPy.mag import MagData
# from JunoPy.waves import WavesData
# from JunoPy.jedi import JediData

# def test_parse_mag():
#     """
#     Test to parse a magnetometer file.
#     """
#     filename = 'data/mag/fgm_jno_l3_2018043pc_r1s_v01.sts'

#     if not os.path.exists(filename):
#         test_index.test_download_list_files()

#     magdata = MagData(filename=filename)
#     magdata.convert_to_spherical()
#     print(magdata.df.head())


# def test_parse_waves_survey():
#     """
#     Test to parse Waves-Survey file.
#     """
#     filename = "data/waves-survey/WAV_2018043T000000_B_V01.CSV"
#     if not os.path.exists(filename):
#         test_index.test_download_list_files()

#     wavesdata = WavesData(filename=filename)
#     print(wavesdata.df.head())


# def test_parse_jedi():
#     """
#     """
#     filename = 'data/jedi/JED_180_LOERSESP_CDR_2018043_V03.TAB'
#     if not os.path.exists(filename):
#         test_index.test_download_list_files()

#     jedidata = JediData(filename=filename)
#     print(jedidata.df.head())
