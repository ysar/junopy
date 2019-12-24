
import os
import test_index

from JunoPy.mag import MagData
from JunoPy.waves import WavesData
from JunoPy.jedi import JediData

def test_parse_mag():
    """
    Test to parse a magnetometer file.
    """
    filename = 'data/mag/fgm_jno_l3_2018043pc_r1s_v01.sts'

    if not os.path.exists(filename):
        test_index.test_download_list_files()

    magdata = MagData(filename=filename)
    magdata.convert_to_spherical()
    print(magdata.df.head())


def test_parse_waves_survey():
    """
    Test to parse Waves-Survey file.
    """
    filename = "data/waves-survey/WAV_2018043T000000_B_V01.CSV"
    if not os.path.exists(filename):
        test_index.test_download_list_files()

    wavesdata = WavesData(filename=filename)
    print(wavesdata.df.head())


def test_parse_jedi():
    """
    """
    filename = 'data/jedi/JED_180_LOERSESP_CDR_2018043_V03.TAB'
    if not os.path.exists(filename):
        test_index.test_download_list_files()

    jedidata = JediData(filename=filename)
    print(jedidata.df.head())
