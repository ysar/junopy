
import os
import datetime as dt

from JunoPy.download import download_file
from JunoPy import index
from JunoPy import config

config.ConfigVars.set_root(os.path.dirname(__file__))

def test_download():
    """
    Unit test to test file download.
    """
    url = 'https://pds-ppi.igpp.ucla.edu/ditdos/download?id=' + \
          'pds://PPI/JNO-E_J_SS-WAV-3-CDR-SRVFULL-V1.0/DATA/WAVES_SURVEY/' + \
          '2016005_JUPITER_APPROACH/2016_021/'

    file1 = 'WAV_2016021T000000_B_V01.CSV'
    file2 = 'WAV_2016021T000000_E_V01.CSV'
    for f in [file1, file2]:
        download_file(url + f, "tmp-" + f)

    assert(os.path.exists("tmp-" + file1))
    assert(os.path.exists("tmp-" + file2))

    os.remove("tmp-" + file1)
    os.remove("tmp-" + file2)


def test_download_indexfiles():
    """
    Downloads many index files and their associated data files for each 
    instrument type.
    """
    magIndex = index.IndexFile(dataset="mag")
    assert(os.path.exists("index/index_mag.csv"))

    jadeIndex = index.IndexFile(dataset="jade")
    assert(os.path.exists("index/index_jade.csv"))

    jediIndex = index.IndexFile(dataset="jedi")
    assert(os.path.exists("index/index_jedi.csv"))

    wavesSurveyIndex = index.IndexFile(dataset="waves-survey")
    assert(os.path.exists("index/index_waves-survey.csv"))


def test_download_mag_data():
    start = dt.datetime(2019, 1, 2)
    stop = dt.datetime(2019, 1, 3)
    
    magIndex = index.IndexFile(dataset="mag", doDownload=True)
    
    magIndex.get_files_between_dates(
        start, stop, 
        sid="SS 1 MINUTE",
        doDownload=True)
    
    filelist = [
        "data/mag/2019/001/fgm_jno_l3_2019001ss_r60s_v01.sts",
        "data/mag/2019/002/fgm_jno_l3_2019002ss_r60s_v01.sts",
        "data/mag/2019/003/fgm_jno_l3_2019003ss_r60s_v01.sts"
    ]

    for fname in filelist:
        assert os.path.exists(fname)
        os.remove(fname)
    
    
def test_download_jade_data():
    start = dt.datetime(2019, 1, 2)
    stop = dt.datetime(2019, 1, 3)
    jadeIndex = index.IndexFile(dataset="jade", doDownload=True)
    jadeIndex.get_files_between_dates(
        start, stop, 
        sid="_ELC_",
        doDownload=True, getLBL=True)
    
    filelist = [
        "data/jade/2019/001/JAD_L30_LRS_ELC_ANY_CNT_2019001_V04.DAT",
        "data/jade/2019/001/JAD_L30_LRS_ELC_ANY_CNT_2019001_V04.LBL",
        "data/jade/2019/002/JAD_L30_LRS_ELC_ANY_CNT_2019002_V04.DAT",
        "data/jade/2019/002/JAD_L30_LRS_ELC_ANY_CNT_2019002_V04.LBL"
    ]
    for fname in filelist:
        assert os.path.exists(fname)
        os.remove(fname)
    
def test_download_jedi_data():
    start = dt.datetime(2019, 1, 2)
    stop = dt.datetime(2019, 1, 3)
    jediIndex = index.IndexFile(dataset="jedi", doDownload=True)
    jediIndex.get_files_between_dates(
        start, stop, 
        sid="ERSESP",
        doDownload=True)
    
    filelist = [
        "data/jedi/2019/001/JED_090_LOERSESP_CDR_2019001_V03.TAB",
        "data/jedi/2019/001/JED_180_LOERSESP_CDR_2019001_V03.TAB",
        "data/jedi/2019/001/JED_270_HIERSESP_CDR_2019001_V03.TAB",
        "data/jedi/2019/002/JED_090_LOERSESP_CDR_2019002_V03.TAB",
        "data/jedi/2019/002/JED_180_LOERSESP_CDR_2019002_V03.TAB",
        "data/jedi/2019/002/JED_270_HIERSESP_CDR_2019002_V03.TAB",
        "data/jedi/2019/003/JED_090_LOERSESP_CDR_2019003_V03.TAB",
        "data/jedi/2019/003/JED_270_HIERSESP_CDR_2019003_V03.TAB",
    ]
    
    for fname in filelist:
        assert os.path.exists(fname)
        os.remove(fname)

# def test_download_waves_survey_data():
#     start = dt.datetime(2019, 1, 2)
#     stop = dt.datetime(2019, 1, 3)
#     wavesSurveyIndex = index.IndexFile(dataset="waves-survey", doDownload=True)
#     wavesSurveyIndex.get_files_between_dates(
#         start, stop, 
#         sid="SURVEY",
#         doDownload=True)
