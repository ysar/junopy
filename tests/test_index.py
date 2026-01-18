
import datetime as dt
from pathlib import Path
from junopy import index
from junopy import config

def test_download_indexfiles():
    """
    Downloads many index files and their associated data files for each
    instrument type.
    """
    # Clear all existing index files.
    for p in config.ConfigVars.get("INDEX_DIR").iterdir():
        p.unlink()

    index_mag = index.IndexFile(dataset="mag", do_download=True)
    assert index_mag.index_file.exists()

    index_jade = index.IndexFile(dataset="jade", do_download=True)
    assert index_jade.index_file.exists()

    index_jedi = index.IndexFile(dataset="jedi", do_download=True)
    assert index_jedi.index_file.exists()

    index_waves_survey = index.IndexFile(dataset="waves-survey", do_download=True)
    assert index_waves_survey.index_file.exists()


def test_parse_mag_indexfile():
    """
    Test to parse magnetometer index file.
    """
    st = dt.datetime(2020, 1, 1)
    et = dt.datetime(2020, 1, 1, 1)

    mag_index = index.IndexFile(dataset='mag')

    received = mag_index.get_files_between_dates(st, et, sid='PC 1 SECOND')

    expected = [
        Path('.junodata/data/mag/2020/001/fgm_jno_l3_2020001pc_r1s_v01.sts').absolute(),
    ]

    assert set(received) == set(expected)


def test_parse_wavessurvey_indexfile():
    """
    Test to parse waves survey index file.
    """
    st = dt.datetime(2020, 1, 1, 2)
    et = dt.datetime(2020, 1, 1, 23)

    waves_index = index.IndexFile(dataset='waves-survey')

    received = waves_index.get_files_between_dates(st, et)

    expected = [
        Path('.junodata/data/waves-survey/2020/001/WAV_2020001T000000_B_V02.CSV').absolute(),
        Path('.junodata/data/waves-survey/2020/001/WAV_2020001T000000_E_V02.CSV').absolute(),
    ]

    assert set(received) == set(expected)

# def test_download_list_files():
#     """
#     Download a list of files between two particular times
#     """
#     st = dt.datetime(2018, 2, 11, 0, 40)
#     et = dt.datetime(2018, 2, 13, 1, 23)


#     for dset in DATASETS.keys():
#         print('\n')
#         index = IndexFile(dataset=dset)
#         flist = index.get_files_between_dates(st, et, DoDownload=True)
