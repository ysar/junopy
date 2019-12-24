
import datetime as dt

from JunoPy.index import IndexFile
from JunoPy.library import DATASETS

def test_parse_indexfile():
    """
    Test to parse index file.
    """
    st = dt.datetime(2018, 1, 11, 0, 40)
    et = dt.datetime(2018, 1, 15, 1, 23)

    magindex = IndexFile(dataset='mag')

    filelist = magindex.get_files_between_dates(
        st, et, sid='PC 1 SECOND')

    filelist_check = [
        'DATA/JUPITER/PC/PERI-10/fgm_jno_l3_2018012pc_r1s_v01.lbl',
        'DATA/JUPITER/PC/PERI-11/fgm_jno_l3_2018013pc_r1s_v01.lbl',
        'DATA/JUPITER/PC/PERI-11/fgm_jno_l3_2018014pc_r1s_v01.lbl',
    ]

    assert(x in filelist_check for x in filelist)


def test_download_list_files():
    """
    Download a list of files between two particular times
    """
    st = dt.datetime(2018, 2, 11, 0, 40)
    et = dt.datetime(2018, 2, 13, 1, 23)


    for dset in DATASETS.keys():
        print('\n')
        index = IndexFile(dataset=dset)
        flist = index.get_files_between_dates(st, et, DoDownload=True)
