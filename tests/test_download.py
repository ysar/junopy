
import os
import datetime as dt

from JunoPy.download import download_file
from JunoPy.index import IndexFile
from JunoPy.library import PDS_URL, DATASETS


def test_download():
    """
    Unit test to test file download.
    """
    url = 'https://pds-ppi.igpp.ucla.edu/ditdos/download?id=' + \
          'pds://PPI/JNO-E_J_SS-WAV-3-CDR-SRVFULL-V1.0/DATA/WAVES_SURVEY/' + \
          '2016005_JUPITER_APPROACH/2016_021/'

    file1 = 'WAV_2016021T000000_B_V01.CSV'
    file2 = 'WAV_2016021T000000_E_V01.CSV'
    for f in [url + file1, url + file2]:
        download_file(f)

    assert(os.path.exists(file1))
    assert(os.path.exists(file2))

    os.remove(file1)
    os.remove(file2)


def test_download_indexfile():
    """
    Download an index file.
    """
    magindex = IndexFile(dataset='mag')
    assert(os.path.exists('index/INDEX_mag.csv'))
