import datetime as dt
import glob

import JunoPy

from JunoPy.index import IndexFile
from JunoPy.mag import MagData

st = dt.datetime(2016, 7, 8, 17, 30)
et = dt.datetime(2019, 5, 31, 16)

# magindex = IndexFile(dataset='mag')
# magfilelist = magindex.get_files_between_dates(
#     st, et, DoDownload=True, sid='SS 1 SECOND')

for magfile in glob.glob('data/mag/*.sts'):
    print('Converting file - {:s}'.format(magfile), end='\r')
    magdata = MagData(filename=magfile, ishdf=False)
    magdata.convert_to_spherical()
    magdata.save_to_hdf5()