import pandas as pd
import glob
import os

from JunoPy import index
from JunoPy import library


JADE_SID_TYPES = {
    'ELECTRONS': '_ELC_',
    'ION_SPECIES': '_ION_ANY_',
    'ION_TOF': '_TOF_',
    'ION_LOGICALS': '_LOG_'
}

def jade_download_data(start, stop, dtype=None, datafolder='data/jade/'):
    """
    Download JADE .DAT and .LBL files corresponding to the given times.
    """
    if not dtype:
        raise Exception("JADE data type not provided.")

    sid = JADE_SID_TYPES[dtype]

    jadeindex = index.IndexFile(dataset='jade')

    for day in pd.date_range(start, stop):

        dataFiles = glob.glob(
            datafolder + f"/{day.strftime('%Y')}/{day.strftime('%Y%j')}/" \
            + f'/{dtype}/' \
            + f"/JAD_*{sid}*{day.strftime('%Y%j')}*.DAT.gz")

        lblFiles = glob.glob(
            datafolder + f"/{day.strftime('%Y')}/{day.strftime('%Y%j')}/" \
            + f'/{dtype}/' \
            + f"/JAD_*{sid}*{day.strftime('%Y%j')}*.LBL.gz")

        if len(dataFiles) == 0 or len(lblFiles) == 0:
            filesExist = False
        else:
            filesExist = True
            print(dataFiles)

        if not filesExist:
            
            outfolder = \
                datafolder \
                + f"/{day.strftime('%Y')}/{day.strftime('%Y%j')}/" \
                + f'/{dtype}/'
            
            downloadedFiles = jadeindex.get_files_between_dates(
                day, day, 
                DoDownload=True, sid=sid, 
                outfolder=outfolder, getLBL=True)
            
            for filename in downloadedFiles:
                library.create_gz_zip_file(os.path.join(
                    outfolder, filename))
    
        