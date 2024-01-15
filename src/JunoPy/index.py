"""
:mod:`JunoPy` downloads data from the PDS by making use of the standard PDS3/PDS4
data formats. Each dataset stored on the PDS contains an `INDEX` directory with 
an ASCII-text `INDEX` file. This file contains a list of the files associated 
with this dataset, where each file typically covers a range of times or dates.
Also given in the index file is the URL or path to the file from the root 
directory of the project. 

"""
import os
import datetime as dt
import pandas as pd
import numpy as np

from pathlib import Path, PosixPath
from JunoPy import download
from JunoPy.library import DATA_ID, PDS_URL, EXTENSIONS, TIMEFORMAT


class IndexFile(object):
    """
    This class operates on the index files, which are present in every
    PDS repository. We need to fetch the index files and identify the URL
    of the relevant data files.
    
    :param str dataset: String describing the dataset. 
                        (mag/waves-survey/jedi, Default - '').
    """

    outfolder = Path('index')

    def __init__(self, dataset='', doDownload=False):

        if dataset:
            self.dataset = dataset
            self.url = PDS_URL + DATA_ID[dataset] + '/INDEX/INDEX.TAB'
        else:
            raise Exception('Specify dataset.')
        
        self.outfolder.mkdir(parents=True, exist_ok=True)
            
        indexfile = self.outfolder / f"index_{self.dataset}.csv"
        
        if not indexfile.exists() or doDownload:
            print('Downloading index file ...')
            self.download_index_file()
            
        self.filename = str(indexfile)
        
        self.parse_index_file()


    def download_index_file(self):
        """
        Downloads the index file from PDS for the specified dataset.
        """
        outfile = self.outfolder / 'INDEX.TAB'

        download.download_file(self.url, str(outfile))

        assert(outfile.exists())

        outfile = outfile.rename(
            self.outfolder / f"index_{self.dataset}.csv") 


    def parse_index_file(self):
        """
        Read in the index file. PDS has conveniently chosen to use comma
        separated values.
        """
        columns = ['VOLUME_ID', 'SID', 'DATA_SET_ID', 'PRODUCT_ID', 'START_TIME',
           'STOP_TIME', 'FILE_SPECIFICATION_NAME', 'CR_DATE',
           'PRODUCT_LABEL_MD5CHECKSUM']

        # if self.dataset == 'waves-survey':
        #     dateparse = pd.api.extensions.no_default
        # else:
        #     dateparse = lambda x: dt.datetime.strptime(x, )

        index = pd.read_csv(
            self.filename,
            names=columns, 
            parse_dates=['START_TIME', 'STOP_TIME'], 
            date_format=TIMEFORMAT[self.dataset],
            skiprows=1)
        
        index['START_TIME'] = pd.to_datetime(index['START_TIME'])
        index['STOP_TIME'] = pd.to_datetime(index['STOP_TIME'])

        self.df = index


    def get_files_between_dates(
        self, startdate, enddate, 
        sid='', 
        doDownload=False, 
        outfolder=None, 
        getLBL=False,
        doMultiprocessing=False):
        """
        Given two dates, find all index files corresponding to time between.
        
        :param startdate: Start time for the needed interval.
        :param enddate: End time for the needed interval. 
        :param str sid: SID is a string that desribes different data in the 
                        the same `INDEX` file. For e.g. magnetometer data with
                        60 s resolution has '60s' appended to the filename
                        in the `SID` column. 
        :param bool DoDownload: Boolean to start download (False by default).
        
        """
        df = self.df
        df_masked = pd.DataFrame()
        
        if sid:
            df_masked = df[df['SID'].str.contains(sid)]
        
        if not outfolder:
            outfolder = Path('data') / self.dataset
        else:
            outfolder = Path(outfolder)

        outfolder.mkdir(parents=True, exist_ok=True)

        # Ignore unwanted files from waves
        if self.dataset == 'waves-survey':
            for s in ['BINT', 'BUNC', 'EINT', 'EUNC']:
                df_masked = df_masked[
                    ~df_masked['PRODUCT_ID'].str.contains(s)]

        # New version with pathlib
        relevantFileUrls = []
        relevantFileOutNames = []

        for _, row in df_masked.iterrows():

            doAppend = False
            
            if (np.abs((startdate - row['START_TIME']).total_seconds()) < \
                np.abs(dt.timedelta(hours=24).total_seconds())):
                doAppend = True
                
            if (np.abs((enddate - row['STOP_TIME']).total_seconds()) < \
                np.abs(dt.timedelta(hours=24).total_seconds())):
                doAppend = True
                
            if (row['START_TIME'] > startdate and row['STOP_TIME'] < enddate):
                doAppend = True

            if doAppend:

                fileurl = PosixPath(
                    row['FILE_SPECIFICATION_NAME'].strip()
                    ).with_suffix(EXTENSIONS[self.dataset])
                
                relevantFileUrls.append(
                    PosixPath(PDS_URL) \
                    / DATA_ID[self.dataset] \
                    / fileurl)
                
                relevantFileOutNames.append(
                    outfolder \
                    / Path(row['START_TIME'].strftime('%Y')) \
                    / Path(row['START_TIME'].strftime('%j')) \
                    / fileurl.name 
                )
        
        if not doDownload:
            return relevantFileOutNames
    
        if not relevantFileUrls:
            print("No relevant files to download from the PDS for the"
                  " selected time.")
            return relevantFileOutNames

        # Create list of files actually to be downloaded
        downloadFileUrls = []
        downloadFileOutNames = []

        for url, fOut in zip(relevantFileUrls, relevantFileOutNames):

            # First we check if file exists
            if Path(fOut).exists(): continue
            if Path(fOut.with_suffix('.hdf5').exists()): continue

            fOut.parent.mkdir(parents=True, exist_ok=True)

            downloadFileOutNames.append(fOut)
            downloadFileUrls.append(url)

            if getLBL:
                downloadFileOutNames.append(fOut.with_suffix('.LBL'))
                downloadFileUrls.append(url.with_suffix('.LBL'))


        if not downloadFileUrls:
            print(f"Files already exist. Nothing to download. ({self.dataset})")
            return relevantFileOutNames

        if doMultiprocessing:
            download.download_multiple_files(
                downloadFileUrls, downloadFileOutNames)
        else:
            for url, fOut in zip(downloadFileUrls, downloadFileUrls):
                print('Downloading to file - {:s}'.format(fOut))
                download.download_file(url, fOut)
        
        return relevantFileOutNames
