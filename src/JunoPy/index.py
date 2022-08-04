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

from JunoPy.library import DATASETS, PDS_URL, EXTENSIONS, TIMEFORMAT
from JunoPy.download import download_file, download_multiple_files

class IndexFile(object):
    """
    This class operates on the index files, which are present in every
    PDS repository. We need to fetch the index files and identify the URL
    of the relevant data files.
    
    :param str dataset: String describing the dataset. 
                        (mag/waves-survey/jedi, Default - '').
    """

    outfolder = 'index'

    def __init__(self, dataset='', dodownload=False):

        if dataset:
            self.dataset = dataset
            self.url = PDS_URL + DATASETS[dataset] + '/INDEX/INDEX.TAB'
        else:
            raise Exception('Specify dataset.')
        
        if not os.path.exists(self.outfolder):
            os.mkdir(self.outfolder)
            
        if not os.path.exists(
            os.path.join(self.outfolder, 
                         'index_{:s}.csv'.format(self.dataset))) or dodownload:
            print('Downloading index file...')
            self.download_index_file()
            
        self.filename = os.path.join(
            self.outfolder, 'index_{:s}.csv'.format(self.dataset))
        
        
        self.parse_index_file()


    def download_index_file(self):
        """
        Downloads the index file from PDS for the specified dataset.
        """
        
        outfile = os.path.join(self.outfolder, 'INDEX.TAB')
        download_file(self.url, os.path.join(self.outfolder, 'INDEX.TAB'))
        assert(os.path.exists(outfile))

        os.rename(
            os.path.join(self.outfolder, 'INDEX.TAB'),
            os.path.join(self.outfolder, 'index_{:s}.csv'.format(self.dataset))) 


    def parse_index_file(self):
        """
        Read in the index file. PDS has conveniently chosen to use comma
        separated values.
        """
        columns = ['VOLUME_ID', 'SID', 'DATA_SET_ID', 'PRODUCT_ID', 'START_TIME',
           'STOP_TIME', 'FILE_SPECIFICATION_NAME', 'CR_DATE',
           'PRODUCT_LABEL_MD5CHECKSUM']

        dateparse = lambda x: dt.datetime.strptime(x, TIMEFORMAT[self.dataset])

        index = pd.read_csv(
            self.filename,
            names=columns, 
            parse_dates=['START_TIME', 'STOP_TIME'], date_parser=dateparse,
            skiprows=1)
        
        self.df = index
        
#         index = pd.read_csv(self.filename)

#         keys = [c.strip() for c in index.columns]
#         index.columns = keys
#         index[index.columns] = index.apply(lambda x: x.str.strip())

#         # --------------------------------------------------------------------
#         # index['START_TIME'] = pd.to_datetime(
#         #     index['START_TIME'], format=TIMEFORMAT[self.dataset])
#         # index['STOP_TIME'] = pd.to_datetime(
#         #     index['STOP_TIME'], format=TIMEFORMAT[self.dataset])
#         # index['CR_DATE'] = pd.to_datetime(
#         #     index['CR_DATE'], format=TIMEFORMAT[self.dataset][:-3])

#         # Unfortunately some PDS data has different time format bugs, so the
#         # above vectorized implementation doesnt work well, specifically for
#         # JADE data.
#         # --------------------------------------------------------------------
        
#         print(index['START_TIME'])
#         raise Exception()
        
#         for i, data in index.iterrows():

#             try:
#                 data['START_TIME'] = pd.to_datetime(
#                     data['START_TIME'], format=TIMEFORMAT[self.dataset])
#             except ValueError:
#                 data['START_TIME'] = pd.to_datetime(
#                     data['START_TIME'], format=TIMEFORMAT[self.dataset][:-3])

#             try:
#                 data['STOP_TIME'] = pd.to_datetime(
#                     data['STOP_TIME'], format=TIMEFORMAT[self.dataset])
#             except ValueError:
#                 data['STOP_TIME'] = pd.to_datetime(
#                     data['STOP_TIME'], format=TIMEFORMAT[self.dataset][:-3])

#             data['CR_DATE'] = pd.to_datetime(
#                 data['CR_DATE'], format=TIMEFORMAT[self.dataset][:-3])
            
        

#         self.df = index


    def get_files_between_dates(
        self, startdate, enddate, sid='', DoDownload=False):
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
        
        # print(df_masked.info())
        filelist = []
  
        for i, row in df_masked.iterrows():
            
            doAppend = False
            
            # if (startdate - row['START_TIME'])
            
            
            
            if (np.abs((startdate - row['START_TIME']).total_seconds()) < \
                np.abs(dt.timedelta(hours=24).total_seconds())):
                doAppend = True
                
            if (np.abs((enddate - row['STOP_TIME']).total_seconds()) < \
                np.abs(dt.timedelta(hours=24).total_seconds())):
                doAppend = True
                
            if (row['START_TIME'] > startdate and row['STOP_TIME'] < enddate):
                doAppend = True
                
            if doAppend:
                filelist.append(row['FILE_SPECIFICATION_NAME'].strip())
                

        # filelist = [x.strip() for x in df_masked['FILE_SPECIFICATION_NAME']]
        # print(startdate, enddate, filelist) # filelist)
        
        if DoDownload and len(filelist) > 0:
            urls = []

            for f in filelist:

                # urls.append(PDS_URL + DATASETS[self.dataset] + '/' + \
                #            f[:-4] + EXTENSIONS[self.dataset])
                
                print('Downloading file - {:s}'.format(f[:-4]))

                downfile = os.path.join(
                    'data', 
                    self.dataset, 
                    f[:-4].split('/')[-1] + EXTENSIONS[self.dataset])

                if not os.path.exists(downfile):

                    download_file(
                        PDS_URL + DATASETS[self.dataset] + '/' + \
                        f[:-4] + EXTENSIONS[self.dataset],
                        downfile)

            # download_multiple_files(
            #     urls, outfolder=os.path.join('data', self.dataset))
            
            

        return filelist
