
import os
import datetime as dt
import pandas as pd

from JunoPy.library import DATASETS, PDS_URL, EXTENSIONS, TIMEFORMAT
from JunoPy.download import download_file, download_multiple_files

class IndexFile(object):

    outfolder = 'index'

    def __init__(self, dataset=''):
        """
        This script operates on the index files, which are present in every
        PDS repository. We need to fetch the index files and identify the URL
        of the relavent data files.
        """
        if dataset:
            self.dataset = dataset
            self.url = PDS_URL + DATASETS[dataset] + '/INDEX/INDEX.TAB'
        else:
            raise Exception('Specify dataset.')

        self.download_index_file()
        self.parse_index_file()


    def download_index_file(self):
        """
        Fetch the index file.
        """
        download_file(self.url, outfolder=self.outfolder)
        assert(os.path.exists(os.path.join(self.outfolder, 'INDEX.TAB')))

        os.rename(
            os.path.join(self.outfolder, 'INDEX.TAB'),
            os.path.join(self.outfolder, 'index_{:s}.csv'.format(self.dataset)))

        self.filename = os.path.join(
            self.outfolder, 'index_{:s}.csv'.format(self.dataset))


    def parse_index_file(self):
        """
        Read in the index file. PDS has conveniently chosen to use comma
        separated values.
        """
        index = pd.read_csv(self.filename)

        keys = [c.strip() for c in index.columns]
        index.columns = keys
        index[index.columns] = index.apply(lambda x: x.str.strip())

        # --------------------------------------------------------------------
        # index['START_TIME'] = pd.to_datetime(
        #     index['START_TIME'], format=TIMEFORMAT[self.dataset])
        # index['STOP_TIME'] = pd.to_datetime(
        #     index['STOP_TIME'], format=TIMEFORMAT[self.dataset])
        # index['CR_DATE'] = pd.to_datetime(
        #     index['CR_DATE'], format=TIMEFORMAT[self.dataset][:-3])

        # Unfortunately some PDS data has different time format bugs, so the
        # above vectorized implementation doesnt work well, specifically for
        # JADE data.
        # --------------------------------------------------------------------
        for i, data in index.iterrows():

            try:
                data['START_TIME'] = pd.to_datetime(
                    data['START_TIME'], format=TIMEFORMAT[self.dataset])
            except ValueError:
                data['START_TIME'] = pd.to_datetime(
                    data['START_TIME'], format=TIMEFORMAT[self.dataset][:-3])

            try:
                data['STOP_TIME'] = pd.to_datetime(
                    data['STOP_TIME'], format=TIMEFORMAT[self.dataset])
            except ValueError:
                data['STOP_TIME'] = pd.to_datetime(
                    data['STOP_TIME'], format=TIMEFORMAT[self.dataset][:-3])

            data['CR_DATE'] = pd.to_datetime(
                data['CR_DATE'], format=TIMEFORMAT[self.dataset][:-3])

        self.df = index


    def get_files_between_dates(
        self, startdate, enddate, sid='', DoDownload=False):
        """
        Given two dates, find all index files corresponding to time between.
        Optional argument - SID, DoDownload
        """
        df = self.df
        df_masked = df.loc[(df['START_TIME'] > startdate) &
                           (df['STOP_TIME'] <= enddate)]

        if df_masked.empty:
            df_masked = df.loc[(df['START_TIME'] < startdate) &
                               (df['STOP_TIME'] >= enddate) &
                               (startdate - df['START_TIME'] < dt.timedelta(hours=24))]

            print('Date range too short - downloading one file.')

        # print(df_masked)

        if sid:
            df_masked = df_masked[df_masked['SID'].str.contains(sid)]

        filelist = [x.strip() for x in df_masked['FILE_SPECIFICATION_NAME']]

        if DoDownload:
            urls = []

            for f in filelist:

                urls.append(PDS_URL + DATASETS[self.dataset] + '/' + \
                           f[:-4] + EXTENSIONS[self.dataset])

            download_multiple_files(
                urls, outfolder=os.path.join('data', self.dataset))

        return filelist
