
from pathlib import Path
from junopy import utils
from junopy import config

CURRENT_DIR = Path(__file__).parent
config.ConfigVars.set_root(CURRENT_DIR)

class TestDownloadUtils:

    def test_download(self):
        """
        Unit test to test file download.
        """
        url = 'https://pds-ppi.igpp.ucla.edu/ditdos/download?id=' + \
            'pds://PPI/JNO-E_J_SS-WAV-3-CDR-SRVFULL-V1.0/DATA/WAVES_SURVEY/' + \
            '2016005_JUPITER_APPROACH/2016_021/'

        file1 = 'WAV_2016021T000000_B_V01.CSV'
        file2 = 'WAV_2016021T000000_E_V01.CSV'

        tmp_dir = config.ConfigVars.get("TMP_DIR")

        for f in [file1, file2]:
            tmp_file = tmp_dir / f
            download_success = utils.download.download_file(url + f, tmp_file)
            assert download_success
            assert tmp_file.exists()
            tmp_file.unlink()


    def test_download_parallel(self):
        """
        Unit test to test file download.
        """
        url = 'https://pds-ppi.igpp.ucla.edu/ditdos/download?id=' + \
            'pds://PPI/JNO-E_J_SS-WAV-3-CDR-SRVFULL-V1.0/DATA/WAVES_SURVEY/' + \
            '2016005_JUPITER_APPROACH/2016_021/'

        file1 = 'WAV_2016021T000000_B_V01.CSV'
        file2 = 'WAV_2016021T000000_E_V01.CSV'

        tmp_dir = config.ConfigVars.get("TMP_DIR")

        download_success = utils.download.download_multiple_files(
            [url + file1, url + file2],
            [tmp_dir / file1, tmp_dir / file2]
        )

        assert download_success

        (tmp_dir / file1).unlink()
        (tmp_dir / file2).unlink()




# def test_download_mag_data():
#     start = dt.datetime(2019, 1, 2)
#     stop = dt.datetime(2019, 1, 3)

#     magIndex = index.IndexFile(dataset="mag", do_download=True)

#     magIndex.get_files_between_dates(
#         start, stop,
#         sid="SS 1 MINUTE",
#         do_download=True)

#     filelist = [
#         "data/mag/2019/001/fgm_jno_l3_2019001ss_r60s_v01.sts",
#         "data/mag/2019/002/fgm_jno_l3_2019002ss_r60s_v01.sts",
#         "data/mag/2019/003/fgm_jno_l3_2019003ss_r60s_v01.sts"
#     ]

#     for fname in filelist:
#         assert os.path.exists(fname)
#         os.remove(fname)


# def test_download_jade_data():
#     start = dt.datetime(2019, 1, 2)
#     stop = dt.datetime(2019, 1, 3)
#     jadeIndex = index.IndexFile(dataset="jade", do_download=True)
#     jadeIndex.get_files_between_dates(
#         start, stop,
#         sid="_ELC_",
#         do_download=True,
#         do_parallel=True,
#         get_label=True)

#     filelist = [
#         "data/jade/2019/001/JAD_L30_LRS_ELC_ANY_CNT_2019001_V04.DAT",
#         "data/jade/2019/001/JAD_L30_LRS_ELC_ANY_CNT_2019001_V04.LBL",
#         "data/jade/2019/002/JAD_L30_LRS_ELC_ANY_CNT_2019002_V04.DAT",
#         "data/jade/2019/002/JAD_L30_LRS_ELC_ANY_CNT_2019002_V04.LBL"
#     ]
#     for fname in filelist:
#         assert os.path.exists(fname)
#         os.remove(fname)

# def test_download_jedi_data():
#     start = dt.datetime(2019, 1, 2)
#     stop = dt.datetime(2019, 1, 3)
#     jediIndex = index.IndexFile(dataset="jedi", do_download=True)
#     jediIndex.get_files_between_dates(
#         start, stop,
#         sid="ERSESP",
#         do_download=True,
#         do_parallel=True,
#     )

#     filelist = [
#         "data/jedi/2019/001/JED_090_LOERSESP_CDR_2019001_V03.TAB",
#         "data/jedi/2019/001/JED_180_LOERSESP_CDR_2019001_V03.TAB",
#         "data/jedi/2019/001/JED_270_HIERSESP_CDR_2019001_V03.TAB",
#         "data/jedi/2019/002/JED_090_LOERSESP_CDR_2019002_V03.TAB",
#         "data/jedi/2019/002/JED_180_LOERSESP_CDR_2019002_V03.TAB",
#         "data/jedi/2019/002/JED_270_HIERSESP_CDR_2019002_V03.TAB",
#         "data/jedi/2019/003/JED_090_LOERSESP_CDR_2019003_V03.TAB",
#         "data/jedi/2019/003/JED_270_HIERSESP_CDR_2019003_V03.TAB",
#     ]

#     for fname in filelist:
#         assert os.path.exists(fname)
#         os.remove(fname)

# def test_download_waves_survey_data():
#     start = dt.datetime(2019, 1, 2)
#     stop = dt.datetime(2019, 1, 3)
#     wavesSurveyIndex = index.IndexFile(dataset="waves-survey", doDownload=True)
#     wavesSurveyIndex.get_files_between_dates(
#         start, stop,
#         sid="SURVEY",
#         doDownload=True)
