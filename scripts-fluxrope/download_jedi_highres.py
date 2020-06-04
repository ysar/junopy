#!/usr/bin/env python3
"""
Script to download only the JEDI data relevant to the MVA analysis.
"""

import datetime as dt

from JunoPy.index import IndexFile
from JunoPy.jedi import JediData


st = dt.datetime(2018, 7, 8, 17, 30)
et = dt.datetime(2019, 5, 31, 16)


jediindex = IndexFile(dataset='jedi')
jedifilelist = jediindex.get_files_between_dates(
    st, et, DoDownload=True, sid='HIERSESP')
    
