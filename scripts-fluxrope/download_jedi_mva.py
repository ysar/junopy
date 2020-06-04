#!/usr/bin/env python3
"""
Script to download only the JEDI data relevant to the MVA analysis.
"""

import datetime as dt

from JunoPy.index import IndexFile
from JunoPy.jedi import JediData


def read_csv_file():
    start_times = []
    end_times = []
    durations = []
    with open('mva_epochs.csv', 'r') as f:
        f.readline()
        for line in f:
            start_times.append(dt.datetime.strptime(line.split(',')[0], '%Y-%jT%H:%M:%S'))
            end_times.append(dt.datetime.strptime(line.split(',')[1], '%Y-%jT%H:%M:%S'))
            durations.append(float(line.split(',')[2]))
    return start_times, end_times, durations
            
    

    
# ================================================================
if __name__ == '__main__':
    
    start_times, end_times, durations = read_csv_file()
    
    jediindex = IndexFile(dataset='jedi')
    
    for start, end in zip(start_times, end_times):
    
        jedifilelist = jediindex.get_files_between_dates(
            start, end, DoDownload=True, sid='ERSESP')
    
    