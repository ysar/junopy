"""
Contains methods to download files using :mod:`wget`.
"""

import os
# import wget
import glob
import time
import urllib.request

def download_file(fileurl, fileoutname):
    """Opens a URL and downloads it's contents onto disk."""
    with urllib.request.urlopen(fileurl) as g:
        with open(fileoutname, 'wb+') as f:
            f.write(g.read()) #.decode('UTF-8'))
            

def download_multiple_files(urls, outfolder=''):
    """
    Use :mod:`multiprocessing` to download multiple files at a time.
    
    :param urls: List of URLs of the files to be downloaded.
    :type urls: list(str)
    :param str outfolder: Output directory. Default is current directory.
    """
    from multiprocessing import Pool

    if outfolder and not os.path.exists(outfolder):
        os.makedirs(outfolder)

    outfiles = []
    print('Fetching data ... ')
    print('Progress bar updates every 5 seconds ...')
    pool = Pool()
    down = []
    numfiles = len(urls)

    for url in urls:

        outfile = os.path.join(outfolder, url.split('/')[-1])
        if os.path.exists(outfile):
            # os.remove(outfile)
            continue

        down.append(
            pool.apply_async(
                download_file,
                args=(url, outfile),
            )
        )

    while True:

        incomplete_count = sum(1 for x in down if not x.ready())
        per_complete = float(numfiles - incomplete_count) / numfiles
        if incomplete_count == 0:
            print('{:6.1f}% |'.format(100) + '#' * 50 + '|')
            break
        else:
            print('{:6.1f}% |'.format(per_complete * 100) +
                  '#' * int(per_complete * 50) +
                  '.' * (50 - int(per_complete * 50)) +
                  '|', end='\r')

        time.sleep(1)

    pool.close()
    pool.join()
