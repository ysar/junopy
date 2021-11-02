"""
Contains methods to download files using :mod:`wget`.
"""

import os
import wget
import glob
import time

def download_file(url, outfolder=''):
    """
    Downloads a file given its url. Replaces any file with the same name in the
    output directory.
    
    :param str url: The URL for the file to be downloaded.
    :param str outfolder: Output directory. Default is current directory.
    """
    if outfolder and not os.path.exists(outfolder):
        os.makedirs(outfolder)

    outfile = os.path.join(outfolder, url.split('/')[-1])
    if os.path.exists(outfile):
        os.remove(outfile)

    wget.download(url, bar=None, out=outfile)


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
                wget.download,
                args=(url,),
                kwds={
                    'bar': None,
                    'out': outfile
                 }
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
