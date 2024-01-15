"""
Contains methods to download files.
"""
import time
import urllib.request

def download_file(fileurl, fileoutname, CHUNK=16*1024*16):
    """Opens a URL and downloads it's contents onto disk."""
    with urllib.request.urlopen(fileurl) as g:
        with open(fileoutname, 'wb+') as f:
            # f.write(g.read())
            for chunk in iter(lambda: g.read(CHUNK), b''):
                f.write(chunk) #.decode('UTF-8'))
            

def download_multiple_files(urls, outFiles):
    """
    Use :mod:`multiprocessing` to download multiple files at a time.
    
    :param urls: List of URLs of the files to be downloaded.
    :type urls: list(str)
    :param str outfolder: Output directory. Default is current directory.
    """
    from multiprocessing import Pool

    assert (len(urls) > 0)
    assert (len(urls) == len(outFiles))


    with Pool(processes=4) as pool:

        processes = []
        for url, fOut in zip(urls, outFiles):
            processes.append(
                pool.apply_async(download_file, (url, fOut,)))
            
        while True:
            incompleteCount = sum(1 for p in processes if not p.ready())
            perComplete = float(len(urls) - incompleteCount) / len(urls)
            if incompleteCount == 0:
                print('{:6.1f}% |'.format(100) + '#' * 50 + '|')
                break
            else:
                print('{:6.1f}% |'.format(perComplete * 100) +
                      '#' * int(perComplete * 50) +
                      '.' * (50 - int(perComplete * 50)) +
                      '|', end='\r')

            time.sleep(1)