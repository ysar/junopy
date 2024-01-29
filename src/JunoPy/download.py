"""
Contains methods to download files.
"""
import os
import time
import uuid
import urllib.request
import urllib.error

from JunoPy import config

def download_file(url, outFile, CHUNK=16*1024*1024):
    """Opens a URL and downloads it's contents onto disk."""
    tmpFile = os.path.join(config.ConfigVars.get("TMP_DIR"), str(uuid.uuid4()))
    
    try:
        with urllib.request.urlopen(url) as g:
            with open(tmpFile, 'wb+') as f:
                # f.write(g.read())
                for chunk in iter(lambda: g.read(CHUNK), b''):
                    f.write(chunk) #.decode('UTF-8'))

        os.replace(tmpFile, outFile)

    except urllib.error.URLError:
        print(80*'')
        print('Could not download - \n{:s}', url)
        print(80*'')
            

def download_multiple_files(urls, outFiles, dchunk=16*1024*1024):
    """
    Use :mod:`multiprocessing` to download multiple files at a time.
    
    :param urls: List of URLs of the files to be downloaded.
    :type urls: list(str)
    :param str outfolder: Output directory. Default is current directory.
    """
    raise Exception("Threaded downloads are currently bugged. Do not use.")
    from threading import Thread

    assert (len(urls) > 0)
    assert (len(urls) == len(outFiles))

    nThreads = 4
    # with Pool(processes=nProcs) as pool:

    nChunks = len(urls) // nThreads

    for i in range(nChunks):

        urlsChunk = urls[i * nThreads : (i + 1) * nThreads]
        outFilesChunk = outFiles[i * nThreads : (i + 1) * nThreads]

        threads = []
        
        for url, fOut in zip(urlsChunk, outFilesChunk):
            threads.append(
                Thread(target=download_file, 
                       args=(url, fOut,), 
                       kwargs={'CHUNK': dchunk}))
            threads[-1].start()
            
        for t in threads:
            t.join(timeout=300)

        print(80*'-')
        print(
            f"Parallel download (nThreads={nThreads}) =" \
            + f"{i/nChunks*100:.1f}% ({i * nThreads} / {len(urls)} files)",
            end='\r')

    print("Download complete.")
    print(80*'-')
        

            
        # while True:
        #     incompleteCount = sum(1 for p in processes if not p.ready())
        #     perComplete = float(len(urls) - incompleteCount) / len(urls)
        #     if incompleteCount == 0:
        #         print('{:6.1f}% |'.format(100) + '#' * 50 + '|')
        #         break
        #     else:
        #         print('{:6.1f}% |'.format(perComplete * 100) +
        #               '#' * int(perComplete * 50) +
        #               '.' * (50 - int(perComplete * 50)) +
        #               '|', end='\r')

        #     time.sleep(1)