# Auxiliary Modules

## Index Files

The `index` module provides an interface to the index (INDEX.TAB) files stored
on the PDS for every PDS3 compliant package. The index file for each repository
is a table with the metadata for various data files that are stored in that
repository. The `index` module makes use of this generalization to identify
relevant files between different datasets. Each repository or dataset may have
different directory structures, but the relative path to each file, as well as
the measurement start and stop times, are stored in the index file for that
dataset.

First, you need to instantiate an object for the `IndexFile` class.

```python
from junopy import index
mag_index = index.IndexFile(dataset='mag', do_download=True)
```

There are a few methods implemented for the `IndexFile` class, but perhaps the
most useful is `get_files_between_dates`, which will download the relevant
datafiles from the repository. You can call this method after instantiation of
the object.

```python
filelist = magIndex.get_files_between_dates(
    start,
    stop,
    sid='SS 1 SECOND'
    do_download=True,
    do_parallel=True,
    outfolder='data/mag',
    get_label=True
)
```

Where `start` and `stop` are expected of type `dt.Datetime`. `sid` is a
unique identifier that is defined for a particular data type, e.g. for the
magnetometer `"SS_1_SECOND"` is an identifier for 1 s cadence magnetometer
data in the Jovian Sun-State (JSS) coordinate system. Different datasets
and repositories use different identifiers to separate different kinds of
data. A table of different `sid` values for each data type is provided in
each of the module sections.

## Downloading files

`junopy` uses the following methods to download files, although in most cases
you will not need to use them directly.  SPICE kernels and PDS3 datasets for
different instruments get be downloaded uses the more appropriate
`download_spice_kernels` method or the instruments' respective `IndexFile`
classes.

```python
from junopy import utils

files_to_download = [url1, url2, url3]
outfolder = "foo/"
outnames = [outfolder + f for f in files_to_download]

utils.download.download_multiple_files(files_to_download, output_files=outnames)
```

`download_multiple_files` uses different Python threads to download and write
multiple files simultaneously to disk.  Please use this with caution to prevent
rate-limit issues with the PDS.
