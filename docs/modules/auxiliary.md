# Auxiliary Modules

## `index`

The `index` module provides an interface to the index (INDEX.TAB) files stored 
on the PDS for every PDS3 compliant package. The index file for each repository
is a table with the metadata for various data files that are stored in that
repository. The `index` module makes use of this generalization to identify 
relevant files between different datasets. Each repository or dataset may have 
different directory structures, but the relative path to each file, as well as
the measurement start and stop times, are stored in the index file for that 
dataset.

### Usage

The `index` module contains a class for reading and operating with INDEX.TAB 
files. First, you need to instantiate an object for the `IndexFile` class.

    from JunoPy import index 

    magIndex = index.IndexFile(dataset='mag', doDownload=True)

There are a few methods implemented for the `IndexFile` class, but perhaps the 
most useful is `get_files_between_dates`, which will download the relevant 
datafiles from the repository. You can call this method after instantiation of
the object.

    filelist = magIndex.get_files_between_dates(
        start,
        stop,
        sid='SS 1 SECOND'
        doDownload=True,
        doMultiprocessing=True,
        outfolder='data/mag',
        getLBL=True
    )

    Where `start` and `stop` are expected of type `dt.Datetime`. `sid` is a 
    unique identifier that is defined for a particular data type, e.g. for the 
    magnetometer `"SS_1_SECOND"` is an identifier for 1 s cadence magnetometer
    data in the Jovian Sun-State (JSS) coordinate system. Different datasets
    and repositories use different identifiers to separate different kinds of 
    data. A table of different `sid` values for each data type is given along 
    with 