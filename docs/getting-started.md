# Getting Started

## Dependencies

junopy uses several external open-source Python packages that are listed in the
`pyproject.toml` file under `dependencies`. If installing using pip, these
dependencies will be installed automatically, or you can choose to install them
individually.

## Installation

I use [`uv`](https://docs.astral.sh/uv/) for installing `junopy` and it is very 
easy to install this way. `uv` will create a virtual environment for this 
project, so you needn't worry about corrupting your global environment. (You
can also explicitly create a virtual environment using `uv venv`)


If you do not forsee making changes to the `junopy` source code, clone the git repository, `cd` into
the `junopy` directory and run - 
```
$ uv pip install .
```

If you do forsee making changes to `junopy`, use the `-e` flag. (for editable install)
```
$ uv pip install -e .
```

`uv` will install all necessary dependencies.

## Usage

After installation, it is useful to first create a project directory that will
store all data files and scripts. I find it useful to separate junopy's core
from the day-to-day analysis code, which I place in a separate directory,
e.g. `project`. I also recommend creating separate directories in the project
directory to store the `data` and `spice` kernels. All this can be done either
manually or by using the Makefile that is provided in the top-level directory.

    $ make project

This will create a new directory called `project` which you can move and
rename as necessary.

Using the Makefile is not necessary, you can also `cd` into a directory of
your choosing and run the following commands. The important thing is that there
needs to be a `data` and `spice` directory to store data.

    $ mkdir -p project/data
    $ mkdir -p project/spice/kernels

### 1. Downloading SPICE kernels

First it is a good idea to download the SPICE kernels relevant to Juno. This
is achieved by running the `download_spice_kernels` method defined in the
`spice` module. E.g. save the following in a new script `download_data.py`.

```python
# Contents of project/download_data.py
from junopy import spice
spice.download_spice_kernels(dodownload=True)
```

Run this script within the `project` folder.
This method uses the `urllib.request` module to open the SPICE NAIF and access
the data via `https`. It downloads only the latest kernels from the `fk`, `pck`,
`lsk`, and `sclk` directories, and the *reconstructed* `spk` trajetory kernels
between years 2016-2025. It creates a list of these kernels and downloads them
to the `spice/kernels` directory in the project folder. It will also create a
meta kernel file called `JunoMetaKernel.txt` in the `spice` directory
with a list of these kernels.

    $ cd project
    $ python download_data.py

And you should get a directory structure like so -

    project
    ├── data
    ├── download_data.py
    └── spice
        ├── JunoMetaKernel.txt
        └── kernels


### 2. Downloading Data

junopy downloads data from the [PDS](https://pds-ppi.igpp.ucla.edu/index.jsp)
by examining the `INDEX.TAB` files that are associated with each dataset. This
index file has a list of all data files organized by time and data ID. Before
downloading any data files for the instruments, we need to download and parse
the corresponding index files. This is achieved by instantiating an `IndexFile`
object. E.g.,

```python
import junopy
from junopy import index

magindex = index.IndexFile(dataset='mag', dodownload=True)
```

The `magindex` object contains the information stored in the magnetometer's
`INDEX.TAB` file. It stored the index file as a .csv in the `index` folder.

    project
    ├── data
    ├── download_data.py
    ├── index
    │   ├── index_mag.csv
    └── spice
        ├── JunoMetaKernel.txt
        └── kernels

The `magindex` object contains a `pd.DataFrame` that has all the tabular data
from `INDEX.TAB`. Now, we can call a method to download data between two dates
stored as `datetime.datetime`.

```python
import datetime as dt

start = dt.datetime(2017, 1, 1)
stop = dt.datetime(2017, 1, 5)

magindex.get_files_between_dates(
    start, stop, DoDownload=True, sid="SS 1 SECOND")
```

Here "SS 1 SECOND" is an identifier for magnetometer data files with 1s cadence
in the JSS coordinate system. Other instruments may have different kinds of
data with their own identifiers. By default, the data files are stored in the
`data/<instrument>/` folder, but you can also specify a folder using the
`outfolder=` keyword argument.

>**Note**: I usually download the entire chunk of magnetometer data and save it on
my local hard drive, but not for the plasma instruments as their
data files are much larger and can easily fill all available disk space. There
are some methods provided in the `jade.py` and `jedi.py` modules that can fetch
relevant data for you on demand, or you can use the above logic to implement
your own workflow.

### 3. Reading Data

Once the files have been downloaded and are present in the appropriate
directory, they can be read using the `<Instrument>Data` classes defined in each
module. This is explained in more detail in the **Usage** tab for each instrument.

In general, data files can be read by instantiating a class object for the
instrument data type. For e.g., for the magnetometer, you can read a file as

```python
from junopy import mag

magfile = 'data/mag/fgm_jno_l3_2017001ss_r1s_v01.sts'
magdata = MagData(filename=magfile, is_hdf=False)
```

`magdata` is an object that has metadata and data associated with that file.

Now, you should be able to plot the data by examining the `magdata.df` DataFrame.

```python
import matplotlib.pyplot as plt

magdata.df.plot(y=['Bx', 'By', 'Bz'])
```

In many cases, you may not need to specify the file.  Some convenience methods
are provided to gather relevant data on demand.  For example, `get_mag_data`
will return a `MagData` instance over a time period specified by the user.

```python
import datetime as dt
from junopy import mag

start = dt.datetime(2017, 1, 1)
stop = dt.datetime(2017, 1, 3)
sid = "SS_1_MINUTE"

magdata = mag.get_mag_data(start, stop, do_download=True, sid=sid)
    # magdata will contain all magnetometer data in JSS coordinate system
    # between start and stop times.
```

Using the convenience methods is often the most effective way to identify and
download data between particular times.  Notice that using `get_mag_data` does
not require an explicit interaction with the `index` module; that is done in
the background.
