# junopy

`junopy` is a Python package for analyzing plasma and fields data collected by the [Juno](https://science.nasa.gov/mission/juno/) 
spacecraft.  It uses data from the Planetary Plasma Interactions ([PPI](https://pds-ppi.igpp.ucla.edu/)) 
node of the Planetary Data System ([PDS](https://pds.nasa.gov/)).  

`junopy` assumes that the PDS data is stored using the PDS3 structure. 
Practically all functionality related to searching for and downloading data relies on `INDEX.TAB` 
files that were present in every PDS3 archive.  Eventually, if there is interest from the community 
and if the PDS or Juno team converts all data to the PDS4 format, and deletes all PDS3 data/metadata, 
this tool will need to be updated.

`junopy` is an unofficial tool for personal use.  There may be bugs and issues
present.  It is open-source and free to use and modify under the terms of the 
MIT license. Please consider citing it.

## Contributors:
Yash Sarkango ([sarkango@umn.edu](mailto:sarkango@umn.edu))

## FAQs

1. **Why only particle and fields instruments?**  
I have mainly worked with the magnetometer, Waves, JEDI, and JADE instruments so those 
are included.  It is possible to add functionality to analyze data from other 
instruments (and contributions are welcome).

2. **How can I cite this tool?**  
Please use the following -  
Sarkango, Y. (2026). junopy: A Python package for processing particles and fields data from the Juno 
spacecraft [Computer software]. Zenodo. [https://doi.org/10.5281/zenodo.23898717](https://doi.org/10.5281/zenodo.23898717)

3. **Was this funded by NASA? Is this sponsored by the Juno mission?**  
No this was not directly sponsored or funded by NASA.  It is not an official tool 
used by the mission team or the instrument teams.

4. **What operating systems are supported?**  
It should work well on MacOS/Linux and Windows via WSL.  In theory it should also 
work natively on Windows.

5. **What library are you using for the magnetic field calculations?**  
[`iupitermag`](https://github.com/ysar/iupitermag) 

6. **How does this compare to MIDL or Autoplot?**  
(To my knowledge) It is much easier to use MIDL or Autoplot if you want to quickly
visualize data from a particular instrument.  `junopy` has some rudimentary plotting
methods but its main purpose is not to plot the data, rather it provides access
to the underlying data.  There is a lot you can do with `junopy`, however (unless
you already have it up and running) it will take a bit more effort to do it.  One
use case where `junopy` really shines is when you need to combine data from 
different instruments.
