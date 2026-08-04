# junopy

`junopy` is an unofficial Python package for analyzing plasma and fields data
collected by the [Juno](https://science.nasa.gov/mission/juno/) 
spacecraft.  It uses data from the Planetary
Plasma Interactions ([PPI](https://pds-ppi.igpp.ucla.edu/)) 
node of the Planetary Data System ([PDS](https://pds.nasa.gov/)).  

Juno data was originally stored on the PDS-PPI node in a PDS3-compliant archive.  
However, PDS has been encouraging space missions to upgrade to the newer
PDS4 archive.  Despite this directive, for now `junopy` continues to
use the PDS3 structure.  Practically all functionality related to searching for
and downloading data relies on `INDEX.TAB` files that were present in every
PDS3 archive.  Eventually, if there is interest from the community and if the PDS 
or Juno team converts all data to the PDS4 format, and deletes all PDS3 data, 
this tool will need to be updated.

`junopy` is an unofficial tool for personal use.  There may be bugs and issues
present.  It is open-source and free to use and modify under the terms of the 
MIT license.

Although unaffiliated, this tool has benefited greatly from 
conversations with many people on the Juno mission and instrument teams.  Thank
you all and hope this is of some use!

## Contributors:
Yash Sarkango ([sarkango@umn.edu](mailto:sarkango@umn.edu))

## FAQs

1. **Why only particle and fields instruments?**  
I mainly worked with the magnetometer, Waves, JEDI, and JADE instruments so those 
are included.  It is possible to add functionality to analyze data from other 
instruments (and contributions are welcome).  It will take some effort.

2. **Have you published this? How can I cite this tool?**  
There is no scientific publication that describes `junopy` but we might use Zenodo
in case we make it public.

3. **Was this funded by NASA? Is this sponsored by the Juno mission?**  
No this was not directly sponsored or funded by NASA.  It is not an official tool 
used by the mission team or the instrument teams.

4. **Can I use this tool for my papers and proposals?**  
If you have access to this repository, then yes! Please don't share it with anyone 
yet though. We plan to make this tool public eventually. 

5. **What operating systems are supported?**  
It should work well on MacOS/Linux and Windows via WSL.  In theory it should also 
work natively on Windows.

6. **What library are you using for the magnetic field calculations?**  
[`iupitermag`](https://github.com/ysar/iupitermag) 
in conjunction with `junopy` to do all sorts of calculations.

7. **How does this compare to MIDL or Autoplot?**  
(To my knowledge) It is much easier to use MIDL or Autoplot if you want to quickly
visualize data from a particular instrument.  `junopy` has some rudimentary plotting
methods but its main purpose is not to plot the data, rather it provides access
to the underlying data.  There is a lot you can do with `junopy`, however (unless
you already have it up and running) it will take a bit more effort to do it.  One
use case where `junopy` really shines is when you need to combine data from 
different instruments.
