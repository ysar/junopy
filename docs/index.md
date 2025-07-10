# junopy

`junopy` is an unofficial Python package for analyzing plasma and fields data
collected by the Juno spacecraft.  It uses data from the Planetary
Plasma Interactions (PPI) node of the Planetary Data System (PDS).  

I started work on `junopy` around 2018, when Juno data was relatively new.  
At the time, Juno data was stored on the PDS-PPI node in a PDS3-compliant format.  
Since then, the PDS has been encouraging space missions to upgrade to the newer
PDS4 format.  Despite this directive, for now `junopy` continues to
use the PDS3 structure.  Practically all functionality related to searching for
and downloading data relies on `INDEX.TAB` files that were present in every
PDS3 archive.

Eventually, if there is interest from the community and if the PDS or Juno team
converts all data to the PDS4 format, this tool may need to be upgraded.

`junopy` is an unofficial tool for personal use.  There may be bugs and issues
present.  It is open-source and free to use and modify under the terms of the 
GNU-GPLv3 license.

Author: Yash Sarkango
