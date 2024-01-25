# Magnetometer Data

## Additional Info

### Magnetometer Data SIDs

Thew following table contains the list of SIDs to be used with the magnetometer
index files.

| SID | Cadence (s) | Coordinate System | 
| --- | ----------- | ----------------- |
| SS_1_SECOND | 1   | Jupiter Sun-State (JSS) |  
| SS_1_min    | 60  | Jupiter Sun-State (JSS) | 
| SUNSTATE    | 1/64 (or native) | Jupiter Sun-State (JSS) | 
| PC_1_SECOND | 1   | Planetocentric |  
| PC_1_min    | 60  | Planetocentric | 
| SUNSTATE    | 1/64 (or native) | Planetocentric | 
| PL_1_SECOND | 1   | Payload |  
| PL_1_min    | 60  | Payload | 
| SUNSTATE    | 1/64 (or native) | Payload | 

For most applications, you should only download the Sun-State data files as you
can always convert the magnetometer data while reading the data files using the 
`coordsystem=` keyword argument (assuming you have the relevant frame 
SPICE kernels.)

### Supported Coordinate Systems

The following coordinate systems are supported by JunoPy for the 

