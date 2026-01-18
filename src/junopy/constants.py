# Copyright (c) 2018 Yash Sarkango
# Licensed under the MIT License

R_JUPITER = 71492  # km
R_GANYMEDE = 2634  # km

PDS_URL = "https://search-pdsppi.igpp.ucla.edu/ditdos/download?id=pds://PPI"

DATA_IDS = {
    "mag": "JNO-J-3-FGM-CAL-V1.0",
    "waves-burst": "JNO-E_J_SS-WAV-3-CDR-BSTFULL-V1.0",
    "waves-survey": "JNO-E_J_SS-WAV-3-CDR-SRVFULL-V2.0",
    "jedi": "JNO-J-JED-3-CDR-V1.0",
    "jade": "JNO-J_SW-JAD-3-CALIBRATED-V1.0",
}

EXTENSIONS = {
    "mag": ".sts",
    "waves-burst": ".DAT",
    "waves-survey": ".CSV",
    "jedi": ".TAB",
    "jade": ".DAT",
}

TIME_FORMATS = {
    "mag": "%Y-%jT%H:%M:%S.%f",
    "waves-burst": "%Y-%m-%dT%H:%M:%S.%f",
    "waves-survey": "ISO8601",
    "jedi": "%Y-%jT%H:%M:%S.%f",
    "jade": "%Y-%jT%H:%M:%S.%f",
}

# These are strings that are used to filter and identify specific data 
# in the PDS INDEX.TAB files.
SID = {
   'mag_1s': "SS 1 SECOND",      # For magnetometer data, junopy defaults to 
   'mag_64s': "SS 1 MINUTE",     # the Juno_JSS coordinate system. 
   'mag_native': "SUNSTATE",
   'waves_survey': "SURVEY",     # Identified both E and B spectra in survey mode.
   'jedi_electrons': "ERSESP",   # prefix 'HI' for high-res, 'LO' for low-res, or leave as is for both.
   'jedi_protons': "ERSTOFXER",  # prefix 'HI' for high-res, 'LO' for low-res, or leave as is for both.
   'jedi_heavies': "NONPTOFXPHR",
   'jedi_proton_pulseheight': "ERSTOFXPHR",
   'jedi_heavies_pulseheight': "NONPTOFXER",   # Verify this!! 
   'jade_electrons': "_ELC_",    # Gets both low-res and high-res electron data.
   'jade_species': 'ION_ANY',    # JADE ion species data.
   'jade_tof': 'ION_TOF',        # JADE ion time-of-flight data. 
} 