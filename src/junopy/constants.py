# junopy - A Python API for the Juno spacecraft's particle and fields datasets.
# Copyright (C) 2018  Yash Sarkango

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

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
