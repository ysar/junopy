# # junopy - A Python API for the Juno spacecraft's particle and fields datasets.
# # Copyright (C) 2018  Yash Sarkango

# # This program is free software: you can redistribute it and/or modify
# # it under the terms of the GNU General Public License as published by
# # the Free Software Foundation, either version 3 of the License, or
# # (at your option) any later version.

# # This program is distributed in the hope that it will be useful,
# # but WITHOUT ANY WARRANTY; without even the implied warranty of
# # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# # GNU General Public License for more details.

# # You should have received a copy of the GNU General Public License
# # along with this program.  If not, see <https://www.gnu.org/licenses/>.

# import datetime as dt
# import pandas as pd
# import glob

# from JunoPy import config

# PATTERNS = {
#     "ELECTRONS": "S_ELC_",
#     "ION_SPECIES": "_ION_ANY_",
#     "ION_LOGICALS": "_ION_LOG_",
#     "ION_TOF": "_ION_TOF_",
# }


# def get_jade_data_files(start, stop, **kwargs):
#     """
#     Return combined JADE low and high-res data.
#     """
#     getHRS = kwargs.pop("getHRS", True)
#     getLRS = kwargs.pop("getLRS", True)
#     datatype = kwargs.pop("datatype", "ELECTRONS")

#     relevantFiles = []
#     for day in pd.date_range(start, stop):

#         relevantFiles += glob.glob(
#             config.ConfigVars.get("JADE_DATA_DIR")\
#             + day.strftime("/%Y/%j/") \
#             + f"*{PATTERNS[datatype]}*.DAT*")

#         if not getLRS:
#             relevantFiles = [
#                 f for f in relevantFiles if "_LRS" not in f and "_LLS" not in f]

#         if not getHRS:
#             relevantFiles = [
#                 f for f in relevantFiles if "_HRS" not in f and "_HLS" not in f]

#     return relevantFiles
