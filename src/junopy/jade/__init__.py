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

# from .__download import JADE_SID_TYPES
# from .__download import jade_download_data
# from .get import *
# from .read import * # read_jade_binary_file, create_dict_from_jade_df
# from .sensor import *
# from .plot import *
# from .calib import *


# from .get import get_jade_data_files
# from .read import read_jade_binary_file, JadeData
# from .calib import get_Geff_electrons

# Let's think a bit about how the JADE processing can be improved. Firstly,
# after reading in the raw data, we should separate it by sensor.
# JadeData
#   |--- JadeSensorData
#   |--- JadeSensorData
#   |--- JadeSensorData
# Another question to ask is if it is necessary to separate high-res vs low-res
# data?

# Jade Processing workflow looks somewhat like this
# 1. We read in the binary file into DataFrames
# 2. We truncate the DataFrame based on what is needed
# 3. We convert the DataFrame into Dictionary with stacked arrays
# 4. We separate the array data into different sensors.
# 5.
