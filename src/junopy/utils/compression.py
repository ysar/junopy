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

import os
import shutil
from isal import igzip


def create_gz_zip_file(filename: str) -> None:
    """
    Compresses a file using igzip into a .gz file.

    Args:
        filename (str):  Name of the file to be compressed.  This method will
            create a new file with the same name with the '.gz' extension.  The
            original file will be removed.
    """
    with open(filename, 'rb') as f_in:
        with igzip.open(filename + '.gz', 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
    os.remove(filename)
