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

import logging
from typing import Dict
from pathlib import Path

class ConfigVars:
    """
    Centralized configuration management for directory paths.

    Maintains a dictionary of configuration variables and automatically
    manages dependent directory paths when the root directory changes.
    """

    __config: Dict[str, Path] = {"ROOT_DIR": Path("./").absolute()}

    # Directory structure template
    __DIR_STRUCTURE = {
        "DATA_DIR": Path("data"),
        "SPICE_DIR": Path("spice"),
        "INDEX_DIR": Path("index"),
        "TMP_DIR": Path("_tmp"),
        # Subdirectories of DATA_DIR
        "MAG_DATA_DIR": Path("data") / "mag",
        "JADE_DATA_DIR": Path("data") / "jade",
        "JEDI_DATA_DIR": Path("data") / "jedi",
        "WAVES_SURVEY_DATA_DIR": Path("data") / "waves_survey",
    }

    @classmethod
    def get(cls, name: str) -> Path:
        """Get a configuration value by name.

        Args:
            name: Name of the configuration variable

        Returns:
            The configuration value

        Raises:
            KeyError: If the configuration variable doesn't exist
        """
        try:
            return cls.__config[name]
        except KeyError:
            logging.error(f"Configuration variable '{name}' not found")
            raise

    @classmethod
    def set_root(cls, value: Path, create_dirs: bool = True) -> None:
        """
        Set the root directory and update all dependent paths.

        Args:
            value: Path to the new root directory
            create_dirs: Whether to create directories automatically
        """
        # Convert to absolute path and normalize
        abs_path = value.absolute()
        cls.__config["ROOT_DIR"] = abs_path
        logging.info(f"Root directory changed to {cls.__config['ROOT_DIR']}")

        cls.update_dependent_vars(create_dirs=create_dirs)

    @classmethod
    def update_dependent_vars(cls, create_dirs: bool = False) -> None:
        """
        Update all directory paths based on the current root directory.

        Args:
            create_dirs (bool): Whether to create a new directory.
        """
        root_dir = cls.__config["ROOT_DIR"]

        # Update all directory paths
        for var_name, rel_path in cls.__DIR_STRUCTURE.items():
            cls.__config[var_name] = root_dir / rel_path
            if create_dirs and "DIR" in var_name:
                cls.__config[var_name].mkdir(parents=True, exist_ok=True)
                logging.info(f"Created directory: {cls.__config[var_name]}")

        # Log important directory changes
        for dir_type in ["DATA", "SPICE", "INDEX"]:
            logging.info(f"{dir_type} directory: {cls.get(f'{dir_type}_DIR')}")

    @classmethod
    def get_all(cls) -> Dict[str, Path]:
        """Get a copy of the complete configuration dictionary."""
        return cls.__config.copy()

# Initialize the configuration
ConfigVars.update_dependent_vars()
