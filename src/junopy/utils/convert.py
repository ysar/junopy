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
"""
Some functions for speed-energy conversions, coordinate conversions, etc.
"""
import numpy as np
from scipy.constants import speed_of_light, m_e


def calc_relativistic_momentum(energy: float, mass: float = m_e) -> float:
    """
    Calculate relativistic momentum from kinetic energy. [SI units]

    Args:
        energy (float): The energy of the particle in J
        mass (float): The mass of the particle in kg

    Returns:
        momentum (float): The relativistic momentum of the particle in SI units.
    """
    return (
        1.0
        / speed_of_light
        * np.sqrt(
            (energy + mass * speed_of_light**2) ** 2 - mass**2 * speed_of_light**4
        )
    )


def calc_speed(energy: float, mass: float = m_e) -> float:
    """
    Calculate non-relativistic speed from kinetic energy. [SI units]

    Args:
        energy (float): The energy of the particle in J
        mass (float): The mass of the particle in kg

    Returns:
        speed (float): The non-relativistic speed of the particle (m/s)
    """
    return np.sqrt(2 * energy / mass)


def calc_relativistic_speed(energy: float, mass: float = m_e) -> float:
    """
    Calculate relativistic speed from kinetic energy. [SI units]

    Args:
        energy (float): The energy of the particle in J
        mass (float): The mass of the particle in kg

    Returns:
        speed (float): The relativistic speed of the particle (m/s)

    """
    return speed_of_light * np.sqrt(
        1.0 - 1.0 / (energy / mass / speed_of_light**2 + 1.0) ** 2
    )


def calc_relativistic_energy(speed: float, mass: float = m_e) -> float:
    """
    Calculate relativistic kinetic energy. [SI units]

    Args:
        speed (float): The speed of the particle (m/s)
        mass (float): The mass of the particle in kg

    Returns:
        energy (float): The kinetic energy of the particle in J, calculated
            using the relativistic equation.

    """
    p = mass * speed / np.sqrt(1 - (speed / speed_of_light) ** 2)

    return (
        np.sqrt((p * speed_of_light) ** 2 + mass**2 * speed_of_light**4)
        - mass * speed_of_light**2
    )


def convert_pos_spherical(pos: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Given cartesian coordinates, return the radius, local-time, and latitude.

    Note - This returns local-time, which is not the same as longitude.

    Args:
        pos (np.ndarray): Array of input positions with shape (3, N).

    Returns:
        radius (np.ndarray): Array of radii with shape (N,)
        local_time (np.ndarray): Local times with shape (N,)
        latitude (np.ndarray): Latitude with shape (N,)
    """
    radius = np.sqrt(pos[0] ** 2 + pos[1] ** 2 + pos[2] ** 2)
    local_time = 6 / 90 * np.arctan2(pos[1], pos[0]) * 180 / np.pi + 12
    latitude = np.arctan2(pos[2], np.sqrt(pos[0] ** 2 + pos[1] ** 2)) * 180 / np.pi

    return radius, local_time, latitude
