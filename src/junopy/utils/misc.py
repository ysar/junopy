import numpy as np


## Charge exchange cross-sections
def cxcs_Hplus_O(energy: float) -> float:
    """
    H+ - O charge-exchange cross-section from Lindsay and Stebbings
    ([2005](https://doi.org/10.1029/2005JA011298)).

    Args:
        energy (float): Energy of the impinging ion in keV

    Returns:
        float: Cross-section of the reaction in cm^2
    """
    a1 = 2.91
    a2 = 0.0886
    a3 = 50.9
    a4 = 4.73
    a5 = -0.862
    a6 = 0.0306
    return (a1 - a2 * np.log(energy)) ** 2 * (1 - np.exp(-a3 / energy)) ** 2 + (
        a4 - a5 * np.log(energy)
    ) * (1 - np.exp(-a6 / energy)) ** 2


def cxcs_Oplus_O(energy: float) -> float:
    """
    O+ - O charge-exchange cross-section from Lindsay and Stebbings
    ([2005](https://doi.org/10.1029/2005JA011298)).

    Args:
        energy (float): Energy of the impinging ion in keV

    Returns:
        float: Cross-section of the reaction in cm^2
    """
    a1 = 4.07
    a2 = 0.269
    a3 = 4125
    return (a1 - a2 * np.log(energy)) ** 2 * (1 - np.exp(-a3 / energy)) ** 0.8


def cxcs_Hplus_H2(energy: float) -> float:
    """
    H+ - H_2 charge-exchange cross-section from Tabata and Shirai
    ([2000](https://doi.org/10.1006/adnd.2000.0835)).

    Args:
        energy (float): Energy of the impinging ion in keV

    Returns:
        float: Cross-section of the reaction in cm^2
    """
    a1 = 2.12e2
    a2 = 1.721
    a3 = 6.7e-4
    a4 = 3.239e-1
    a5 = 4.34e-3
    a6 = 1.296
    a7 = 1.42e-1
    a8 = 9.34
    a9 = 2.997

    rydberg = 1.361e-2  # keV
    threshold_energy = 2.5e-3  # keV

    e1 = energy - threshold_energy

    f1 = a1 * (e1 / rydberg) ** a2
    f2 = f1 / (1 + (e1 / a3) ** (a2 + a4))

    f1_2 = a5 * (e1 / rydberg) ** a2
    f3_2 = f1_2 / (1 + (e1 / a6) ** (a2 + a7) + (e1 / a8) ** (a2 + a9))

    return f2 + f3_2


def cxcs_Hplus_Oplus(energy: float) -> float:
    """
    O+ - O charge-exchange cross-section from Lindsay and Stebbings
    ([2005](https://doi.org/10.1029/2005JA011298)).

    Args:
        energy (float): Energy of the impinging ion in keV

    Returns:
        float: Cross-section of the reaction in cm^2
    """
    return 6e-18 / 1e-16


def cxcs_Hplus_Oplusplus(energy: float) -> float:
    """
    H+ - O++ charge-exchange cross-section

    Args:
        energy (float): Energy of the impinging ion in keV

    Returns:
        float: Cross-section of the reaction in cm^2
    """
    return 1e-18 / 1e-16


def cxcs_Hplus_Splus(energy: float) -> float:
    """
    H+ - S+ cross-section

    Args:
        energy (float): Energy of the impinging ion in keV

    Returns:
        float: Cross-section of the reaction in cm^2
    """
    return 4e-18 / 1e-16


def cxcs_Hplus_Splusplus(energy: float) -> float:
    """
    H+ - S++ cross-section

    Args:
        energy (float): Energy of the impinging ion in keV

    Returns:
        float: Cross-section of the reaction in cm^2
    """
    return 6e-18 / 1e-16


def get_percentile(arr: np.ndarray, p: float):
    """
    Returns the log-percentile of an input array after discounting the zero
    values and NaNs.

    Args:
        arr (np.ndarray): The input array of any shape.
        p (float): The percentile between 0-100.

    Returns:
        result (float): The log-percentile of the input array.
    """
    return np.exp(np.percentile(np.log(arr[arr > 0 & ~np.isnan(arr)]), p))


def bagenal_delamere_density(radius: np.ndarray):
    """
    Get the density using equation from Bagenal and Delamere (2011).

    Args:
        radius (np.ndarray): Numpy array with different radii.

    Returns:
        result (np.ndarray): Number densities in cm^-3 .
    """
    a1 = 1987.0
    b1 = -8.2
    a2 = 14.0
    b2 = -3.2
    a3 = 0.05
    b3 = -0.65

    return (
        a1 * (radius / 6) ** (-b1)
        + a2 * (radius / 6) ** (-b2)
        + a3 * (radius / 6) ** (-b3)
    )
