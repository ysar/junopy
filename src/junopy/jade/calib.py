"""
This submodule contains code to calculate the geometric factor for the different
JADE sensors, following the supporting information of Allegrini et al.
([2020](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2019JA027693)).
The way to calculate the geometric factor varies for high-rate electron
measurements, low-rate electron measurements, and high/low-rate ion measurements.

The general equation for the electron sensors is -

```
G_eff (delta, E, Anode, time) = G_simulations \n
    * max_epsilon \n
    * [1 - exp(-k * Eratio**(1-n) * (1 - exp(-rho * Eratio**n))] \n
    * tau(time) \n
    * D(delta) \n
    * F(anode)
```

Where `Eratio = E / max_energy`

Here `(max_epsilon, k, max_energy, rho, n)` depend on the JADE sensor. `tau` is
a function of time/date. `D(delta)` depends on deflection angle and sensor.

"""

import datetime as dt

import numpy as np

TIME_DEPENDENCE_EPOCH1 = dt.datetime.strptime("2016-186T00:00", "%Y-%jT%H:%M")
TIME_DEPENDENCE_EPOCH2 = dt.datetime.strptime("2019-142T00:00", "%Y-%jT%H:%M")
TIME_DEPENDENCE_EPOCH3 = dt.datetime.strptime("2022-180T18:35", "%Y-%jT%H:%M")
TIME_DEPENDENCE_EPOCH4 = dt.datetime.strptime("2023-062T00:00", "%Y-%jT%H:%M")
TIME_DEPENDENCE_EPOCH5 = dt.datetime.strptime("2025-090T19:42", "%Y-%jT%H:%M")

MCP_EFFICIENCY_EPOCH1 = TIME_DEPENDENCE_EPOCH2
MCP_EFFICIENCY_EPOCH2 = dt.datetime.strptime("2022-180T18:15", "%Y-%jT%H:%M")
MCP_EFFICIENCY_EPOCH3 = dt.datetime.strptime("2025-090T19:42", "%Y-%jT%H:%M")


def _calc_time_dependent_sensitivity_e060(time) -> np.ndarray:
    """
    The sensitivity of E060 is assumed to be time-independent.
    """
    # fmt: off
    f0 = np.array([0.718, 1.069, 1.111, 1.085, 1.076, 1.069, 1.027, 1.074, 1.058, 1.053, 1.022, 1.007, 1.023, 1.002, 0.971, 0.636])
    f1 = np.array([3.03, -3.95, -2.73, -0.728, 0.834, 2.09, 2.58, 1.29, 1.84, 1.17, 1.32, -0.698, -2.37, -4.40, -5.10, 5.82]) * 1.0e-2
    f2 = np.array([-7.27, 10.4, 7.45, 2.87, -1.48, -5.20, -7.04, -3.41, -4.21, -2.74, -2.16, 3.02, 7.04, 11.4, 11.5, -20.2]) * 1.0e-3
    f3 = np.array([8.66, -6.41, -5.41, -2.99, -0.29, 2.28, 3.82, 0.79, 1.40, 0.75, 0.39, -3.18, -5.48, -7.75, -5.27, 18.7]) * 1.0e-4
    # fmt: on

    if time < TIME_DEPENDENCE_EPOCH1:
        y = 0.0
    else:
        y = (time - TIME_DEPENDENCE_EPOCH1).total_seconds() / (24 * 3600 * 365)

    return f0 + f1 * y + f2 * y**2 + f3 * y**3


def _calc_time_dependent_sensitivity_e180(time) -> np.ndarray:
    """
    The sensitivity of E180 is time-independent.
    """
    # fmt: off
    if time < dt.datetime.strptime("2022-181", "%Y-%j"):
        f0 = np.array([0.583, 0.994, 1.052, 1.055, 1.035, 1.023, 1.083, 1.107, 1.112, 1.126, 1.065, 1.072, 1.065, 1.030, 0.930, 0.668])
        f1 = np.array([11.5, -0.7, 0.5, 3.33, 3.74, 4.91, 2.86, 1.13, 0.908, -1.65, 0.146, -5.19, -6.55, -7.99, -6.72, -0.2]) * 1.0e-2
        f2 = np.array([-27.5, 0.95, -1.45, -9.34, -8.56, -12.2, -8.80, -6.39, -6.51, -0.19, -3.71, 15.0, 18.9, 24.2, 20.0, 5.61]) * 1.0e-3
        f3 = np.array([29.1, 7.94, 7.58, 13.6, 9.47, 11.4, 8.51, 7.10, 7.72, 2.32, 4.09, -17.1, -22.2, -29.5, -27.9, -12.1]) * 1.0e-4
    else:
        f0 = np.array([10.35, 4.395, -0.799, -6.440, -3.735, -1.617, -2.116, -2.526, -2.268, -1.288, 0.436, -1.511, 0.005, 3.169, 9.545, 10.40])
        f1 = np.array([-4.33, -1.61, 0.732, 3.27, 2.12, 1.18, 1.41, 1.57, 1.47, 1.03, 0.276, 1.15, 0.515, -0.890, -3.61, -4.27])
        f2 = np.array([6.43, 2.52, -0.92, -4.69, -3.07, -1.73, -2.05, -2.28, -2.14, -1.50, -0.42, -1.72, -0.84, 1.22, 5.03, 6.17]) * 1.0e-1
        f3 = np.array([-31.1, -12.9, 3.79, 22.3, 14.7, 8.28, 9.83, 10.9, 10.3, 7.17, 2.10, 8.42, 4.39, -5.71, -23.3, -29.1]) * 1e-3
    # fmt: on

    if time < TIME_DEPENDENCE_EPOCH1:
        y = 0.0
    else:
        y = (time - TIME_DEPENDENCE_EPOCH1).total_seconds() / (24 * 3600 * 365)

    return f0 + f1 * y + f2 * y**2 + f3 * y**3


def _calc_simulation_geom_factor(energies) -> np.ndarray:
    return 0.8 * 0.9 * (2.37e-4 + 9.17e-11 * energies)


def _calc_deflection_correction(delta) -> np.ndarray:
    d = delta * np.pi / 180.0
    return 1.0 + (-5.30e-8 / 3.12e-5) * d + (-6.92e-9 / 3.12e-5) * d**2


def _calc_time_dependence_detector_efficiency_e060(time) -> float:
    t = (time - dt.datetime(2016, 7, 4)).total_seconds() / (24 * 3600 * 365)

    if time < TIME_DEPENDENCE_EPOCH1:
        tau0 = 1.06
        tau1 = 0.0
    elif time < TIME_DEPENDENCE_EPOCH2:
        tau0 = 1.06
        tau1 = -0.064
    elif time < TIME_DEPENDENCE_EPOCH3:
        tau0 = 1.27
        tau1 = -0.086
    elif time < TIME_DEPENDENCE_EPOCH4:
        tau0 = 1.99
        tau1 = -0.164
    elif time < TIME_DEPENDENCE_EPOCH5:
        tau0 = 1.43
        tau1 = -0.0801
    else:
        tau0 = 1.74
        tau1 = -0.14

    return tau0 + tau1 * t


def _calc_time_dependence_detector_efficiency_e180(time) -> float:
    t = (time - dt.datetime(2016, 7, 4)).total_seconds() / (24 * 3600 * 365)

    if time < TIME_DEPENDENCE_EPOCH1:
        tau0 = 1.18
        tau1 = 0.0
    elif time < TIME_DEPENDENCE_EPOCH2:
        tau0 = 1.18
        tau1 = -0.091
    elif time < TIME_DEPENDENCE_EPOCH3:
        tau0 = 1.31
        tau1 = -0.115
    elif time < TIME_DEPENDENCE_EPOCH4:
        tau0 = 2.55
        tau1 = -0.251
    elif time < TIME_DEPENDENCE_EPOCH5:
        tau0 = 1.55
        tau1 = -0.101
    else:
        tau0 = 1.95
        tau1 = -0.101

    return tau0 + tau1 * t


def _calc_mcp_efficiency_e060(time, energies) -> np.ndarray:
    # E060 MCPs assumed to be time-independent.
    energy_max = 0.98e3
    rho = 1.58
    alpha = 1.67
    k = 3.58
    epsilon_max = 0.96

    inner_exp = 1 - np.exp(-rho * (energies / energy_max) ** alpha)
    return epsilon_max * (
        1 - np.exp(-k * (energies / energy_max) ** (1 - alpha) * inner_exp)
    )


def _calc_mcp_efficiency_e180(time, energies) -> np.ndarray:
    if time < MCP_EFFICIENCY_EPOCH1:
        energy_max = 0.45e3
        rho = 0.64
        alpha = 1.67
        k = 4.38
        epsilon_max = 1.19
    elif time < MCP_EFFICIENCY_EPOCH2:
        energy_max = 0.68e3
        rho = 1.09
        alpha = 1.67
        k = 3.06
        epsilon_max = 1.41
    else:
        energy_max = 0.64e3
        rho = 1.17
        alpha = 1.67
        k = 3.48
        epsilon_max = 1.26

    inner_exp = 1 - np.exp(-rho * (energies / energy_max) ** alpha)
    return epsilon_max * (
        1 - np.exp(-k * (energies / energy_max) ** (1 - alpha) * inner_exp)
    )
