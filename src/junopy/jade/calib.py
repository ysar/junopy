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


class ElectronSensorConstants:
    def __init__(
        self,
        max_epsilon: float = np.nan,
        k: float = np.nan,
        max_energy: float = np.nan,
        rho: float = np.nan,
        n: float = np.nan,
        tau0_t0: float = np.nan,
        tau0_t1: float = np.nan,
        tau0_t2: float = np.nan,
        tau1_t1: float = np.nan,
        tau1_t2: float = np.nan,
        c0: float = np.nan,
        c1: float = np.nan,
        c2: float = np.nan,
        b0: float = np.nan,
        b1: float = np.nan,
        f: float = np.nan,
    ):
        """
        Initialize the constants for a JADE electron sensor.  These constants are
        related to the sensitivity and calibration of the sensor and can be
        used to calculate the effective geometric factor.

        Attributes:
            max_epsilon (float): Maximum epsilon value for the sensor.
            k (float): Sensitivity factor.
            max_energy (float): Maximum energy for the sensor.
            rho (float): Density factor.
            n (float): Exponent for energy dependence.
            tau0_t0 (float): Time correction for reference time 1.
            tau0_t1 (float): Time correction for reference time 2.
            tau0_t2 (float): Time correction for reference time 3.
            tau1_t1 (float): Time correction slope for reference time 1.
            tau1_t2 (float): Time correction slope for reference time 2.
            c0 (float): Coefficient for deflection correction.
            c1 (float): Coefficient for deflection correction.
            c2 (float): Coefficient for deflection correction.
            b0 (float): Coefficient for anode correction.
            b1 (float): Coefficient for anode correction.
            f (float): Anode factor.
        """
        self.max_epsilon = max_epsilon
        self.k = k
        self.max_energy = max_energy
        self.rho = rho
        self.n = n
        self.tau0_t0 = tau0_t0
        self.tau0_t1 = tau0_t1
        self.tau0_t2 = tau0_t2
        self.tau1_t1 = tau1_t1
        self.tau1_t2 = tau1_t2
        self.c0 = c0
        self.c1 = c1
        self.c2 = c2
        self.b0 = b0
        self.b1 = b1
        self.f = f


def _get_sensor_constants(sensor_id) -> ElectronSensorConstants:
    if sensor_id == 104:  # E-060
        sensor_const = ElectronSensorConstants(
            max_epsilon=0.96,
            k=3.58,
            max_energy=0.98e3,
            rho=1.58,
            n=1.67,
            tau0_t0=1.04,
            tau0_t1=1.04,
            tau0_t2=1.23,
            tau1_t1=-0.062,
            tau1_t2=-0.076,
            c0=3.12e-5,
            c1=-5.3e-8,
            c2=-6.92e-9,
            b0=1.088,
            b1=-9.766e-3,
            f=0.745,
        )

    elif sensor_id == 106:  # E-180
        sensor_const = ElectronSensorConstants(
            max_epsilon=1.22,
            k=2.78,
            max_energy=0.91e3,
            rho=2.44,
            n=1.67,
            tau0_t0=1.16,
            tau0_t1=1.16,
            tau0_t2=1.20,
            tau1_t1=-0.079,
            tau1_t2=-0.074,
            c0=3.12e-5,
            c1=-5.3e-8,
            c2=-6.92e-9,
            b0=1.048,
            b1=-5.372e-3,
            f=0.745,
        )

    return sensor_const


def _calc_time_correction_tau(t, sensor_const):
    c = sensor_const
    ref_time1 = dt.datetime(2016, 7, 4)
    ref_time2 = dt.datetime(2019, 5, 22)
    deltat1 = (t - ref_time1).total_seconds() / (3600 * 24 * 365)
    deltat2 = (t - ref_time2).total_seconds() / (3600 * 24 * 365)

    if deltat1 < 0:
        return c.tau0_t0
    elif deltat2 < 0:
        return c.tau0_t1 + c.tau1_t1 * deltat1
    else:
        return c.tau0_t2 + c.tau1_t2 * deltat1


def _calc_deflection_correction_D(delta, sensor_const):
    # delta is 'DIM2_ELEVATION'
    d = delta * np.pi / 180
    c = sensor_const
    return 1.0 + c.c1 / c.c0 * d + c.c2 / c.c0 * d**2


def _calc_anode_correction_F(sensor_const):
    c = sensor_const
    F = np.zeros(16)
    F[0] = (c.b0 + c.b1 * 16) * c.f
    F[-1] = (c.b0 + c.b1 * 1) * c.f
    for i in range(1, 15):
        F[i] = c.b0 + c.b1 * (16 - i)

    return F


def _calc_simulation_geom_factor(energies):
    return 0.8 * 0.9 * (2.37e-4 + 9.17e-11 * energies)


def _calc_sensor_dependence_exp(sensor_const, energies):
    c = sensor_const
    inner_exp = 1 - np.exp(-c.rho * (energies / c.max_energy) ** c.n)
    return c.max_epsilon * (
        1 - np.exp(-c.k * (energies / c.max_energy) ** (1 - c.n) * inner_exp)
    )
