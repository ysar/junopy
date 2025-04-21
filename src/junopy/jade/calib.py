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
This submodule contains code to calculate the geometric factor for the different
JADE sensors. The way to calculate the geometric factor varies for
    - High Rate electrons Geff(Anode, Energy, time)
    - Low Rate electrons Geff()
    - High / Low Rate Ions

The general equation is -

G(delta, E, Anode, time) = G_simulations \
    * maxEpsilon \
    * [1 - exp(-k * Eratio**(1-n) * (1 - exp(-rho * Eratio**n))] \
    * tau(time) \
    * D(delta) \
    * F(anode)

Eratio = E / maxE

Here (maxEpsilon, k, maxE, rho, n) depend on sensor only.
     (tau) depends only on time of measurement
     D(delta) depends on deflection angle and sensor

"""

import numpy as np
import datetime as dt
import pandas as pd


def get_sensor_constants(sensorId):

    if sensorId == 104:    # E-060
        sensorConst = {
            'maxEpsilon': 0.96,
            'k': 3.58,
            'maxE': 0.98e3,
            'rho': 1.58,
            'n': 1.67,
            'tau0_t0': 1.04,
            'tau0_t1': 1.04,
            'tau0_t2': 1.23,
            'tau1_t1': -0.062,
            'tau1_t2': -0.076,
            'c0': 3.12e-5,
            'c1': -5.3e-8,
            'c2': -6.92e-9,
            'b0': 1.088,
            'b1': -9.766e-3,
            'f': 0.745
        }

    elif sensorId == 106:    # E-180
        sensorConst = {
            'maxEpsilon': 1.22,
            'k': 2.78,
            'maxE': 0.91e3,
            'rho': 2.44,
            'n': 1.67,
            'tau0_t0': 1.16,
            'tau0_t1': 1.16,
            'tau0_t2': 1.20,
            'tau1_t1': -0.079,
            'tau1_t2': -0.074,
            'c0': 3.12e-5,
            'c1': -5.3e-8,
            'c2': -6.92e-9,
            'b0': 1.048,
            'b1': -5.372e-3,
            'f': 0.745
        }

    return sensorConst


def calc_time_correction_tau(t, sensorConst):
    c = sensorConst
    refTime1 = dt.datetime(2016, 7, 4)
    refTime2 = dt.datetime(2019, 5, 22)
    deltat1 = (t - refTime1).total_seconds() / (3600 * 24 * 365)
    deltat2 = (t - refTime2).total_seconds() / (3600 * 24 * 365)

    if deltat1 < 0:
        return c['tau0_t0']
    elif deltat2 < 0:
        return c['tau0_t1'] + c['tau1_t1'] * deltat1
    else:
        return c['tau0_t2'] + c['tau1_t2'] * deltat1


def calc_deflection_correction_D(delta, sensorConst):
    # delta is 'DIM2_ELEVATION'
    d = delta * np.pi / 180
    c = sensorConst
    return (1. + c['c1'] / c['c0'] * d + c['c2'] / c['c0'] * d**2)

def calc_anode_correction_F(sensorConst):
    c = sensorConst
    F = np.zeros(16)
    F[0] = (c['b0'] + c['b1'] * 16) * c['f']
    F[-1] = (c['b0'] + c['b1'] * 1) * c['f']
    for i in range(1, 15):
        F[i] = (c['b0'] + c['b1'] * (16 - i))

    return F


def calc_simulation_geom_factor(energies):
    return 0.8 * 0.9 * (2.37e-4 + 9.17e-11 * energies)

def calc_sensor_dependence_exp(sensorConst, energies):
    c = sensorConst
    innerExp = 1 - np.exp(-c['rho'] * (energies / c['maxE'])**c['n'])
    return c['maxEpsilon'] \
        * (1 - np.exp(-c['k'] * (energies / c['maxE'])**(1 - c['n']) * innerExp))
