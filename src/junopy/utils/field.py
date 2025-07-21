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
Methods to calculate field quantities like the JRM magnetic field
model and the CON2020 current sheet model.
"""

import os
import numpy as np
import con2020
import scipy
import JupiterMag


def calc_assoc_legendre_poly(theta, N):
    """
    Calculates Gaussian normalized Legendre polynomials
    Assumes that the input is of the form P_nm(cos(theta)) NOT P_nm(x)
    """
    sintheta = np.sin(theta)
    costheta = np.cos(theta)

    P = np.zeros((N + 1, N + 1))
    K = np.zeros((N + 1, N + 1))
    diffP = np.zeros((N + 1, N + 1))

    P[0, 0] = 2
    diffP[0, 0] = 0

    n = 1
    P[1, 1] = sintheta * P[n - 1, n - 1]
    diffP[1, 1] = sintheta * diffP[n - 1, n - 1] + costheta * P[n - 1, n - 1]
    P[1, 0] = costheta * P[0, 0]
    diffP[1, 0] = costheta * diffP[0, 0] - sintheta * P[0, 0]

    for n in range(2, N + 1):
        P[n, n] = sintheta * P[n - 1, n - 1]
        diffP[n, n] = sintheta * diffP[n - 1, n - 1] + costheta * P[n - 1, n - 1]

        for m in range(0, n):
            K[n, m] = ((n - 1) ** 2 - m**2) / ((2 * n - 1) * (2 * n - 3))
            P[n, m] = costheta * P[n - 1, m] - K[n, m] * P[n - 2, m]
            diffP[n, m] = (
                costheta * diffP[n - 1, m]
                - sintheta * P[n - 1, m]
                - K[n, m] * diffP[n - 2, m]
            )

    return P, diffP


def calc_internal_field(r_a, theta, phi, g, h, N):
    """
    Calculates Psi, Br, Btheta, Bphi based on the spherical harmonics
    representation.
    """
    Psi = 0
    Br = 0
    Btheta = 0
    Bphi = 0

    a_r = 1.0 / r_a

    P, diffP = calc_assoc_legendre_poly(theta, N)

    sintheta = np.sin(theta)
    costheta = np.cos(theta)
    sinphi = np.sin(phi)
    cosphi = np.cos(phi)
    inv_sintheta = 0.0
    if np.abs(sintheta) > 1.0e-3:
        inv_sintheta = 1.0 / sintheta

    for n in range(0, N + 1):
        m = 0
        sinphi_prev = 0.0
        cosphi_prev = 1.0
        # Psi += a * (a_r)**(n + 1) * P[n, m] * g[n, m]
        Br += (a_r) ** (n + 2) * (n + 1) * P[n, m] * g[n, m]
        Btheta -= (a_r) ** (n + 2) * diffP[n, m] * g[n, m]
        Bphi += inv_sintheta * (a_r) ** (n + 2) * P[n, m] * m * (-h[n, m])

        for m in range(1, n + 1):
            sinmphi = sinphi_prev * cosphi + cosphi_prev * sinphi
            cosmphi = cosphi_prev * cosphi - sinphi_prev * sinphi

            Br += (
                (a_r) ** (n + 2)
                * (n + 1)
                * P[n, m]
                * (g[n, m] * cosmphi + h[n, m] * sinmphi)
            )

            Btheta -= (
                (a_r) ** (n + 2) * diffP[n, m] * (g[n, m] * cosmphi + h[n, m] * sinmphi)
            )

            Bphi += (
                inv_sintheta
                * (a_r) ** (n + 2)
                * P[n, m]
                * m
                * (g[n, m] * sinmphi - h[n, m] * cosmphi)
            )

            sinphi_prev = sinmphi
            cosphi_prev = cosmphi

        B = np.zeros(3)
        B[0] = (Br * sintheta * cosphi) + (Btheta * costheta * cosphi) - Bphi * sinphi
        B[1] = (Br * sintheta * sinphi) + (Btheta * costheta * sinphi) + Bphi * cosphi
        B[2] = (Br * costheta) - (Btheta * sintheta)

    return B


def read_coefficients_file(filename, degree=None, doNormalize=False):
    """Reads the coefficients and return after normalization."""

    N = degree if degree else 30
    g = np.zeros((N + 1, N + 1))
    h = np.zeros((N + 1, N + 1))

    with open(filename, "r") as f:
        f.readline()
        for n in range(1, N + 1):
            for m in range(0, n + 1):
                g[n, m] = float(f.readline().split()[1])
                # print('g', n, m, g[n,m])

            for m in range(1, n + 1):
                h[n, m] = float(f.readline().split()[1])
                # print('h', n, m, h[n,m])

    # Normalize the Gauss coefficients
    if doNormalize:
        S = np.zeros((N + 1, N + 1))
        S[0, 0] = 1
        for n in range(1, N + 1):
            S[n, 0] = S[n - 1, 0] * (2 * n - 1) / n
            S[n, 1] = S[n, 0] * np.sqrt(n * 2.0 / (n + 1))

            g[n, 0] = g[n, 0] * S[n, 0]
            h[n, 0] = h[n, 0] * S[n, 0]
            g[n, 1] = g[n, 1] * S[n, 1]
            h[n, 1] = h[n, 1] * S[n, 1]

            for m in range(2, n + 1):
                S[n, m] = S[n, m - 1] * np.sqrt((n - m + 1) / (n + m))
                g[n, m] = g[n, m] * S[n, m]
                h[n, m] = h[n, m] * S[n, m]

    return g, h, N


def calc_currentsheet_field(pos, equation_type="analytic", docurve=True):
    if not docurve:
        tilt = 0.001

    def_model = con2020.Model(equation_type=equation_type, xt__cs_tilt_degs=tilt)

    return def_model.Field(pos[0], pos[1], pos[2])


def trace_to_planet(start, g, h, N, deltaS=0.2, scheme="Euler"):
    streamline = np.zeros((10000, 3)) * np.nan

    streamline[0] = start

    pos = start * 1.0

    #     for n in range(1, N + 1):
    #         for m in range(0, n + 1):
    #             if not (n == 1 and m == 0):
    #                 g[n, m] = 0.
    #             # print('g', n, m, g[n,m])

    #         for m in range(1, n + 1):
    #             h[n, m] = 0.

    i = 0
    while np.sqrt(np.sum(pos**2)) >= 0.99 and np.sqrt(np.sum(pos**2)) <= 100:
        r_a = np.sqrt(np.sum(pos**2))
        theta = np.arccos(pos[2] / r_a)
        phi = np.arctan2(pos[1], pos[0])

        Bx_int, By_int, Bz_int = calc_internal_field(r_a, theta, phi, g, h, N)

        # Bxyz = jovian_jrm33_order13_internal_xyz(pos[0], pos[1], pos[2])
        # Bx_int = Bxyz[0]
        # By_int = Bxyz[1]
        # Bz_int = Bxyz[2]

        B_cs = calc_currentsheet_field(pos)

        Bx = Bx_int + B_cs[0, 0]
        By = By_int + B_cs[0, 1]
        Bz = Bz_int + B_cs[0, 2]

        Bmag = np.sqrt(Bx**2 + By**2 + Bz**2)

        bunit = np.array([Bx / Bmag, By / Bmag, Bz / Bmag])

        if "Euler" in scheme:
            pos += deltaS * bunit
            streamline[i] = pos * 1.0
            i += 1
            # print(np.sqrt(np.sum(pos**2)), end='\r')

    return streamline


def trace_field_bequator(start, g, h, N, deltaS=0.2, scheme="RK4"):
    pos = start * 1.0
    sign_prev = 1
    sign_next = 1

    while np.sqrt(np.sum(pos**2)) >= 0.99 and np.sqrt(np.sum(pos**2)) <= 100:
        r_a = np.sqrt(np.sum(pos**2))
        theta = np.arccos(pos[2] / r_a)
        phi = np.arctan2(pos[1], pos[0])

        Bx_int, By_int, Bz_int = JupiterMag.Internal.Field(pos[0], pos[1], pos[2])
        # calc_internal_field(r_a, theta, phi, g, h, N)

        # Bxyz = jovian_jrm33_order13_internal_xyz(pos[0], pos[1], pos[2])
        # Bx_int = Bxyz[0]
        # By_int = Bxyz[1]
        # Bz_int = Bxyz[2]

        B_cs = calc_currentsheet_field(pos)

        Bx = Bx_int + B_cs[0, 0]
        By = By_int + B_cs[0, 1]
        Bz = Bz_int + B_cs[0, 2]

        Bmag = np.sqrt(Bx**2 + By**2 + Bz**2)

        bunit = np.array([Bx / Bmag, By / Bmag, Bz / Bmag])

        if "Euler" in scheme:
            pos += deltaS * bunit
            print(np.sqrt(np.sum(pos**2)), end="\r")

        return pos


def get_B(pos, g, h, N, equation_type="analytic", normalize=True, docurve=True):
    """Get the magnetic field unit vector."""
    r_a = np.sqrt(np.sum(pos**2))
    theta = np.arccos(pos[2] / r_a)
    phi = np.arctan2(pos[1], pos[0])

    B_int = calc_internal_field(r_a, theta, phi, g, h, N)

    # Bx_int, By_int, Bz_int = JupiterMag.Internal.Field(pos[0],pos[1],pos[2])
    #     # calc_internal_field(r_a, theta, phi, g, h, N)

    # B_int = np.array([Bx_int[0], By_int[0], Bz_int[0]])
    B_cs = calc_currentsheet_field(pos, equation_type=equation_type, docurve=docurve)
    B = B_int + B_cs

    if normalize:
        return B / np.sqrt(np.sum(B**2))
    else:
        return B


def trace_field_dormandprince(
    start, gCoeff, hCoeff, N, h=0.2, Atol=1.0, equation_type="analytic", docurve=True
):
    streamline = np.zeros((10000, 3)) * np.nan
    streamline[0, :] = start * 1.0

    pos = start * 1.0

    # Define the Dormand-Prince coefficients
    c = np.array([0.0, 1.0 / 5, 3.0 / 10, 4.0 / 5, 8.0 / 9, 1.0, 1.0])
    a = np.array(
        [
            [0, 0, 0, 0, 0, 0],
            [1.0 / 5, 0, 0, 0, 0, 0],
            [3.0 / 40, 9.0 / 40, 0, 0, 0, 0],
            [44.0 / 45, -56.0 / 15, 32.0 / 9, 0, 0, 0],
            [19372.0 / 6561, -25360.0 / 2187, 64448.0 / 6561, -212.0 / 729, 0.0, 0.0],
            [
                9017.0 / 3168,
                -355.0 / 33,
                46732.0 / 5247,
                49.0 / 176,
                -5103.0 / 18656,
                0.0,
            ],
            [35.0 / 384, 0, 500.0 / 1113, 125.0 / 192, -2187.0 / 6784, 11.0 / 84],
        ]
    )

    b = [35.0 / 384, 0, 500.0 / 1113, 125.0 / 192, -2187.0 / 6784, 11.0 / 84, 0]
    bcap = [
        5179.0 / 57600,
        0,
        7571.0 / 16695,
        393.0 / 640,
        -92097.0 / 339200,
        187.0 / 2100,
        1.0 / 40,
    ]

    k = np.zeros((7, 3))

    Rtol = 0.0

    facmax = 4
    facmin = 0.125

    count = 1
    while np.sqrt(np.sum(pos**2)) >= 0.999 and np.sqrt(np.sum(pos**2)) <= 100:
        pos_cap = pos * 1.0
        pos_new = pos * 1.0

        err = 2.0

        while err > 1.0:
            k[0, :] = get_B(
                pos_new,
                gCoeff,
                hCoeff,
                N,
                equation_type=equation_type,
                normalize=True,
                docurve=docurve,
            )
            ds = 0

            for i in range(1, 7):
                posStep = pos_new * 1.0
                for j in range(i):
                    posStep += h * a[i, j] * k[i - 1, :]

                k[i, :] = get_B(
                    posStep,
                    gCoeff,
                    hCoeff,
                    N,
                    equation_type=equation_type,
                    normalize=True,
                    docurve=docurve,
                )

                pos_new += h * b[i] * k[i, :]
                pos_cap += h * bcap[i] * k[i, :]

                ds += h * b[i] * k[i, :]

            sc = Atol + np.maximum(np.abs(pos), np.abs(pos_new)) * Rtol
            err = np.sqrt(1 / 3 * np.sum(((pos - pos_new) / sc) ** 2))

            h = h * np.minimum(facmax, np.maximum(facmin, 0.8 * (1 / err) ** (1.0 / 4)))
            # print(err, end='\r')
            # print(np.sqrt(np.sum(pos**2)), h, err, end='\r')

        pos = pos_new * 1.0
        streamline[count, :] = pos * 1.0
        count += 1
        # print(np.sqrt(np.sum(pos**2)), h, end='\r')

    return streamline  # [~np.all(streamline == 0, axis=1)]


class TraceMShell:
    G_COEFF = np.array([[0.0, 0.0], [410993.4, 0.0]])
    H_COEFF = np.array([[0.0, 0.0], [0.0, 0.0]])
    N_COEFF = 1

    def __init__(self, M):
        """
        Class to help trace a field line at a certain M-shell using a dipole +
        current sheet approximation.
        """
        self.pos0 = np.array([M, 0.0, 0.0])
        self.get_fieldline()
        self.calc_B0()
        self.calc_dBds()

    def _trace_field_to_planet(self, h=0.0002, Atol=0.01):
        """
        Trace a field line from some starting position "start" to the planet.
        """
        start = self.pos0

        streamline = np.zeros((10000, 3)) * np.nan
        s = np.zeros(10000) * np.nan

        streamline[0, :] = start * 1.0
        s[0] = 0.0

        pos = start * 1.0

        # Define the Dormand-Prince coefficients
        c = np.array([0.0, 1.0 / 5, 3.0 / 10, 4.0 / 5, 8.0 / 9, 1.0, 1.0])
        a = np.array(
            [
                [0, 0, 0, 0, 0, 0],
                [1.0 / 5, 0, 0, 0, 0, 0],
                [3.0 / 40, 9.0 / 40, 0, 0, 0, 0],
                [44.0 / 45, -56.0 / 15, 32.0 / 9, 0, 0, 0],
                [
                    19372.0 / 6561,
                    -25360.0 / 2187,
                    64448.0 / 6561,
                    -212.0 / 729,
                    0.0,
                    0.0,
                ],
                [
                    9017.0 / 3168,
                    -355.0 / 33,
                    46732.0 / 5247,
                    49.0 / 176,
                    -5103.0 / 18656,
                    0.0,
                ],
                [35.0 / 384, 0, 500.0 / 1113, 125.0 / 192, -2187.0 / 6784, 11.0 / 84],
            ]
        )

        b = [35.0 / 384, 0, 500.0 / 1113, 125.0 / 192, -2187.0 / 6784, 11.0 / 84, 0]
        bcap = [
            5179.0 / 57600,
            0,
            7571.0 / 16695,
            393.0 / 640,
            -92097.0 / 339200,
            187.0 / 2100,
            1.0 / 40,
        ]

        k = np.zeros((7, 3))

        Rtol = 0.0

        facmax = 4
        facmin = 0.125

        count = 1

        while np.sqrt(np.sum(pos**2)) >= 0.999 and np.sqrt(np.sum(pos**2)) <= 100:
            pos_cap = pos * 1.0
            pos_new = pos * 1.0

            err = 2.0
            while err > 1.0:
                k[0, :] = field.get_B(
                    pos_new,
                    self.G_COEFF,
                    self.H_COEFF,
                    self.N_COEFF,
                    equation_type="analytic",
                    normalize=True,
                    docurve=False,
                )

                ds = 0

                for i in range(1, 7):
                    posStep = pos_new * 1.0
                    for j in range(i):
                        posStep += h * a[i, j] * k[i - 1, :]

                    k[i, :] = field.get_B(
                        posStep,
                        self.G_COEFF,
                        self.H_COEFF,
                        self.N_COEFF,
                        equation_type="analytic",
                        normalize=True,
                        docurve=False,
                    )

                    pos_new += h * b[i] * k[i, :]
                    pos_cap += h * bcap[i] * k[i, :]

                    # ds += np.sqrt(np.sum((h * b[i] * k[i, :])**2))

                sc = Atol + np.maximum(np.abs(pos), np.abs(pos_new)) * Rtol
                err = np.sqrt(1 / 3 * np.sum(((pos - pos_new) / sc) ** 2))

                h = h * np.minimum(
                    facmax, np.maximum(facmin, 0.8 * (1 / err) ** (1.0 / 4))
                )

            ds = np.sqrt(np.sum((pos_new - pos) ** 2))
            pos = pos_new * 1.0

            streamline[count, :] = pos * 1.0
            s[count] = s[count - 1] + ds

            count += 1

        return streamline, s

    def get_fieldline(self):
        """
        Trace in both directions to the planet and get the coordinates of
        the entire field line.
        """
        traceS, sS = self._trace_field_to_planet(h=0.0002)
        traceN, sN = self._trace_field_to_planet(h=-0.0002)

        fieldline = np.concatenate(
            [np.flip(truncate_nan_rows(traceN), axis=0), truncate_nan_rows(traceS)[1:]]
        )

        sN = np.flip(sN[~np.isnan(sN)])
        sS = sS[~np.isnan(sS)]

        sN = np.abs(sN - sN[0])
        sS = sS[1:] + sN[-1]

        s = np.concatenate([sN, sS])

        self.trace = pd.DataFrame(
            data={
                "x": fieldline[:, 0],
                "y": fieldline[:, 1],
                "z": fieldline[:, 2],
                "s": s,
            }
        )

        self.trace["B"] = self.trace.apply(
            lambda x: np.sqrt(
                np.sum(
                    field.get_B(
                        np.array([x["x"], x["y"], x["z"]]),
                        self.G_COEFF,
                        self.H_COEFF,
                        1,
                        normalize=False,
                        docurve=False,
                    )
                    ** 2
                )
            ),
            axis=1,
        )

        # Find location of minimum B
        i = np.argmin(self.trace["B"])

        # Shift s so that s=0 is the magnetic equator
        self.trace.loc[:, "s"] = self.trace.loc[:, "s"] - self.trace.loc[i, "s"]

        # Find s0
        i0 = np.argmin(
            (self.trace["x"] - self.pos0[0]) ** 2
            + (self.trace["y"] - self.pos0[1]) ** 2
            + (self.trace["z"] - self.pos0[2]) ** 2
        )

        self.s0 = self.trace.loc[i0, "s"]

        # Find B0
        B0 = self.trace.loc[i0, "B"]

        # Sort dataframe so s is increasing
        self.trace.sort_values("s", axis=0, inplace=True)

    def calc_B0(self):
        """
        For axissymmetric field, the minimum field strength has to be at the
        equator.
        """
        magVector = field.get_B(
            self.pos0,
            self.G_COEFF,
            self.H_COEFF,
            self.N_COEFF,
            normalize=False,
            docurve=False,
        )

        self.B0 = np.sqrt(np.sum(magVector**2))

    def calc_dBds(self):
        self.trace["dBds"] = np.gradient(self.trace["B"], self.trace["s"])
        self.calc_dBds = CubicSpline(self.trace["s"], self.trace["dBds"])


def truncate_nan_rows(a):
    mask = np.all(np.isnan(a) | np.equal(a, 0), axis=1)
    return a[~mask]


def trace_particle(E, m, alpha0, pos0, tmax):
    """
    Trace the particle with energy E, mass m, initial pitch angle alpha0,
    position pos0 to its location along the field line at tmax.
    If tmax < 0, it will trace the particle backwards in time to tmax.
    """
    v0 = calc_relativistic_v(E / e, m=m)

    trace = TraceMShell(M=pos0[0])

    # Calculate mu that remains constant over the motion

    p = calc_p_from_kineticenergy(E / e, m=m)
    mu = p**2 / (2 * m) * np.sin(alpha0) ** 2 / trace.B0 / 1e-9

    # solve the initial value problem
    def func(_, F, mu, m):
        return np.array(
            [
                -mu
                / (m / np.sqrt(1 - (v0 / speed_of_light) ** 2))
                * trace.calc_dBds(F[1] / R_JUPITER)
                * 1e-9
                / R_JUPITER,
                F[0],
            ]
        )

    y0 = np.array([v0 * np.cos(alpha0), trace.s0 * R_JUPITER])

    return (
        solve_ivp(
            func,
            (0, tmax),
            y0,
            args=(mu, m),
            atol=(0.1, 0.0001 * R_JUPITER),
            dense_output=True,
        ),
        trace,
    )


def calc_bounce_period(E0, alpha0, trace=None, mShell=None, m=m_e):
    """
    Calculate the bounce period for a particle. You can pass either a `trace`
    object of type `TraceMShell` or pass an M-shell value. Tracing is expensive,
    so if calculating bounce periods for particles at different energies and
    pitch angles, it is faster to pass in a trace object instead.

    Args:
        E0: energy [eV]
        alpha0: equatorial pitch angle [radians]
        trace: TraceMShell object
        mShell: M-shell
        m: particle mass [kg]

    Returns:
        Bounce Period in [s]
    """
    v0 = calc_relativistic_v(E0 * e, m=m)

    if not trace and mShell:
        trace = TraceMShell(M=mShell)

    elif trace and not mShell:
        # If trace is already passed, do nothing.
        pass

    else:
        raise Exception("Can only pass `trace` or `mShell` but not both.")

    trace.trace["u"] = np.arccos(
        np.sqrt(1 - np.sin(alpha0) ** 2 * trace.trace["B"] / trace.B0)
    )

    traceFiltered = trace.trace[(~np.isnan(trace.trace["u"])) & (trace.trace["s"] > 0)]

    calc_dBds_u = CubicSpline(traceFiltered["u"], traceFiltered["dBds"])

    def eval_func(u):
        return 2 * mpmath.sin(u) / calc_dBds_u(u)

    return (
        4
        * trace.B0
        / np.sin(alpha0) ** 2
        / (v0 / 71492e3)
        * mpmath.quad(eval_func, [alpha0, np.pi / 2], method="tanh-sinh")
    )


def get_M_shell(x: float, y: float, z: float) -> tuple[float, float, float]:
    """Get the M shell location of Juno using Ray-Tracing."""
    JupiterMag.Con2020.Config(equation_type="analytic")
    T_analytic = JupiterMag.TraceField(x, y, z, IntModel="jrm33", ExtModel="Con2020")

    trace = T_analytic.TraceDict()
    Bmag = np.sqrt(trace["Bx"][0] ** 2 + trace["By"][0] ** 2 + trace["Bz"][0] ** 2)
    rad = np.sqrt(trace["x"][0] ** 2 + trace["y"][0] ** 2 + trace["z"][0] ** 2)
    lat = np.arctan2(
        trace["z"][0] ** 2, np.sqrt(trace["x"][0] ** 2 + trace["y"][0] ** 2)
    )
    lon = np.arctan2(trace["y"][0], trace["x"][0])

    M = rad[np.argmin(Bmag)]
    lat = lat[np.argmin(Bmag)]
    lon = lon[np.argmin(Bmag)]

    return M, lat * 180 / np.pi, lon * 180 / np.pi


def get_M_shell_time(t):
    """Get the M-shell at ONE time."""
    pos, _ = get_juno_position([t], coordsystem="IAU_JUPITER")
    pos = pos[0]
    return get_M_shell(pos[0], pos[1], pos[2])


def get_M_shell_times(times):
    position, _ = get_juno_position(times, coordsystem="IAU_JUPITER")
    resultMshells = np.zeros(position.shape[0])
    for i, pos in enumerate(position):
        resultMshells[i], _, _ = get_M_shell(pos[0], pos[1], pos[2])
    return resultMshells


def calc_delta_alfven_juno(times, satellite="IO"):
    """
    Use the magnetic field model to calculate the longitudinal separation
    between Juno's mapped location on Io's orbital plane and Io.
    This is incorrect as it does not account for the propagation time between
    Juno and Io.
    """
    posJuno, _ = spice.get_juno_position(times, coordsystem="IAU_JUPITER")

    deltaAlfven = np.zeros(len(times))
    print(posJuno.shape)
    print(posJuno)
    with withspiceloaded() as context:
        for i, pos in enumerate(posJuno):
            print(pos)
            # Trace field in IAU_JUPITER
            JupiterMag.Con2020.Config(equation_type="analytic")
            T_analytic = JupiterMag.TraceField(
                pos[0], pos[1], pos[2], IntModel="jrm33", ExtModel="Con2020"
            )
            trace = T_analytic.TraceDict()
            fieldtrace = np.array([trace["x"][0], trace["y"][0], trace["z"][0]])

            print("Fieldtrace shape - ", fieldtrace.shape)

            # Find rotation matrix to PHI-O system
            print(times[i].strftime("%Y-%m-%dT%H:%M:%S"))
            time = spice.str2et(times[i].strftime("%Y-%m-%dT%H:%M:%S"))

            R = spice.pxform("IAU_JUPITER", "JUNO_{:s}_PHI_O".format(satellite), time)
            print("R shape - ", R.shape)

            # Convert fieldline position array to PHI-O
            fieldtracePhiO = np.matmul(R, fieldtrace).T

            izero = np.where(np.diff(np.signbit(fieldtracePhiO[:, 2])))[0]

            rEqJuno = 0.5 * (fieldtrace[:, izero] + fieldtrace[:, izero + 1])
            posSat, _ = spice.spkpos(satellite, time, "IAU_JUPITER", "NONE", "JUPITER")

            posSat = np.asarray(posSat) / 71492

            print("posSat", posSat)
            print("rEqJuno", rEqJuno)

            deltaAlfven[i] = (
                np.arccos(
                    np.dot(rEqJuno.T, posSat)
                    / np.sqrt(np.sum(rEqJuno**2))
                    / np.sqrt(np.sum(posSat**2))
                )
                * 180
                / np.pi
            )

            xJunoPhiO = 0.5 * (fieldtracePhiO[izero, 0] + fieldtracePhiO[izero + 1, 0])

            if xJunoPhiO < 0:  # Wake
                deltaAlfven[i] = 180 - deltaAlfven[i]

    return deltaAlfven
