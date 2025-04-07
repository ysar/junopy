import numpy as np
import pandas as pd
import mpmath

from . import field
from .library import R_JUPITER

from scipy.constants import m_e, e, speed_of_light
from scipy.interpolate import CubicSpline
from scipy.integrate import solve_ivp


def calc_p_from_kineticenergy(energy, m=m_e):
    """
    Calculate relativistic momentum from kinetic energy. [SI units]
    """
    return 1. / speed_of_light * np.sqrt(
        (energy + m * speed_of_light**2)**2 - m**2 * speed_of_light**4)

def calc_v(energy, m=m_e):
    """
    Calculate non-relativistic speed from kinetic energy. [SI units]"""
    return np.sqrt(2 * energy / m)


def calc_relativistic_v(energy, m=m_e):
    """Calculate relativistic speed from kinetic energy. [SI units]"""
    return speed_of_light * np.sqrt(
        1. - 1. / (energy / m / speed_of_light**2 + 1.)**2)


def calc_relativistic_E(velocity, m=m_e):
    """Calculate relativistic kinetic energy. [SI units]"""

    p = m * velocity / np.sqrt(1 - (velocity / speed_of_light)**2)

    kineticEnergy = \
        np.sqrt((p * speed_of_light)**2 + m**2 * speed_of_light**4)\
        - m * speed_of_light**2
    
    return kineticEnergy


class TraceMShell:
    
    G_COEFF = np.array([[0., 0.], [410993.4, 0.]])
    H_COEFF = np.array([[0., 0.], [0., 0.]])
    N_COEFF = 1

    def __init__(self, M):
        """
        Class to help trace a field line at a certain M-shell using a dipole +
        current sheet approximation.
        """
        self.pos0 = np.array([M, 0., 0.])
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
        s[0] = 0.
    
        pos = start * 1.0
    
        # Define the Dormand-Prince coefficients
        c = np.array([0., 1./5, 3./10, 4./5, 8./9, 1., 1.])
        a = np.array([
            [0, 0, 0, 0, 0, 0],
            [1./5, 0, 0, 0, 0, 0],
            [3./40, 9./40, 0, 0, 0, 0],
            [44./45, -56./15, 32./9, 0, 0, 0],
            [19372./6561, -25360./2187, 64448./6561, -212./729, 0., 0.],
            [9017./3168, -355./33, 46732./5247, 49./176, -5103./18656, 0.],
            [35./384, 0, 500./1113, 125./192, -2187./6784, 11./84]])

        b = [35./384, 0, 500./1113, 125./192, -2187./6784, 11./84, 0]
        bcap = [5179./57600, 0, 7571./16695, 393./640, -92097./339200, 187./2100, 1./40]
    
        k = np.zeros((7, 3))
    
        Rtol = 0.
    
        facmax = 4
        facmin = 0.125
    
        count = 1
     
        while np.sqrt(np.sum(pos**2)) >= 0.999 and np.sqrt(np.sum(pos**2)) <= 100:
            
            pos_cap = pos * 1.0
            pos_new = pos * 1.0
        
            err = 2.
            while err > 1.0:

                k[0, :] = field.get_B(
                    pos_new, self.G_COEFF, self.H_COEFF, self.N_COEFF, 
                    equation_type='analytic', normalize=True, docurve=False)
            
                ds = 0
            
                for i in range(1, 7):

                    posStep = pos_new * 1.0
                    for j in range(i):
                        posStep += h * a[i, j] * k[i - 1, :]

                    k[i, :] = field.get_B(
                        posStep, self.G_COEFF, self.H_COEFF, self.N_COEFF, 
                        equation_type='analytic', normalize=True, docurve=False)

                    pos_new += h * b[i] * k[i, :]
                    pos_cap += h * bcap[i] * k[i, :]

                    # ds += np.sqrt(np.sum((h * b[i] * k[i, :])**2))

                sc = Atol + np.maximum(np.abs(pos), np.abs(pos_new)) * Rtol
                err = np.sqrt(1 / 3 * np.sum(((pos - pos_new) / sc)**2))
            
                h = h * np.minimum(facmax, np.maximum(facmin, 0.8 * (1 / err)**(1./4)))
            
            ds = np.sqrt(np.sum((pos_new - pos)**2))
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
        
        fieldline = np.concatenate([
            np.flip(truncate_nan_rows(traceN), axis=0),
            truncate_nan_rows(traceS)[1:]])
        
        sN = np.flip(sN[~np.isnan(sN)])
        sS = sS[~np.isnan(sS)]
        
        sN = np.abs(sN - sN[0])
        sS = sS[1:] + sN[-1]

        s = np.concatenate([sN, sS])

        self.trace = pd.DataFrame(
            data = {
                'x': fieldline[:, 0],
                'y': fieldline[:, 1],
                'z': fieldline[:, 2],
                's': s
            }
        )

        self.trace['B'] = self.trace.apply(
            lambda x: np.sqrt(np.sum(field.get_B(
                np.array([x['x'], x['y'], x['z']]), 
                self.G_COEFF, self.H_COEFF, 1, normalize=False, docurve=False)**2)), 
            axis=1)


        # Find location of minimum B
        i = np.argmin(self.trace['B'])
        
        # Shift s so that s=0 is the magnetic equator
        self.trace.loc[:, 's'] = (self.trace.loc[:, 's'] - self.trace.loc[i, 's'])
        
        # Find s0
        i0 = np.argmin(
            (self.trace['x'] - self.pos0[0])**2 + \
            (self.trace['y'] - self.pos0[1])**2 + \
            (self.trace['z'] - self.pos0[2])**2)

        self.s0 = self.trace.loc[i0, 's']

        # Find B0
        B0 = self.trace.loc[i0, 'B']

        # Sort dataframe so s is increasing
        self.trace.sort_values('s', axis=0, inplace=True)

    
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
        
        self.trace['dBds'] = np.gradient(self.trace['B'], self.trace['s'])
        self.calc_dBds = CubicSpline(self.trace['s'], self.trace['dBds'])



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
    mu = p**2 / (2 * m) * np.sin(alpha0)**2 / trace.B0 / 1e-9


    # solve the initial value problem
    def func(_, F, mu, m):
        return np.array([
            -mu / (m / np.sqrt(1 - (v0 / speed_of_light)**2)) 
                * trace.calc_dBds(F[1] / R_JUPITER) * 1e-9 / R_JUPITER, 
            F[0]])

    y0 = np.array([v0 * np.cos(alpha0), trace.s0 * R_JUPITER])

    return solve_ivp(
        func, (0, tmax), y0, args=(mu, m),
        atol=(0.1, 0.0001 * R_JUPITER), dense_output=True), trace


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

        
    trace.trace['u'] = np.arccos(
        np.sqrt(1 - np.sin(alpha0)**2 * trace.trace['B'] / trace.B0))

    traceFiltered = trace.trace[
        (~np.isnan(trace.trace['u'])) & (trace.trace['s'] > 0)]
    
    calc_dBds_u = CubicSpline(
        traceFiltered['u'], traceFiltered['dBds'])

    def eval_func(u):
        return 2 * mpmath.sin(u) / calc_dBds_u(u)

    return 4\
        * trace.B0\
        / np.sin(alpha0)**2\
        / (v0 / 71492e3)\
        * mpmath.quad(eval_func, [alpha0, np.pi / 2], method='tanh-sinh')
