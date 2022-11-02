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

    P[0, 0] = 1
    diffP[0, 0] = 0
    
    n = 1
    P[1, 1] = sintheta * P[n - 1, n - 1]
    diffP[1, 1] = sintheta * diffP[n - 1, n - 1] + \
        costheta * P[n - 1, n - 1]
    P[1, 0] = costheta * P[0, 0]
    diffP[1, 0] = costheta * diffP[0, 0] - sintheta * P[0, 0]
    
    for n in range(2, N + 1):
        P[n, n] = sintheta * P[n - 1, n - 1]
        diffP[n, n] = sintheta * diffP[n - 1, n - 1] + \
            costheta * P[n - 1, n - 1]

        for m in range(0, n):
            K[n, m] = ((n - 1)**2 - m**2) / ((2 * n - 1) * (2 * n - 3))
            P[n, m] = costheta * P[n - 1, m] - K[n, m] * P[n - 2, m]
            diffP[n, m] = costheta * diffP[n - 1, m] - \
                sintheta * P[n - 1, m] - K[n, m] * diffP[n - 2, m]

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

    a_r = 1. / r_a

    P, diffP = calc_assoc_legendre_poly(theta, N)

    sintheta = np.sin(theta)
    costheta = np.cos(theta)
    sinphi = np.sin(phi)
    cosphi = np.cos(phi)
    inv_sintheta = 0.0
    if(np.abs(sintheta) > 1.0e-3):
        inv_sintheta = 1. / sintheta

    for n in range(0, N + 1):

        m = 0
        sinphi_prev = 0.0
        cosphi_prev = 1.0
        # Psi += a * (a_r)**(n + 1) * P[n, m] * g[n, m]
        Br += (a_r)**(n + 2) * (n + 1) * P[n, m] * g[n, m]
        Btheta -= (a_r)**(n + 2) * diffP[n, m] * g[n, m]
        Bphi += inv_sintheta * (a_r)**(n + 2) * P[n, m] * m * (- h[n, m])

        for m in range(1, n + 1):

            sinmphi = sinphi_prev * cosphi + cosphi_prev * sinphi
            cosmphi = cosphi_prev * cosphi - sinphi_prev * sinphi

            Br += (a_r)**(n + 2) * (n + 1) * P[n, m] * (
                g[n, m] * cosmphi + h[n, m] * sinmphi)

            Btheta -= (a_r)**(n + 2) * diffP[n, m] * (
                g[n, m] * cosmphi + h[n, m] * sinmphi)

            Bphi += inv_sintheta * (a_r)**(n + 2) * P[n, m] * m * (
                g[n, m] * sinmphi - h[n, m] * cosmphi)

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
    
    with open(filename, 'r') as f:
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
            S[n, 1] = S[n, 0] * np.sqrt(n * 2. / (n + 1))

            g[n, 0] = g[n, 0] * S[n, 0]
            h[n, 0] = h[n, 0] * S[n, 0]
            g[n, 1] = g[n, 1] * S[n, 1]
            h[n, 1] = h[n, 1] * S[n, 1]

            for m in range(2, n + 1):
                S[n, m] = S[n, m - 1] * np.sqrt((n - m + 1) / (n + m))
                g[n, m] = g[n, m] * S[n, m]
                h[n, m] = h[n, m] * S[n, m]
            
    return g, h, N


def calc_currentsheet_field(pos, equation_type='analytic'):
    def_model = con2020.Model(equation_type=equation_type)
    return def_model.Field(pos[0], pos[1], pos[2])


def trace_to_planet(start, g, h, N, deltaS=0.2, scheme='Euler'):
    
    streamline = np.zeros((10000, 3)) * np.NaN
    
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
        
        if 'Euler' in scheme:
            pos += deltaS * bunit
            streamline[i] = pos * 1.0 
            i += 1
            # print(np.sqrt(np.sum(pos**2)), end='\r')
            
    return streamline
         
    
def trace_field_bequator(start, g, h, N, deltaS=0.2, scheme='RK4'):
    
    pos = start * 1.0
    sign_prev = 1
    sign_next = 1
    
    while np.sqrt(np.sum(pos**2)) >= 0.99 and np.sqrt(np.sum(pos**2)) <= 100:
        
        r_a = np.sqrt(np.sum(pos**2))
        theta = np.arccos(pos[2] / r_a)
        phi = np.arctan2(pos[1], pos[0])
      
        Bx_int, By_int, Bz_int = JupiterMag.Internal.Field(pos[0],pos[1],pos[2]) 
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
        
        if 'Euler' in scheme:
            pos += deltaS * bunit
            print(np.sqrt(np.sum(pos**2)), end='\r')
            
        return pos
    
            
def get_bunit(pos, g, h, N, equation_type='analytic'):
    """Get the magnetic field unit vector."""
    r_a = np.sqrt(np.sum(pos**2))
    theta = np.arccos(pos[2] / r_a)
    phi = np.arctan2(pos[1], pos[0])
    
    # B_int = calc_internal_field(r_a, theta, phi, g, h, N)
    
    Bx_int, By_int, Bz_int = JupiterMag.Internal.Field(pos[0],pos[1],pos[2]) 
        # calc_internal_field(r_a, theta, phi, g, h, N)
        
    B_int = np.array([Bx_int[0], By_int[0], Bz_int[0]])
    B_cs = calc_currentsheet_field(pos, equation_type=equation_type)
    B = B_int + B_cs

    return B / np.sqrt(np.sum(B**2))
    
    
        
def trace_field_dormandprince(start, gCoeff, hCoeff, N, h=0.2, Atol=1.0, equation_type='analytic'):
    
    streamline = np.zeros((10000, 3)) * np.NaN
    streamline[0, :] = start * 1.0
    
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
    while np.sqrt(np.sum(pos**2)) >= 0.99 and np.sqrt(np.sum(pos**2)) <= 100:
        
        pos_cap = pos * 1.0
        pos_new = pos * 1.0
        
        err = 2.
        
        while err > 1.0:

            k[0, :] = get_bunit(pos_new, gCoeff, hCoeff, N, equation_type=equation_type)

            for i in range(1, 7):

                posStep = pos_new * 1.0
                for j in range(i):
                    posStep += h * a[i, j] * k[i - 1, :]

                k[i, :] = get_bunit(posStep, gCoeff, hCoeff, N)

                pos_new += h * b[i] * k[i, :]
                pos_cap += h * bcap[i] * k[i, :]

            sc = Atol + np.maximum(np.abs(pos), np.abs(pos_new)) * Rtol
            err = np.sqrt(1 / 3 * np.sum(((pos - pos_new) / sc)**2))
            
            h = h * np.minimum(facmax, np.maximum(facmin, 0.8 * (1 / err)**(1./4)))
            # print(err, end='\r')
            #print(np.sqrt(np.sum(pos**2)), h, err, end='\r')
            
        pos = pos_new * 1.0
        streamline[count, :] = pos * 1.0
        count += 1
        # print(np.sqrt(np.sum(pos**2)), h, end='\r')
            
    return streamline # [~np.all(streamline == 0, axis=1)]
    
    
    
def trace_field_scipy(start, gCoeff, hCoeff, N):
    
    pass
    