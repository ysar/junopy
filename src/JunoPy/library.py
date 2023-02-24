"""
Library
-------
Contains constants and hard-links, also contains some general functions.

Constants:
:param int R_JUPITER: Radius of Jupiter in km
:param dict DATASETS: 
    Dictionary containing technical names of the datasets.
    Keys - mag, waves-burst, waves-survey, jedi, jade 
:param dict EXTENSIONS: 
    Dictionary containing file extensions for each dataset.
:param dict TIMEFORMAT:
    Dictionary containing time formats for different instruments.
    
"""
import os
import numpy as np
import datetime as dt
import pandas as pd
import JupiterMag

from JunoPy.spice import get_juno_position

R_JUPITER = 71492

PDS_URL = 'https://pds-ppi.igpp.ucla.edu/ditdos/download?id=pds://PPI/'

DATASETS = {
    'mag': 'JNO-J-3-FGM-CAL-V1.0',
    'waves-burst': 'JNO-E_J_SS-WAV-3-CDR-BSTFULL-V1.0',
    'waves-survey': 'JNO-E_J_SS-WAV-3-CDR-SRVFULL-V2.0',
    'jedi': 'JNO-J-JED-3-CDR-V1.0',
    'jade': 'JNO-J_SW-JAD-3-CALIBRATED-V1.0'
}

EXTENSIONS = {
    'mag': '.sts',
    'waves-burst': '.DAT',
    'waves-survey': '.CSV',
    'jedi': '.TAB',
    'jade': '.DAT'
}

TIMEFORMAT = {
    'mag': '%Y-%jT%H:%M:%S.%f',
    'waves-burst': '%Y-%m-%dT%H:%M:%S.%f',
    'waves-survey': '%Y-%m-%dT%H:%M:%S.%f',
    'jedi': '%Y-%jT%H:%M:%S.%f',
    'jade': '%Y-%jT%H:%M:%S.%f'
}


def convert_pos_spherical(pos):
    """
    Given coordinates (x,y,z) in pos[0], pos[1], pos[2], return R, LT, Latitude
    
    :param list pos: 3-value list containing spacecraft position in cartesian 
                     coordinates
    ::
    """
    R = np.sqrt(pos[0]**2 + pos[1]**2 + pos[2]**2)
    LT = 6 / 90 * np.arctan2(pos[1], pos[0]) * 180 / np.pi + 12
    Latitude = np.arctan2(pos[2], np.sqrt(pos[0]**2 + pos[1]**2)) * 180 / np.pi

    return R, LT, Latitude


def find_ticks_between(start, end, nticks=6):
    """
    Given an axes object make the ticks appropriately for subplots.
    """
    time_array = pd.date_range(
        start=start, end=end, periods=nticks).to_pydatetime().tolist()

    position, _ = get_juno_position(time_array)


    tickstring = []
    for i, pos in enumerate(position):

        R, LT, Latitude = convert_pos_spherical(pos)
        
        if i == 0:
            tickstring.append(
            'Time - {:s}'.format(time_array[i].strftime('%H:%M:%S')) + '\n' +
            'R  - {:7.2f}'.format(R) + '\n' +
            'LT - {:7.2f}'.format(LT) + '\n' +
            'Latitude - {:7.2f}'.format(Latitude))
        else:
            tickstring.append(
                time_array[i].strftime('%H:%M:%S') + '\n' +
                '{:7.2f}'.format(R) + '\n' +
                '{:7.2f}'.format(LT) + '\n' +
                '{:7.2f}'.format(Latitude))

    return time_array, tickstring


def find_ticks_between_GPhiO(start, end, nticks=6):
    """
    Given an axes object make the ticks appropriately for subplots.
    """
    time_array = pd.date_range(
        start=start, end=end, periods=nticks).to_pydatetime().tolist()

    position, _ = get_juno_position(time_array,
                                    coordsystem='JUNO_GANYMEDE_PHI_O',
                                    center='503')


    tickstring = []
    for i, pos in enumerate(position):
        
        if i == 0:
            tickstring.append(
                'Time - {:s}'.format(time_array[i].strftime('%H:%M:%S')) + '\n' +
                'X (R) - {:7.2f}'.format(pos[0]) + '\n' +
                'Y (R) - {:7.2f}'.format(pos[1]) + '\n' +
                'Z (R) - {:7.2f}'.format(pos[2]))
        else:
            tickstring.append(
                '{:s}'.format(time_array[i].strftime('%H:%M:%S')) + '\n' +
                '{:7.2f}'.format(pos[0]) + '\n' +
                '{:7.2f}'.format(pos[1]) + '\n' +
                '{:7.2f}'.format(pos[2]))

    return time_array, tickstring

def get_M_shell(x, y, z):
        """Get the M shell location of Juno using Ray-Tracing."""
        JupiterMag.Con2020.Config(equation_type='analytic')
        T_analytic = JupiterMag.TraceField(x, y, z,
                                IntModel='jrm33', ExtModel='Con2020')
        
        trace = T_analytic.TraceDict()
        Bmag = np.sqrt(trace['Bx'][0]**2 + trace['By'][0]**2 + trace['Bz'][0]**2)
        rad = np.sqrt(trace['x'][0]**2 + trace['y'][0]**2 + trace['z'][0]**2)
        lat = np.arctan2(trace['z'][0]**2, np.sqrt(trace['x'][0]**2 + trace['y'][0]**2))
        lon = np.arctan2(trace['y'][0], trace['x'][0])

        M = rad[np.argmin(Bmag)]
        lat = lat[np.argmin(Bmag)]
        lon = lon[np.argmin(Bmag)]

        return M, lat * 180 / np.pi, lon * 180 / np.pi

def find_ticks_between_Mshell(start, end, nticks=6, timeformat='%H:%M:%S'):
    """
    Given an axes object make the ticks appropriately for subplots.
    We would like to have :
    M shell, Latitude, SIII-Longitude
    """
    time_array = pd.date_range(
        start=start, end=end, periods=nticks).to_pydatetime().tolist()

    position, _ = get_juno_position(time_array, coordsystem='IAU_JUPITER')
    position_mag, _ = get_juno_position(time_array, coordsystem='JUNO_MAG_VIP4')

    tickstring = []
    for i, pos in enumerate(position):

        M, _, lon = get_M_shell(pos[0], pos[1], pos[2])

        mlat = np.arctan2(
            position_mag[i, 2], 
            np.sqrt(position_mag[i, 0]**2 + position_mag[i, 1]**2)) \
                * 180 / np.pi

        if i == 0:
            tickstring.append(
            'Time   {:s}'.format(time_array[i].strftime(timeformat)) + '\n' +
            'M     {:7.2f}'.format(M) + '\n' +
            'MLat   {:7.2f}'.format(mlat) + '\n' +
            'Lon   {:7.2f}'.format(lon))
        else:
            tickstring.append(
                time_array[i].strftime(timeformat) + '\n' +
                '{:7.2f}'.format(M) + '\n' +
                '{:7.2f}'.format(mlat) + '\n' +
                '{:7.2f}'.format(lon))

    return time_array, tickstring


## Charge exchange cross-sections
def cxcs_Hplus_O(E):
    """
    Lindsay and Stebbings (2005)
    """
    a1 = 2.91
    a2 = 0.0886
    a3 = 50.9
    a4 = 4.73
    a5 = -0.862
    a6 = 0.0306
    return (a1 - a2 * np.log(E))**2 * (1 - np.exp(-a3 / E))**2 \
         + (a4 - a5 * np.log(E)) * (1 - np.exp(-a6 / E))**2

def cxcs_Oplus_O(E):
    """
    Lindsay and Stebbings (2005)
    """
    a1 = 4.07
    a2 = 0.269
    a3 = 4125
    return (a1 - a2 * np.log(E))**2 * (1 - np.exp(-a3 / E))**0.8

def cxcs_Hplus_H2(Ein):
    """
    Tabata and Shirai (2000)
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
    
    E = Ein - 2.5e-3
    ER = 1.361e-2
    
    f1 = 1e-16 * a1 * (E / ER)**a2
    f2 = f1 / (1 + (E / a3)**(a2 + a4))
    
    f1_2 = 1e-16 * a5 * (E / ER)**a2
    f3_2 = f1_2 / (1 + (E / a6)**(a2 + a7) + (E / a8)**(a2 + a9))
    
    return (f2 + f3_2) / 1e-16  #E + 2.5e-3, 