#
#   Contains constants, URLs etc.
#
import os
import numpy as np
import datetime as dt


PDS_URL = 'https://pds-ppi.igpp.ucla.edu/ditdos/download?id=pds://PPI/'

DATASETS = {
    'mag': 'JNO-J-3-FGM-CAL-V1.0',
    'waves-burst': 'JNO-E_J_SS-WAV-3-CDR-BSTFULL-V1.0',
    'waves-survey': 'JNO-E_J_SS-WAV-3-CDR-SRVFULL-V1.0',
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
    """
    R = np.sqrt(pos[0]**2 + pos[1]**2 + pos[2]**2)
    LT = 6 / 90 * np.arctan2(pos[1], pos[0]) * 180 / np.pi + 12
    Latitude = np.arctan2(pos[2], np.sqrt(pos[0]**2 + pos[1]**2)) * 180 / np.pi

    return R, LT, Latitude


def find_ticks(time, duration='day', nticks=6):
    """
    Given an axes object make the ticks appropriately for subplots.
    """
    from JunoPy.spice import InitSpiceKernels, get_juno_position

    if duration == 'day':
        starttime = dt.datetime(time.year, time.month, time.day)
        endtime = dt.datetime(time.year, time.month, time.day) + \
                  dt.timedelta(days=1)
        spacing = dt.timedelta(seconds =  24 * 3600 / nticks)

    elif duration == '4hrs':
        starttime = time + dt.timedelta(hours=0)
        endtime = starttime + dt.timedelta(hours=4)
        spacing = dt.timedelta(seconds =  4 * 3600 / nticks)

    time_array = np.arange(
        starttime,
        endtime + spacing, spacing).astype(dt.datetime)

    InitSpiceKernels()
    position, _ = get_juno_position(time_array)

    tickstring = []
    for i, pos in enumerate(position):

        R, LT, Latitude = convert_pos_spherical(pos)

        # if i == 0:
        #     tickstring.append(
        #         'Time\nRadial Distance ($R_J$) \nLocal Time\nLatitude ($^\\circ$)')
        # else:
        tickstring.append(
            time_array[i].strftime('%H:%M') + '\n' +
            '{:7.2f}'.format(R) + '\n' +
            '{:7.2f}'.format(LT) + '\n' +
            '{:7.2f}'.format(Latitude))

    return time_array, tickstring
