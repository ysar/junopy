
import os
import numpy as np
import datetime as dt
import matplotlib.pyplot as plt

from JunoPy.mag import MagData
from JunoPy.waves import WavesData
from JunoPy.library import find_ticks


def test_plot():
    """
    Plot all quantities.
    """
    fig, ax = plt.subplots(6, figsize=(8.5, 11),
                           gridspec_kw={'hspace': 0, 'right':0.9, 'top':0.9})

    ax[:4] = t_plot_mag_data(ax[:4])
    ax[4], ax[5] = t_plot_waves_data(ax[4:6])

    date = dt.datetime.strptime('2018-043', '%Y-%j')
    timearray, tickstring = find_ticks(date)


    for a in ax:

        a.set_xlim([timearray[0], timearray[-1]])
        a.set_xticks(timearray)
        a.set_xticklabels(tickstring)

        a.label_outer()

    a.text(-0.05, -0.5,
           'Time\nr($R_J$) \nLT\nLat($^\\circ$)',
           transform=a.transAxes,
           horizontalalignment='right', weight='bold')

    plt.suptitle(dt.datetime.strftime(date, '%Y-%m-%d') + \
                 ', DOY={:s}'.format(dt.datetime.strftime(date, '%j')), y=0.92)

    plt.savefig('test_plot.png', dpi=300)



def t_plot_mag_data(axs):
    """
    Test plotting for magnetometer data.
    """
    filename = 'data/mag/fgm_jno_l3_2018043pc_r1s_v01.sts'

    magdata = MagData(filename=filename)
    magdata.convert_to_spherical()

    return magdata.plot(vars=['Br', 'Btheta', 'Bphi', 'Bmag'], ax=axs)


def t_plot_waves_data(axs):
    """Test plot waves survey data."""
    filename1 = "data/waves-survey/WAV_2018043T000000_B_V01.CSV"
    filename2 = "data/waves-survey/WAV_2018043T000000_E_V01.CSV"
    wavesdata1 = WavesData(filename=filename1)
    wavesdata2 = WavesData(filename=filename2)

    return wavesdata1.plot(ax=axs[0]), wavesdata2.plot(ax=axs[1])


def t_plot_jedi_data(axs):
    """Test plot the JEDI data.
    Things to plot -
        - Proton Pitch angle
        - Proton E v/s time
        - Heavy pitch angle
        - Heavy Energy v/s time
    """
    filename = "data/jedi/JED_180_LOERSESP_CDR_2018043_V03.TAB"
    jedidata = JediData(filename)

    return jedidata.plot(ax=axs)
