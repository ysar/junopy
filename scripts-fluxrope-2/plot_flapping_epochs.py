import os
import re
import glob
import numpy as np
import datetime as dt
import matplotlib.pyplot as plt

from matplotlib.ticker import AutoMinorLocator

from JunoPy.mag import MagData, mva
from JunoPy.library import find_ticks_between


def plot_magdata(magdata, start, end):
    """Plot data between start and end."""
    
    # Slice the magnetometer data.
    # magdata.df = magdata.df.loc[start:end]
    
    # Find the ticks.
    timearray, tickstring = find_ticks_between(start, end)
    
    # Create the figure.
    fig, axs = plt.subplots(4, 1, 
                            gridspec_kw={'hspace': 0, 'right':0.9, 'top':0.9}, 
                            figsize=(11, 8.5), 
                            dpi=300)
    
    axs[:] = magdata.plot(vars=['Br', 'Btheta', 'Bphi', 'Bmag'], axs=axs)
    
    for a in axs:
        a.set_xlim([timearray[0], timearray[-1]])
        a.set_xticks(timearray)
        a.set_xticklabels(tickstring)
        a.label_outer()
    a.tick_params(axis='x', width=2, rotation=0)
    a.xaxis.label.set_visible(False)
    
    plt.suptitle(dt.datetime.strftime(start, '%Y-%m-%d') + \
                 ', DOY={:s}'.format(dt.datetime.strftime(start, '%j')), y=0.92)
    
    fig.savefig('out_flapping_examples/' + dt.datetime.strftime(start, '%Y-%jT%H:%M') + '.png')
    plt.close(fig)
    
    
# =================================================================
if __name__ == '__main__':
    
    dates = []
    with open('events_flapping.csv', 'r') as f:
        f.readline()
        for line in f:
            dates.append(line)
        
    for date in dates:
        print(date.split(',')[0].strip())
        
        datesplit = date.split(',')
        start, end = datesplit[0].strip(), datesplit[1].strip()
        deltaT = dt.timedelta(hours=2)
        start = dt.datetime.strptime(start, '%Y-%jT%H:%M:%S') - deltaT
        end = dt.datetime.strptime(end, '%Y-%jT%H:%M:%S') + deltaT
        magfile = glob.glob('data/mag/*{:s}ss*.hdf5'.format(start.strftime('%Y%j')))[0]
        magdata = MagData(filename=magfile)
        plot_magdata(magdata, start, end)