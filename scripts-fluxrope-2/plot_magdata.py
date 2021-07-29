import os
import re
import glob
import time
import numpy as np
import datetime as dt
import matplotlib.pyplot as plt

from matplotlib.ticker import AutoMinorLocator

from JunoPy.mag import MagData
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
    
    fig.savefig('out_magdata_4hr/' + dt.datetime.strftime(start, '%Y-%jT%H:%M') + '.png')
    plt.close(fig)
    

# =================================================================================    
if __name__ == '__main__':
    
    magfiles = glob.glob('data/mag/*.hdf5')
    
    if not os.path.exists('out_magdata_4hr'):
        os.mkdir('out_magdata_4hr')
    
    nfiles = len(magfiles)
    for i, magfile in enumerate(magfiles):

        stime = time.time()

        # Identify the date.
        match = re.search('(\d{7})', magfile).group(0)
        year, doy = int(match[:4]), int(match[-3:])

        # Read in the file. 
        magdata = MagData(filename=magfile)

        # Assign initial Start and End dates.
        start = dt.datetime.strptime('{:d}-{:d}'.format(year, doy), '%Y-%j')
        end = start + dt.timedelta(hours=4)

        for j in range(6):
            plot_magdata(magdata, start, end)
            start = end + dt.timedelta(hours=0)
            end = start + dt.timedelta(hours=4)

        etime = time.time()
        duration = etime - stime
        print('Completed - {:d}/{:d}, Time - {:.2f}'.format(i + 1, nfiles, duration), end='\r')
        