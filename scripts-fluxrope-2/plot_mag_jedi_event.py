import glob
import numpy as np
import pandas as pd
import datetime as dt

import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt 
# import matplotlib.colors as colors
# from matplotlib.patches import Rectangle
# from matplotlib.collections import PatchCollection

from JunoPy.jedi import get_jedi_data
from JunoPy.index import IndexFile
from JunoPy.mag import MagData
from JunoPy.library import find_ticks_between


def read_csv(filename):
    
    start_times = []
    end_times = []
    with open(filename, 'r') as f:
        
        f.readline()
        for line in f:
            start_time = dt.datetime.strptime(line.split(',')[0].strip(), '%Y-%jT%H:%M:%S')
            end_time = dt.datetime.strptime(line.split(',')[1].strip(), '%Y-%jT%H:%M:%S')
            
            yield start_time, end_time
            

def plot_data(magdata, electron_data, proton_data, nonproton_data, 
              start=None, end=None, outfolder=None, sevent=None, eevent=None):
    """Plot data between start and end."""
    
    # Find the ticks.
    timearray, tickstring = find_ticks_between(start, end)
    
    cmap = 'turbo'
    
    ndata = 1
    if electron_data: ndata += 2
    if proton_data: ndata += 2
    if nonproton_data[0]: ndata += 6
    
    # Create the figure.
    fig, axs = plt.subplots(ndata, 1, 
                            gridspec_kw={'hspace': 0.2, 'right':0.9, 'top':0.9}, 
                            figsize=(6.5, 11), 
                            dpi=300)
    
    axs[0:1] = magdata.plot(vars=['Btheta'], axs=axs[0:1])
    
    i = 1
    if electron_data:
        electron_data.plot_energy_spectra(ax=axs[i], cmap=cmap)
        electron_data.plot_pitchangle_spectra(ax=axs[i + 1], cmap=cmap)
        i += 2
    
    if proton_data:
        proton_data.plot_energy_spectra(ax=axs[i], cmap=cmap)
        proton_data.plot_pitchangle_spectra(ax=axs[i + 1], cmap=cmap)
        i += 2
        
    if nonproton_data[0]:
        for npdata in nonproton_data:
            # npdata.plot_energy_spectra(ax=axs[i], cmap='jet')
            npdata.plot_flux_timeseries(ax=axs[i])
            npdata.plot_pitchangle_spectra(ax=axs[i + 1], cmap=cmap)
            i += 2
        
    
#     nonproton_data[0].plot_pitchangle_spectra(
#         ax=axs[i], combine=[nonproton_data[1], nonproton_data[2]])
#     i += 1
    
    for a in axs:
        a.set_xlim([timearray[0], timearray[-1]])
        a.set_xticks(timearray)
        a.set_xticklabels(tickstring)
        a.label_outer()
        a.axvline(sevent, ls='dotted', color='r', lw=2.) # [0.2, 0.2, 0.2])
        a.axvline(eevent, ls='dotted', color='r', lw=2.)
        
        for item in ([a.yaxis.label] + a.get_xticklabels() + a.get_yticklabels()):
            item.set_fontsize(10)
            
            
    annotations = ['a)', 'b)', 'c)', 'd)', 'e)', 'f)', 'g)', 'h)', 'i)', 'j)', 'k)']
    for i, a in enumerate(axs):
        a.text(0.96, 0.8, annotations[i], transform=a.transAxes, weight='bold',
                    bbox=dict(facecolor='w', edgecolor='k', pad=2.))

    a.tick_params(axis='x', width=2, rotation=0)
    a.xaxis.label.set_visible(False)
    
    plt.suptitle(dt.datetime.strftime(start, '%Y-%m-%d') + \
                 ', DOY={:s}'.format(dt.datetime.strftime(start, '%j')), y=0.92)
   
    
    # plt.show()
    fig.savefig('20210717/jedi/' + dt.datetime.strftime(start, '%Y-%jT%H:%M') + '.png', 
                bbox_inches='tight',
                facecolor='white')
    
    plt.close(fig)

# --------------------------------------------------------------------
if __name__ == '__main__':
    
    for startevent, endevent in read_csv('auto_plasmoid_list_20210525_window10.csv'):

        print('Working on - {:}'.format(startevent))
        start = startevent - dt.timedelta(minutes=10)
        end   = endevent + dt.timedelta(minutes=10)

        if startevent.day != start.day:
            start = dt.datetime(startevent.year, startevent.month, startevent.day, 0, 1, 30)

        if endevent.day != end.day:
            end = dt.datetime(endevent.year, endevent.month, endevent.day, 23, 58, 30)

        magfile = glob.glob('data/mag/*{:s}ss*.hdf5'.format(startevent.strftime('%Y%j')))[0]
        magdata = MagData(filename=magfile)
        magdata.df = magdata.df.loc[start:end]
        
        try:
            electron_data, proton_data, nonproton_data = get_jedi_data(start, end, dtbin='30S')
        
            
            plot_data(magdata, electron_data, proton_data, nonproton_data, 
                      start=start, end=end, sevent=startevent, eevent=endevent)
        except:
            continue
 