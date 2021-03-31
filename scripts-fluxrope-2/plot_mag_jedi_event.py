import glob
import numpy as np
import pandas as pd
import datetime as dt

import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt 
import matplotlib.colors as colors
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection

from JunoPy.jedi import JediData, merge_data_species
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
            

def plot_data(magdata, ion_data=None, electron_data=None, start=None, end=None, outfolder=None, sevent=None, eevent=None):
    """Plot data between start and end."""
    
    # Find the ticks.
    timearray, tickstring = find_ticks_between(start, end)
    
    # Create the figure.
    fig, axs = plt.subplots(4 + 4, 1, 
                            gridspec_kw={'hspace': 0.1, 'right':0.9, 'top':0.9}, 
                            figsize=(9, 11), 
                            dpi=150)
    
    axs[:4] = magdata.plot(vars=['Br', 'Btheta', 'Bphi', 'Bmag'], axs=axs[:4])
    
    if ion_data:
        ion_data.plot_flux_spectra(ax=axs[4])
        xticks = ion_data.plot_pitchangle_spectra(ax=axs[6], vmin=1e5, vmax=1e7, cmap='jet')
    
    if electron_data:
        electron_data.plot_flux_spectra(ax=axs[5])   
        xticks = electron_data.plot_pitchangle_spectra(ax=axs[7], vmin=1e5, vmax=1e7, cmap='jet')
        
        
    
    for a in axs[:-2]:
        a.set_xlim([timearray[0], timearray[-1]])
        a.set_xticks(timearray)
        a.set_xticklabels(tickstring)
        a.label_outer()
        a.axvline(sevent, ls='dotted', color='r') # [0.2, 0.2, 0.2])
        a.axvline(eevent, ls='dotted', color='r')
        
    for a in axs[-2:]:
        
        start_int = np.datetime64(start).astype('int64') * 1000
        end_int = np.datetime64(end).astype('int64') * 1000
        a.set_xlim([start_int, end_int])
        a.set_xticks(xticks)
        a.set_xticklabels(tickstring)
        a.label_outer()

        
    for a in axs[-2:]:
        start_int = np.datetime64(sevent).astype('int64') * 1000
        end_int = np.datetime64(eevent).astype('int64') * 1000
        a.axvline(start_int, ls='dotted', color='r', zorder=3)
        a.axvline(end_int, ls='dotted', color='r', zorder=3)
            
    annotations = ['a)', 'b)', 'c)', 'd)', 'e)', 'f)', 'g)', 'h)']
    for i, a in enumerate(axs):
        if i < 4:
             a.text(0.96, 0.8, annotations[i], transform=a.transAxes, weight='bold')
        else:
            a.text(0.96, 0.8, annotations[i], transform=a.transAxes, color='r', weight='bold')
            
    a.tick_params(axis='x', width=2, rotation=0)
    a.xaxis.label.set_visible(False)
    
    plt.suptitle(dt.datetime.strftime(start, '%Y-%m-%d') + \
                 ', DOY={:s}'.format(dt.datetime.strftime(start, '%j')), y=0.92)
    
    if outfolder:
        fig.savefig(outfolder + dt.datetime.strftime(start, '%Y-%jT%H:%M') + '.png')
        plt.close(fig)
    
    fig.savefig('20210329/jedi/' + dt.datetime.strftime(start, '%Y-%jT%H:%M') + '.png')
    plt.close(fig)

# --------------------------------------------------------------------
if __name__ == '__main__':
    
    for startevent, endevent in read_csv('auto_plasmoid_list_20210329_window10.csv'):
        
        try:
            print('Working on - {:}'.format(startevent))
            start = startevent - dt.timedelta(minutes=20)
            end   = endevent + dt.timedelta(minutes=20)

            if startevent.day != start.day:
                start = dt.datetime(startevent.year, startevent.month, startevent.day, 0, 0, 1)

            if endevent.day != end.day:
                end = dt.datetime(endevent.year, endevent.month, endevent.day, 23, 59, 59)

            magfile = glob.glob('data/mag/*{:s}ss*.hdf5'.format(startevent.strftime('%Y%j')))[0]
            magdata = MagData(filename=magfile)
            magdata.df = magdata.df.loc[start:end]

            jedifiles = glob.glob('data/jedi/*_CDR_{:s}_*.TAB'.format(startevent.strftime('%Y%j')))
            if not jedifiles:
                jediindex = IndexFile(dataset='jedi')
                jediindex.get_files_between_dates(
                    startevent, endevent, DoDownload=True, sid='ERSISP')
                jediindex.get_files_between_dates(
                    startevent, endevent, DoDownload=True, sid='ERSESP')

            jedifiles = glob.glob('data/jedi/*_CDR_{:s}_*.TAB'.format(startevent.strftime('%Y%j')))
            jedidatas = []
            for jedifile in jedifiles:
                jedidatas.append(JediData(filename=jedifile))
                print(jedifile)

            ion_data = []
            electron_data = []
            for jdata in jedidatas:
                jdata.df = jdata.df.loc[start:end]
                jdata.calc_broadband_flux()

                if jdata.metadata['particletype'].str.contains('ELECTRON').any():
                    electron_data.append(jdata)
                else:
                    ion_data.append(jdata)

            try:
                ion_data = merge_data_species(ion_data, species='ion')
            except IndexError:
                pass

            try:
                electron_data = merge_data_species(electron_data, species='electron')
            except IndexError:
                pass

            try:
                plot_data(magdata, ion_data=ion_data, electron_data=electron_data, start=start, end=end, sevent=startevent, eevent=endevent)
            except:
                print('Couldnt plot file - {:}'.format(startevent))
                pass

            print('--------')
        
        except:
            continue