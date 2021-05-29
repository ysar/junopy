"""
We will use this script to detect flux ropes. Let's start simple. 

Given a file - filename.hdf5

1. Read in the file.
2. Detect all Btheta crossings. (no just flux ropes but also O-lines)
3. Use a sliding 10 minute window to detect flux ropes. 

"""
import os
import re
import glob
import time
import numpy as np
import datetime as dt
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import pandas as pd
import copy

from matplotlib.ticker import AutoMinorLocator
from scipy.signal import argrelextrema

from JunoPy.mag import MagData, AutoFluxRope
from JunoPy.library import find_ticks_between


def plot_magdata_line(magdata, ts, te, saveplot=False):
    """Plot data between start and end."""
    
    # Find the ticks.
    timearray, tickstring = find_ticks_between(magdata.df.index[0], magdata.df.index[-1])
    
    # Create the figure.
    fig, axs = plt.subplots(4, 1, 
                            gridspec_kw={'hspace': 0, 'right':0.9, 'top':0.9}, 
                            figsize=(11, 8.5), 
                            dpi=300)
    
    axs[:] = magdata.plot(vars=['Br', 'Btheta', 'Bphi', 'Bmag'], axs=axs)
        
    for tstart, tend in zip(ts, te):
        df_line = magdata.df.loc[tstart:tend].copy()
        if df_line.loc[tstart, 'Btheta'] < 0:
            axs[1].plot(df_line.index, df_line['Btheta'].values, color='b')
        else:
            axs[1].plot(df_line.index, df_line['Btheta'].values, color='r')
    
    for a in axs:
        a.set_xticks(timearray)
        a.set_xticklabels(tickstring)
        a.set_xlim([timearray[0], timearray[-1]])
        a.label_outer()
        
    a.tick_params(axis='x', width=2, rotation=0)
    a.xaxis.label.set_visible(False)
    
    plt.suptitle(dt.datetime.strftime(ts[0], '%Y-%m-%d') + \
                ', DOY={:s}'.format(dt.datetime.strftime(ts[0], '%j')), y=0.92)
    
    if saveplot:
        fig.savefig('20210525/out_auto_fluxrope/' + dt.datetime.strftime(ts[0], '%Y-%j') + '_lines.png')
        plt.close(fig)



def plot_mva_composite(filename, start, end, saveplot=False):
    """Plot the MVA composite"""
    start_cut = start - dt.timedelta(minutes=7)
    end_cut   = end   + dt.timedelta(minutes=7)
    magdata = MagData(filename=filename)
    magdata.df = magdata.df.loc[start_cut:end_cut].copy()

    magdata_mva = MagData(filename=filename)
    magdata_mva.df = magdata_mva.df.loc[start:end].copy()
    try:
        magdata_mva.perform_mva()
    except np.linalg.LinAlgError:
        print(start, end)
        print(magdata_mva.df)
        raise Exception('Stopped.')
        
        
    fig = plt.figure(figsize=(11, 8.5), constrained_layout=False, dpi=300)

    gs1 = fig.add_gridspec(
        nrows=4, ncols=1, top=0.99, bottom=0.55, left=0.15, right=0.7, hspace=0)
    ax_magdata = list()
    for i in range(4):
        ax_magdata.append(fig.add_subplot(gs1[i, 0]))

    gs2 = fig.add_gridspec(
        nrows=3, ncols=1, top=0.45, bottom=0.1, left=0.15, right=0.66, hspace=0)
    ax_mvaB = list()
    for i in range(3):
        ax_mvaB.append(fig.add_subplot(gs2[i, 0]))

    gs3 = fig.add_gridspec(
        nrows=2, ncols=1, top=0.45, bottom=0.1, left=0.67, right=0.9, hspace=0)
    ax_hodo = list()
    for i in range(2):
        ax_hodo.append(fig.add_subplot(gs3[i, 0]))

    gs4 = fig.add_gridspec(
        nrows=1, ncols=1, top=0.9, bottom=0.55, left=0.71, right=0.9)
    ax_text = fig.add_subplot(gs4[0, 0])

    # 1. Plot the mag data +/- 10 min from epoch. ------------------------------------
    ax_magdata = magdata.plot(
        vars=['Br', 'Btheta', 'Bphi', 'Bmag'], axs=ax_magdata)
    timearray, tickstring = find_ticks_between(magdata.df.index[0], magdata.df.index[-1])
    ylabels = ['B$_r$', 'B$_\\theta$', 'B$_\phi$', '|B|']
    for i, a in enumerate(ax_magdata):
        a.set_xticks([])
        a.set_xticks(timearray)
        a.set_xticklabels(tickstring)
        a.set_ylabel(ylabels[i] + ' (nT)')
        a.set_xlim(magdata.df.index[0], magdata.df.index[-1])
        a.label_outer()
        if start and end:
            a.axvspan(start, end, alpha=0.2, color=[0.2, 0.2, 0.2])
    a.tick_params(axis='x', width=2, rotation=0)
    a.xaxis.label.set_visible(False)


    # 2. Plot the MVA B vectors. -------------------------------------------------------
    ax_mvaB = magdata_mva.plot(vars=['B1', 'B2', 'B3'], axs=ax_mvaB)
    timearray2, tickstring2 = find_ticks_between(magdata_mva.df.index[0], magdata_mva.df.index[-1])
    ylabels = ['B$_L$', 'B$_M$', 'B$_N$']
    for i, a in enumerate(ax_mvaB):
        a.set_xticks([])
        a.set_xticks(timearray2)
        a.set_xticklabels(tickstring2, fontdict={'horizontalalignment':'right'})
        a.set_xlim(magdata_mva.df.index[0], magdata_mva.df.index[-1])
        a.set_ylabel(ylabels[i] + ' (nT)')
        a.label_outer()

    a.tick_params(axis='x', width=2, rotation=0)
    a.xaxis.label.set_visible(False)

    dy = (ax_mvaB[-1].get_ylim()[1] - ax_mvaB[-1].get_ylim()[0]) / 2
    
    ax_mvaB[0].set_ylim([magdata_mva.df['B1'].mean() - dy,
                         magdata_mva.df['B1'].mean() + dy])
    ax_mvaB[1].set_ylim([magdata_mva.df['B2'].mean() - dy,
                         magdata_mva.df['B2'].mean() + dy])
    
    # 3. Plot the hodograms. -----------------------------------------------------------
    combo = [['B3', 'B1'],
             ['B3', 'B2']]
    Bmax = np.max(np.array(
            [magdata_mva.df['B1'].abs().max(), 
             magdata_mva.df['B2'].abs().max(), 
             magdata_mva.df['B3'].abs().max()]))  * 1.01
    
    
    for i, a in enumerate(ax_hodo):
        
        a.scatter(magdata_mva.df[combo[i][0]][0], magdata_mva.df[combo[i][1]][0], marker='o', color='b', s=50, label='Start')
        a.scatter(magdata_mva.df[combo[i][0]][-1], magdata_mva.df[combo[i][1]][-1], marker='x', color='r', s=50, label='End')
        a.plot(magdata_mva.df[combo[i][0]], magdata_mva.df[combo[i][1]],
               color='k',lw=0.4)
        
        a.set_xlabel(combo[i][0] + ' (nT)')
        a.set_ylabel(combo[i][1] + ' (nT)')
        a.set_aspect('equal', adjustable='datalim')
        a.xaxis.set_minor_locator(AutoMinorLocator())
        a.yaxis.set_minor_locator(AutoMinorLocator())
        a.grid('on', which='minor', lw=0.5, ls='dotted')
        a.grid(color='k', ls='dashed', lw=0.25)
        a.yaxis.tick_right()
        a.yaxis.set_label_position("right")
        a.label_outer()
    
    ax_hodo[0].legend(loc='upper right')
    ax_hodo[0].set_xlabel('B$_N$ (nT)')
    ax_hodo[0].set_ylabel('B$_L$ (nT)')
    ax_hodo[1].set_xlabel('B$_N$ (nT)')
    ax_hodo[1].set_ylabel('B$_M$ (nT)')
        
    # 4. Add Text ----------------------------------------------------------------------
    ax_text.axis('off')

    date = dt.datetime.strptime(start.strftime('%Y-%j'), '%Y-%j')

    text = []
    text.append(dt.datetime.strftime(date, '%Y-%m-%d') + ', DOY={:s}'.format(dt.datetime.strftime(date, '%j')))
    text.append('\n\n')
    text.append('MVA Output:\n--------\n\nCartesian Eigenvectors:\n')
    text.append('x$_L$ = ({:.2f}, {:.2f}, {:.2f})\n'.format(magdata_mva.x1[0], magdata_mva.x1[1], magdata_mva.x1[2]))
    text.append('x$_M$ = ({:.2f}, {:.2f}, {:.2f})\n'.format(magdata_mva.x2[0], magdata_mva.x2[1], magdata_mva.x2[2]))
    text.append('x$_N$ = ({:.2f}, {:.2f}, {:.2f})\n\n'.format(magdata_mva.x3[0], magdata_mva.x3[1], magdata_mva.x3[2]))
    text.append('--------\nRatio of eigenvalues:\n')
    text.append('$\lambda_N$ / $\lambda_M$ = {:.2f}\n'.format(magdata_mva.eigval[2] / magdata_mva.eigval[1]))
    text.append('$\lambda_M$ / $\lambda_L$ = {:.2f}\n'.format(magdata_mva.eigval[1] / magdata_mva.eigval[0]))
    text.append('$P = {:.2f}$\n'.format(magdata_mva.P))
    text = ''.join(text)
    ax_text.text(0.5, 0.5, text, horizontalalignment='center', verticalalignment='center',
                 bbox=dict(boxstyle="round",
                   ec=(0.9, 0.9, 0.9),
                   fc=(0.9, 0.9, 0.9),
                   ))
    
    if saveplot:
        plotfile = '20210525/mva/' + dt.datetime.strftime(start, '%Y-%jT%H:%M:%S')
        fig.savefig(plotfile)
        plt.close(fig)
        
        

# ========================================================================
if __name__ == '__main__':
    
    allfiles = glob.glob('data/mag/*.hdf5')
    nfiles = len(allfiles)
    
    with open('auto_plasmoid_list_20210525_window10.csv', 'w') as f:
        
        f.write('Start, End, Motion, PlasmoidType\n')
        
        for ifile, filename in enumerate(allfiles):
            
            print('{:d} / {:d} - {:s}'.format(ifile, nfiles, filename))
            
            magdata = MagData(filename=filename)
            
            # 1. Find all Btheta reversals.
            crossings = AutoFluxRope.find_btheta_crossings(magdata)
            
            if len(crossings) == 0:
                continue
                
            # 2. Find all plasmoid / flux rope candidates based on filters.
            ts, te = AutoFluxRope.find_plasmoid_candidates(magdata, crossings, window=10)
            
            if not ts and not te:
                continue
                
            # 3. Plot the plasmoids for the entire day.
            plot_magdata_line(magdata, ts, te, saveplot=True)
                
            # 3. Add times to the database.
            for tstart, tend in zip(ts, te):
                
                if tstart and tend:
                    
                    if magdata.df.loc[tstart, 'Btheta'] < 0:
                        motion = 'planetward'
                    else:
                        motion = 'tailward'
                        
                    magdata_mva = MagData(filename=filename)
                    magdata_mva.df = magdata_mva.df.loc[tstart:tend].copy()
                    magdata_mva.perform_mva()
                        
                    if np.abs(magdata_mva.df['B2'].median()) > np.abs(magdata_mva.df.loc[tstart, 'B2']) and \
                       np.abs(magdata_mva.df['B2'].median()) > np.abs(magdata_mva.df.loc[tend, 'B2']):
                        plasmoidtype = 'flux-rope'
                    else:
                        plasmoidtype = 'o-line'
                    
                    print(motion, plasmoidtype)
                    
                    f.write('{:s}, {:s}, {:s}, {:s}\n'.format(
                        tstart.strftime('%Y-%jT%H:%M:%S'), 
                        tend.strftime('%Y-%jT%H:%M:%S'),
                        motion,
                        plasmoidtype))

                    # 4. For each plasmoid event, perform MVA and save plot.
                    plot_mva_composite(filename, tstart, tend, saveplot=True)

                # raise Exception('Check plot.')