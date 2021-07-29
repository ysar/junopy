import os
import re
import glob
import numpy as np
import datetime as dt
import pandas as pd
import matplotlib.pyplot as plt

from matplotlib.ticker import AutoMinorLocator
from scipy.special import jv
from lmfit import Parameters, minimize, report_fit

from JunoPy.mag import MagData, ForceFreeFluxRope
from JunoPy.library import find_ticks_between
from JunoPy.spice import InitSpiceKernels, get_juno_position




def plot_mva_composite(FFFR, start, end, savefile=False):

    fig = plt.figure(figsize=(11, 8.5), constrained_layout=False, dpi=300)
    plt.rcParams['figure.facecolor'] = 'white'
    plt.rcParams['axes.facecolor'] = 'white'

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
    ax_magdata[:] = FFFR.magdata.plot(
        vars=['Br', 'Btheta', 'Bphi', 'Bmag'], axs=ax_magdata[:])
    
    timearray, tickstring = find_ticks_between(FFFR.magdata.df.index[0], FFFR.magdata.df.index[-1])
    ylabels = ['B$_r$ (nT)', 'B$_\\theta$ (nT)', 'B$_\phi$ (nT)', '|B| (nT)', 'Freq (Hz)']
    for i, a in enumerate(ax_magdata):
        a.set_xticks([])
        a.set_xticks(timearray)
        a.set_xticklabels(tickstring)
        a.set_ylabel(ylabels[i])
        a.set_xlim(FFFR.magdata.df.index[0], FFFR.magdata.df.index[-1])
        a.label_outer()
        
            # a.axvspan(start, end, alpha=0.2, color=[0.2, 0.2, 0.2])
    a.tick_params(axis='x', width=2, rotation=0)
    a.xaxis.label.set_visible(False)

    for a in ax_magdata:
        if start and end:
            a.axvline(start, color='magenta', ls='dashed')
            a.axvline(end, color='magenta', ls='dashed')
            
    # 2. Plot the MVA B vectors. -------------------------------------------------------
    ax_mvaB = FFFR.magdata_mva.plot(vars=['B1', 'B2', 'B3'], Bfit=True, axs=ax_mvaB)
    timearray2, tickstring2 = find_ticks_between(FFFR.magdata_mva.df.index[0], 
                                                 FFFR.magdata_mva.df.index[-1])
    
    ylabels = ['B$_N$', 'B$_M$', 'B$_L$']
    for i, a in enumerate(ax_mvaB):
        a.set_xticks([])
        a.set_xticks(timearray2)
        a.set_xticklabels(tickstring2, fontdict={'horizontalalignment':'right'})
        a.set_xlim(FFFR.magdata_mva.df.index[0], FFFR.magdata_mva.df.index[-1])
        a.set_ylabel(ylabels[i] + ' (nT)')
        a.label_outer()
        
    ax_mvaB[0].legend(['Data', 'Model'], loc='upper left')
    a.tick_params(axis='x', width=2, rotation=0)
    a.xaxis.label.set_visible(False)

    dy = (ax_mvaB[-1].get_ylim()[1] - ax_mvaB[-1].get_ylim()[0]) / 2
    
    ax_mvaB[0].set_ylim([FFFR.magdata_mva.df['B1'].mean() - dy,
                         FFFR.magdata_mva.df['B1'].mean() + dy])
    ax_mvaB[1].set_ylim([FFFR.magdata_mva.df['B2'].mean() - dy,
                         FFFR.magdata_mva.df['B2'].mean() + dy])
    
    # 3. Plot the hodograms. -----------------------------------------------------------
    combo = [['B3', 'B1'],
             ['B3', 'B2']]
    Bmax = np.max(np.array(
            [FFFR.magdata_mva.df['B1'].abs().max(), 
             FFFR.magdata_mva.df['B2'].abs().max(), 
             FFFR.magdata_mva.df['B3'].abs().max()]))  * 1.01
    
    
    for i, a in enumerate(ax_hodo):
        
        a.scatter(FFFR.magdata_mva.df[combo[i][0]][0], FFFR.magdata_mva.df[combo[i][1]][0], 
                  marker='o', color='b', s=50, label='Start')
        
        a.scatter(FFFR.magdata_mva.df[combo[i][0]][-1], FFFR.magdata_mva.df[combo[i][1]][-1], 
                  marker='x', color='r', s=50, label='End')
        
        a.plot(FFFR.magdata_mva.df[combo[i][0]], FFFR.magdata_mva.df[combo[i][1]],
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
    
    ax_hodo[0].legend(loc='upper left')
    ax_hodo[0].set_xlabel('B$_L$ (nT)')
    ax_hodo[0].set_ylabel('B$_N$ (nT)')
    ax_hodo[1].set_xlabel('B$_L$ (nT)')
    ax_hodo[1].set_ylabel('B$_M$ (nT)')
        
    # 4. Add Text ----------------------------------------------------------------------
    ax_text.axis('off')

    date = dt.datetime.strptime(start.strftime('%Y-%j'), '%Y-%j')

    text = []
    text.append(dt.datetime.strftime(date, '%Y-%m-%d') + 
                ', DOY={:s}'.format(dt.datetime.strftime(date, '%j')))
    
    text.append('\n\n')
    text.append('MVA Output:\n--------\n\nCartesian Eigenvectors:\n')
    text.append('N = ({:.2f}, {:.2f}, {:.2f})\n'.format(FFFR.magdata_mva.x1[0], 
                                                        FFFR.magdata_mva.x1[1], 
                                                        FFFR.magdata_mva.x1[2]))
    
    text.append('M = ({:.2f}, {:.2f}, {:.2f})\n'.format(FFFR.magdata_mva.x2[0], 
                                                        FFFR.magdata_mva.x2[1], 
                                                        FFFR.magdata_mva.x2[2]))
    
    text.append('L = ({:.2f}, {:.2f}, {:.2f})\n\n'.format(FFFR.magdata_mva.x3[0], 
                                                          FFFR.magdata_mva.x3[1], 
                                                          FFFR.magdata_mva.x3[2]))
    text.append('--------\nRatio of eigenvalues:\n')
    text.append('$\lambda_N$ / $\lambda_M$ = {:.2f}\n'.format(FFFR.magdata_mva.eigval[2] / 
                                                              FFFR.magdata_mva.eigval[1]))
    
    text.append('$\lambda_M$ / $\lambda_L$ = {:.2f}\n'.format(FFFR.magdata_mva.eigval[1] / 
                                                              FFFR.magdata_mva.eigval[0]))
    
    text.append('\n')
    text.append('--------\n')
    text.append('Force-free flux-rope fit:\n')
    text.append('\n')
    text.append('Step 1 $\chi^2_r$ = {:.2f} nT$^2$\n'.format(FFFR.fitreport_1.redchi))
    text.append('Step 2 $\chi^2_r$ = {:.2f} nT$^2$\n'.format(FFFR.fitreport_2.redchi))
    text.append('$\\alpha$ = {:.2f}, '.format(FFFR.fitreport_1.params['alpha'].value))
    text.append('IF = {:.2f}\n'.format(FFFR.fitreport_1.params['IF'].value))
    text.append('$B_0$ = {:.2f} nT\n'.format(FFFR.fitreport_2.params['B0'].value))
    
    thetaA = FFFR.fitreport_1.params['thetaA'].value
    phiA = FFFR.fitreport_1.params['phiA'].value
    
    axial_dir = np.zeros(3)
    axial_dir[0] = np.sin(thetaA) * np.cos(phiA)
    axial_dir[1] = np.sin(thetaA) * np.sin(phiA)
    axial_dir[2] = np.cos(thetaA)
    
    if np.dot(axial_dir, FFFR.magdata_mva.x2) < 0:
        axial_dir *= -1
        
    text.append('$x_A$ = ({:.2f}, {:.2f}, {:.2f})'.format(
        axial_dir[0],
        axial_dir[1],
        axial_dir[2]))
    
    text = ''.join(text)
    
    ax_text.text(0.65, 0.6, text, horizontalalignment='center', verticalalignment='center',
                 bbox=dict(boxstyle="round",
                   ec=(0.9, 0.9, 0.9),
                   fc=(0.9, 0.9, 0.9),
                   ))
    
    ax_magdata[0].text(0.96, 0.8, 'a)', transform=ax_magdata[0].transAxes, weight='bold')
    ax_magdata[1].text(0.96, 0.8, 'b)', transform=ax_magdata[1].transAxes, weight='bold')
    ax_magdata[2].text(0.96, 0.8, 'c)', transform=ax_magdata[2].transAxes, weight='bold')
    ax_magdata[3].text(0.96, 0.8, 'd)', transform=ax_magdata[3].transAxes, weight='bold')
    ax_mvaB[0].text(0.96, 0.8, 'e)', transform=ax_mvaB[0].transAxes, weight='bold')
    ax_mvaB[1].text(0.96, 0.8, 'f)', transform=ax_mvaB[1].transAxes, weight='bold')
    ax_mvaB[2].text(0.96, 0.8, 'g)', transform=ax_mvaB[2].transAxes, weight='bold')
    ax_hodo[0].text(0.9, 0.8, 'h)', transform=ax_hodo[0].transAxes, weight='bold')
    ax_hodo[1].text(0.9, 0.8, 'i)', transform=ax_hodo[1].transAxes, weight='bold')
    
    if savefile:
        plotfile = 'out_mva_forcefree/' + dt.datetime.strftime(start, '%Y-%jT%H:%M:%S') + '_fitted'
        fig.savefig(plotfile + '.png', transparent=False, bbox_inches='tight')
        plt.close(fig)
        
        
def read_csv(filename):
    
    start_times = []
    end_times = []
    with open(filename, 'r') as f:
        
        f.readline()
        for line in f:
            linesplit = line.split(',')
            start_time = dt.datetime.strptime(linesplit[0].strip(), '%Y-%jT%H:%M:%S')
            end_time = dt.datetime.strptime(linesplit[1].strip(), '%Y-%jT%H:%M:%S')
            direction = linesplit[2]
            plasmoidtype = linesplit[3]
            
            yield start_time, end_time, direction, plasmoidtype
        
        
# ---------------------------------------------------------
if __name__ == '__main__':
    
    filename = 'auto_plasmoid_list_20210525_window10.csv'
    
    
    df_fit = pd.DataFrame(columns=['redchi1', 'redchi2', 'alpha', 'IF', 'B0', 'dtheta'])
    
    
    
    for start, end, direction, plasmoidtype in read_csv(filename):
    
    # if 'o-line' in plasmoidtype:
    #    continue
        
        magfile = glob.glob('data/mag/*{:s}ss*.hdf5'.format(start.strftime('%Y%j')))[0]
        # magdata = MagData(filename=magfile)
        
        FFFR = ForceFreeFluxRope(magfile, start, end)

        if 'flux-rope' in plasmoidtype:
            
            print('Plotting file - {:s}'.format(start.strftime('%Y-%jT%H:%M')))
            
            plot_mva_composite(FFFR, start, end, savefile=True)
            
            
            thetaA = FFFR.fitreport_1.params['thetaA'].value
            phiA = FFFR.fitreport_1.params['phiA'].value
            axial_dir = np.zeros(3)
            axial_dir[0] = np.sin(thetaA) * np.cos(phiA)
            axial_dir[1] = np.sin(thetaA) * np.sin(phiA)
            axial_dir[2] = np.cos(thetaA)

            if np.dot(axial_dir, FFFR.magdata_mva.x2) < 0:
                axial_dir *= -1
            
            df_fit = df_fit.append({
                'redchi1': FFFR.fitreport_1.redchi,
                'redchi2': FFFR.fitreport_2.redchi,
                'alpha'  : FFFR.fitreport_1.params['alpha'].value,
                'IF'     : FFFR.fitreport_1.params['IF'].value,
                'B0'     : FFFR.fitreport_2.params['B0'].value,
                'dtheta' : np.arccos(np.dot(axial_dir, FFFR.magdata_mva.x2)) * 180 / np.pi
            }, ignore_index=True)

            df_fit.to_csv('fffr_20210608.csv')