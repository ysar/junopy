import os
import numpy as np
import pandas as pd
import datetime as dt
import h5py

from matplotlib.ticker import AutoMinorLocator

from JunoPy.spice import InitSpiceKernels, get_juno_position


class MagData(object):

    columns = [
        "year", "DOY", "hour", "minute", "second", "msec", "DecD",
        "Bx", "By", "Bz", "Brange", "X", "Y", "Z"
    ]

    def __init__(self, filename=None, ishdf=True):
        """
        Class to parse magnetometer data files. Stores data into a pandas
        dataframe.
        """
        if filename:
            if ishdf:
                self.filename = filename
                self.df = pd.read_hdf(filename)
            else:
                self.filename = filename
                self.parse_mag_file()

    def __repr__(self):
        """String representation for class. Writes the dataframe head"""
        return self.df.head().__repr__()
    
    def count_header_lines(self):
        """Count header lines in the magnetometer file."""
        with open(self.filename, 'r') as f:

            n_start = 0
            n_end = 0
            for i, line in enumerate(f):
                if "OBJECT" in line:
                    if "END_OBJECT" in line:
                        n_end += 1
                    else:
                        n_start += 1

                if n_end == n_start:
                    break

        return i + 1


    def parse_mag_file(self):
        """
        Parse magnetometer file into a pandas dataframe.
        """
        nheader = self.count_header_lines()

        self.df = pd.read_csv(self.filename,
                              names=self.columns,
                              delim_whitespace=True,
                              skiprows=nheader)

        self.df['time'] = self.df.apply(
            lambda x: dt.datetime.strptime(
                '{:04.0f}-{:03.0f}T{:02.0f}:{:02.0f}:{:02.0f}.{:03.0f}'.format(
                    x['year'],
                    x['DOY'],
                    x['hour'],
                    x['minute'],
                    x['second'],
                    x['msec']
                ), '%Y-%jT%H:%M:%S.%f'
            ), axis=1)

        self.df.set_index(['time'], inplace=True)
        self.df.drop(
            inplace=True,
            columns=['year', 'DOY', 'hour', 'minute', 'second', 'msec']
            )

        self.df['Bmag'] = np.sqrt(
            self.df['Bx']**2 + self.df['By']**2 + self.df['Bz']**2)


    def convert_to_spherical(self):
        """
        Calculate magnetic field values in spherical coordinates based
        on the XYZ location.
        """
        # Calculate R and LT for axes
        self.df['R'] = np.sqrt(
            self.df['X']**2 + self.df['Y']**2 + self.df['Z']**2)

        self.df['LT'] = 6 / 90 * \
            np.arctan2(self.df['Y'], self.df['X']) * 180 / np.pi + 12

        self.df['Latitude'] = np.arctan2(
            self.df['Z'],
            np.sqrt(self.df['X']**2 + self.df['Y']**2)
            ) * 180 / np.pi

        phi = np.arctan2(self.df['Y'], self.df['X'])
        theta = np.pi / 2 - self.df['Latitude'] * np.pi / 180

        self.df['Br'] = np.sin(theta) * np.cos(phi) * self.df['Bx'] + \
            np.sin(theta) * np.sin(phi) * self.df['By'] + \
            np.cos(theta) * self.df['Bz']

        self.df['Btheta'] = np.cos(theta) * np.cos(phi) * self.df['Bx'] + \
            np.cos(theta) * np.sin(phi) * self.df['By'] - \
            np.sin(theta) * self.df['Bz']

        self.df['Bphi'] = - np.sin(phi) * self.df['Bx'] + \
            np.cos(phi) * self.df['By']


    def plot(self, vars=None, axs=None, Bfit=None):
        """
        Plot magnetometer data. Returns an axes object.
        """
        for i, ax in enumerate(axs):

#             self.df.plot(
#                 y=vars[i],
#                 ax=ax,
#                 grid=True,
#                 color='k',
#                 legend=False,
#                 lw=0.75
#             )
            
            ax.plot(self.df.index, self.df[vars[i]], color='k', lw=0.75)
            ax.grid()
            
            ax.xaxis.set_minor_locator(AutoMinorLocator())
            ax.grid('on', which='minor', lw='0.5', ls='dotted')

            xmin, xmax, ymin, ymax = ax.axis()
            if ymin < 0 < ymax:
                ax.axhline(y=0, color='k', ls='dashed', lw=1)
                
            if Bfit is not None:
                ax.plot(self.df.index, Bfit[i], color='b', lw=0.7)

            ax.set_ylabel(vars[i] + ' (nT)')
            ax.set_xlim([self.df.index[0], self.df.index[-1]])

        return axs
    
    def plot_between(self, start, end, **kwargs):
        """
        Plot magnetometer data. Returns an axes object.
        """
        vars = kwargs.get('vars')
        axs = kwargs.get('axs')
        Bfit = kwargs.get('Bfit')
        
        for i, ax in enumerate(axs):

#             self.df.plot(
#                 y=vars[i],
#                 ax=ax,
#                 grid=True,
#                 color='k',
#                 legend=False,
#                 lw=0.75
#             )
            
            ax.plot(self.df.loc[start:end].index, self.df.loc[start:end][vars[i]], 
                    color=kwargs.get('color', 'k'), lw=0.75)
            ax.grid()
            
            ax.xaxis.set_minor_locator(AutoMinorLocator())
            ax.grid('on', which='minor', lw='0.5', ls='dotted')

            xmin, xmax, ymin, ymax = ax.axis()
            if ymin < 0 < ymax:
                ax.axhline(y=0, color='k', ls='dashed', lw=1)
                
            if Bfit is not None:
                ax.plot(self.df.index, Bfit[i], color='b', lw=0.7)

            ax.set_ylabel(vars[i] + ' (nT)')
            ax.set_xlim([self.df.index[0], self.df.index[-1]])

        return axs
    
    
    def plot_fft(self, start=None, end=None, **kwargs):
        """
        Plot the Fast fourier transform spectra.
        """
        import numpy.fft as fft
        
        vars = kwargs.get('vars')
        axs = kwargs.get('axs')
        xlim = kwargs.get('xlim')
        
        for i, ax in enumerate(axs):
            
            # Perform fft.
            N = len(self.df[vars[i]])
            rfft = fft.fft(self.df[vars[i]])
            rfreq = fft.fftfreq(len(self.df[vars[i]]))
            
            # Plot fft.
            ax.loglog(rfreq, rfft)
            # ax.set_xlim(xlim[0], xlim[1])
            
        return axs
        
            

    def plot_electron_cyc_freq(self, ax=None):
        """
        Calculate and plot the electron cyclotron frequency.
        """
        self.df['omega_ce'] = 1.60217e-19 * self.df['Bmag'] * 1e-9 \
                              / 9.1e-31 / 2 / np.pi
        self.df.plot(
            y='omega_ce',
            ax=ax,
            color='k',
            legend=False,
            lw=0.4
        )
        return ax
    
    
    def save_to_hdf5(self):
        """
        Save the current instance to an hdf5 file.
        """
        hdf5_filename = self.filename[:-4] + '.hdf5'
        self.df.to_hdf(hdf5_filename, 
                       os.path.split(self.filename)[-1][:-4], 
                       mode='w')
        os.remove(self.filename)
    
    
def mva(magdata):
    """
    Perform minimum variance analysis (MVA).
    Input - MagData object containing pandas dataframe in self.df
    Ouput - MagData object with additional variables
            B1, B2, B3, eigenval, x1, x2, x3
    """
    magdata.x1 = np.zeros(3)
    magdata.x2 = np.zeros(3)
    magdata.x3 = np.zeros(3)
    magdata.eigval = np.zeros(3)
    magdata.df['B1'] = 0
    magdata.df['B2'] = 0
    magdata.df['B3'] = 0

    M = np.zeros((3, 3))
    for m, mcom in enumerate(['Bx', 'By', 'Bz']):
        for n, ncom in enumerate(['Bx', 'By', 'Bz']):
            M[m, n] = (magdata.df[mcom] * magdata.df[ncom]).mean() - \
                      magdata.df[mcom].mean() * magdata.df[ncom].mean()

    eigval, eigvec = np.linalg.eig(M)

    imax = np.argmax(eigval)
    imin = np.argmin(eigval)
    imid = list(set([0, 1, 2]) - set([imin, imax]))[0]

    magdata.eigval[0] = eigval[imin]
    magdata.eigval[1] = eigval[imid]
    magdata.eigval[2] = eigval[imax]
    magdata.x1 = eigvec[:, imin]
    magdata.x2 = eigvec[:, imid]
    magdata.x3 = eigvec[:, imax]
    
    if magdata.x3[-1] < 0:
        magdata.x1 *= -1
        magdata.x2 *= -1
        magdata.x3 *= -1
        
    magdata.df['B1'] = magdata.df.apply(
        lambda x: x['Bx'] * magdata.x1[0]
                  + x['By'] * magdata.x1[1]
                  + x['Bz'] * magdata.x1[2], axis=1)
    magdata.df['B2'] = magdata.df.apply(
        lambda x: x['Bx'] * magdata.x2[0]
                  + x['By'] * magdata.x2[1]
                  + x['Bz'] * magdata.x2[2], axis=1)
    magdata.df['B3'] = magdata.df.apply(
        lambda x: x['Bx'] * magdata.x3[0]
                  + x['By'] * magdata.x3[1]
                  + x['Bz'] * magdata.x3[2], axis=1)

    return magdata


def force_free_fit(magfile, start, end):
    """
    Given a magdata where MVA has already been performed,
    fit a force free flux rope. 
    
    Steps for fitting.
    ------------------
    1. Unit normalize the magnetic field raw data.
    2. Perform MVA analysis.
    3. Given the following parameters - 
        
        r0 - radius of flux rope
        phiA, thetaA - orientation of fluxrope axis.
        
        Fit a force free model with the initial parameters being
        rinit - distance between start and end times
        phiMVA, thetaMVA - orientation of the x2 eigenvector from MVA.
        
    4. Calculate the modeled magnetic field in flux rope coords. 
    5. Convert the modeled field to MVA coordinates.
    """
    from scipy.special import jv
    from scipy.optimize import least_squares
    
    
    # Read the magnetometer file.
    magdata = MagData(filename=magfile)
    start_cut = start - dt.timedelta(minutes=5)
    end_cut   = end   + dt.timedelta(minutes=5)
    magdata.df = magdata.df.loc[start_cut:end_cut]
    
    # Unit normalization
    magdata.df['Bx'] = magdata.df['Bx'] / magdata.df['Bmag']
    magdata.df['By'] = magdata.df['By'] / magdata.df['Bmag']
    magdata.df['Bz'] = magdata.df['Bz'] / magdata.df['Bmag']
    magdata.df['Bmag_saved'] = magdata.df['Bmag'] 
    magdata.df['Bmag'] = 1.
    
    # Re-read into new object for MVA analysis.
    magdata_mva = MagData(filename=magfile)
    magdata_mva.df = magdata_mva.df.loc[start:end]
    magdata_mva.df['Bx'] = magdata_mva.df['Bx'] / magdata_mva.df['Bmag']
    magdata_mva.df['By'] = magdata_mva.df['By'] / magdata_mva.df['Bmag']
    magdata_mva.df['Bz'] = magdata_mva.df['Bz'] / magdata_mva.df['Bmag']
    magdata_mva.df['Bmag_saved'] = magdata.df['Bmag'] 
    magdata_mva.df['Bmag'] = 1.
    magdata_mva = mva(magdata_mva)
    
    
    # Calculate initial orientation for the flux rope axis.
    epoch = magdata.df.index[0]
    init = np.zeros(6)
    init[0] = np.arccos(magdata_mva.x2[2] / np.sqrt(magdata_mva.x2[0]**2 + magdata_mva.x2[1]**2))
    init[1] = np.arctan2(magdata_mva.x2[1], magdata_mva.x2[0])
    init[2] = (start + (end - start) / 2 - epoch).total_seconds()
    init[3] = ((end - start) / 2).total_seconds()
    init[4] = 0.2
    init[5] = 0.4
    
    shat = np.zeros(3)
    index_end = magdata.df.index.get_loc(end, method='nearest')
    index_start = magdata.df.index.get_loc(start, method='nearest')
        
    shat[0] = magdata.df['X'][index_end] - magdata.df['X'][index_start]
    shat[1] = magdata.df['Y'][index_end] - magdata.df['Y'][index_start]
    shat[2] = magdata.df['Z'][index_end] - magdata.df['Z'][index_start]
    shat = shat / np.sqrt(np.sum(shat**2))
    
    shat_xy = shat
    shat_xy[2] = 0
    shat_xy = shat_xy / np.sqrt(np.sum(shat_xy**2))
     
    # Force-free flux rope fitting.
    def calc_cost_function(x, magdata, start, end):
        """x contains the following parameters
            
            x[0] - thetaA
            x[1] - phiA
            x[2] - t0 (in units of seconds from epoch)
            x[3] - T0 (seconds)
            x[4] - Y0 (seconds)
            x[5] - alpha (units?)
        """
        thetaA = x[0]
        phiA = x[1]
        zhat = np.zeros(3)
        zhat[0] = np.sin(thetaA) * np.cos(phiA)
        zhat[1] = np.sin(thetaA) * np.sin(phiA)
        zhat[2] = np.cos(thetaA)
        
        yhat = np.cross(zhat, shat)
        xhat = np.cross(yhat, zhat)
        
        tstart = x[2] - x[3]
        tend = x[2] + x[3]
        nPoints = 0
        xisqr = 0
        for i, dpoint in magdata.df.loc[epoch + dt.timedelta(seconds=tstart):
                                        epoch + dt.timedelta(seconds=tend)].iterrows():
                
            ti = (i - epoch).total_seconds()
            Ti = np.sqrt(ti**2 + x[4]**2)
            BAi = jv(0, x[5] * Ti)
            BTi = jv(1, x[5] * Ti)

            BAi_xyz = BAi * zhat
            BTi_xyz = BTi * shat_xy
            Bi_xyz = BAi_xyz + BTi_xyz
            
            xisqr += (Bi_xyz[0] - dpoint['Bx'])**2 + \
                     (Bi_xyz[1] - dpoint['By'])**2 + \
                     (Bi_xyz[2] - dpoint['Bz'])**2
                
        
            nPoints += 1

        xisqr = xisqr / nPoints
        
        return xisqr
            

    # Force free flux rope fitting
    print(init)
    res = least_squares(calc_cost_function, init, args=(magdata, start, end))
        
    # Once fitted, we create a new data set with the 'time, Bx, By, Bz' for 
    # the fitted flux rope.
    
        
    tstart = res.x[2] - res.x[3]
    tend = res.x[2] + res.x[3]
    thetaA = res.x[0]
    phiA = res.x[1]
    zhat = np.zeros(3)
    zhat[0] = np.sin(thetaA) * np.cos(phiA)
    zhat[1] = np.sin(thetaA) * np.sin(phiA)
    zhat[2] = np.cos(thetaA)
    yhat = np.cross(zhat, shat)
    xhat = np.cross(yhat, zhat)
        
    
    nPoints = len(magdata.df.loc[epoch + dt.timedelta(seconds=tstart):
                                 epoch + dt.timedelta(seconds=tend)])
    magdata_mva.df_fr = pd.DataFrame(columns=['time', 'Bx', 'By', 'Bz'], index=range(nPoints))
    
    ipoint = 0
    for i, dpoint in magdata.df.loc[epoch + dt.timedelta(seconds=tstart):
                                    epoch + dt.timedelta(seconds=tend)].iterrows():
        time = i
        ti = (i - epoch).total_seconds()
        Ti = np.sqrt(ti**2 + res.x[4]**2)
        BAi = jv(0, res.x[5] * Ti)
        BTi = jv(1, res.x[5] * Ti)
        BAi_xyz = BAi * zhat
        BTi_xyz = BTi * shat_xy
        Bi_xyz = BAi_xyz + BTi_xyz
        
        magdata_mva.df_fr['time'][ipoint] = time 
        magdata_mva.df_fr['Bx'][ipoint] = Bi_xyz[0]
        magdata_mva.df_fr['By'][ipoint] = Bi_xyz[1]
        magdata_mva.df_fr['Bz'][ipoint] = Bi_xyz[2]
        ipoint += 1
    
    
    magdata_mva.df_fr.set_index(['time'], inplace=True)
    
    return res, magdata, magdata_mva
    