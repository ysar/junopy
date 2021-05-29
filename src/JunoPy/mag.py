import os
import numpy as np
import pandas as pd
import datetime as dt
import h5py
import copy
import itertools


from scipy.special import jv
from lmfit import Parameters, minimize, report_fit
from matplotlib.ticker import AutoMinorLocator
from scipy.signal import argrelextrema

from JunoPy.library import R_JUPITER
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
    
    
    def calculate_bendback(self):
        """
        Calculate bendback angle alpha.
        Note : tan(alpha) / R = (B_phi / (R * B_R)) is often plotted.
        alpha is the angle between the magnetic field in the (R, phi) plane and the radial direction.
        """
        self.df['alpha'] = np.arctan2(self.df['Bphi'], self.df['Br']) * 180 / np.pi
    
    
    def plot(self, vars=None, axs=None, Bfit=None):
        """
        Plot magnetometer data. Returns an axes object.
        """
        
        if 'Btheta' in vars:
            var_string = {
                'Btheta': r'$B_\theta$',
                'Bphi': r'$B_\phi$',
                'Br': r'$B_r$',
                'Bmag': r'|B|',
            }
        
        if 'alpha' in vars:
            var_string.update({'alpha': r'$\alpha'})
        # print(var_string)
        
        for i, ax in enumerate(axs):
            
            ax.plot(self.df.index, self.df[vars[i]], color='k', lw=0.75)
            ax.grid()
            
            # ax.xaxis.set_minor_locator(AutoMinorLocator())
            ax.grid('on', which='minor', lw='0.5', ls='dotted')

            xmin, xmax, ymin, ymax = ax.axis()
            if ymin < 0 < ymax:
                ax.axhline(y=0, color='k', ls='dashed', lw=1)
                
            if Bfit:
                ax.plot(self.df.index, self.df['Bfit_' + chr(i + 120)], color='b', lw=0.7)
            
            if 'Btheta' in vars:
                ax.set_ylabel(var_string[vars[i]] + ' (nT)')
            elif 'alpha' in vars:
                ax.set_ylabel(var_string[vars[i]] + '$^\circ$')
            else:
                ax.set_ylabel(vars[i] + ' (nT)')
            # ax.set_xlim([self.df.index[0], self.df.index[-1]])

        return axs
    
    def plot_between(self, start, end, **kwargs):
        """
        Plot magnetometer data. Returns an axes object.
        """
        vars = kwargs.get('vars')
        axs = kwargs.get('axs')
        Bfit = kwargs.get('Bfit')
        
        for i, ax in enumerate(axs):
            
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
        
        if start and end:
            df = self.df.loc[start:end].copy()
            
        for i, ax in enumerate(axs):
            
            # Perform fft.
            N = len(df[vars[i]])
            rfft = fft.fft(df[vars[i]])
            rfreq = fft.fftfreq(len(df[vars[i]]))
            
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
    
    
    def perform_mva(self):
        """
        Perform minimum variance analysis (MVA).
        Input - MagData object containing pandas dataframe in self.df
        Ouput - MagData object with additional variables
                B1, B2, B3, eigenval, x1, x2, x3
        """
        self.x1 = np.zeros(3)
        self.x2 = np.zeros(3)
        self.x3 = np.zeros(3)
        self.eigval = np.zeros(3)
        self.df['B1'] = 0
        self.df['B2'] = 0
        self.df['B3'] = 0

        M = np.zeros((3, 3))
        for m, mcom in enumerate(['Bx', 'By', 'Bz']):
            for n, ncom in enumerate(['Bx', 'By', 'Bz']):
                M[m, n] = (self.df[mcom] * self.df[ncom]).mean() - \
                           self.df[mcom].mean() * self.df[ncom].mean()

        eigval, eigvec = np.linalg.eig(M)

        imax = np.argmax(eigval)
        imin = np.argmin(eigval)
        imid = list(set([0, 1, 2]) - set([imin, imax]))[0]

        self.eigval[0] = eigval[imin]
        self.eigval[1] = eigval[imid]
        self.eigval[2] = eigval[imax]
        self.x1 = eigvec[:, imin]
        self.x2 = eigvec[:, imid]
        self.x3 = eigvec[:, imax]

        if self.x3[-1] < 0:
            self.x1 *= -1
            self.x2 *= -1
            self.x3 *= -1

        self.df['B1'] = self.df.apply(
            lambda x: x['Bx'] * self.x1[0]
                      + x['By'] * self.x1[1]
                      + x['Bz'] * self.x1[2], axis=1)
        self.df['B2'] = self.df.apply(
            lambda x: x['Bx'] * self.x2[0]
                      + x['By'] * self.x2[1]
                      + x['Bz'] * self.x2[2], axis=1)
        self.df['B3'] = self.df.apply(
            lambda x: x['Bx'] * self.x3[0]
                      + x['By'] * self.x3[1]
                      + x['Bz'] * self.x3[2], axis=1)
        
        self.P = 100 * ((np.sqrt(self.eigval[2]) - np.sqrt(self.eigval[0]))
                     * (np.sqrt(self.eigval[2]) - np.sqrt(self.eigval[1]))
                     * (np.sqrt(self.eigval[1]) - np.sqrt(self.eigval[0]))
                     / np.sqrt(self.eigval[2]**3))


# ==================================================================
class ForceFreeFluxRope(object):
    
    def __init__(self, magfile, start, end):
        """Main function for force free fit."""
        
        magdata = MagData(filename=magfile)
        
        start_cut = start - dt.timedelta(minutes=5)
        end_cut   = end   + dt.timedelta(minutes=5)
        magdata.df = magdata.df.loc[start_cut:end_cut]
        
        magdata_mva = MagData(filename=magfile)
        magdata_mva.df = magdata_mva.df.loc[start:end]
        magdata_mva.perform_mva()

        magdata_mva.df['Bmag_saved'] = magdata_mva.df['Bmag']
        magdata_mva.df['Bx'] /= magdata_mva.df['Bmag']
        magdata_mva.df['By'] /= magdata_mva.df['Bmag']
        magdata_mva.df['Bz'] /= magdata_mva.df['Bmag']

        Ti, initParams, T0, mid = self.init_params(magdata_mva)

        eTraj = np.zeros(3)
        eTraj[0] = magdata_mva.df['X'][-1] - magdata_mva.df['X'][0]
        eTraj[1] = magdata_mva.df['Y'][-1] - magdata_mva.df['Y'][0]
        eTraj[2] = magdata_mva.df['Z'][-1] - magdata_mva.df['Z'][0]
        eTraj = eTraj / np.linalg.norm(eTraj)

        out = minimize(self.objective, initParams, args=(Ti, magdata_mva, eTraj), method='ampgo')

        fitParams = Parameters()
        fitParams.add('alpha', value=out.params['alpha'].value, vary=False)
        fitParams.add('B0', value=1.0, vary=True)
        fitParams.add('IF', value=out.params['IF'].value, vary=False)
        fitParams.add('thetaA', value=out.params['thetaA'].value, vary=False)
        fitParams.add('phiA', value=out.params['phiA'].value, vary=False)

        magdata_mva.df['Bx'] *= magdata_mva.df['Bmag']
        magdata_mva.df['By'] *= magdata_mva.df['Bmag']
        magdata_mva.df['Bz'] *= magdata_mva.df['Bmag']

        out2 = minimize(self.objective, fitParams, args=(Ti, magdata_mva, eTraj), method='cg')


        Bx_fit, By_fit, Bz_fit = self.model_flux_rope(out2.params, Ti, eTraj)
        B_fit = np.asarray([Bx_fit, By_fit, Bz_fit]).T
        B1_fit = np.dot(B_fit, np.asarray([magdata_mva.x1, magdata_mva.x2, magdata_mva.x3]).T)
        
        print(B1_fit.shape)
        
        self.magdata = magdata
        self.magdata_mva = magdata_mva
        self.magdata_mva.df['Bfit_x'] = B1_fit[:, 0]
        self.magdata_mva.df['Bfit_y'] = B1_fit[:, 1] 
        self.magdata_mva.df['Bfit_z'] = B1_fit[:, 2]
        self.fitreport_1 = out
        self.fitreport_2 = out2
    
    
    def model_flux_rope(self, params, t, eTraj):
        """Creates a flux rope based on parameters."""
        alpha = params['alpha']
        IF = params['IF']
        thetaA = params['thetaA']
        phiA = params['phiA']
        B0 = params['B0']

        x = 2 * t * np.sqrt(1 - IF**2)
        r = np.sqrt(IF**2 + x**2)

        eZ = np.zeros(3)
        eZ[0] = np.sin(thetaA) * np.cos(phiA)
        eZ[1] = np.sin(thetaA) * np.sin(phiA)
        eZ[2] = np.cos(thetaA)

        eY = np.cross(eZ, eTraj)
        eX = np.cross(eY, eZ)

        eTrajTangent = eTraj - np.dot(eTraj, eZ) * eZ

        B = np.zeros((len(x), 3))
        for i, xi in enumerate(x):
            gamma = np.arctan2(xi, IF)
            eTangent = eTrajTangent * np.cos(gamma) + \
                       np.cross(eZ, eTrajTangent) * np.sin(gamma) + \
                       eZ * np.dot(eZ, eTrajTangent) * (1 - np.cos(gamma))
                       # Rodrigues' Formula

            B_a = B0 * jv(0, alpha * r[i])
            B_t = B0 * jv(1, alpha * r[i])

            B[i, :] = B_a * eZ + B_t * eTangent

        return B[:, 0], B[:, 1], B[:, 2]

    
    def init_params(self, mvadata):
        """
        Initialize flux rope parameters.
        """
        thetaA = np.arccos(mvadata.x2[2])
        phiA = np.arctan2(mvadata.x2[1], mvadata.x2[0])

        eTraj = np.zeros(3)
        eTraj[0] = mvadata.df['X'][-1] - mvadata.df['X'][0]
        eTraj[1] = mvadata.df['Y'][-1] - mvadata.df['Y'][0]
        eTraj[2] = mvadata.df['Z'][-1] - mvadata.df['Z'][0]

        start = mvadata.df.index[0]
        end = mvadata.df.index[-1]
        mid = start + (end - start) / 2
        T0 = (end - start).total_seconds()

        Ti = list()
        for time, row in mvadata.df.iterrows():
            Ti.append((time - mid).total_seconds() / T0)

        Ti = np.asarray(Ti)

        initParams = Parameters()
        initParams.add('alpha', value=2.4048, min=-2.4048, max=2.4048)
        initParams.add('B0', value=1.0, vary=False)
        initParams.add('IF', value=0.5, min=0., max=1.)
        initParams.add('thetaA', value=thetaA, min=0., max=np.pi)
        initParams.add('phiA', value=phiA, min=-np.pi, max=np.pi)

        return Ti, initParams, T0, mid

    
    def objective(self, params, t, magdata, eTraj):
        """Objective function to minimize."""
        Bx, By, Bz = self.model_flux_rope(params, t, eTraj)

        residuals = np.zeros((3, len(Bx)))
        residuals[0, :] = Bx - magdata.df['Bx'].values
        residuals[1, :] = By - magdata.df['By'].values
        residuals[2, :] = Bz - magdata.df['Bz'].values

        # return np.abs(residuals).flatten()
        return np.sum(np.sum(residuals**2, axis=0)) / len(Bx)
        
        
# ==================================================================
class AutoFluxRope(object):
    
    def __init__(self):
        """This is a namespace class which contains methods necessary
        for the automated detection of flux ropes. It does not take 
        any arguments and should not be instantiated."""
        raise Exception('Namespace class, should not be instantiated.')
    
    
    @staticmethod
    def find_btheta_crossings(magdata, start=None, end=None):
        """
        Give a magnetometer dataset, find all locations where Btheta changes sign.
        At this point we don't care about whether the crossing is from North-to-south
        or vice-versa as there can be large changes superimposed onto the smaller ones.
        We also don't care atm about how large of a change there is in Btheta.
        """
        if start and end:
            df = magdata.df.loc[start:end].copy()
        else:
            df = magdata.df.copy()

        df['sign_btheta'] = np.sign(df['Btheta'].values)

        return df[np.abs(df['sign_btheta'].diff()) > 0].index 

    
    @staticmethod
    def find_plasmoid_start_end(magdata, crossing, window=5, n=10):
        """
        This is done is multiple steps
        
        For each indentified Btheta crossing, find the nearest maxima or minima. How?
        1. Cut away +/- 5 min from the Btheta crossing. 
        2. Identify all points where d(Btheta)/dt = 0 (or where it changes sign.)
        3. Store maxima and minima points separately. 
        
        Compare each set of maxima and minima for a particular crossing and choose the best fit.
        1. Create max, min pairs.
        2. Create a model dataset by interpolating the Btheta values between the two points.
        3. Calculate LS error for the model and data for each max, min pair.
        4. Select the pair with the least LS error. 
        
        Returns the best start, end times for the Btheta crossing event.
        """
        
        # Cut a chunk of data near the crossing, according to the window specified.
        
        # FILTER:
        # To save time and get the best data. We concentrate on certain radial ranges and local times
        if magdata.df.loc[crossing, 'R'] > 90 * R_JUPITER or magdata.df.loc[crossing, 'LT'] > 5:
            return None
        
        # FILTER 1: Select based on window size < 10 min.
        df = magdata.df.loc[crossing - dt.timedelta(minutes=window): 
                            crossing + dt.timedelta(minutes=window)].copy()
    
        # Identify all maxima and minima inside the analysis window and create pairs.
        maxima = df.iloc[argrelextrema(df['Btheta'].values, np.greater_equal, order=n)[0]]
        maxima = maxima[maxima['Btheta'] > 0].index
        maxima = list(set(maxima))
        minima = df.iloc[argrelextrema(df['Btheta'].values, np.less_equal, order=n)[0]]
        minima = minima[minima['Btheta'] < 0].index
        minima = list(set(minima))
        
        pairs = itertools.product(maxima, minima)
        
        # Obtain the start, stop pair which gives the minimum least-square error when
        # the data is approximated as a straight line.
        
        
        min_rsqr = 10
        rsqr = 100
        
        tstart_min = False
        tend_min = False
        
        for positive, negative in pairs:
            
            if not ((positive < crossing < negative) or (negative < crossing < positive)):   
                # We want the crossing to be inside the start, end pair interval.
                continue 
                
            
            # FILTER 3: based on the extreme values of Btheta (ignore all less than std. dev.)  
            
            if np.abs(df.loc[positive, 'Btheta'] - df.loc[negative, 'Btheta']) < 2:
                continue
                
            df_stddev = magdata.df.loc[crossing - dt.timedelta(minutes=5*window): 
                                       crossing + dt.timedelta(minutes=5*window)].copy()
            
            if np.abs(df.loc[positive, 'Btheta'] - df.loc[negative, 'Btheta']) < df_stddev['Btheta'].std():
                continue
            
            
            if positive < negative:
                start = positive
                end = negative
            elif negative < positive:
                start = negative
                end = positive
                
            # FILTER : We don't want any events which are 0.5 times or larger than the total window.
            if (end - start) > dt.timedelta(minutes=0.5*window) or (end - start) < dt.timedelta(seconds=10):
                continue
                
            # FILTER : We don't want a minor excursion into the negative values
            bthetapos = np.abs(df.loc[crossing, 'Btheta'] - df.loc[start, 'Btheta'])
            bthetaneg = np.abs(df.loc[end, 'Btheta'] - df.loc[crossing, 'Btheta'])
            if (bthetapos / bthetaneg > 6 or bthetaneg / bthetapos > 6):
                continue
            
                
            df_data = magdata.df.loc[start:end].copy()    # Create a model dataset.
            df_data['Btheta_model'] = np.NaN
            df_data.loc[start,'Btheta_model'] = df_data.loc[start, 'Btheta'] 
            df_data.loc[end,'Btheta_model'] = df_data.loc[end, 'Btheta']
            df_data['Btheta_model'] = df_data['Btheta_model'].interpolate() 
            
            
            avg_btheta = df_data['Btheta'].mean()
            df_data['ls_square_model'] = (df_data['Btheta'] - df_data['Btheta_model'])**2
            df_data['ls_square_mean'] = (df_data['Btheta'] - avg_btheta)**2
            
            rsqr = 1 - (df_data['ls_square_model'].sum() / df_data['ls_square_mean'].sum())
                # Calculate the least square
            
            
            # FILTER 2: r-value must be larger than 0.85 for good linear fit.
            if rsqr < 0.85:
                continue
                
            # FILTER 4: MVA based filters. 
            magdata_mva = copy.copy(magdata)
            magdata_mva.df = magdata_mva.df.loc[start:end].copy()
            if (magdata_mva.df['Br'].abs().mean() > 2 or magdata_mva.df['Br'].abs().median() > 2):
                continue
            magdata_mva.perform_mva()
            
            if np.abs(magdata_mva.x3[2]) < 0.8:
                continue
            
            if magdata_mva.P < 4.5:    # Degeneracy test based on Oleviera et al. 2020
                continue
                
            if (np.abs(magdata_mva.eigval[2] / magdata_mva.eigval[1]) < 3) or (np.abs(magdata_mva.eigval[1] / magdata_mva.eigval[0]) < 3):
                continue
                
            if (magdata_mva.df['B1'].abs().mean() > 2 or magdata_mva.df['B1'].abs().median() > 2):
                continue
                
                
            
                
            # Check for minima.    
            if rsqr < min_rsqr:
                tstart_min = start
                tend_min = end
                
        if tstart_min and tend_min:
            return tstart_min, tend_min
    
    
    @staticmethod
    def find_plasmoid_candidates(magdata, crossings, **kwargs):
        """Given a list of crossings, find all plasmoid candidates. Typically used in scripts."""
        
        starttimes = []
        endtimes = []
        for crossing in crossings:
            times = AutoFluxRope.find_plasmoid_start_end(magdata, crossing, **kwargs)
            if times:
                starttimes.append(times[0])
                endtimes.append(times[1])
            
        pairs = set(zip(starttimes, endtimes))
        
        starttimes = []
        endtimes = []
        for p in pairs:
            starttimes.append(p[0])
            endtimes.append(p[1])
            
        return starttimes, endtimes