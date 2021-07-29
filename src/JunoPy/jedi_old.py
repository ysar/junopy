import os
import numpy as np
import pandas as pd
import datetime as dt

from matplotlib import ticker
import matplotlib.colors as colors

def create_jedi_metadata(columns, lowE, highE, particletype):
    """Return metadata for the datafile."""

    df = pd.DataFrame()
    df['channels'] = columns
    df['indices'] = df.index.values
    
    # Filter only the FLUX columns
    df = df[df['channels'].str.contains("FLUX")]
    
    # Get the indices for the filtered columns
    indices = [int(x) for x in df['indices'].values]
    
    # Get other metadata for the filtered columns
    df['particletype'] = [particletype[x] for x in indices]
    df['lowEbins'] = [float(lowE[x]) for x in indices]
    df['highEbins'] = [float(highE[x]) for x in indices]
    
    df = df[~df['particletype'].str.contains("NONPARTICLE")]
    
    # Filter out the lowest energy channels for the electrons. 
    

    return df
    

class JediData(object):
    
    nSSD = 6
    isMerged = False
    nfiles = 1
    species = 'ion'
    
    def __init__(self, filename=None):
        """
        Class to read in JEDI data.
        """
        if filename:
            self.filename = filename
            self.parse_jedi_file()


    def parse_jedi_file(self):
        """
        Method to parse a JEDI data file (.csv) format
        """
        with open(self.filename, 'r') as f:

            columns = f.readline()
            columns = [x.strip() for x in columns.split(',')]

            low_E_limit = f.readline()
            low_E_limit = [x.strip() for x in low_E_limit.split(',')]

            high_E_limit = f.readline()
            high_E_limit = np.array([
                x.strip() for x in high_E_limit.split(',')]) # if x.strip()])

            particle_type = f.readline()
            particle_type = [
                x.strip() for x in particle_type.split(',')]#  if x.strip()]

            assert (len(low_E_limit) == \
                    len(high_E_limit) == \
                    len(particle_type))

            f.readline()
        
        self.metadata = create_jedi_metadata(columns, low_E_limit, high_E_limit, particle_type)
        
        self.df = pd.read_csv(self.filename, names=columns, skiprows=5)

        self.df['time'] = self.df.apply(
            lambda x: dt.datetime.strptime(x['UTC'], '%Y-%m-%dT%H:%M:%S'),
            axis=1)
        
        self.df = self.df.drop('UTC', 1)
        self.df.set_index(['time'], inplace=True)





    def plot_flux_spectra(self, ax, vmin=100, vmax=1e4, cmap='jet'):
        """
        Type Plot -
                P_PitchAngle, P_EvsTime, Heavy_PitchAngle, Heavy_EvsTime
        """
        import matplotlib.colors as colors

        self.df['endtime'] = self.df.apply(lambda x: x.name + dt.timedelta(seconds=x['Duration']), axis=1)

        x = np.zeros(len(self.df) + 1, dtype=object)
        x[:-1] = self.df.index
        x[-1] = self.df['endtime'][-1]
        

        columns = [x for x in self.metadata['channels']]

        freq_lower, freq_upper = self.create_new_energy_bins()
            

        y = np.zeros(len(freq_lower) * 2)
        j = 0
        for k in range(len(freq_lower)):
            y[j] = float(freq_lower[k])
            y[j + 1] = float(freq_upper[k])
            j += 2

        c0 = self.df[columns].values.T
        c = np.zeros((y.shape[0] - 1, x.shape[0] - 1))
        for i in range(c.shape[0]):
            for j in range(c.shape[1]):
                if i % 2 != 0:
                    c[i, j] = np.NaN
                else:
                    c[i, j] = c0[int(i/2), j]

        mp = ax.pcolormesh(x, y, c, cmap=cmap, norm=colors.LogNorm(vmin=vmin, vmax=vmax))

        ax.set_yscale('log')
        ax.set_ylim([20, y[-1]])
        ax.set_ylabel('Energy [keV]')

        size = 0.02
        bbox = ax.get_position()
        cax = ax.figure.add_axes([bbox.x1 + 0.01, bbox.y0, 0.02, 0.8 * (bbox.y1- bbox.y0)])
        cbar = ax.figure.colorbar(mp, cax=cax)
        cbar.set_ticks(ticker.LogLocator(numticks=4, subs='all'))
        # cbar.set_label('counts/s/cm$^2$/str/keV', rotation=270, labelpad=10.0)
        
        ax.set_ylabel('JEDI Flux\n{:s}\n'.format(self.species) + '\nkeV')
        
        return ax
    
    
    def plot_flux_spectra_corrected(self, ax, vmin=100, vmax=1e4, cmap='jet'):
        """
        Type Plot -
                P_PitchAngle, P_EvsTime, Heavy_PitchAngle, Heavy_EvsTime
        """
        import matplotlib.colors as colors

        self.df['endtime'] = self.df.apply(lambda x: x.name + dt.timedelta(seconds=x['Duration']), axis=1)

        x = np.zeros(len(self.df) + 1, dtype=object)
        x[:-1] = self.df.index
        x[-1] = self.df['endtime'][-1]
        

        columns = [x for x in self.metadata['channels']]

        freq_lower, freq_upper = self.create_new_energy_bins()
            

        y = np.zeros(len(freq_lower) * 2)
        j = 0
        for k in range(len(freq_lower)):
            y[j] = float(freq_lower[k])
            y[j + 1] = float(freq_upper[k])
            j += 2

        c0 = self.df[columns].values.T
        print(c0)
        c = np.zeros((y.shape[0] - 1, x.shape[0] - 1))
        for i in range(c.shape[0]):
            for j in range(c.shape[1]):
                if i % 2 != 0:
                    c[i, j] = np.NaN
                else:
                    c[i, j] = c0[int(i/2), j]

        mp = ax.pcolormesh(x, y, c, cmap=cmap, norm=colors.LogNorm(vmin=vmin, vmax=vmax))

        ax.set_yscale('log')
        ax.set_ylim([20, y[-1]])
        ax.set_ylabel('Energy [keV]')

        size = 0.02
        bbox = ax.get_position()
        cax = ax.figure.add_axes([bbox.x1 + 0.01, bbox.y0, 0.02, 0.8 * (bbox.y1- bbox.y0)])
        cbar = ax.figure.colorbar(mp, cax=cax)
        cbar.set_ticks(ticker.LogLocator(numticks=4, subs='all'))
        # cbar.set_label('counts/s/cm$^2$/str/keV', rotation=270, labelpad=10.0)
        
        ax.set_ylabel('JEDI Flux\n{:s}\n'.format(self.species) + '\nkeV')
        
        return ax


    def plot_low_res_line_plots_unmerged(self, ax, iTele=0):
        """
        Plot different lines showing the flux at different energies.
        """
        import matplotlib.pyplot as plt
        from matplotlib import cm
        from matplotlib.ticker import AutoMinorLocator
        
        for iTele in range(self.nSSD):
            columns = [x for x in self.FluxChannels.columns if 'T{:d}'.format(iTele) in x]

            start = 0.0
            stop = 1.0
            nlines = len(columns)
            cm_subsection = np.linspace(start, stop, nlines)
            colors = [ cm.jet(x) for x in cm_subsection ]

            for i, col in enumerate(columns):
                ax.plot(self.df.index, self.df[col], '.', lw=0.1, markersize=2.0, color=colors[i])
        
        ax.set_yscale('log')
        ax.xaxis.set_minor_locator(AutoMinorLocator())
        ax.grid('on', which='minor', lw='0.5', ls='dotted')
        ax.set_xlim([self.df.index[0], self.df.index[-1]])
        
        ax.set_ylabel('JEDI e$^{-}$ Flux\n' + 'Sensor - ' + 
                      os.path.split(self.filename)[-1][4:7] + '\ncounts/s\n/cm$^2$/str/keV')
        
        # ax.set_ylabel('Differential Intensity')

        # legendtext = ['{:s}-{:s} keV'.format(lower, upper) for lower, upper in zip(
        #     self.FluxChannels.lowE,
        #     self.FluxChannels.highE)]
        #
        # for i, text in enumerate(legendtext):
        #     ax.text(0, i / len(legendtext))
        #
        # ax.legend(, loc='lower left', borderaxespad=0.0, framealpha=1)

        return ax
    
        
    def calc_broadband_flux(self):
        """
        Given the JediData dataframe, calculate the integrated flux from all energy channels.
        Broadband flux is separate for each telescope/SSD.
        """
        for i in range(self.nSSD):
            
            # Get all channels containing the SSD ID. 
            channels = self.metadata['channels'].str.contains('T{:d}'.format(i))

            columns = self.metadata['channels'][channels]
            lowEbins = self.metadata['lowEbins'][channels]
            highEbins = self.metadata['highEbins'][channels]
            diffE = (highEbins - lowEbins).to_numpy()

            self.df['BroadbandFlux_T{:d}'.format(i)] = self.df.apply(
                lambda x: (x[columns].values * diffE).sum(), axis=1)
            
    
    def calc_narrowband_flux(self, low, high):
        """
        Given the JediData dataframe, calculate the integrated flux from all energy channels.
        Broadband flux is separate for each telescope/SSD.
        """
        for i in range(self.nSSD):
            
            # Get all channels containing the SSD ID. 
            channels = self.metadata['channels'].str.contains('T{:d}'.format(i))

            columns = self.metadata['channels'][channels]
            lowEbins = self.metadata['lowEbins'][channels]
            highEbins = self.metadata['highEbins'][channels]
            diffE = (highEbins - lowEbins).to_numpy()
            
            ifilter = (lowEbins > low & highEbins < high)
            print(ifilter)
            raise Exception()
            
            self.df['NarrowbandFlux_T{:d}_{:.2f}_{:.2f}'.format(i, low, high)] = self.df.apply(
                lambda x: (x[columns].values * diffE).sum(), axis=1)
            
            
    def create_new_energy_bins(self):
        """Create new energy bins for merged data."""
        
        lowbin_list = []
        highbin_list = []
        
        if self.isMerged:            
            chars = ['x', 'y', 'z']
            
            n_bins = int(len(self.metadata['lowEbins_x']) / self.nSSD)

            # Save the first few bins
            energy_bins_lower  = self.metadata['lowEbins_x'].to_numpy()[:n_bins]
            energy_bins_higher = self.metadata['highEbins_x'].to_numpy()[:n_bins]
        
            # Continuously check for min/max
            for i in range(self.nfiles):
                for j in range(1, self.nSSD):

                    lowE  = self.metadata['lowEbins_{:s}'.format(chars[i])].to_numpy() 
                    highE = self.metadata['highEbins_{:s}'.format(chars[i])].to_numpy()

                    energy_bins_lower  = np.minimum(energy_bins_lower,   lowE[j * n_bins: (j + 1) * n_bins])
                    energy_bins_higher = np.maximum(energy_bins_higher, highE[j * n_bins: (j + 1) * n_bins])
                    
        else:
            
            n_bins = int(len(self.metadata['lowEbins']) / self.nSSD)

            # Save the first few bins
            energy_bins_lower  = self.metadata['lowEbins'].to_numpy()[:n_bins]
            energy_bins_higher = self.metadata['highEbins'].to_numpy()[:n_bins]
        
            # Continuously check for min/max
            for j in range(1, self.nSSD):

                lowE  = self.metadata['lowEbins'].to_numpy() 
                highE = self.metadata['highEbins'].to_numpy()

                energy_bins_lower  = np.minimum(energy_bins_lower,   lowE[j * n_bins: (j + 1) * n_bins])
                energy_bins_higher = np.maximum(energy_bins_higher, highE[j * n_bins: (j + 1) * n_bins])
                    

        return energy_bins_lower, energy_bins_higher
    
    
    
    
    def plot_pitchangle_scatter(self, ax, size=0.5):
        """
        Plot the pitchangle scatter for the current dataset on the ax.
        """
        if self.nfiles > 1:
            chars = ['_x', '_y', '_z']
        else:
            chars = ['']

        for j in range(self.nfiles):
            for i in range(self.nSSD):
                mp = ax.scatter(
                    self.df.index, 
                    self.df['Pitch_Angle_T{:d}{:s}'.format(i, chars[j])],
                    c=self.df['BroadbandFlux_T{:d}{:s}'.format(i, chars[j])],
                    s=size,
                    marker='s',
                    cmap='jet', 
                    norm=colors.LogNorm(vmin=10**4, vmax=10**6))
        
        ax.set_ylim([0, 180])
        ax.set_ylabel('{:s}\n'.format(self.species) + 'PA $^\circ$')
        ax.set_yticks([0, 90, 180])

        size = 0.02
        bbox = ax.get_position()
        cax = ax.figure.add_axes([bbox.x1 + 0.01, bbox.y0, 0.02, 0.8 * (bbox.y1- bbox.y0)])
        cbar = ax.figure.colorbar(mp, cax=cax)
        # cbar.set_label('counts/s/cm$^2$/str/keV', rotation=270, labelpad=10.0)

    
    def plot_pitchangle_spectra(self, ax, start=None, end=None, nticks=6, vmin=1e5, vmax=5e6, cmap='coolwarm'):
        """
        Plot the pitch angle spectra .
        
        We need to use pcolormesh for better control.
        pcolormesh uses three main inputs:
            X[i, j] - X[i+1, j], array of N + 1, M + 1 points
            Y[i, j] - Y[i, j+1], array of N + 1, M + 1 points
            C[i, j], array of N, M points
            
        We can contruct X and Y on our own using the timedeltas and specific pitch angles.
        
        For C , we need to know the number of measurements within a given time period and pitch angle interval.  
        
        """
        if self.nfiles > 1:
            chars = ['_x', '_y', '_z']
        else:
            chars = ['']
    
        if start and end:
            n_times = int((end - start).total_seconds() / 60) 
            time_array = pd.date_range(
                start=start, end=end, periods=n_times + 1).to_pydatetime().tolist()
        else:
            n_times = int((self.df.index[-1] - self.df.index[0]).total_seconds() / 60) 
            time_array = pd.date_range(
                start=self.df.index[0], end=self.df.index[-1], periods=n_times + 1).to_pydatetime().tolist()    
                # Time array returns list
        
        n_angles = 18
        pitch_angle_array = np.linspace(0, 180, n_angles + 1)
        
        X, Y = np.meshgrid(time_array, pitch_angle_array)   # The size for the meshgrid is (N + 1) x (M + 1)
        
        C = np.zeros((n_angles, n_times))                   # There is a weird transpose that happens after meshgrid. 
        
        """
        Now we need to construct C. For this we need to loop over each time bin, and each pitch angle bin.
        We need to collect all data points from all files and SSDs which lie inside this 2D bin. 
        We will collect the number of data points (npoints) within the bin and the Broadband Flux (broadbandflux).
        """
        times_all = np.array([], dtype=object)
        pa_all = np.array([])
        flux_all = np.array([])

        for j in range(self.nfiles):
            for iax in range(self.nSSD):
                times_all = np.concatenate((times_all, self.df.index.values))
                pa_all = np.concatenate((pa_all, self.df['Pitch_Angle_T{:d}{:s}'.format(iax, chars[j])].to_numpy()))
                flux_all = np.concatenate((flux_all, self.df['BroadbandFlux_T{:d}{:s}'.format(iax, chars[j])].to_numpy()))
        
        
        times = np.zeros(len(times_all), dtype=object)                     # Quick fix for numpy / datetime incompatibilities.
        for i in range(len(times_all)):
            times[i] = dt.datetime.utcfromtimestamp(times_all[i] / 1e9)
            
            
        npoints_array = np.zeros((n_angles, n_times))
        bbandflux_array = np.zeros((n_angles, n_times))
        
        for i in range(n_times):
            
            pa_filter = pa_all[(times > time_array[i]) & (times <= time_array[i + 1])]
            flux_filter = flux_all[(times > time_array[i]) & (times <= time_array[i + 1])]
            
            for j in range(n_angles):
                
                flux = flux_filter[(pa_filter > pitch_angle_array[j]) & (pa_filter <= pitch_angle_array[j + 1])] 
                
                npoints_array[j, i] = len(flux)
                bbandflux_array[j, i] = np.sum(flux)
                
        
        npoints_array /= np.sum(npoints_array)
        
        C = bbandflux_array * 1.0 # / npoints_array
        
        vmin = 10**(np.log10(np.median(C)) - 1.01)
        vmax = 10**(np.log10(np.median(C)) + 1.01)
            
        
        mp = ax.pcolormesh(X, Y, C, cmap=cmap, norm=colors.LogNorm()) # vmin=vmin, vmax=vmax))

         # c = pd.to_datetime(x_pa, unit='ns')
            
        ax.set_ylim([0, 180])
        if 'ion' in self.species:
            ax.set_ylabel('H$^+$\n' + 'PA $^\circ$')
        else:
            ax.set_ylabel('{:s}\n'.format(self.species) + 'PA $^\circ$')
        ax.set_yticks([0, 90, 180])


        size = 0.02
        bbox = ax.get_position()
        cax = ax.figure.add_axes([bbox.x1 + 0.01, bbox.y0, 0.02, 0.8 * (bbox.y1- bbox.y0)])
        
        # cticks = np.array([vmin, vmax])
        cbar = ax.figure.colorbar(mp, cax=cax)# ,
                                  # format=ticker.LogFormatter(10, labelOnlyBase=False))
        
        cbar.set_ticks(ticker.LogLocator(numticks=3, subs='all'))
        
    
    def plot_pitchangle_spectra_narrowband(self, low, high, ax, start=None, end=None, nticks=6, vmin=1e5, vmax=5e6, cmap='coolwarm'):
        """
        Plot the pitch angle spectra .
        
        We need to use pcolormesh for better control.
        pcolormesh uses three main inputs:
            X[i, j] - X[i+1, j], array of N + 1, M + 1 points
            Y[i, j] - Y[i, j+1], array of N + 1, M + 1 points
            C[i, j], array of N, M points
            
        We can contruct X and Y on our own using the timedeltas and specific pitch angles.
        
        For C , we need to know the number of measurements within a given time period and pitch angle interval.  
        
        """
        
        # First we get the data corresponding to that channel
        # self.calc_narrowband_flux(low, high)
        
        self.df.apply(
                lambda x: (x[columns].values * diffE).sum(), axis=1)
        
    
        if self.nfiles > 1:
            chars = ['_x', '_y', '_z']
        else:
            chars = ['']
    
        if start and end:
            n_times = int((end - start).total_seconds() / 60) 
            time_array = pd.date_range(
                start=start, end=end, periods=n_times + 1).to_pydatetime().tolist()
        else:
            n_times = int((self.df.index[-1] - self.df.index[0]).total_seconds() / 60) 
            time_array = pd.date_range(
                start=self.df.index[0], end=self.df.index[-1], periods=n_times + 1).to_pydatetime().tolist()    
                # Time array returns list
        
        n_angles = 18
        pitch_angle_array = np.linspace(0, 180, n_angles + 1)
        
        X, Y = np.meshgrid(time_array, pitch_angle_array)   # The size for the meshgrid is (N + 1) x (M + 1)
        
        C = np.zeros((n_angles, n_times))                   # There is a weird transpose that happens after meshgrid. 
        
        """
        Now we need to construct C. For this we need to loop over each time bin, and each pitch angle bin.
        We need to collect all data points from all files and SSDs which lie inside this 2D bin. 
        We will collect the number of data points (npoints) within the bin and the Broadband Flux (broadbandflux).
        """
        times_all = np.array([], dtype=object)
        pa_all = np.array([])
        flux_all = np.array([])

        for j in range(self.nfiles):
            for iax in range(self.nSSD):
                times_all = np.concatenate((times_all, self.df.index.values))
                pa_all = np.concatenate((pa_all, self.df['Pitch_Angle_T{:d}{:s}'.format(iax, chars[j])].to_numpy()))
                flux_all = np.concatenate((flux_all, self.df['BroadbandFlux_T{:d}{:s}'.format(iax, chars[j])].to_numpy()))
        
        
        times = np.zeros(len(times_all), dtype=object)                     # Quick fix for numpy / datetime incompatibilities.
        for i in range(len(times_all)):
            times[i] = dt.datetime.utcfromtimestamp(times_all[i] / 1e9)
            
            
        npoints_array = np.zeros((n_angles, n_times))
        bbandflux_array = np.zeros((n_angles, n_times))
        
        for i in range(n_times):
            
            pa_filter = pa_all[(times > time_array[i]) & (times <= time_array[i + 1])]
            flux_filter = flux_all[(times > time_array[i]) & (times <= time_array[i + 1])]
            
            for j in range(n_angles):
                
                flux = flux_filter[(pa_filter > pitch_angle_array[j]) & (pa_filter <= pitch_angle_array[j + 1])] 
                
                npoints_array[j, i] = len(flux)
                bbandflux_array[j, i] = np.sum(flux)
                
        
        npoints_array /= np.sum(npoints_array)
        
        C = bbandflux_array * 1.0 # / npoints_array
        
        vmin = 10**(np.log10(np.median(C)) - 1.01)
        vmax = 10**(np.log10(np.median(C)) + 1.01)
            
        
        mp = ax.pcolormesh(X, Y, C, cmap=cmap, norm=colors.LogNorm(vmin=vmin, vmax=vmax))

         # c = pd.to_datetime(x_pa, unit='ns')
            
        ax.set_ylim([0, 180])
        if 'ion' in self.species:
            ax.set_ylabel('H$^+$\n' + 'PA $^\circ$')
        else:
            ax.set_ylabel('{:s}\n'.format(self.species) + 'PA $^\circ$')
        ax.set_yticks([0, 90, 180])


        size = 0.02
        bbox = ax.get_position()
        cax = ax.figure.add_axes([bbox.x1 + 0.01, bbox.y0, 0.02, 0.8 * (bbox.y1- bbox.y0)])
        
        cticks = np.array([vmin, vmax])
        cbar = ax.figure.colorbar(mp, cax=cax)# ,
                                  # format=ticker.LogFormatter(10, labelOnlyBase=False))
        
        cbar.set_ticks(ticker.LogLocator(numticks=3, subs='all'))
        

    
        
def merge_data_species(data, species):
    """Given various JediData objects for the same species, merge them."""

    data_merge = JediData()
    
    if len(data) == 1:
        data[0].nfiles = len(data)
        data[0].species = species
        return data[0]
    else:
        
        data_merge.df = pd.merge_ordered(data[0].df, data[1].df, on='time')
        data_merge.metadata = pd.merge_ordered(data[0].metadata, data[1].metadata, on='channels')

        # Merge fluxes from same energy channels. 
        for channel in data_merge.metadata['channels']:
            data_merge.df[channel] = data_merge.df[channel + '_x'] + data_merge.df[channel + '_y']
        
    # Set time as index
    data_merge.df.set_index(['time'], inplace=True)
    data_merge.df['Duration'] = np.maximum(data_merge.df['Duration_x'], data_merge.df['Duration_y'])
    
    # Modify in case we have three files
    if len(data) == 3:
        data_merge.df = pd.merge_ordered(data_merge.df, data[2].df, on='time', suffixes=(None, "_z"))
        data_merge.metadata = pd.merge_ordered(data_merge.metadata, data[2].metadata, on='channels', suffixes=(None, "_z"))
        
        # Merge fluxes from same energy channels. 
        for channel in data_merge.metadata['channels']:
            data_merge.df[channel] = data_merge.df[channel] + data_merge.df[channel +  '_z']

        data_merge.df['Duration'] = np.maximum(data_merge.df['Duration'], data_merge.df['Duration_z'])
        
    # Remove NaNs
    data_merge.df = data_merge.df.dropna()
    
    data_merge.isMerged = True
    data_merge.nfiles = len(data)
    data_merge.species = species
    
    return data_merge