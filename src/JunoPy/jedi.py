import os
import numpy as np
import pandas as pd
import datetime as dt
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
        cax = ax.figure.add_axes([bbox.x1, bbox.y0, 0.02, bbox.y1- bbox.y0])
        cbar = ax.figure.colorbar(mp, cax=cax)
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
            
            
    def create_new_energy_bins(self):
        """Create new energy bins for merged data."""
        
        if not self.isMerged:
            raise Exception('Should only be done for merged data.')
            
        chars = ['x', 'y', 'z']

        lowbin_list = []
        highbin_list = []

        n_bins = int(len(self.metadata['lowEbins_x']) / self.nSSD)

        # Save the first few bins
        energy_bins_lower  = self.metadata['lowEbins_x'].to_numpy()[:n_bins]
        energy_bins_higher = self.metadata['highEbins_x'].to_numpy()[:n_bins]
        
        # Continuously check for min/max
        for i in range(self.nfiles):
            for j in range(1, self.nSSD):

                lowE  = self.metadata['lowEbins_{:s}'.format(chars[i])].to_numpy() 
                highE = self.metadata['lowEbins_{:s}'.format(chars[i])].to_numpy()

                energy_bins_lower  = np.minimum(energy_bins_lower,   lowE[j * n_bins: (j + 1) * n_bins])
                energy_bins_higher = np.maximum(energy_bins_higher, highE[j * n_bins: (j + 1) * n_bins])

        return energy_bins_lower, energy_bins_higher
    
    
    
    
    def plot_pitchangle_scatter(self, ax):
        """
        Plot the pitchangle scatter for the current dataset on the ax.
        """
        chars = ['x', 'y', 'z']

        for j in range(self.nfiles):
            for i in range(self.nSSD):
                mp = ax.scatter(
                    self.df.index, 
                    self.df['Pitch_Angle_T{:d}_{:s}'.format(i, chars[j])],
                    c=self.df['BroadbandFlux_T{:d}_{:s}'.format(i, chars[j])],
                    s=0.5,
                    marker='s',
                    cmap='jet', 
                    norm=colors.LogNorm(vmin=10**4, vmax=10**6))

        ax.set_ylim([0, 180])
        ax.set_ylabel('PA ($^\circ$)')

        size = 0.02
        bbox = ax.get_position()
        cax = ax.figure.add_axes([bbox.x1, bbox.y0, 0.02, bbox.y1- bbox.y0])
        cbar = ax.figure.colorbar(mp, cax=cax)
        # cbar.set_label('counts/s/cm$^2$/str/keV', rotation=270, labelpad=10.0)

    
    def plot_pitchangle_spectra(self, ax, start=None, end=None, nticks=6, vmin=1e5, vmax=5e6, cmap='coolwarm'):
        """
        Plot the pitch angle spectra using hist2D . 
        """
        chars = ['x', 'y', 'z']
    
        if start and end:
            nbins = int((end - start).total_seconds() / 60)
        else:
            nbins = int((self.df.index[-1]- self.df.index[0]).total_seconds() / 60)

        x_pa = np.array([], dtype=object)
        y_pa = np.array([])
        weights = np.array([])

        for j in range(self.nfiles):
            for iax in range(self.nSSD):
                x_pa = np.concatenate((x_pa, self.df.index.values))
                y_pa = np.concatenate((y_pa, self.df['Pitch_Angle_T{:d}_{:s}'.format(iax, chars[j])].to_numpy()))
                weights = np.concatenate((weights, self.df['BroadbandFlux_T{:d}_{:s}'.format(iax, chars[j])].to_numpy())) 

        # c = pd.to_datetime(x_pa, unit='ns')

        mp = ax.hist2d(x_pa, y_pa, bins=[nbins, 18], weights=weights, cmin=0.5, cmap=cmap, vmin=vmin, vmax=vmax) # norm=colors.LogNorm(vmin=vmin, vmax=vmax))

        xticks = np.linspace(x_pa[0], x_pa[-1], 6)
        
        ax.set_ylim([0, 180])
        ax.set_ylabel('PA $^\circ$')
        ax.set_yticks([0, 90, 180])


        size = 0.02
        bbox = ax.get_position()
        cax = ax.figure.add_axes([bbox.x1, bbox.y0, 0.02, bbox.y1- bbox.y0])
        cbar = ax.figure.colorbar(mp[3], cax=cax)
        # cbar.set_label('counts/s/cm$^2$/str/keV', rotation=270, labelpad=10.0)

        return xticks
    
        
def merge_data_species(data, species):
    """Given various JediData objects for the same species, merge them."""
    
    data_merge = JediData()
    
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