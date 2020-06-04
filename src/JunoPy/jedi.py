import os
import numpy as np
import pandas as pd
import datetime as dt


class JediChannel(object):

    def __init__(self, type, columns, lowE, highE, typeparticle, iEnergy):
        """Class to store infromation about the different energy channels"""
        self.type = type
        self.indices = [
            i for i, x in enumerate(columns) if (i in iEnergy and type in x)]
        self.columns = [columns[i] for i in self.indices]
        self.typeparticle = [typeparticle[i] for i in self.indices]
        self.lowE = np.array([lowE[i] for i in self.indices])
        self.highE = np.array([highE[i] for i in self.indices])

        # print(self.type, self.columns, self.typeparticle, self.lowE, self.highE)


class JediData(object):
    
    nSSD = 6
    
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

            iEnergies = [i for i, x in enumerate(low_E_limit) if x]
            # low_E_limit = np.array([float(x) for x in low_E_limit if x])

            high_E_limit = f.readline()
            high_E_limit = np.array([
                x.strip() for x in high_E_limit.split(',')]) # if x.strip()])

            particle_type = f.readline()
            particle_type = [
                x.strip() for x in particle_type.split(',')]#  if x.strip()]


            # assert (len(iEnergies) == \
            #         len(low_E_limit) == \
            #         len(high_E_limit) == \
            #         len(particle_type))

            f.readline()

        self.FluxChannels = JediChannel(
            'FLUX', columns, low_E_limit, high_E_limit, particle_type, iEnergies)

        self.df = pd.read_csv(self.filename, names=columns, skiprows=5)

        self.df['time'] = self.df.apply(
            lambda x: dt.datetime.strptime(x['UTC'], '%Y-%m-%dT%H:%M:%S'),
            axis=1)

        self.df.set_index(['time'], inplace=True)





    def plot_flux_spectra(self, ax, iTele=0):
        """
        Type Plot -
                P_PitchAngle, P_EvsTime, Heavy_PitchAngle, Heavy_EvsTime
        """
        import matplotlib.colors as colors

        self.df['endtime'] = self.df.apply(lambda x: x.name + dt.timedelta(seconds=x['Duration']), axis=1)

        x = np.zeros(len(self.df) + 1, dtype=object)
        x[:-1] = self.df.index
        x[-1] = self.df['endtime'][-1]

        columns = [x for x in self.FluxChannels.columns if 'T{:d}'.format(iTele) in x]
        freq_lower = self.FluxChannels.lowE[iTele * len(columns): (iTele + 1) * len(columns)]
        freq_upper = self.FluxChannels.highE[iTele * len(columns): (iTele + 1) * len(columns)]

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

        mp = ax.pcolormesh(x, y, c, cmap='jet', norm=colors.LogNorm(), vmin=100, vmax=10**6)
        ax.set_yscale('log')
        ax.set_ylim([20, y[-1]])
        ax.set_ylabel('Energy [keV]')

        size = 0.02
        bbox = ax.get_position()
        cax = ax.figure.add_axes([bbox.x1, bbox.y0, 0.02, bbox.y1- bbox.y0])
        cbar = ax.figure.colorbar(mp, cax=cax)
        cbar.set_label('counts/s/cm$^2$/str/keV', rotation=270, labelpad=10.0)
        
        ax.set_ylabel('JEDI e$^{-}$ Flux\n' + 'Sensor - ' + 
                      os.path.split(self.filename)[-1][4:7] + '\nkeV')
        
        return ax


    def plot_low_res_line_plots(self, ax, iTele=0):
        """
        Plot different lines showing the flux at different energies.
        """
        import matplotlib.pyplot as plt
        import matplotlib.pylab as pl
        from matplotlib.ticker import AutoMinorLocator

        columns = [x for x in self.FluxChannels.columns if 'T{:d}'.format(iTele) in x]
#         colors = pl.cm.jet(np.linspace(0, 1, len(columns)))

#         mp = self.df.plot(
#             y=columns,
#             ax=ax,
#             grid=True,
#             legend=False,
#             lw=0.4,
#             color=colors
#         )
        
        self.df.plot(
            y=columns,
            ax=ax,
            grid=True,
            style='.',
            legend=False,
            lw=0.1,
            colormap='jet',
            markersize=2.0
        )
        
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

    
    @staticmethod
    def calc_broadband_flux_static(x, lowE, highE):
        sum = 0
        for i, val in enumerate(x.values):
            if val > 0:
                sum += val * (float(upE[i]) - float(lowE[i]))
        return sum
    
        
    def calc_broadband_flux():
        """
        Given the JediData dataframe, calculate the integrated flux from all energy channels.
        """
        for iSSD in range(self.nSSD):
            
            columns = [x for x in self.FluxChannels.columns if 'T{:d}'.format(iSSD) in x]
            freqlower = self.FluxChannels.lowE[iSSD * len(columns) : (iSSD + 1) * len(columns)]
            frequpper = self.FluxChannels.highE[iSSD * len(columns) : (iSSD + 1) * len(columns)]
            
            self.df['BroadbandFlux_T{:d}'.format(iSSD)] = self.df.apply(
                lambda x: calc_broadband_flux(
                    x[columns], freqlower, frequpper), axis=1)
            
            
    
    def plot_pitch_angle_spectra(self, ax):
        """
        Plot the pitchangle spectra for the current dataset on the ax.
        """
        