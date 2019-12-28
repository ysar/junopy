
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

        print(self.type, self.columns, self.typeparticle, self.lowE, self.highE)


class JediData(object):

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


            assert (len(iEnergies) == \
                    len(low_E_limit) == \
                    len(high_E_limit) == \
                    len(particle_type))

            f.readline()

        self.FluxChannels = JediChannel(
            'FLUX', columns, low_E_limit, high_E_limit, particle_type, iEnergies)

        self.df = pd.read_csv(self.filename, names=columns, skiprows=5)

        self.df['time'] = self.df.apply(
            lambda x: dt.datetime.strptime(x['UTC'], '%Y-%m-%dT%H:%M:%S'),
            axis=1)

        self.df.set_index(['time'], inplace=True)





    def plot(self, ax, typeplot):
        """
        Type Plot -
                P_PitchAngle, P_EvsTime, Heavy_PitchAngle, Heavy_EvsTime
        """
        import matplotlib.colors as colors

        x = self.df.index
        y = self.channels['freq']
        X, Y = np.meshgrid(x, y)
        C = self.df[self.channels['names']].values.T

        mp = ax.pcolormesh(X, Y, C,
                           cmap='cividis',
                           norm=colors.LogNorm())

        ax.set_yscale('log')
        ax.set_ylabel(self.yvar + self.units)

        size = 0.02
        bbox = ax.get_position()
        cax = ax.figure.add_axes([bbox.x1, bbox.y0, 0.02, bbox.y1- bbox.y0])
        cbar = ax.figure.colorbar(mp, cax=cax)
        cbar.set_label(self.cbarlabel, rotation=270, labelpad=10.0)

        return ax
