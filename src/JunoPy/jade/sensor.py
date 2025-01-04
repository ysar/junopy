
import numpy as np
import pandas as pd
import datetime as dt
from copy import deepcopy
from scipy.constants import m_p, m_e, e

from JunoPy import mag
from .read import JadeBinaryData
from .calib import get_sensor_constants
from .calib import calc_anode_correction_F
from .calib import calc_deflection_correction_D
from .calib import calc_time_correction_tau
from .calib import calc_simulation_geom_factor
from .calib import calc_sensor_dependence_exp


def _create_energy_bins(energyArr):
    energyBins = np.zeros(len(energyArr) + 1)
    energyBins[1:-1] = np.sqrt(energyArr[:-1] * energyArr[1:])
    energyBins[0] = energyBins[1]**2 / energyBins[2]
    energyBins[-1] = energyBins[-2]**2 / energyBins[-3]
    return energyBins


class SensorData(object):

    def __init__(self, df, id):

        self.id = id

        if self.id == 104:
            self.name = "E060"
        elif self.id == 106:
            self.name = "E180"
        elif self.id == 131:
            self.name = "Species-3"
        elif self.id == 132:
            self.name = "Species-4"
        elif self.id == 133:
            self.name = "Species-5"

        self.times = df.index.to_pydatetime()

        self.data_raw = JadeBinaryData.create_dict_from_df(df)


    def calc_geometric_factor_electrons(self, isHighRes=True):

        delta = self.data_raw['DIM2_ELEVATION']
        energies = self.data_raw['DIM1_E']
        times = self.times

        sensorConst = get_sensor_constants(sensorId=self.id)

        if isHighRes:
            F = calc_anode_correction_F(sensorConst)[:, np.newaxis, np.newaxis] # shape (16, nE, nTime)
            D = calc_deflection_correction_D(delta, sensorConst)                # shape (16, nE, nTime)
        else:
            F = D = 1.

        tau = np.zeros(len(times))
        for i in range(len(times)):
            tau[i] = calc_time_correction_tau(times[i], sensorConst)
        
        tau = tau[np.newaxis, np.newaxis, :]                        # shape (16, nE, nTime)

        gSim = calc_simulation_geom_factor(energies)
        x = calc_sensor_dependence_exp(sensorConst, energies)

        # Energy Geometric Factor in units of cm2 sr eV/eV
        self.geomFactorEff = gSim * x * tau * D * F 
    

    def calc_geometric_factor_ion(self):
        energies = self.data_raw['DIM1_E']
        self.geomFactorEff = 6.84e-5 \
                + (2.52e-6 + np.log10(energies)) \
                + (6.84e-5 * energies**(-0.28))
        
        if self.id == 105 or self.id == 121:
            self.geomFactorEff *= 12
    

    def convert_units(self, inUnits='countrate', outUnits='dnf', m=None, q=1):

        toConvert = ['DATA', 'DATA_SIGMA', 'BACKGROUND', 'BACKGROUND_SIGMA']

        if 'countrate' in inUnits:
            assert((self.data_raw['DATA_UNITS'] == 2).all())
            if 'def' in outUnits:
                for k in toConvert:
                    self.data_raw[k] /= self.geomFactorEff    # cm-2 s-1 sr-1
                self.data_raw['DATA_UNITS'][:] = 6

            elif 'dnf' in outUnits:
                for k in toConvert:
                    # Differential Intensity [cm-2 s-1 sr-1 keV-1]
                    self.data_raw[k] /= (self.geomFactorEff * q * self.data_raw['DIM1_E'] * 1e-3)
                self.data_raw['DATA_UNITS'][:] = 7

            elif 'psd' in outUnits:
                for k in toConvert:
                    # Phase space density [s3 cm-6]
                    speedSquared = 2 * self.data_raw['DIM1_E'] * q * e / m  # m2 s-2 
                    self.data_raw[k] /= (self.geomFactorEff * 0.5e8 * speedSquared**2)
                

    def calc_pitchangles(self, useInternal=False):

        elevation = self.data_raw['DIM2_ELEVATION'] * np.pi / 180
        azimuth = self.data_raw['DIM2_AZIMUTH_DESPUN'] * np.pi / 180
        
        if useInternal:
            magvector = self.data_raw['MAG_VECTOR']
            Bmag = np.sqrt(np.sum(magvector**2, axis=0))
            Bx = magvector[0, :] / Bmag
            By = magvector[1, :] / Bmag
            Bz = magvector[2, :] / Bmag

        else:
            magdata = mag.get_mag_data(
                self.times[0], self.times[-1],
                coordsystem='Juno_JADE_despun'
                ).df.loc[:, ['Bx', 'By', 'Bz', 'Bmag']]
        
            magdata = magdata.reindex(self.times, method='nearest')
            
            for c in ['Bx', 'By', 'Bz']:
                magdata[c] /= magdata['Bmag']

            Bx = magdata['Bx'].values
            By = magdata['By'].values
            Bz = magdata['Bz'].values

        self.data_raw['PITCH_ANGLE'] = np.arccos(-(
            np.cos(elevation) * np.cos(azimuth) * Bx + \
            np.cos(elevation) * np.sin(azimuth) * By + \
            np.sin(elevation) * Bz
            )) * 180 / np.pi
        

    def rebin_sensor_data(
            self, 
            minE=None, maxE=None, minPA=None, maxPA=None, nPAbins=12,
            addBackground=False):
        """
        Rebin the sensor data from (nLookDir, nE, nTime) -> (nPA, nE, nTime).
        This is useful for averaging later on.
        """
        d = self.data_raw

        energies = d['DIM1_E']
        data = d['DATA']
        pitchangles = d['PITCH_ANGLE']
        uncertainty = d['DATA_SIGMA']

        if addBackground:
            data += d['BACKGROUND']
            
        # Set invalid values to zero or NaN
        data[data < -999990] = np.NaN
        data[(data < 0) & (data > -999990)] = 0.
        
        uncertainty[uncertainty < -999990] = np.NaN
        uncertainty[(uncertainty < 0) & (uncertainty > -999990)] = 0.


        # Filter out the appropriate energies
        _energies = energies[0, :, 0]
        condition = np.full(_energies.shape, True)
        if minE:
            condition *= (_energies >= minE)
        
        if maxE:
            condition *= (_energies <= maxE)

        energies = energies[:, condition, :]
        data = data[:, condition, :]
        uncertainty = uncertainty[:, condition, :]
        pitchangles = pitchangles[:, condition, :]
        
        energyBins = _create_energy_bins(_energies)
        nEbins = len(energyBins) - 1

        paBins = np.linspace(0, 180, nPAbins + 1)

        # Create masks - 
        # 1. Mask invalid pitch angle values (fortunately this is only on one axis)
        # 2. Mask any NaN values prior to averaging
        
        dataRebinned = np.zeros((nPAbins, nEbins, data.shape[2]))
        uncertaintyRebinned = np.zeros((nPAbins, nEbins, data.shape[2]))

        nanMask = ~np.isnan(data)
        sigmaMask = ~np.isnan(uncertainty)

        for i in range(nPAbins):
            paMask = (paBins[i] <= pitchangles) & (paBins[i + 1] > pitchangles)
            mask = paMask & nanMask
            with np.errstate(divide='ignore', invalid='ignore'):

                dataRebinned[i, :, :] = np.sum(data * mask, axis=0) / np.sum(mask, axis=0)

                uncertaintyRebinned[i, :, :] = \
                    np.sum(uncertainty**2 * mask * sigmaMask, axis=0) / np.sum(mask * sigmaMask, axis=0)

        dataRebinned[dataRebinned == 0.] = np.NaN

        uncertaintyRebinned = np.sqrt(uncertaintyRebinned)
        uncertaintyRebinned[uncertaintyRebinned == 0.] = np.NaN

        self.dataR = dataRebinned
        self.uncertaintyR = uncertaintyRebinned

        self.energyBins = energyBins
        self.paBins = paBins


    def calc_energy_spectra(self, weighSolidAngle=False):
        """
        Given some data that has already been rebinned into (nPA, nE, nTime),
        create the energy spectra (nE, nTime) by averaging over the different 
        pitch angles.
        """
        paBins = self.paBins

        if weighSolidAngle:
            deltaPitchAngle = (paBins[1] - paBins[0]) * np.pi / 180
            midPitchAngle = 0.5 * (paBins[1:] + paBins[:-1]) * np.pi / 180
            solidAngles = np.sin(midPitchAngle) \
                            * np.sin(deltaPitchAngle * 0.5)\
                            * 4 * np.pi
    
            weights = solidAngles
        else:
            weights = np.zeros(len(paBins) - 1) + 1.

        # We do this for numpy broadcasting magic
        weights = weights[:, np.newaxis, np.newaxis]

        nanmask = ~np.isnan(self.dataR)
        self.dataR[~nanmask] = 0.

        # We use this context as there is no guarantee that nanmask has 
        # atleast one True value, which would make the denominator 0. 
        # 0 in the denominator would lead to a NaN, which is perfectly 
        # fine in our case. Physically, it means that no measurements
        # were made for that particular bin.

        with np.errstate(divide='ignore', invalid='ignore'):
            energySpectra = np.sum(self.dataR * weights * nanmask, axis=0)\
                            / np.sum(weights * nanmask, axis=0)

        energySpectra[energySpectra == 0.] = np.NaN

        self.energySpectra = pd.DataFrame(
            index=self.times,
            columns=[f'Ebin_{x}' for x in range(len(self.energyBins) - 1)],
            data=energySpectra.T)


    def calc_pitch_angle_spectra(self, weighEnergy=True):
        """
        Given some data that has already been rebinned into (nPA, nE, nTime),
        create the pitch angle spectra. Default weighing along the energy axis is 
        the energy bandwidth. 
        """
        energyBins = self.energyBins

        if weighEnergy:
            weights = np.sqrt(energyBins[1:] * energyBins[:-1])
        else:
            weights = np.zeros(len(energyBins) - 1) + 1.
        
        # We do this for numpy broadcasting magic
        weights = weights[np.newaxis, :, np.newaxis]

        nanmask = ~np.isnan(self.dataR)
        self.dataR[~nanmask] = 0.

        with np.errstate(divide='ignore', invalid='ignore'):
            pitchAngleSpectra = np.sum(self.dataR * weights * nanmask, axis=1)\
                                / np.sum(weights * nanmask, axis=1)
        
        pitchAngleSpectra[pitchAngleSpectra == 0.] = np.NaN
        
        self.pitchAngleSpectra = pd.DataFrame(
            index=self.times,
            columns=['PAbin_{:d}'.format(x) for x in range(len(self.paBins) - 1)],
            data=pitchAngleSpectra.T)
        

    def rebin_to_perpparallel(self, m=m_e):
        """
        Using the pitchangle to rebin the data from 
            (PA, E, time) to (vPerp, vParallel, time)
        """
        # Reprocessed data has dimensions (nPA, nE, time)

        pitchAnglesMid = 0.5 * (self.paBins[1:] + self.paBins[:-1])
        energiesMid = np.sqrt(self.energyBins[1:] * self.energyBins[:-1])
        speeds = np.sqrt(2 * energiesMid * e / m)
        
        paMesh, vMesh = np.meshgrid(pitchAnglesMid, speeds, indexing='ij')
        paMesh *= np.pi / 180.

        vPerp = vMesh * np.sin(paMesh)
        vParallel = vMesh * np.cos(paMesh)
        
        pass





    def calc_accumulation_time_seconds(self):
        fmt = "%Y-%jT%H:%M:%S.%f"
        zipped = zip(self.data_raw['DIM0_UTC_UPPER'], 
                     self.data_raw['DIM0_UTC_LOWER'])
        
        self.accumulationTime = np.array([
            (dt.datetime.strptime(x, fmt) - dt.datetime.strptime(y, fmt))\
                .total_seconds() for x, y in zipped])
    
    def valid_times_copy(self, times):
        """
        Data present in sensors that will need to be truncated.
        """
        x = deepcopy(self)

        indices = (self.times != self.times)
        for t in times:
            indices += (t == self.times)

        x.times = x.times[indices]

        for key in ['DATA', 'DATA_SIGMA', 'BACKGROUND', 'BACKGROUND_SIGMA']:
            x.data_raw[key] = x.data_raw[key][:, :, indices]

        x.accumulationTime = x.accumulationTime[indices]

        x.geomFactorEff = x.geomFactorEff[:, :, indices]

        try:
            x.energySpectra = x.energySpectra.iloc[indices].copy()
            x.pitchAngleSpectra = x.pitchAngleSpectra.iloc[indices, :].copy()
        except AttributeError:
            pass
        
        try:
            x.dataR = x.dataR[:, :, indices]
            x.uncertaintyR = x.uncertaintyR[:, :, indices]
        except AttributeError:
            pass

        return x
        



