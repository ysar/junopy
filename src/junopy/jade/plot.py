# """
# junopy - A Python API for the Juno spacecraft's particle and fields datasets.
# Copyright (C) 2018  Yash Sarkango

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
# """
# import numpy as np
# import datetime as dt

# import matplotlib as mpl
# from matplotlib import cm
# from matplotlib import colors


# def plot_timeseries_spectra(
#         dfSpectra, bins, ax,
#         vmin=None, vmax=None, cmap=None, zorder=1,
#         plotcolorbar=True, cbarlabel=None):

#     # x data is time bins
#     times = dfSpectra.index.to_pydatetime()
#     deltaT = times[1] - times[0]
#     timeBins = (times - deltaT)
#     timeBins = np.append(
#         timeBins, timeBins[-1] + (timeBins[1] - timeBins[0]))

#     # y data is just passed bins

#     # z data are spectral values
#     z = dfSpectra.values.T

#     # Correct for non-contigous data.
#     for i in range(len(times)):
#         if (timeBins[i + 1] - timeBins[i]).total_seconds() > 301:
#             z[:, i] = np.nan

#     if not vmin and not vmax:
#         def get_percentile(arr, p):
#             return np.exp(np.percentile(np.log(arr[arr > 0 & ~np.isnan(arr)]), p))

#         vmin = get_percentile(z, 1)
#         vmax = get_percentile(z, 99)

#     norm = colors.LogNorm(vmin=vmin, vmax=vmax)

#     if not cmap:
#         cmap = "turbo"

#     cmap = mpl.colormaps.get_cmap(cmap)
#     cmap.set_bad(color='0.5')


#     mp = ax.pcolormesh(
#         timeBins, bins, z, norm=norm, cmap=cmap, zorder=zorder)

#     ax.set_ylim([bins[0], bins[-1]])

#     if plotcolorbar:
#         bbox = ax.get_position()
#         cax = ax.figure.add_axes(
#             [bbox.x1 + 0.01, bbox.y0, 0.02, 0.8 * (bbox.y1- bbox.y0)])
#         cbar = ax.figure.colorbar(mp, cax=cax)
#         cbar.set_label(cbarlabel)

#     return ax, mp


# def plot_energy_pa_spectra(
#         ax, sensorData, start=None, stop=None,
#         doNormalizeEnergyBins=False,
#         cmap='turbo', plotcolorbar=True, vmin=None, vmax=None):

#     if start and stop:
#         indices = (sensorData.times >= start) & (sensorData.times <= stop)
#     else:
#         indices = np.full((sensorData.dataR.shape[-1],), True)

#     # z data are spectral values
#     z = np.nanmean(sensorData.dataR[:, :, indices], axis=2).T

#     if doNormalizeEnergyBins:

#         z = np.apply_along_axis(
#             lambda x: (x - x.min()) / (x.max() - x.min()), 1, z)
#         norm = colors.Normalize
#     else:
#         norm = colors.LogNorm

#     if not vmin and not vmax:
#         def get_percentile(arr, p):
#             return np.exp(np.percentile(np.log(arr[arr > 0 & ~np.isnan(arr)]), p))

#         vmin = get_percentile(z, 1)
#         vmax = get_percentile(z, 99)

#     mp = ax.pcolormesh(
#         sensorData.paBins, sensorData.energyBins, z,
#         norm=norm(vmin=vmin, vmax=vmax), cmap=cmap)

#     ax.set_ylim([sensorData.energyBins[0], sensorData.energyBins[-1]])
#     ax.set_xlim([0, 180])
#     ax.set_yscale('log')
#     ax.set_xticks([0, 45, 90, 135, 180])

#     if plotcolorbar:
#         bbox = ax.get_position()
#         cax = ax.figure.add_axes(
#             [bbox.x1 + 0.01, bbox.y0, 0.02, 0.8 * (bbox.y1- bbox.y0)])
#         cbar = ax.figure.colorbar(mp, cax=cax)

#     return ax, mp


# def plot_anode_energy_time_spectra(
#         axs, sensorData, data='data',
#         vmin=None, vmax=None, cmap=None, zorder=1,
#         plotcolorbar=True, locator=None):

#     # x data is time bins
#     times = sensorData.times
#     deltaT = times[1] - times[0]
#     timeBins = (times - deltaT)
#     timeBins = np.append(
#         timeBins, timeBins[-1] + (timeBins[1] - timeBins[0]))

#     # y data are the energy bins
#     y = sensorData.energyBins

#     # z data are spectral values
#     if data == 'data':
#         d = sensorData.data_raw['DATA']
#     elif data == 'data+bg':
#         sensorData.data_raw['DATA'][np.isnan(sensorData.data_raw['DATA'])] = 0.
#         d = sensorData.data_raw['DATA'] + sensorData.data_raw['BACKGROUND']
#     elif data == 'pitchangle':
#         d = sensorData.data_raw['PITCH_ANGLE']
#     elif data == 's/n':
#         sensorData.data_raw['DATA'][np.isnan(sensorData.data_raw['DATA'])] = 0.
#         d = sensorData.data_raw['DATA'] / sensorData.data_raw['BACKGROUND']

#     if not data == 'pitchangle':
#         if not vmin and not vmax:
#             def get_percentile(arr, p):
#                 return np.exp(np.percentile(np.log(arr[arr > 0 & ~np.isnan(arr)]), p))

#             vmin = get_percentile(d, 1)
#             vmax = get_percentile(d, 99)

#         norm = colors.LogNorm(vmin=vmin, vmax=vmax)
#     else:
#         vmin = 0
#         vmax = 180
#         norm = colors.Normalize(vmin=vmin, vmax=vmax)

#     for i in range(4):
#         for j in range(4):

#             iAnode = 4 * i + j
#             axs[i, j].set_facecolor('0.7')
#             axs[i, j].pcolormesh(timeBins, y, d[iAnode, :, :], norm=norm, cmap=cmap)
#             axs[i, j].set_yscale('log')
#             axs[i, j].set_xlabel("Time")
#             axs[i, j].set_ylabel("Energy [eV]")
#             axs[i, j].text(
#                 0.1, 0.8, f"Anode={iAnode}", backgroundcolor='w',
#                 transform=axs[i, j].transAxes, fontsize=7)
#             if locator:
#                 axs[i, j].xaxis.set_major_locator(locator)
#             axs[i, j].label_outer()


#     mp = cm.ScalarMappable(norm=norm, cmap=cmap)

#     ax = axs[3, 3]
#     bbox = ax.get_position()
#     cax = ax.figure.add_axes(
#         [bbox.x1 + 0.01, bbox.y0, 0.02, 2.5 * (bbox.y1- bbox.y0)])
#     cbar = ax.figure.colorbar(mp, cax=cax)
#     # cax.text(0.9, 1.1, "Counts s$^{-1}$", transform=cax.transAxes)

#     return axs, mp
