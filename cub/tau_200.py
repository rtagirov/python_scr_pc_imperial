import numpy as np
import matplotlib.pyplot as plt

from netCDF4 import Dataset

import sys
import auxplt
import random
import os

from matplotlib.ticker import AutoMinorLocator
from matplotlib.ticker import MultipleLocator

reload(auxplt)

from tqdm import tqdm

tau_vero = np.loadtxt('5top_points.dat')

nx = 512
ny = 512
nz = 324

#nrr = 1000
nrr = 50

rn =  np.arange(0, 512)
apn = np.arange(0, 25)
#apn = np.arange(0, 1)

maxtau0 = np.zeros(len(apn))
maxtau1 = np.zeros(len(apn))
maxtau2 = np.zeros(len(apn))
maxtau3 = np.zeros(len(apn))
maxtau4 = np.zeros(len(apn))

base_pdf_name = 'ext_'

pdfs = ''

for m in tqdm(range(len(apn))):

    name = base_pdf_name + str(m)

    dp = np.arange(nz + m)

#    nc_file = '/mnt/SSD/sim/tau_200/bin/result_t_tau.123000.ext.' + str(apn[m]) + '.nc'
    nc_file = '/mnt/SSD/sim/tau_200/bin/result_t_tau.118000.ext.' + str(apn[m]) + '.nc'

    data = Dataset(nc_file, mode = 'r')

    tau = data.variables['tau']

    maxtau0[m] = np.max(tau[0, :, :])
    maxtau1[m] = np.max(tau[1, :, :])
    maxtau2[m] = np.max(tau[2, :, :])
    maxtau3[m] = np.max(tau[3, :, :])
    maxtau4[m] = np.max(tau[4, :, :])

    fig, ax = plt.subplots(nrows = 1, ncols = 1, figsize = (10, 5))

    for k in range(nrr):

        i = random.choice(rn)
        j = random.choice(rn)

        od = tau[:, i, j]

        ax.plot(dp, od, color = 'gray')

#        if (m == 0):

#            ax.plot(np.arange(5), tau_vero[i * j, :], color = 'red')

    ax.set_yscale('log')

    ax.set_xlim(0.000, 49.00)
    ax.set_ylim(1e-10, 1.0e+3)

#    if (m == 0):

#        ax.set_ylim(1e-10, 1e+3)
#    ax.set_ylim(1e-10, 1e-1)

    ax.axvline(m, linestyle = '--', color = 'k')

    ax.set_xlabel('Depth Index')
    ax.set_ylabel('Optical depth')

    ax.xaxis.set_minor_locator(AutoMinorLocator(10))

    auxplt.savepdf(name, './')

    plt.close('all')

    pdfs += ' ' + name + '.pdf'

os.system('pdftk ' + pdfs + ' output overall.pdf')

fig, ax = plt.subplots(nrows = 1, ncols = 1, figsize = (10, 5))

#ax.plot(apn, maxtau0, color = 'black',  label = '0th depth point')
ax.plot(apn, maxtau1, color = 'red',    label = '2st layer')
ax.plot(apn, maxtau2, color = 'blue',   label = '3nd layer')
ax.plot(apn, maxtau3, color = 'orange', label = '4rd layer')
ax.plot(apn, maxtau4, color = 'green',  label = '5th layer')

ax.axhline(y = 0.2, linestyle = '--', color = 'k')

#ax.scatter(np.arange(1), np.max(tau_vero[:, 0]), color = 'black')
#ax.scatter(np.arange(1), np.max(tau_vero[:, 1]), color = 'red')
#ax.scatter(np.arange(1), np.max(tau_vero[:, 2]), color = 'blue')
#ax.scatter(np.arange(1), np.max(tau_vero[:, 3]), color = 'orange')
#ax.scatter(np.arange(1), np.max(tau_vero[:, 4]), color = 'green')

ax.set_yscale('log')

ax.set_xlabel('Number of layers added to the cube')
ax.set_ylabel('Maximum of optical depth over the layer')

ax.set_ylim(1e-1, 30)
#ax.set_ylim(8e-8, 1e-4)
ax.set_xlim(0, 24)

ax.xaxis.set_minor_locator(AutoMinorLocator(5))

leg = ax.legend(framealpha = 1, loc = 1, handletextpad = 1, prop = {'size': 7.5})

for obj in leg.legendHandles: obj.set_linewidth(3.0)

auxplt.savepdf('tau_200_max_300G', './')
