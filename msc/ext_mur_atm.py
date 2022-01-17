import numpy as np
import netCDF4
import random

import matplotlib.pyplot as plt

import sys

from tqdm import tqdm

path = '/mnt/HDD/extr_atms/'

pres_nc = netCDF4.Dataset(path + 'result_pres.123000.nc')['pres']
temp_nc = netCDF4.Dataset(path + 'result_temp.123000.nc')['temp']
dens_nc = netCDF4.Dataset(path + 'result_rhox.123000.nc')['rhox']

rn = np.arange(0, 512)

#plt.close('all')

for k in tqdm(range(8, 100)):

    i = random.choice(rn)
    j = random.choice(rn)

    pres = np.array(list(pres_nc[:, i, j]))
    temp = np.array(list(temp_nc[:, i, j]))
    dens = np.array(list(dens_nc[:, i, j]))

    p = np.array([x for m, x in enumerate(pres) if m % 4 == 0])
    T = np.array([x for m, x in enumerate(temp) if m % 4 == 0])
    d = np.array([x for m, x in enumerate(dens) if m % 4 == 0])

    ne =       np.ones(len(p))
    vt = 1.5 * np.ones(len(p))

    f1 = np.zeros(len(p))
    f2 = np.zeros(len(p))

    np.savetxt(path + 'muram_hd_' + str(k), \
    np.transpose((d, T, p, ne, f1, f2, vt)), \
    fmt = ('%10.4E', '%10.4E', '%10.4E', '%10.4E', '%10.4E', '%10.4E', '%10.4E'), delimiter = '  ')

#    plt.plot(d, T)

#plt.show()
