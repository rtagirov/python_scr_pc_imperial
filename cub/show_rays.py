import netCDF4

import matplotlib.pyplot as plt
import numpy as np

import random
import itertools

import sys
import os

z = netCDF4.Dataset('Z_onTau.123000.nc')['Z']
T = netCDF4.Dataset('T_onTau.123000.nc')['T']
p = netCDF4.Dataset('P_onTau.123000.nc')['P']
d = netCDF4.Dataset('rho_onTau.123000.nc')['R']

z = np.array(z) / 1e+5
T = np.array(T)
p = np.array(p)
d = np.array(d)

plt.close('all')

fig, ax = plt.subplots(nrows = 1, ncols = 1, figsize = (5, 5))

r1 = random.sample(range(512), 10)
r2 = random.sample(range(512), 10)

for i, j in itertools.product(r1, r2):

    ax.plot(z[:, i, j], T[:, i, j])

ax.set_xlabel('Height, km')
ax.set_ylabel('Temperature, K')

fig.savefig('temp.pdf', bbox_inches = 'tight')

r1 = random.sample(range(512), 100)
r2 = random.sample(range(512), 10)

if not os.path.isdir('./atms'): os.mkdir('./atms')

n = 1

for i, j in itertools.product(r1, r2):

    np.savetxt('./atms/atm.' + str(n), \
               np.column_stack([z[:, i, j], T[:, i, j], p[:, i, j], d[:, i, j]]), \
               fmt = ('%6.1f', '%7.1f', '%7.5e', '%7.5e'), delimiter = '  ')

    n += 1
