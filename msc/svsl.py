import numpy as np
import matplotlib.pyplot as plt

import sys
import paths
import pltaux

mod = sys.argv[1]

fmt = sys.argv[2]

hs =  np.loadtxt(paths.it0h + 'redmur/' + fmt + '/' + mod + '/ATM_STR', skiprows = 2, usecols = [1])
hl =  np.loadtxt(paths.it0h + 'murkur/' + fmt + '/' + mod + '/ATM_STR', skiprows = 2, usecols = [1])

Ts =  np.loadtxt(paths.it0h + 'redmur/' + fmt + '/' + mod + '/ATM_STR', skiprows = 2, usecols = [2])
Tl =  np.loadtxt(paths.it0h + 'murkur/' + fmt + '/' + mod + '/ATM_STR', skiprows = 2, usecols = [2])

ns =  np.loadtxt(paths.it0h + 'redmur/' + fmt + '/' + mod + '/ATM_STR', skiprows = 2, usecols = [3])
nl =  np.loadtxt(paths.it0h + 'murkur/' + fmt + '/' + mod + '/ATM_STR', skiprows = 2, usecols = [3])

dhs = np.loadtxt(paths.it0h + 'redmur/' + fmt + '/' + mod + '/ATM_STR', skiprows = 2, usecols = [4])
dhl = np.loadtxt(paths.it0h + 'murkur/' + fmt + '/' + mod + '/ATM_STR', skiprows = 2, usecols = [4])

rs =  np.loadtxt(paths.it0h + 'redmur/' + fmt + '/' + mod + '/ATM_STR', skiprows = 2, usecols = [5])
rl =  np.loadtxt(paths.it0h + 'murkur/' + fmt + '/' + mod + '/ATM_STR', skiprows = 2, usecols = [5])

Tgs = np.loadtxt(paths.it0h + 'redmur/' + fmt + '/' + mod + '/ATM_STR', skiprows = 2, usecols = [6])
Tgl = np.loadtxt(paths.it0h + 'murkur/' + fmt + '/' + mod + '/ATM_STR', skiprows = 2, usecols = [6])

ngs = np.loadtxt(paths.it0h + 'redmur/' + fmt + '/' + mod + '/ATM_STR', skiprows = 2, usecols = [7])
ngl = np.loadtxt(paths.it0h + 'murkur/' + fmt + '/' + mod + '/ATM_STR', skiprows = 2, usecols = [7])

#plt.close('all')

fig, ax = plt.subplots(nrows = 2, ncols = 2, figsize = (15, 10))

ax[0, 0].plot(hs, Ts)
ax[0, 0].plot(hl, Tl)

#ax[0, 0].set_xlim(0, 1000)
#ax[0, 0].set_ylim(3000, 12000)

ax[0, 1].plot(rs, ns)
ax[0, 1].plot(rl, nl)

#ax[0, 1].set_ylim(1e+13, 1e+18)

ax[1, 0].plot(hs, Tgs)
ax[1, 0].plot(hl, Tgl)

#ax[1, 0].set_ylim(-210, 50)
#ax[1, 0].set_xlim(0, 1000)

ax[1, 1].plot(rs, ngs)
ax[1, 1].plot(rl, ngl)

#ax[1, 1].set_xlim(1.000, 1.0015)
#ax[0, 1].set_xlim(1.000, 1.0015)

ax[0, 1].set_yscale('log')

pltaux.savepdf(paths.figdir, 'grad')

#plt.show()

#fig, ax = plt.subplots(nrows = 1, ncols = 1, figsize = (10, 5))

#ax.plot(hs, dhs)
#ax.plot(hl, dhl)

#ax.set_xlim(0, 1000)
#ax.set_ylim(0, 60)

#plt.show()
