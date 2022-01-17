import numpy as np

import matplotlib.pyplot as plt

import paths
import pltaux

plt.close('all')

h1_kur = np.loadtxt(paths.it0h + 'murkur/kur/1/ATM_STR', skiprows = 2, usecols = [1])
T1_kur = np.loadtxt(paths.it0h + 'murkur/kur/1/ATM_STR', skiprows = 2, usecols = [2])
n1_kur = np.loadtxt(paths.it0h + 'murkur/kur/1/ATM_STR', skiprows = 2, usecols = [3])

h2_kur = np.loadtxt(paths.it0h + 'murkur/kur/2/ATM_STR', skiprows = 2, usecols = [1])
T2_kur = np.loadtxt(paths.it0h + 'murkur/kur/2/ATM_STR', skiprows = 2, usecols = [2])
n2_kur = np.loadtxt(paths.it0h + 'murkur/kur/2/ATM_STR', skiprows = 2, usecols = [3])

h3_kur = np.loadtxt(paths.it0h + 'murkur/kur/3/ATM_STR', skiprows = 2, usecols = [1])
T3_kur = np.loadtxt(paths.it0h + 'murkur/kur/3/ATM_STR', skiprows = 2, usecols = [2])
n3_kur = np.loadtxt(paths.it0h + 'murkur/kur/3/ATM_STR', skiprows = 2, usecols = [3])

h4_kur = np.loadtxt(paths.it0h + 'murkur/kur/4/ATM_STR', skiprows = 2, usecols = [1])
T4_kur = np.loadtxt(paths.it0h + 'murkur/kur/4/ATM_STR', skiprows = 2, usecols = [2])
n4_kur = np.loadtxt(paths.it0h + 'murkur/kur/4/ATM_STR', skiprows = 2, usecols = [3])

h5_kur = np.loadtxt(paths.it0h + 'murkur/kur/5/ATM_STR', skiprows = 2, usecols = [1])
T5_kur = np.loadtxt(paths.it0h + 'murkur/kur/5/ATM_STR', skiprows = 2, usecols = [2])
n5_kur = np.loadtxt(paths.it0h + 'murkur/kur/5/ATM_STR', skiprows = 2, usecols = [3])

h1_mur = np.loadtxt(paths.it0h + 'murkur/mur/1/ATM_STR', skiprows = 2, usecols = [1])
T1_mur = np.loadtxt(paths.it0h + 'murkur/mur/1/ATM_STR', skiprows = 2, usecols = [2])
n1_mur = np.loadtxt(paths.it0h + 'murkur/mur/1/ATM_STR', skiprows = 2, usecols = [3])

h2_mur = np.loadtxt(paths.it0h + 'murkur/mur/2/ATM_STR', skiprows = 2, usecols = [1])
T2_mur = np.loadtxt(paths.it0h + 'murkur/mur/2/ATM_STR', skiprows = 2, usecols = [2])
n2_mur = np.loadtxt(paths.it0h + 'murkur/mur/2/ATM_STR', skiprows = 2, usecols = [3])

h3_mur = np.loadtxt(paths.it0h + 'murkur/mur/3/ATM_STR', skiprows = 2, usecols = [1])
T3_mur = np.loadtxt(paths.it0h + 'murkur/mur/3/ATM_STR', skiprows = 2, usecols = [2])
n3_mur = np.loadtxt(paths.it0h + 'murkur/mur/3/ATM_STR', skiprows = 2, usecols = [3])

h4_mur = np.loadtxt(paths.it0h + 'murkur/mur/4/ATM_STR', skiprows = 2, usecols = [1])
T4_mur = np.loadtxt(paths.it0h + 'murkur/mur/4/ATM_STR', skiprows = 2, usecols = [2])
n4_mur = np.loadtxt(paths.it0h + 'murkur/mur/4/ATM_STR', skiprows = 2, usecols = [3])

h5_mur = np.loadtxt(paths.it0h + 'murkur/mur/5/ATM_STR', skiprows = 2, usecols = [1])
T5_mur = np.loadtxt(paths.it0h + 'murkur/mur/5/ATM_STR', skiprows = 2, usecols = [2])
n5_mur = np.loadtxt(paths.it0h + 'murkur/mur/5/ATM_STR', skiprows = 2, usecols = [3])

fig, ax = plt.subplots(nrows = 1, ncols = 2, figsize = (15, 5))

ax[0].plot(h1_kur, T1_kur)
ax[0].plot(h2_kur, T2_kur)
ax[0].plot(h3_kur, T3_kur)
ax[0].plot(h4_kur, T4_kur)
ax[0].plot(h5_kur, T5_kur)

ax[0].plot(h1_mur, T1_mur)
ax[0].plot(h2_mur, T2_mur)
ax[0].plot(h3_mur, T3_mur)
ax[0].plot(h4_mur, T4_mur)
ax[0].plot(h5_mur, T5_mur)

ax[0].set_xlabel('Height, [km]')
ax[0].set_ylabel('Temperature, [K]')

ax[1].plot(h1_kur, n1_kur)
ax[1].plot(h2_kur, n2_kur)
ax[1].plot(h3_kur, n3_kur)
ax[1].plot(h4_kur, n4_kur)
ax[1].plot(h5_kur, n5_kur)

ax[1].plot(h1_mur, n1_mur)
ax[1].plot(h2_mur, n2_mur)
ax[1].plot(h3_mur, n3_mur)
ax[1].plot(h4_mur, n4_mur)
ax[1].plot(h5_mur, n5_mur)

ax[1].set_yscale('log')

ax[1].set_xlabel('Height, [km]')
ax[1].set_ylabel('Number density, [cm$^{-3}$]')

pltaux.savepdf(paths.figdir, 'murvskur')
