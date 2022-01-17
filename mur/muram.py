import numpy as np
import matplotlib.pyplot as plt

import paths

plt.close('all')

h_mur = np.loadtxt(paths.it0h + 'murtst/mur/ATM_STR', skiprows = 2, usecols = [1])
h_kur = np.loadtxt(paths.it0h + 'murtst/kur/ATM_STR', skiprows = 2, usecols = [1])
h_utm = np.loadtxt(paths.it0h + 'murtst/utm/ATM_STR', skiprows = 2, usecols = [1])
h_fal = np.loadtxt(paths.it0h + 'murtst/fal/ATM_STR', skiprows = 2, usecols = [1])

r_mur = np.loadtxt(paths.it0h + 'murtst/mur/ATM_STR', skiprows = 2, usecols = [5])
r_kur = np.loadtxt(paths.it0h + 'murtst/kur/ATM_STR', skiprows = 2, usecols = [5])
r_utm = np.loadtxt(paths.it0h + 'murtst/utm/ATM_STR', skiprows = 2, usecols = [5])
r_fal = np.loadtxt(paths.it0h + 'murtst/fal/ATM_STR', skiprows = 2, usecols = [5])

T_mur = np.loadtxt(paths.it0h + 'murtst/mur/ATM_STR', skiprows = 2, usecols = [2])
T_kur = np.loadtxt(paths.it0h + 'murtst/kur/ATM_STR', skiprows = 2, usecols = [2])
T_utm = np.loadtxt(paths.it0h + 'murtst/utm/ATM_STR', skiprows = 2, usecols = [2])
T_fal = np.loadtxt(paths.it0h + 'murtst/fal/ATM_STR', skiprows = 2, usecols = [2])

n_mur = np.loadtxt(paths.it0h + 'murtst/mur/ATM_STR', skiprows = 2, usecols = [3])
n_kur = np.loadtxt(paths.it0h + 'murtst/kur/ATM_STR', skiprows = 2, usecols = [3])
n_utm = np.loadtxt(paths.it0h + 'murtst/utm/ATM_STR', skiprows = 2, usecols = [3])
n_fal = np.loadtxt(paths.it0h + 'murtst/fal/ATM_STR', skiprows = 2, usecols = [3])

Tg_mur = np.loadtxt(paths.it0h + 'murtst/mur/ATM_STR', skiprows = 2, usecols = [6])
Tg_kur = np.loadtxt(paths.it0h + 'murtst/kur/ATM_STR', skiprows = 2, usecols = [6])
Tg_utm = np.loadtxt(paths.it0h + 'murtst/utm/ATM_STR', skiprows = 2, usecols = [6])
Tg_fal = np.loadtxt(paths.it0h + 'murtst/fal/ATM_STR', skiprows = 2, usecols = [6])

ng_mur = np.loadtxt(paths.it0h + 'murtst/mur/ATM_STR', skiprows = 2, usecols = [7])
ng_kur = np.loadtxt(paths.it0h + 'murtst/kur/ATM_STR', skiprows = 2, usecols = [7])
ng_utm = np.loadtxt(paths.it0h + 'murtst/utm/ATM_STR', skiprows = 2, usecols = [7])
ng_fal = np.loadtxt(paths.it0h + 'murtst/fal/ATM_STR', skiprows = 2, usecols = [7])

dh_mur = np.loadtxt(paths.it0h + 'murtst/mur/ATM_STR', skiprows = 2, usecols = [4])
dh_kur = np.loadtxt(paths.it0h + 'murtst/kur/ATM_STR', skiprows = 2, usecols = [4])
dh_utm = np.loadtxt(paths.it0h + 'murtst/utm/ATM_STR', skiprows = 2, usecols = [4])
dh_fal = np.loadtxt(paths.it0h + 'murtst/fal/ATM_STR', skiprows = 2, usecols = [4])

fig, ax = plt.subplots(nrows = 2, ncols = 2, figsize = (15, 10))

ax[0, 0].plot(h_fal, T_fal)
ax[0, 0].plot(h_mur, T_mur)
ax[0, 0].plot(h_kur, T_kur)
ax[0, 0].plot(h_utm, T_utm)

ax[0, 0].set_xlim(0, 1000)
ax[0, 0].set_ylim(3000, 12000)

ax[0, 1].plot(r_fal, n_fal)
ax[0, 1].plot(r_mur, n_mur)
ax[0, 1].plot(r_kur, n_kur)
ax[0, 1].plot(r_utm, n_utm)

ax[0, 1].set_ylim(1e+13, 1e+18)

ax[1, 0].plot(h_fal, Tg_fal)
ax[1, 0].plot(h_mur, Tg_mur)
ax[1, 0].plot(h_kur, Tg_kur)
ax[1, 0].plot(h_utm, Tg_utm)

ax[1, 0].set_ylim(-210, 50)
ax[1, 0].set_xlim(0, 1000)

ax[1, 1].plot(r_fal, ng_fal)
ax[1, 1].plot(r_mur, ng_mur)
ax[1, 1].plot(r_kur, ng_kur)
ax[1, 1].plot(r_utm, ng_utm)

ax[1, 1].set_xlim(1.000, 1.0015)
ax[0, 1].set_xlim(1.000, 1.0015)

ax[0, 1].set_yscale('log')

plt.show()

fig, ax = plt.subplots(nrows = 1, ncols = 1, figsize = (10, 5))

ax.plot(h_fal, dh_fal)
ax.plot(h_mur, dh_mur)
ax.plot(h_kur, dh_kur)
ax.plot(h_utm, dh_utm)

ax.set_xlim(0, 1000)
ax.set_ylim(0, 60)

plt.show()
