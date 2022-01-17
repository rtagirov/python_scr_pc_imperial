import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as const

from matplotlib.ticker import AutoMinorLocator

import sys

if not '../aux/' in sys.path: sys.path.append('../aux/')

import importlib
import auxplt
import auxsys
import paths

importlib.reload(auxplt)
importlib.reload(auxsys)
importlib.reload(paths)

h_q_fal, T_q_fal, x, x, x = np.loadtxt(paths.it1f + 'var/Q/fal/atm.inp', unpack = True)
h_f_fal, T_f_fal, x, x, x = np.loadtxt(paths.it1f + 'var/F/fal/atm.inp', unpack = True)

wvl_q_fal, fh_q_fal = np.loadtxt(paths.it1f + 'var/Q/fal/form.height', unpack = True)
wvl_f_fal, fh_f_fal = np.loadtxt(paths.it1f + 'var/F/fal/form.height', unpack = True)

fT_q_fal = np.zeros(len(fh_q_fal))
fT_f_fal = np.zeros(len(fh_f_fal))

for i, w in enumerate(fh_q_fal):

    idx_q = np.abs(fh_q_fal[i] - h_q_fal).argmin()

    fT_q_fal[i] = T_q_fal[idx_q]

for i, w in enumerate(fh_f_fal):

    idx_f = np.abs(fh_f_fal[i] - h_f_fal).argmin()

    fT_f_fal[i] = T_f_fal[idx_f]

plt.close('all')

fig, ax = plt.subplots(nrows = 1, ncols = 1, figsize = (const.golden * 5, 5))

col = ['g', 'r', 'c', 'm', 'b', 'k', 'y']

auxplt.figpar(3, 3, 15)

fig.tight_layout()

plt.plot(wvl_q_fal / 10.0, fT_q_fal, label = 'Quiet Sun')

plt.xlim(325, 350)
plt.ylim(4000, 10000)

plt.xlabel(r'$Wavelength, [nm]$')
plt.ylabel(r'$T_f - T_q$')

leg = plt.legend(framealpha = 1, loc = 1, handletextpad = 1, prop = {'size': 17.0})

for obj in leg.legendHandles: obj.set_linewidth(3.0)

auxplt.savepdf('var/ftemp')
#plt.show()
