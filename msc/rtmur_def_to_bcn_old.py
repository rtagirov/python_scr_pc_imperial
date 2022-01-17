import numpy as np
import matplotlib.pyplot as plt

from matplotlib.ticker import AutoMinorLocator
from matplotlib.ticker import MultipleLocator

import importlib
import sys

if not '../aux/' in sys.path: sys.path.append('../aux/')

import paths;  importlib.reload(paths)
import spec;   importlib.reload(spec)
import nessy;  importlib.reload(nessy)
import auxsys; importlib.reload(auxsys)
import auxplt; importlib.reload(auxplt)

N = int(sys.argv[1])

fd = np.zeros((N, 901))
fb = np.zeros((N, 901))

frmh = np.zeros((N, 901))

fhi = np.zeros((N, 1802000))

height = np.zeros((N, 81))

for i in range(N):

    wvl, flu_def = nessy.read_spec(paths.it1f + 'rtmur/def/' + str(i), wvl1 = 1005, wvl2 = 10000)
    wvl, flu_bcn = nessy.read_spec(paths.it1f + 'rtmur/bcn/' + str(i), wvl1 = 1005, wvl2 = 10000)

    wvl = wvl / 10.0

    wvls, fd[i, :] = spec.mean_within_delta(wvl, flu_def, 1.0)
    wvls, fb[i, :] = spec.mean_within_delta(wvl, flu_bcn, 1.0)

    outname = 'rtmur_def_to_cbn_' + str(i)

    np.savez(paths.npz + outname, wv = wvls, fd = fd[i, :], fb = fb[i, :])

    specs = np.load(paths.npz + outname + '.npz')

    wv =       specs['wv']
    fd[i, :] = specs['fd']
    fb[i, :] = specs['fb']

#    wvlh, frmh[i, :] = nessy.weighted_form_height(paths.it1f + 'rtmur/def/' + str(i), wvl1 = 1005, wvl2 = 10000)

#    wvlh /= 10.0

    outname = 'rtmur_def_frmh_' + str(i)

#    np.savez(paths.npz + outname, wv = wvlh, fh = frmh[i, :])

    frmhs = np.load(paths.npz + outname + '.npz')

    wvlh =       frmhs['wv']
    frmh[i, :] = frmhs['fh']

    height[i, :] = np.loadtxt(paths.it1h + 'rtmur/def/' + str(i) + '/strat.out', skiprows = 1, usecols = [1]) 

#    wvlh_hr, fhi[i, :], x = nessy.read_tau(paths.it1f + 'rtmur/def/' + str(i) + '/', wvl1 = 1005, wvl2 = 10000)

#    wvlh_hr /= 10.0

    outname = 'rtmur_def_frmh_hr_' + str(i)

#    np.savez(paths.npz + outname, wv = wvlh_hr, fh = fhi[i, :])

    frmh_hs = np.load(paths.npz + outname + '.npz')

    wvlh_hr =   frmh_hs['wv']
    fhi[i, :] = frmh_hs['fh']

plt.close('all')

fig, ax = plt.subplots(nrows = 1, ncols = 1, figsize = (6.0, 6.0))

fig.tight_layout()

for i in range(N):

    ax.plot(wv, fd[i, :] / fb[i, :], color = 'gray')
#    ax[1].plot(wvlh, frmh[i, :],          color = 'gray')

    ax.set_xlim(100, 1000)
#    ax[1].set_xlim(100, 1000)

#    ax[j].xaxis.set_minor_locator(AutoMinorLocator(10))
#    ax[j].xaxis.set_major_locator(MultipleLocator(100))

    ax.set_ylabel('Ratio')
#    ax[1].set_ylabel('Formation height, km')

    ax.axhline(y = 1.0, color = 'r', linestyle = '-', linewidth = 0.5)

ax.set_xlabel('Wavelength, [nm]')

#ax[0].set_yscale('log')

#leg0 = ax.legend(framealpha = 1, loc = 1, handletextpad = 1, prop = {'size':7.5})

#for obj in leg0.legendHandles: obj.set_linewidth(3.0)

auxplt.savepdf('rtmur_def_to_bcn_1p')

fig, ax = plt.subplots(nrows = N, ncols = 1, figsize = (6.0, 12.0))

fig.tight_layout()

for i in range(N):

#    ax[i].plot(wvlh_hr, height[i, fhi[i, :].astype(int)], color = 'gray')
    ax[i].plot(wvlh_hr, fhi[i, :], color = 'gray')

    ax[i].set_xlim(100, 1000)
    ax[i].set_ylim(10, 0)

    ax[i].yaxis.set_minor_locator(AutoMinorLocator(5))
#    ax[j].xaxis.set_major_locator(MultipleLocator(100))

    ax[i].set_ylabel('Formation index')

    ax[i].axhline(y = max(height[i, :]), color = 'r', linestyle = '-', linewidth = 0.5)

ax[N - 1].set_xlabel('Wavelength, [nm]')

#ax[0].set_yscale('log')

#leg0 = ax.legend(framealpha = 1, loc = 1, handletextpad = 1, prop = {'size':7.5})

#for obj in leg0.legendHandles: obj.set_linewidth(3.0)

auxplt.savepdf('rtmur_frmh')
