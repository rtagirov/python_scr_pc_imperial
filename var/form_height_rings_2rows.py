import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as const
import more_itertools as mit
import math as m

from matplotlib.ticker import AutoMinorLocator
from matplotlib.ticker import MultipleLocator

from tqdm import tqdm

import sys

if not '../aux/' in sys.path: sys.path.append('../aux/')

import importlib
import auxplt
import auxsys
import auxfunc
import paths
import nessy
import spec

importlib.reload(auxplt)
importlib.reload(auxsys)
importlib.reload(auxfunc)
importlib.reload(paths)
importlib.reload(nessy)
importlib.reload(spec)

def plot_gr(wvl, contr, indices, col, label = ''):

    groups = [list(gr) for gr in mit.consecutive_groups(indices[0].tolist())]
    
    for i, g in enumerate(groups):

        idx = (np.array(g), )

        if len(label) != 0:

            if i == 0: ax[0].plot(wvl[idx], contr[idx], color = col, label = label)
            if i != 0: ax[0].plot(wvl[idx], contr[idx], color = col)

        else:

            ax[0].plot(wvl[idx], -contr[idx], color = col, linestyle = '--')

def formation_height(path, wvl1, wvl2, mode, num):

    nf =  2000 # number of lines in each .tau and .mdisp file

    ivl = 10   # length of the spectral interval of each .tau and .mdisp file (angstroems)
    mid = 5    # distance to the middle of each .tau and .mdisp file (angstroems)

    nint = int(ivl)

    wmin = int(wvl1)
    wmax = int(wvl2)

    imin = wmin - ((wmin - mid) % nint)
    imax = wmax - ((wmax - mid) % nint) + nint

    na = (imax - imin) / nint + 1 # number of arrays

    nw = na * nf

    fhw = []

    for i in tqdm(range(int(na)), ncols = auxfunc.term_width(), desc = 'mode ' + str(mode) + '; run ' + str(num) + '; ' + path):

        idx = str(imin + i * nint)

        f1 = path + '/tau/' + idx + '.tau'
        f2 = path + '/mdisp/' + idx + '.mdisp'

        fh = np.loadtxt(f1, usecols = [1])
        I =  np.loadtxt(f2, usecols = [1])

        fh[np.where(np.isnan(fh))] = 0.0

        fhw.append(sum(I * fh) / sum(I))

    return np.array(fhw)

def ring_weights(mu, mode):

    p = np.sqrt(1 - mu**2)

    w = np.zeros(len(p))

    if mode >= 0 and mode <= 10:

        w[mode] = 1.0

    if mode == 11:

        for i in range(1, len(p) - 1):

            w[i] = p[i + 1]**2 - p[i - 1]**2

        w[0] = p[1]**2

        w[len(p) - 1] = 1 - p[len(p) - 2]**2

    if mode == 12:

        for i in range(len(p) - 1):

            w[i] = p[i + 1]**2 - p[i]**2

        w[len(p) - 1] = 1 - p[len(p) - 1]**2

    if mode == 13:

#        w = np.ones(len(p)) * 11**(-1)
        w = np.ones(len(p)) * 0.1

    return w / sum(w)

#cases = ['Q kur 0', 'F kur 0', 'Q kur 1', 'F kur 1', 'Q fal 1', 'F fal 1']
cases = ['Q kur 0', 'F kur 0', 'Q fal 1', 'F fal 1']

#mu = np.array([1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.05])
#mu = np.array([0.87])

p = np.sqrt(np.arange(0, 11)) / np.sqrt(10)

p_mid = np.zeros(len(p) - 1)

for i in range(len(p_mid)):

    p_mid[i] = (p[i + 1] + p[i]) / 2

mu = np.sqrt(1 - p_mid**2)

mode = int(sys.argv[1])

w = ring_weights(mu, mode)

wvl = np.linspace(1005, 11005, 1001)

fh = np.zeros((len(cases), len(wvl)))

delta = 5

nws = int(m.ceil((wvl[len(wvl) - 1] - wvl[0]) / delta / 10))

fhs = np.zeros((len(cases), nws))

num = 1

for i, case in enumerate(cases):

    ar = case.split()[0]
    mo = case.split()[1]
    it = case.split()[2]

    path = ''

    if it == '0': path = paths.it0f
    if it == '1': path = paths.it1f

    path += 'var_od/' + ar + '/rings/' + mo

#    for j in range(len(mu)):

##        if mode >= 0 and mode <= 10 and j != mode: continue

#        fhr = formation_height(path + '/' + str(mu[j])[:4], wvl[0], wvl[len(wvl) - 1] - 5, mode, num)

#        fh[i, :] += fhr * mu[j] * w[j]

#        num += 1
    
#    np.savez(paths.npz + ar + '_' + mo + '_' + it + '.fhei_od.npz', fh = fh[i, :])

    fheight = np.load(paths.npz + ar + '_' + mo + '_' + it + '.fhei_od.npz')['fh']

    wvls_h, heis = spec.mean_within_delta(wvl / 10, fheight, delta)

    fhs[i, :] = heis

plt.close('all')

fig, ax = plt.subplots(nrows = 1, ncols = 1, figsize = (6.0, 3.375))
#fig, ax = plt.subplots(nrows = 2, ncols = 1, figsize = (6.0, 6.75))

auxplt.figpar(3, 3, 15)

#fig.tight_layout()

#ax.set_xlim(100, 1100)
ax.set_xlim(100, 500)
#ax[1].set_xlim(100, 1100)
ax.set_ylim(30, 2e+3)
#ax[1].set_ylim(-1e+3, 1e+2)
#ax[1].set_ylim(0.5, 1.35)

ax.set_yscale('log')
#ax[1].set_yscale('symlog')

#ax.xaxis.set_major_locator(MultipleLocator(200))
ax.xaxis.set_major_locator(MultipleLocator(50))
#ax.xaxis.set_minor_locator(AutoMinorLocator(4))
ax.xaxis.set_minor_locator(AutoMinorLocator(5))
#ax[1].xaxis.set_major_locator(MultipleLocator(200))
#ax[1].xaxis.set_minor_locator(AutoMinorLocator(4))
#ax[1].yaxis.set_major_locator(MultipleLocator(0.1))
#ax[1].yaxis.set_minor_locator(AutoMinorLocator(4))

ax.fill_between(wvls_h, 140 * np.ones(len(wvls_h)), 400 * np.ones(len(wvls_h)), color = 'gray', alpha = 0.5)

ax.plot(wvls_h, fhs[0, :] - 134.416, color = 'm')
ax.plot(wvls_h, fhs[1, :] - 104.806, color = 'm', linestyle = '--')
ax.plot(wvls_h, fhs[2, :] - 114.388, 'r', alpha = 0.5)
ax.plot(wvls_h, fhs[3, :] - 106.802, 'r', alpha = 0.5, linestyle = '--')

ax.plot(wvls_h, 1e+6 * fhs[2, :], color = 'k',   label = 'Quiet Sun')
ax.plot(wvls_h, 1e+6 * fhs[3, :], color = 'k', label = 'Facula', linestyle = '--')

#ax[1].plot(wvls_h, (fhs[1, :] - 104.806) / (fhs[0, :] - 134.416), color = 'm')
#ax[1].plot(wvls_h, (fhs[3, :] - 106.802) / (fhs[2, :] - 114.388), 'r', alpha = 0.5)

ax.set_ylabel('Formation height, [km]')
#ax[1].set_ylabel('Facula to quiet Sun ratio')
ax.set_xlabel('Wavelength, [nm]')

ax.axvline(x = 210, linestyle = ':', color = 'k')
ax.axvline(x = 450, linestyle = ':', color = 'k')
#ax[1].axvline(x = 210, linestyle = ':', color = 'k')
#ax[1].axvline(x = 400, linestyle = ':', color = 'k')
#ax[1].axhline(y = 1,   linestyle = '--', color = 'k')

leg1 = ax.legend(framealpha = 1, loc = 1, handletextpad = 1, prop = {'size': 12.0})

for obj in leg1.legendHandles: obj.set_linewidth(3.0)

#ax.text(465, 1.3e+3, 'NESSY, LTE, U99',    color = 'm', fontsize = 11)
#ax.text(465, 8.0e+2, 'NESSY, NLTE, FAL99', color = 'r', fontsize = 11)
ax.text(220, 1.3e+3, 'NESSY, LTE, U99',    color = 'm', fontsize = 11)
ax.text(220, 8.0e+2, 'NESSY, NLTE, FAL99', color = 'r', fontsize = 11)

auxplt.savepdf('var/fheight_weighted_rings_2rows_od')
