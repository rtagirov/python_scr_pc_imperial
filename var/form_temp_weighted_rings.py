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

            if i == 0: plt.plot(wvl[idx], contr[idx], color = col, label = label)
            if i != 0: plt.plot(wvl[idx], contr[idx], color = col)

        else:

            plt.plot(wvl[idx], -contr[idx], color = col, linestyle = '--')

def formation_temperature(path, mu, wvl1, wvl2, mode, num):

    h = np.loadtxt(path + '/ATM_STR', usecols = [1], skiprows = 2) / mu
    T = np.loadtxt(path + '/ATM_STR', usecols = [2], skiprows = 2)

    h = np.flip(h, axis = 0)
    T = np.flip(T, axis = 0)

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

    fTw = []

    for i in tqdm(range(int(na)), ncols = auxfunc.term_width(), desc = 'mode ' + str(mode) + '; run ' + str(num) + '; ' + path):

        idx = str(imin + i * nint)

        f1 = path + '/tau/' + idx + '.tau'
        f2 = path + '/mdisp/' + idx + '.mdisp'

        fh = np.loadtxt(f1, usecols = [1])
        I =  np.loadtxt(f2, usecols = [1])

        fT = np.zeros(len(fh))

        for k in range(len(fh)):

            if np.isnan(fh[k]): fh[k] = 0.0

            j = np.searchsorted(h, fh[k], side = 'right')

            fT[k] = T[j - 1] + (fh[k] - h[j - 1]) * (T[j] - T[j - 1]) / (h[j] - h[j - 1])

        fTw.append(sum(I * fT) / sum(I))

    return np.array(fTw)

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

cases = ['Q kur 0', 'F kur 0', 'Q kur 1', 'F kur 1', 'Q fal 1', 'F fal 1']

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

fT = np.zeros((len(cases), len(wvl)))

nws = int(m.ceil((wvl[len(wvl) - 1] - wvl[0]) / 20))

fTs = np.zeros((len(cases), nws))

num = 1

for i, case in enumerate(cases):

    ar = case.split()[0]
    mo = case.split()[1]
    it = case.split()[2]

    path = ''

    if it == '0': path = paths.it0f
    if it == '1': path = paths.it1f

    path += 'var/' + ar + '/rings/' + mo

#    for j in range(len(mu)):

#        if mode >= 0 and mode <= 10 and j != mode: continue

#        fTr = formation_temperature(path + '/' + str(mu[j])[:4], mu[j], wvl[0], wvl[len(wvl) - 1] - 5, mode, num)

#        fT[i, :] += fTr * w[j]

#        num += 1
    
#    np.savez(paths.npz + ar + '_' + mo + '_' + it + '.ftemp.npz', fT = fT[i, :])

    fT[i, :] = np.load(paths.npz + ar + '_' + mo + '_' + it + '.ftemp.npz')['fT']

    wvls, tmps = spec.mean_within_delta(wvl / 10, fT[i, :], 2)

    fTs[i, :] = tmps

kur0 = (fTs[1, :] - fTs[0, :]) / fTs[0, :]
kur1 = (fTs[3, :] - fTs[2, :]) / fTs[2, :]
fal1 = (fTs[5, :] - fTs[4, :]) / fTs[4, :]

plt.close('all')

fig, ax = plt.subplots(nrows = 1, ncols = 1, figsize = (18, 6))

auxplt.figpar(3, 3, 15)

fig.tight_layout()

#if mode >= 0 and mode <= 10: 

#    fig.suptitle(r'$\mu = $ ' + str(mu[mode]) + '$, p = $ ' + str(np.sqrt(1 - mu[mode]**2)), y = 1.01)

#else:

#    fig.suptitle('mode = ' + str(mode), y = 1.01)

idx0_pos = np.where(kur0 >= 0.0)
idx1_pos = np.where(kur1 >= 0.0)
idx2_pos = np.where(fal1 >= 0.0)

idx0_neg = np.where(kur0 < 0.0)
idx1_neg = np.where(kur1 < 0.0)
idx2_neg = np.where(fal1 < 0.0)

plot_gr(wvls, kur0, idx0_pos, 'm', 'NESSY, LTE, U99')
plot_gr(wvls, kur1, idx1_pos, 'g', 'NESSY, NLTE, U99')
plot_gr(wvls, fal1, idx2_pos, 'r', 'NESSY, NLTE, FAL99')
                        
plot_gr(wvls, kur0, idx0_neg, 'm')
plot_gr(wvls, kur1, idx1_neg, 'g')
plot_gr(wvls, fal1, idx2_neg, 'r')

#plt.plot(wvl / 10, kur0, color = 'm', label = 'NESSY, LTE, U99')
#plt.plot(wvl / 10, kur1, color = 'g', label = 'NESSY, NLTE, U99')
#plt.plot(wvl / 10, fal1, color = 'r', label = 'NESSY, NLTE, FAL99')

plt.xlim(100, 1100)

plt.yscale('log')

#plt.ylim(top = 1.1e+3)
plt.ylim(top = 3e-1)

ax.xaxis.set_major_locator(MultipleLocator(100))
ax.xaxis.set_minor_locator(AutoMinorLocator(4))

plt.xlabel('Wavelength, [nm]')
#plt.ylabel(r'$T^\mathrm{form}_f - T^\mathrm{form}_q, [K]$')
plt.ylabel(r'$(T^\mathrm{form}_f - T^\mathrm{form}_q) / T^\mathrm{form}_q, [K]$')

leg = plt.legend(framealpha = 1, loc = 4, handletextpad = 1, prop = {'size': 17.0})

for obj in leg.legendHandles: obj.set_linewidth(3.0)

auxplt.savepdf('var/ftemp_weighted_rings_mode' + str(mode))
