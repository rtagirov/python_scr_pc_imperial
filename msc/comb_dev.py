import numpy as np
import matplotlib.pyplot as plt

from matplotlib.ticker import AutoMinorLocator
from matplotlib.ticker import MultipleLocator

import importlib
import sys
import math

if not '../aux/' in sys.path: sys.path.append('../aux/')

import paths;  importlib.reload(paths)
import spec;   importlib.reload(spec)
import nessy;  importlib.reload(nessy)
import auxsys; importlib.reload(auxsys)
import auxplt; importlib.reload(auxplt)

lambda1 = 1005
#lambda2 = 10000
lambda2 = 3000

Nw = int(lambda2 / 10.0) - math.floor(lambda1 / 10.0) + 1

#wvl, faldef = nessy.read_spec(paths.it1f + 'runtime/def',                             wvl1 = lambda1, wvl2 = lambda2)
#wvl, falbcn = nessy.read_spec(paths.it1f + 'runtime/nlte11_fgrid_red_5_20_eps_2',     wvl1 = lambda1, wvl2 = lambda2)
#wvl, kurdef = nessy.read_spec(paths.it1f + 'runtime/def_kur',                         wvl1 = lambda1, wvl2 = lambda2)
#wvl, kurbcn = nessy.read_spec(paths.it1f + 'runtime/nlte11_fgrid_red_5_20_eps_2_kur', wvl1 = lambda1, wvl2 = lambda2)

#wvl = wvl / 10.0

#wvls, faldef_s = spec.mean_within_delta(wvl, faldef, 1.0)
#wvls, falbcn_s = spec.mean_within_delta(wvl, falbcn, 1.0)
#wvls, kurdef_s = spec.mean_within_delta(wvl, kurdef, 1.0)
#wvls, kurbcn_s = spec.mean_within_delta(wvl, kurbcn, 1.0)

outname = 'rtmur_def_to_bcn_standard'

#np.savez(paths.npz + outname, w = wvls, f1 = faldef_s,
#                                        f2 = falbcn_s,
#                                        f3 = kurdef_s,
#                                        f4 = kurbcn_s)

specs = np.load(paths.npz + outname + '.npz')

faldef = specs['f1']
falbcn = specs['f2']
kurdef = specs['f3']
kurbcn = specs['f4']

h_fal = np.loadtxt(paths.it0h + 'runtime/def/strat.out',     skiprows = 1, usecols = [1])
T_fal = np.loadtxt(paths.it0h + 'runtime/def/strat.out',     skiprows = 1, usecols = [2])
n_fal = np.loadtxt(paths.it0h + 'runtime/def/strat.out',     skiprows = 1, usecols = [3])
h_kur = np.loadtxt(paths.it0h + 'runtime/def_kur/strat.out', skiprows = 1, usecols = [1])
T_kur = np.loadtxt(paths.it0h + 'runtime/def_kur/strat.out', skiprows = 1, usecols = [2])
n_kur = np.loadtxt(paths.it0h + 'runtime/def_kur/strat.out', skiprows = 1, usecols = [3])

N = 41

fd = np.zeros((N, Nw))
fb = np.zeros((N, Nw))

dev = np.zeros((N, Nw))

h = np.zeros((N, 81))
T = np.zeros((N, 81))
n = np.zeros((N, 81))

for i in range(N):

#    wvl, flu_def = nessy.read_spec(paths.it1f + 'rtmur/def/' + str(i), wvl1 = lambda1, wvl2 = lambda2)
#    wvl, flu_bcn = nessy.read_spec(paths.it1f + 'rtmur/bcn/' + str(i), wvl1 = lambda1, wvl2 = lambda2)

#    wvl = wvl / 10.0

#    wvls, fd[i, :] = spec.mean_within_delta(wvl, flu_def, 1.0)
#    wvls, fb[i, :] = spec.mean_within_delta(wvl, flu_bcn, 1.0)

    outname = 'rtmur_def_to_cbn_' + str(i)

#    np.savez(paths.npz + outname, wv = wvls, fd = fd[i, :], fb = fb[i, :])

    specs = np.load(paths.npz + outname + '.npz')

    wv =       specs['wv']
    fd[i, :] = specs['fd']
    fb[i, :] = specs['fb']

    h[i, :] = np.loadtxt('/mnt/HDD/sim/runs/hminus/IT1/rtmur/def/' + str(i) + '/strat.out', skiprows = 1, usecols = [1])
    T[i, :] = np.loadtxt('/mnt/HDD/sim/runs/hminus/IT1/rtmur/def/' + str(i) + '/strat.out', skiprows = 1, usecols = [2])
    n[i, :] = np.loadtxt('/mnt/HDD/sim/runs/hminus/IT1/rtmur/def/' + str(i) + '/strat.out', skiprows = 1, usecols = [3])

plt.close('all')

fig, ax = plt.subplots(nrows = 3, ncols = 1, figsize = (6.0, 6.0))

fig.tight_layout()

for i in range(N):

    dev[i, :] = fd[i, :] / fb[i, :]

    ax[0].plot(wv, dev[i, :], color = 'gray')

    if i == 0: ax[1].plot(h[i, :], T[i, :], color = 'gray', label = 'MuRaM')
    if i != 0: ax[1].plot(h[i, :], T[i, :], color = 'gray')

    ax[2].plot(h[i, :], n[i, :], color = 'gray')

mean_dev = np.mean(dev[:, :], axis = 0)
mean_h =   np.mean(h[:, :],   axis = 0)
mean_T =   np.mean(T[:, :],   axis = 0)
mean_n =   np.mean(n[:, :],   axis = 0)

ax[0].plot(wv,    faldef / falbcn, color = 'b')
ax[0].plot(wv,    kurdef / kurbcn, color = 'r')
ax[0].plot(wv,    mean_dev,        color = 'k')

ax[1].plot(mean_h, mean_T, color = 'k', label = 'MuRaM average')
ax[1].plot(h_fal,  T_fal,  color = 'b', label = 'FAL99-C')
ax[1].plot(h_kur,  T_kur,  color = 'r', label = 'Kurucz')

ax[2].plot(mean_h, mean_n, color = 'k')
ax[2].plot(h_fal,  n_fal,  color = 'b')
ax[2].plot(h_kur,  n_kur,  color = 'r')

ax[0].axvline(x = 170,  color = 'k', linestyle = '--', linewidth = 1.0)
ax[0].axhline(y = 1.00, color = 'k', linestyle = '--', linewidth = 1.0)

ax[0].axhline(y = 1.07, color = 'k', linestyle = '--', linewidth = 1.0)
ax[0].axhline(y = 0.93, color = 'k', linestyle = '--', linewidth = 1.0)

ax[0].set_xlim(100, 300)
ax[1].set_xlim(0.0, 1050)
ax[2].set_xlim(0.0, 1050)

ax[1].set_ylim(3000, 12000)
ax[2].set_ylim(5e+13, 3e+17)

ax[0].set_ylabel('Ratio')
ax[0].set_xlabel('Wavelength, nm')
ax[1].set_ylabel('Temperature, K')
ax[2].set_ylabel(r'Number density, cm$^{-3}$')
ax[2].set_xlabel('Height, km')

ax[2].set_yscale('log')

ax[0].xaxis.set_minor_locator(AutoMinorLocator(4))
ax[0].xaxis.set_major_locator(MultipleLocator(20))
ax[1].xaxis.set_minor_locator(AutoMinorLocator(4))
ax[1].xaxis.set_major_locator(MultipleLocator(100))
ax[2].xaxis.set_minor_locator(AutoMinorLocator(4))
ax[2].xaxis.set_major_locator(MultipleLocator(100))

leg0 = ax[1].legend(framealpha = 1, loc = 1, handletextpad = 1, prop = {'size':7.5})

for obj in leg0.legendHandles: obj.set_linewidth(3.0)

auxplt.savepdf('rtmur_def_to_bcn')
