import numpy as np
import matplotlib.pyplot as plt

from matplotlib.ticker import AutoMinorLocator
from matplotlib.ticker import MultipleLocator
from matplotlib.ticker import LogLocator

import importlib
import math
import sys

if not '../aux/' in sys.path: sys.path.append('../aux/')

import paths;  importlib.reload(paths)
import spec;   importlib.reload(spec)
import phys;   importlib.reload(phys)
import nessy;  importlib.reload(nessy)
import auxsys; importlib.reload(auxsys)
import auxplt; importlib.reload(auxplt)

prefix0 = paths.it0f
prefix1 = paths.it1f

#3 - ATLAS
#2 - NESSY NLTE
#5 - NESSY NLTE FAL

waf = np.loadtxt(paths.atlruns + 'var_m/Q/spec.out', skiprows = 2, usecols = [1])
Q3f = np.loadtxt(paths.atlruns + 'var_m/Q/spec.out', skiprows = 2, usecols = [3])
F3f = np.loadtxt(paths.atlruns + 'var_m/F/spec.out', skiprows = 2, usecols = [3])
U3f = np.loadtxt(paths.atlruns + 'var_m/U/spec.out', skiprows = 2, usecols = [3])
P3f = np.loadtxt(paths.atlruns + 'var_m/P/spec.out', skiprows = 2, usecols = [3])

Q3f = Q3f * phys.c / (waf * 1.0e-7)**2.0 * 1.0e-7 * (phys.r_sun / phys.au)**2.0 * 1.0e-3 * math.pi
F3f = F3f * phys.c / (waf * 1.0e-7)**2.0 * 1.0e-7 * (phys.r_sun / phys.au)**2.0 * 1.0e-3 * math.pi
P3f = P3f * phys.c / (waf * 1.0e-7)**2.0 * 1.0e-7 * (phys.r_sun / phys.au)**2.0 * 1.0e-3 * math.pi
U3f = U3f * phys.c / (waf * 1.0e-7)**2.0 * 1.0e-7 * (phys.r_sun / phys.au)**2.0 * 1.0e-3 * math.pi

idx = np.where((waf >= 100.5) & (waf <= 1300.0))

wa = waf[idx]

Q3 = Q3f[idx]
F3 = F3f[idx]
U3 = U3f[idx]
P3 = P3f[idx]

#wn, Q2h = nessy.read_spec(prefix1 + 'var/Q/kur/', wvl1 = 1005, wvl2 = 13000)
#wn, F2h = nessy.read_spec(prefix1 + 'var/F/kur/', wvl1 = 1005, wvl2 = 13000)
#wn, U2h = nessy.read_spec(prefix1 + 'var/U/kur/', wvl1 = 1005, wvl2 = 13000)
#wn, P2h = nessy.read_spec(prefix1 + 'var/P/kur/', wvl1 = 1005, wvl2 = 13000)

#wn, Q5h = nessy.read_spec(prefix1 + 'var/Q/fal/', wvl1 = 1005, wvl2 = 13000)
#wn, F5h = nessy.read_spec(prefix1 + 'var/F/fal/', wvl1 = 1005, wvl2 = 13000)

#wn = wn / 10.0

#Q2 = spec.mean_over_grid(Q2h, wn, wa)
#F2 = spec.mean_over_grid(F2h, wn, wa)
#U2 = spec.mean_over_grid(U2h, wn, wa)
#P2 = spec.mean_over_grid(P2h, wn, wa)

#Q5 = spec.mean_over_grid(Q5h, wn, wa)
#F5 = spec.mean_over_grid(F5h, wn, wa)

#np.savez(paths.npz + 'contr_var_2', w = wa,
#                                    q2 = Q2,\
#                                    f2 = F2,\
#                                    u2 = U2,\
#                                    p2 = P2,\
#                                    q5 = Q5,\
#                                    f5 = F5,)

contr = np.load(paths.npz + 'contr_var_2.npz')

w =  contr['w']

Q2 = contr['q2']
F2 = contr['f2']
U2 = contr['u2']
P2 = contr['p2']

Q5 = contr['q5']
F5 = contr['f5']

plt.close('all')

#fig, ax = plt.subplots(nrows = 1, ncols = 1, figsize = (6.0, 6.75))
fig, ax = plt.subplots(nrows = 1, ncols = 3, figsize = (18.0, 6.0))

bbox = dict(boxstyle = 'round', ec = (1.0, 0.5, 0.5), fc = (1.0, 0.8, 0.8),)

auxplt.figpar(3, 3, 15)

fig.tight_layout()

plt.subplots_adjust(hspace = 0.15)

ls = ':'
lw = 1.0

for i in range(0, 3):

    ax[i].plot(w, Q3, color = 'g', linewidth = lw, label = 'ATLAS9-U99')
#    ax[i].plot(w, Q2, color = 'c', linewidth = lw, label = 'NESSY-U99')
    ax[i].plot(w, Q5, color = 'r', linewidth = lw, label = 'NESSY-FAL99')

    ax[i].plot(w, F3, color = 'g', linewidth = lw, linestyle = '--')
#    ax[i].plot(w, F2, color = 'c', linewidth = lw, linestyle = '--')
    ax[i].plot(w, F5, color = 'r', linewidth = lw, linestyle = '--')

ax[0].set_yscale('log')
ax[1].set_yscale('log')

for i in range(0,  3):

    ax[i].set_xlabel('Wavelength, [nm]', fontsize = 15)

ax[0].set_xlim(100, 210)
ax[1].set_xlim(210, 400)
ax[2].set_xlim(400, 600)

ax[0].set_ylim(1e-10, 1e-1)
ax[1].set_ylim(1e-2, 2e+0)
ax[2].set_ylim(1.1, 2.3)

ax[0].xaxis.set_major_locator(MultipleLocator(20))
ax[0].xaxis.set_minor_locator(AutoMinorLocator(4))

ax[0].set_ylabel(r'Flux, [W / m$^2$ / nm]', fontsize = 15)

leg = ax[1].legend(framealpha = 1, loc = 4, handletextpad = 1, prop = {'size': 12.0})

for obj in leg.legendHandles: obj.set_linewidth(3.0)

auxplt.savepdf('var/spec_log')
