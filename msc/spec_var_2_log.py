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

#1 - NESSY LTE
#2 - NESSY NLTE
#3 - ATLAS

#4 - NESSY LTE  FAL
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

F3[np.where(wa < 122.5)] = 0.0
U3[np.where(wa < 168.5)] = 0.0

waQ = wa[np.where(Q3 > 0.0)]
Q30 = Q3[np.where(Q3 > 0.0)]

waF = wa[np.where(F3 > 0.0)]
F30 = F3[np.where(F3 > 0.0)]

waU = wa[np.where(U3 > 0.0)]
U30 = U3[np.where(U3 > 0.0)]

waP = wa[np.where(P3 > 0.0)]
P30 = P3[np.where(P3 > 0.0)]

wn, Q1h = nessy.read_spec(prefix0 + 'var/Q/kur/', wvl1 = 1005, wvl2 = 13000)
wn, F1h = nessy.read_spec(prefix0 + 'var/F/kur/', wvl1 = 1005, wvl2 = 13000)
wn, U1h = nessy.read_spec(prefix0 + 'var/U/kur/', wvl1 = 1005, wvl2 = 13000)
wn, P1h = nessy.read_spec(prefix0 + 'var/P/kur/', wvl1 = 1005, wvl2 = 13000)

wn, Q2h = nessy.read_spec(prefix1 + 'var/Q/kur/', wvl1 = 1005, wvl2 = 13000)
wn, F2h = nessy.read_spec(prefix1 + 'var/F/kur/', wvl1 = 1005, wvl2 = 13000)
wn, U2h = nessy.read_spec(prefix1 + 'var/U/kur/', wvl1 = 1005, wvl2 = 13000)
wn, P2h = nessy.read_spec(prefix1 + 'var/P/kur/', wvl1 = 1005, wvl2 = 13000)

wn, Q5h = nessy.read_spec(prefix1 + 'var/Q/fal/', wvl1 = 1005, wvl2 = 13000)
wn, F5h = nessy.read_spec(prefix1 + 'var/F/fal/', wvl1 = 1005, wvl2 = 13000)

wn = wn / 10.0

Q1 = spec.mean_over_grid(Q1h, wn, wa)
F1 = spec.mean_over_grid(F1h, wn, wa)
U1 = spec.mean_over_grid(U1h, wn, wa)
P1 = spec.mean_over_grid(P1h, wn, wa)

Q2 = spec.mean_over_grid(Q2h, wn, wa)
F2 = spec.mean_over_grid(F2h, wn, wa)
U2 = spec.mean_over_grid(U2h, wn, wa)
P2 = spec.mean_over_grid(P2h, wn, wa)

Q5 = spec.mean_over_grid(Q5h, wn, wa)
F5 = spec.mean_over_grid(F5h, wn, wa)

np.savez(paths.npz + 'contr_var_2', w = wa,
                                    q1 = Q1,\
                                    f1 = F1,\
                                    u1 = U1,\
                                    p1 = P1,\
                                    q2 = Q2,\
                                    f2 = F2,\
                                    u2 = U2,\
                                    p2 = P2,\
                                    q5 = Q5,\
                                    f5 = F5,)

contr = np.load(paths.npz + 'contr_var_2.npz')

w =  contr['w']

Q1 = contr['q1']
F1 = contr['f1']
U1 = contr['u1']
P1 = contr['p1']

Q2 = contr['q2']
F2 = contr['f2']
U2 = contr['u2']
P2 = contr['p2']

Q5 = contr['q5']
F5 = contr['f5']

plt.close('all')

fig, ax = plt.subplots(nrows = 2, ncols = 1, figsize = (18.0, 12.00))

bbox = dict(boxstyle = 'round', ec = (1.0, 0.5, 0.5), fc = (1.0, 0.8, 0.8),)

auxplt.figpar(3, 3, 20)

fig.tight_layout()

plt.subplots_adjust(hspace = 0.00)

ls = ':'
lw = 1.0

ax[0].plot(w, Q3, color = 'k', linewidth = lw * 2, label = 'ATLAS9 (LTE, U99)')
ax[0].plot(w, Q1, color = 'm', linewidth = lw, label = 'NESSY (LTE, U99)')
ax[0].plot(w, Q2, color = 'g', linewidth = lw, label = 'NESSY (NLTE, U99)')
ax[0].plot(w, Q5, color = 'r', linewidth = lw, label = 'NESSY (NLTE, FAL99)')

ax[0].xaxis.set_tick_params(direction = 'in', which = 'both')

in1 = fig.add_axes([0.055, 0.722, 0.16, 0.231])

in1.plot(waQ, Q30, color = 'k', linewidth = lw * 2)
in1.plot(w,   Q1,  color = 'm', linewidth = lw)
in1.plot(w,   Q2,  color = 'g', linewidth = lw)
in1.plot(w,   Q5,  color = 'r', linewidth = lw)

in1.set_xlim(100, 300)
in1.set_ylim(1e-10, 1e+0)

in1.set_yscale('log')

in1.tick_params(labelsize = 15)

in1.yaxis.tick_right()

ax[0].text(720, 1.5, 'Quiet Sun', bbox = bbox)

ax[1].plot(w, F3, color = 'k', linewidth = lw * 2, label = 'ATLAS9 (LTE, U99)')
ax[1].plot(w, F1, color = 'm', linewidth = lw, label = 'NESSY (LTE, U99)')
ax[1].plot(w, F2, color = 'g', linewidth = lw, label = 'NESSY (NLTE, U99)')
ax[1].plot(w, F5, color = 'r', linewidth = lw, label = 'NESSY (NLTE, FAL99)')

ax[1].plot(w, P3, color = 'k', linewidth = lw * 2)
ax[1].plot(w, P1, color = 'm', linewidth = lw)
ax[1].plot(w, P2, color = 'g', linewidth = lw)

ax[1].plot(w, U3, color = 'k', linewidth = lw * 2)
ax[1].plot(w, U1, color = 'm', linewidth = lw)
ax[1].plot(w, U2, color = 'g', linewidth = lw)

in2 = fig.add_axes([0.055, 0.268, 0.16, 0.231])

in2.plot(waF, F30, color = 'k', linewidth = lw * 2)
in2.plot(w,   F1,  color = 'm', linewidth = lw)
in2.plot(w,   F2,  color = 'g', linewidth = lw)
in2.plot(w,   F5,  color = 'r', linewidth = lw)

in2.plot(waP, P30, color = 'k', linewidth = lw * 2)
in2.plot(w,   P1,  color = 'm', linewidth = lw)
in2.plot(w,   P2,  color = 'g', linewidth = lw)

in2.plot(waU, U30, color = 'k', linewidth = lw * 2)
in2.plot(w,   U1,  color = 'm', linewidth = lw)
in2.plot(w,   U2,  color = 'g', linewidth = lw)

in2.set_xlim(100, 300)
in2.set_ylim(1e-17, 1e+0)

in2.set_yscale('log')

in2.tick_params(labelsize = 15)

in2.yaxis.tick_right()

ax[1].text(700, 1.65, 'Facula', bbox = bbox)
ax[1].text(600, 0.95, 'Penumbra', bbox = bbox)
ax[1].text(530, 0.22, 'Umbra',    bbox = bbox)

for i in range(0, 2):

#    ax[i].set_xlim(125.0, 800)
    ax[i].set_xlim(180, 800)

    ax[i].xaxis.set_major_locator(MultipleLocator(50))
    ax[i].xaxis.set_minor_locator(AutoMinorLocator(5))

    ax[i].yaxis.set_minor_locator(AutoMinorLocator(5))

    ax[i].set_ylabel(r'Flux, [W / m$^2$ / nm]', fontsize = 20)

ax[0].xaxis.set_tick_params(direction = 'in', which = 'both')

ax[1].set_xlabel('Wavelength, [nm]', fontsize = 20)
ax[0].tick_params(labelbottom = 'off')

#leg0 = ax[0].legend(framealpha = 1, loc = 3, bbox_to_anchor = (0.379, 0.03), handletextpad = 1, prop = {'size': 25.0})
leg0 = ax[0].legend(framealpha = 1, loc = 4, handletextpad = 1, prop = {'size': 25.0})

for obj in leg0.legendHandles: obj.set_linewidth(3.0)

auxplt.savepdf('var/spec_var_2_log')
