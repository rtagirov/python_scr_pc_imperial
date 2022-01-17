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

idx = np.where((waf >= 100.5) & (waf <= 1100.0))

wa = waf[idx]

Q3 = Q3f[idx]
F3 = F3f[idx]
U3 = U3f[idx]
P3 = P3f[idx]

#wn, Q1h = nessy.read_spec(prefix0 + 'var_od/Q/kur/', wvl1 = 1005, wvl2 = 11000)
#wn, F1h = nessy.read_spec(prefix0 + 'var_od/F/kur/', wvl1 = 1005, wvl2 = 11000)
#wn, U1h = nessy.read_spec(prefix0 + 'var_od/U/kur/', wvl1 = 1005, wvl2 = 11000)
#wn, P1h = nessy.read_spec(prefix0 + 'var_od/P/kur/', wvl1 = 1005, wvl2 = 11000)
#
#wn, Q2h = nessy.read_spec(prefix1 + 'var_od/Q/kur/', wvl1 = 1005, wvl2 = 11000)
#wn, F2h = nessy.read_spec(prefix1 + 'var_od/F/kur/', wvl1 = 1005, wvl2 = 11000)
#wn, U2h = nessy.read_spec(prefix1 + 'var_od/U/kur/', wvl1 = 1005, wvl2 = 11000)
#wn, P2h = nessy.read_spec(prefix1 + 'var_od/P/kur/', wvl1 = 1005, wvl2 = 11000)
#
#wn, Q5h = nessy.read_spec(prefix1 + 'var_od/Q/fal/', wvl1 = 1005, wvl2 = 11000)
#wn, F5h = nessy.read_spec(prefix1 + 'var_od/F/fal/', wvl1 = 1005, wvl2 = 11000)
#
#wn = wn / 10.0
#
#Q1 = spec.mean_over_grid(Q1h, wn, wa)
#F1 = spec.mean_over_grid(F1h, wn, wa)
#U1 = spec.mean_over_grid(U1h, wn, wa)
#P1 = spec.mean_over_grid(P1h, wn, wa)
#
#Q2 = spec.mean_over_grid(Q2h, wn, wa)
#F2 = spec.mean_over_grid(F2h, wn, wa)
#U2 = spec.mean_over_grid(U2h, wn, wa)
#P2 = spec.mean_over_grid(P2h, wn, wa)
#
#Q5 = spec.mean_over_grid(Q5h, wn, wa)
#F5 = spec.mean_over_grid(F5h, wn, wa)
#
#np.savez(paths.npz + 'contr_var', w = wa,
#                                  q1 = Q1,\
#                                  f1 = F1,\
#                                  u1 = U1,\
#                                  p1 = P1,\
#                                  q2 = Q2,\
#                                  f2 = F2,\
#                                  u2 = U2,\
#                                  p2 = P2,\
#                                  q5 = Q5,\
#                                  f5 = F5,)

contr = np.load(paths.npz + 'contr_var.npz')

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

FQ1 = (F1 - Q1) / Q1
FQ2 = (F2 - Q2) / Q2
FQ3 = (F3 - Q3) / Q3
FQ5 = (F5 - Q5) / Q5

UQ1 = (U1 - Q1) / Q1
UQ2 = (U2 - Q2) / Q2
UQ3 = (U3 - Q3) / Q3

PQ1 = (P1 - Q1) / Q1
PQ2 = (P2 - Q2) / Q2
PQ3 = (P3 - Q3) / Q3

#RDQ13 = (Q1 - Q3) * 100.0 / Q3
#RDQ23 = (Q2 - Q3) * 100.0 / Q3

#RDF13 = (F1 - F3) * 100.0 / F3
#RDF23 = (F2 - F3) * 100.0 / F3

#RDP13 = (P1 - P3) * 100.0 / P3
#RDP23 = (P2 - P3) * 100.0 / P3

#RDU13 = (U1 - U3) * 100.0 / U3
#RDU23 = (U2 - U3) * 100.0 / U3

#RDFQ13 = (FQ1 - FQ3) * 100.0 / FQ3
#RDFQ23 = (FQ2 - FQ3) * 100.0 / FQ3
#RDFQ43 = (FQ4 - FQ3) * 100.0 / FQ3
#RDFQ53 = (FQ5 - FQ3) * 100.0 / FQ3

#RDUQ13 = (UQ1 - UQ3) * 100.0 / UQ3
#RDUQ23 = (UQ2 - UQ3) * 100.0 / UQ3

#RDPQ13 = (PQ1 - PQ3) * 100.0 / PQ3
#RDPQ23 = (PQ2 - PQ3) * 100.0 / PQ3

plt.close('all')

fig, ax = plt.subplots(nrows = 2, ncols = 1, figsize = (6.0, 6.75))

bbox = dict(boxstyle = 'round', ec = (1.0, 0.5, 0.5), fc = (1.0, 0.8, 0.8),)

auxplt.figpar(3, 3, 15)

fig.tight_layout()

plt.subplots_adjust(hspace = 0.15)

ls = ':'
lw = 1.0

ax[0].axvline(x = 210, linestyle = '--', color = 'k')
ax[0].axvline(x = 450, linestyle = '--', color = 'k')

ax[0].plot(w, FQ3, color = 'k', linewidth = lw * 1.5, label = 'ATLAS9 (LTE, U99)')
ax[0].plot(w, FQ1, color = 'm', linewidth = lw, label = 'NESSY (LTE, U99)')
ax[0].plot(w, FQ2, color = 'g', linewidth = lw, label = 'NESSY (NLTE, U99)')
ax[0].plot(w, FQ5, color = 'r', linewidth = lw, label = 'NESSY (NLTE, FAL99)')

#ax[0].text(140, 4.7e-3, 'Facula', bbox = bbox)
ax[0].text(140, 5e-2, 'Facula', bbox = bbox)

ax[1].plot(w, PQ3, color = 'k', linewidth = lw * 1.5, label = 'ATLAS9 (LTE, U99)')
ax[1].plot(w, PQ1, color = 'm', linewidth = lw, label = 'NESSY (LTE, U99)')
ax[1].plot(w, PQ2, color = 'g', linewidth = lw, label = 'NESSY (NLTE, U99)')

ax[1].plot(w, UQ3, color = 'k', linewidth = lw * 1.5)
ax[1].plot(w, UQ1, color = 'm', linewidth = lw)
ax[1].plot(w, UQ2, color = 'g', linewidth = lw)

ax[1].text(200, -0.3, 'Penumbra', bbox = bbox)
ax[1].text(600, -0.55, 'Umbra',    bbox = bbox)

ax[0].set_yscale('log')

#ax[0].set_ylim(1e-3, 1e+3)
#ax[0].set_ylim(1e-2, 1e+3)
ax[0].set_ylim(1e-2, 5e+2)
ax[1].set_ylim(-1.0, -0.1)

ax[1].yaxis.set_minor_locator(AutoMinorLocator(4))

ax[1].set_xlabel('Wavelength, [nm]', fontsize = 15)

for i in range(0, 2):

    ax[i].set_xlim(100.0, 1100)

    ax[i].xaxis.set_major_locator(MultipleLocator(200))
    ax[i].xaxis.set_minor_locator(AutoMinorLocator(4))

#ax[0].tick_params(labelbottom = 'off')

ax[0].set_ylabel(r'$(S_f - S_q) / S_q$', fontsize = 15)
ax[1].set_ylabel(r'$(S_{\{p,\ u\}} - S_q) / S_q$', fontsize = 15)

leg0 = ax[0].legend(framealpha = 1, loc = 1, handletextpad = 1, prop = {'size': 12.0})
leg1 = ax[1].legend(framealpha = 1, loc = 4, handletextpad = 1, prop = {'size': 12.0})

for obj in leg0.legendHandles: obj.set_linewidth(3.0)
for obj in leg1.legendHandles: obj.set_linewidth(3.0)

auxplt.savepdf('var/contr_log')

sys.exit()

fig, ax = plt.subplots(nrows = 2, ncols = 1, figsize = (12.0, 10.0))

#pltaux.figpar()

fig.tight_layout()

fig.suptitle('Facular contrast, Kurucz models: NESSY vs. ATLAS', y = 1.01)

ax[0].plot(w, np.zeros(len(FQ1)), 'k--')

ax[0].set_xlim(110.5, 1000)

ls = ':'; lw = 1.5

ax[0].plot(w, FQ1, color = 'b', linewidth = lw, label = 'NESSY (LTE)')
ax[0].plot(w, FQ2, color = 'r', linewidth = lw, label = 'NESSY (NLTE)')
ax[0].plot(w, FQ3, color = 'g', linewidth = lw, label = 'ATLAS')

ax[1].plot(w, abs(RDFQ13), color = 'b',                 label = 'NESSY (LTE)  vs ATLAS')
ax[1].plot(w, abs(RDFQ23), color = 'r', linewidth = lw, label = 'NESSY (NLTE) vs ATLAS')

ax[1].plot(w, np.ones(len(RDFQ13)), 'k')

ax[1].set_yscale('log')

ax[1].set_xlim(110.5, 1000)
ax[1].set_ylim(1e-2, 1e+3)

for i in range(0, 2):

    ax[i].xaxis.set_major_locator(MultipleLocator(100))
    ax[i].xaxis.set_minor_locator(AutoMinorLocator(10))

ax[0].yaxis.set_minor_locator(AutoMinorLocator(5))

ax[1].yaxis.set_major_locator(LogLocator(10))
ax[1].yaxis.set_minor_locator(LogLocator(base = 10.0, subs = (2, 3, 4, 5, 6, 7, 8, 9)))

#ax[1].yaxis.tick_right()
ax[1].yaxis.set_ticks_position('both')

ax[0].set_ylabel('Facular Contrast, [W / m$^2$ / nm]', fontsize = 12.5)
ax[1].set_ylabel('(NESSY - ATLAS) / ATLAS, [%]',       fontsize = 12.5)
ax[1].set_xlabel('Wavelength, [nm]',                   fontsize = 12.5)

leg0 = ax[0].legend(framealpha = 1, loc = 1, handletextpad = 1, prop = {'size': 20.5})
leg1 = ax[1].legend(framealpha = 1, loc = 4, handletextpad = 1, prop = {'size': 20.5})

for obj in leg0.legendHandles: obj.set_linewidth(3.0)
for obj in leg1.legendHandles: obj.set_linewidth(3.0)

auxplt.savepdf('var/fcontr_kur_nesatl')

fig, ax = plt.subplots(nrows = 2, ncols = 1, figsize = (12.0, 10.0))

#pltaux.figpar()

fig.tight_layout()

fig.suptitle('Facular contrast: NESSY (using FAL99 models) vs. ATLAS (using Kurucz models)', y = 1.01)

ax[0].plot(w, np.zeros(len(FQ4)), 'k--')

ax[0].set_xlim(110.5, 1000)

ls = ':'; lw = 1.5

#ax[0].plot(w, FQ4, color = 'b', linewidth = lw, label = 'NESSY (LTE)')
#ax[0].plot(w, FQ5, color = 'r', linewidth = lw, label = 'NESSY (NLTE)')
ax[0].plot(w, FQ5, color = 'r', linewidth = lw, label = 'NESSY')
ax[0].plot(w, FQ3, color = 'g', linewidth = lw, label = 'ATLAS')

#ax[1].plot(w, abs(RDFQ43), color = 'b',                 label = 'NESSY (LTE)  vs ATLAS')
#ax[1].plot(w, abs(RDFQ53), color = 'r', linewidth = lw, label = 'NESSY (NLTE) vs ATLAS')
ax[1].plot(w, abs(RDFQ53), color = 'r', linewidth = lw, label = 'NESSY vs ATLAS')

ax[1].plot(w, np.ones(len(RDFQ43)), 'k')

ax[1].set_yscale('log')

ax[1].set_xlim(110.5, 1000)
ax[1].set_ylim(1e-2, 1e+3)

for i in range(0, 2):

    ax[i].xaxis.set_major_locator(MultipleLocator(100))
    ax[i].xaxis.set_minor_locator(AutoMinorLocator(10))

ax[0].yaxis.set_minor_locator(AutoMinorLocator(5))

ax[1].yaxis.set_major_locator(LogLocator(10))
ax[1].yaxis.set_minor_locator(LogLocator(base = 10.0, subs = (2, 3, 4, 5, 6, 7, 8, 9)))

#ax[1].yaxis.tick_right()
ax[1].yaxis.set_ticks_position('both')

ax[0].set_ylabel('Facular Contrast, [W / m$^2$ / nm]', fontsize = 12.5)
ax[1].set_ylabel('|(NESSY - ATLAS) / ATLAS|, [%]',       fontsize = 12.5)
ax[1].set_xlabel('Wavelength, [nm]',                   fontsize = 12.5)

leg0 = ax[0].legend(framealpha = 1, loc = 1, handletextpad = 1, prop = {'size': 20.5})
#leg1 = ax[1].legend(framealpha = 1, loc = 4, handletextpad = 1, prop = {'size': 20.5})

for obj in leg0.legendHandles: obj.set_linewidth(3.0)
#for obj in leg1.legendHandles: obj.set_linewidth(3.0)

auxplt.savepdf('var/fcontr_fal_nesatl')

fig, ax = plt.subplots(nrows = 2, ncols = 1, figsize = (12.0, 10.0))

#pltaux.figpar()

fig.tight_layout()

fig.suptitle('Umbral contrast: NESSY (using Kurucz models) vs. ATLAS (using Kurucz models)', y = 1.01)

ax[0].plot(w, np.zeros(len(UQ1)), 'k--')

ax[0].set_xlim(110.5, 1000)

ls = ':'; lw = 1.5

#ax[0].plot(w, UQ1, color = 'b', linewidth = lw, label = 'NESSY (LTE)')
#ax[0].plot(w, UQ2, color = 'r', linewidth = lw, label = 'NESSY (NLTE)')
ax[0].plot(w, UQ2, color = 'r', linewidth = lw, label = 'NESSY')
ax[0].plot(w, UQ3, color = 'g', linewidth = lw, label = 'ATLAS')

#ax[1].plot(w, abs(RDUQ13), color = 'b',                 label = 'NESSY (LTE)  vs ATLAS')
#ax[1].plot(w, abs(RDUQ23), color = 'r', linewidth = lw, label = 'NESSY (NLTE) vs ATLAS')
ax[1].plot(w, abs(RDUQ23), color = 'r', linewidth = lw, label = 'NESSY vs ATLAS')

ax[1].plot(w, np.ones(len(RDUQ13)), 'k')

ax[1].set_yscale('log')

ax[1].set_xlim(110.5, 1000)
ax[1].set_ylim(1e-2, 1e+3)

for i in range(0, 2):

    ax[i].xaxis.set_major_locator(MultipleLocator(100))
    ax[i].xaxis.set_minor_locator(AutoMinorLocator(10))

ax[0].yaxis.set_minor_locator(AutoMinorLocator(5))

ax[1].yaxis.set_major_locator(LogLocator(10))
ax[1].yaxis.set_minor_locator(LogLocator(base = 10.0, subs = (2, 3, 4, 5, 6, 7, 8, 9)))

ax[1].yaxis.set_ticks_position('both')

ax[0].set_ylabel('Umbral Contrast, [W / m$^2$ / nm]', fontsize = 12.5)
ax[1].set_ylabel('|(NESSY - ATLAS) / ATLAS|, [%]',    fontsize = 12.5)
ax[1].set_xlabel('Wavelength, [nm]',                  fontsize = 12.5)

leg0 = ax[0].legend(framealpha = 1, loc = 4, handletextpad = 1, prop = {'size': 20.5})
#leg1 = ax[1].legend(framealpha = 1, loc = 1, handletextpad = 1, prop = {'size': 20.5})

for obj in leg0.legendHandles: obj.set_linewidth(3.0)
#for obj in leg1.legendHandles: obj.set_linewidth(3.0)

auxplt.savepdf('var/ucontr_kur_nesatl')

fig, ax = plt.subplots(nrows = 2, ncols = 1, figsize = (12.0, 10.0))

#pltaux.figpar()

fig.tight_layout()

fig.suptitle('Penumbral contrast: NESSY (using Kurucz models) vs. ATLAS (using Kurucz models)', y = 1.01)

ax[0].plot(w, np.zeros(len(PQ1)), 'k--')

ax[0].set_xlim(110.5, 1000)

ls = ':'; lw = 1.5

#ax[0].plot(w, PQ1, color = 'b', linewidth = lw, label = 'NESSY (LTE)')
#ax[0].plot(w, PQ2, color = 'r', linewidth = lw, label = 'NESSY (NLTE)')
ax[0].plot(w, PQ2, color = 'r', linewidth = lw, label = 'NESSY')
ax[0].plot(w, PQ3, color = 'g', linewidth = lw, label = 'ATLAS')

#ax[1].plot(w, abs(RDPQ13), color = 'b',                 label = 'NESSY (LTE)  vs ATLAS')
#ax[1].plot(w, abs(RDPQ23), color = 'r', linewidth = lw, label = 'NESSY (NLTE) vs ATLAS')
ax[1].plot(w, abs(RDPQ23), color = 'r', linewidth = lw, label = 'NESSY vs ATLAS')

ax[1].plot(w, np.ones(len(RDPQ13)), 'k')

ax[1].set_yscale('log')

ax[1].set_xlim(110.5, 1000)
ax[1].set_ylim(1e-2, 1e+3)

for i in range(0, 2):

    ax[i].xaxis.set_major_locator(MultipleLocator(100))
    ax[i].xaxis.set_minor_locator(AutoMinorLocator(10))

ax[0].yaxis.set_minor_locator(AutoMinorLocator(5))

ax[1].yaxis.set_major_locator(LogLocator(10))
ax[1].yaxis.set_minor_locator(LogLocator(base = 10.0, subs = (2, 3, 4, 5, 6, 7, 8, 9)))

ax[1].yaxis.set_ticks_position('both')

ax[0].set_ylabel('Penumbral Contrast, [W / m$^2$ / nm]', fontsize = 12.5)
ax[1].set_ylabel('|(NESSY - ATLAS) / ATLAS|, [%]',       fontsize = 12.5)
ax[1].set_xlabel('Wavelength, [nm]',                     fontsize = 12.5)

leg0 = ax[0].legend(framealpha = 1, loc = 4, handletextpad = 1, prop = {'size': 20.5})
#leg1 = ax[1].legend(framealpha = 1, loc = 1, handletextpad = 1, prop = {'size': 20.5})

for obj in leg0.legendHandles: obj.set_linewidth(3.0)
#for obj in leg1.legendHandles: obj.set_linewidth(3.0)

auxplt.savepdf('var/pcontr_kur_nesatl')

fig, ax = plt.subplots(nrows = 2, ncols = 1, figsize = (12.0, 10.0))

#pltaux.figpar()

fig.tight_layout()

fig.suptitle('Kurucz quiet sun model: NESSY vs. ATLAS', y = 1.01)

ax[0].plot(w, np.zeros(len(Q1)), 'k--')

ax[0].set_xlim(110.5, 1000)

ls = ':'; lw = 1.5

ax[0].plot(w, Q1, color = 'b', linewidth = lw, label = 'NESSY (LTE)')
ax[0].plot(w, Q2, color = 'r', linewidth = lw, label = 'NESSY (NLTE)')
ax[0].plot(w, Q3, color = 'g', linewidth = lw, label = 'ATLAS')

ax[1].plot(w, abs(RDQ13), color = 'b',                 label = 'NESSY (LTE)  vs ATLAS')
ax[1].plot(w, abs(RDQ23), color = 'r', linewidth = lw, label = 'NESSY (NLTE) vs ATLAS')

ax[1].plot(w, np.ones(len(RDQ13)), 'k')

ax[1].set_yscale('log')

ax[1].set_xlim(110.5, 1000)
ax[1].set_ylim(1e-2, 1e+3)

for i in range(0, 2):

    ax[i].xaxis.set_major_locator(MultipleLocator(100))
    ax[i].xaxis.set_minor_locator(AutoMinorLocator(10))

ax[0].yaxis.set_minor_locator(AutoMinorLocator(5))

ax[1].yaxis.set_major_locator(LogLocator(10))
ax[1].yaxis.set_minor_locator(LogLocator(base = 10.0, subs = (2, 3, 4, 5, 6, 7, 8, 9)))

#ax[1].yaxis.tick_right()
ax[1].yaxis.set_ticks_position('both')

ax[0].set_ylabel('Flux, [W / m$^2$ / nm]',       fontsize = 12.5)
ax[1].set_ylabel('(NESSY - ATLAS) / ATLAS, [%]', fontsize = 12.5)
ax[1].set_xlabel('Wavelength, [nm]',             fontsize = 12.5)

leg0 = ax[0].legend(framealpha = 1, loc = 1, handletextpad = 1, prop = {'size': 20.5})
leg1 = ax[1].legend(framealpha = 1, loc = 1, handletextpad = 1, prop = {'size': 20.5})

for obj in leg0.legendHandles: obj.set_linewidth(3.0)
for obj in leg1.legendHandles: obj.set_linewidth(3.0)

auxplt.savepdf('var/Q_kur_nesatl')

fig, ax = plt.subplots(nrows = 2, ncols = 1, figsize = (12.0, 10.0))

#pltaux.figpar()

fig.tight_layout()

fig.suptitle('Kurucz facula model: NESSY vs. ATLAS', y = 1.01)

ax[0].plot(w, np.zeros(len(F1)), 'k--')

ax[0].set_xlim(110.5, 1000)

ls = ':'; lw = 1.5

ax[0].plot(w, F1, color = 'b', linewidth = lw, label = 'NESSY (LTE)')
ax[0].plot(w, F2, color = 'r', linewidth = lw, label = 'NESSY (NLTE)')
ax[0].plot(w, F3, color = 'g', linewidth = lw, label = 'ATLAS')

ax[1].plot(w, abs(RDF13), color = 'b',                 label = 'NESSY (LTE)  vs ATLAS')
ax[1].plot(w, abs(RDF23), color = 'r', linewidth = lw, label = 'NESSY (NLTE) vs ATLAS')

ax[1].plot(w, np.ones(len(RDF13)), 'k')

ax[1].set_yscale('log')

ax[1].set_xlim(110.5, 1000)
ax[1].set_ylim(1e-2,  1e+3)

for i in range(0, 2):

    ax[i].xaxis.set_major_locator(MultipleLocator(100))
    ax[i].xaxis.set_minor_locator(AutoMinorLocator(10))

ax[0].yaxis.set_minor_locator(AutoMinorLocator(5))

ax[1].yaxis.set_major_locator(LogLocator(10))
ax[1].yaxis.set_minor_locator(LogLocator(base = 10.0, subs = (2, 3, 4, 5, 6, 7, 8, 9)))

ax[1].yaxis.set_ticks_position('both')

ax[0].set_ylabel('Flux, [W / m$^2$ / nm]',       fontsize = 12.5)
ax[1].set_ylabel('(NESSY - ATLAS) / ATLAS, [%]', fontsize = 12.5)
ax[1].set_xlabel('Wavelength, [nm]',             fontsize = 12.5)

leg0 = ax[0].legend(framealpha = 1, loc = 1, handletextpad = 1, prop = {'size': 20.5})
leg1 = ax[1].legend(framealpha = 1, loc = 1, handletextpad = 1, prop = {'size': 20.5})

for obj in leg0.legendHandles: obj.set_linewidth(3.0)
for obj in leg1.legendHandles: obj.set_linewidth(3.0)

auxplt.savepdf('var/F_kur_nesatl')

fig, ax = plt.subplots(nrows = 2, ncols = 1, figsize = (12.0, 10.0))

#pltaux.figpar()

fig.tight_layout()

fig.suptitle('Kurucz penumbra model: NESSY vs. ATLAS', y = 1.01)

ax[0].plot(w, np.zeros(len(P1)), 'k--')

ax[0].set_xlim(110.5, 1000)

ls = ':'; lw = 1.5

ax[0].plot(w, P1, color = 'b', linewidth = lw, label = 'NESSY (LTE)')
ax[0].plot(w, P2, color = 'r', linewidth = lw, label = 'NESSY (NLTE)')
ax[0].plot(w, P3, color = 'g', linewidth = lw, label = 'ATLAS')

ax[1].plot(w, abs(RDP13), color = 'b',                 label = 'NESSY (LTE)  vs ATLAS')
ax[1].plot(w, abs(RDP23), color = 'r', linewidth = lw, label = 'NESSY (NLTE) vs ATLAS')

ax[1].plot(w, np.ones(len(RDP13)), 'k')

ax[1].set_yscale('log')

ax[1].set_xlim(110.5, 1000)
ax[1].set_ylim(1e-2,  1e+3)

for i in range(0, 2):

    ax[i].xaxis.set_major_locator(MultipleLocator(100))
    ax[i].xaxis.set_minor_locator(AutoMinorLocator(10))

ax[0].yaxis.set_minor_locator(AutoMinorLocator(5))

ax[1].yaxis.set_major_locator(LogLocator(10))
ax[1].yaxis.set_minor_locator(LogLocator(base = 10.0, subs = (2, 3, 4, 5, 6, 7, 8, 9)))

ax[1].yaxis.set_ticks_position('both')

ax[0].set_ylabel('Flux, [W / m$^2$ / nm]',       fontsize = 12.5)
ax[1].set_ylabel('(NESSY - ATLAS) / ATLAS, [%]', fontsize = 12.5)
ax[1].set_xlabel('Wavelength, [nm]',             fontsize = 12.5)

leg0 = ax[0].legend(framealpha = 1, loc = 1, handletextpad = 1, prop = {'size': 20.5})
leg1 = ax[1].legend(framealpha = 1, loc = 1, handletextpad = 1, prop = {'size': 20.5})

for obj in leg0.legendHandles: obj.set_linewidth(3.0)
for obj in leg1.legendHandles: obj.set_linewidth(3.0)

auxplt.savepdf('var/P_kur_nesatl')
