import numpy             as np
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

hkurq = np.loadtxt(paths.it0h + 'var/Q/kur/ATM_STR', usecols = [1], skiprows = 2) - 134.416 # 134.237
tkurq = np.loadtxt(paths.it0h + 'var/Q/kur/ATM_STR', usecols = [2], skiprows = 2)

rkurq = 10.0**np.loadtxt(paths.it1h + 'var/Q/kur/tau_ross.dat', usecols = [7], skiprows = 1)

hkurf = np.loadtxt(paths.it0h + 'var/F/kur/ATM_STR', usecols = [1], skiprows = 2) - 104.806  #104.736
tkurf = np.loadtxt(paths.it0h + 'var/F/kur/ATM_STR', usecols = [2], skiprows = 2)

rkurf = 10.0**np.loadtxt(paths.it1h + 'var/F/kur/tau_ross.dat', usecols = [7], skiprows = 1)

hkuru = np.loadtxt(paths.it0h + 'var/U/kur/ATM_STR', usecols = [1], skiprows = 2) - 455.054  #452.537
tkuru = np.loadtxt(paths.it0h + 'var/U/kur/ATM_STR', usecols = [2], skiprows = 2)

rkuru = 10.0**np.loadtxt(paths.it1h + 'var/U/kur/tau_ross.dat', usecols = [7], skiprows = 1)

hkurp = np.loadtxt(paths.it0h + 'var/P/kur/ATM_STR', usecols = [1], skiprows = 2) - 428.559  # 427.693
tkurp = np.loadtxt(paths.it0h + 'var/P/kur/ATM_STR', usecols = [2], skiprows = 2)

rkurp = 10.0**np.loadtxt(paths.it1h + 'var/P/kur/tau_ross.dat', usecols = [7], skiprows = 1)

hfalq = np.loadtxt(paths.it0h + 'var/Q/fal/ATM_STR', usecols = [1], skiprows = 2) - 114.388  # 114.227
tfalq = np.loadtxt(paths.it0h + 'var/Q/fal/ATM_STR', usecols = [2], skiprows = 2)

rfalq = 10.0**np.loadtxt(paths.it1h + 'var/Q/fal/tau_ross.dat', usecols = [7], skiprows = 1)

hfalf = np.loadtxt(paths.it0h + 'var/F/fal/ATM_STR', usecols = [1], skiprows = 2) - 106.802 # 106.725
tfalf = np.loadtxt(paths.it0h + 'var/F/fal/ATM_STR', usecols = [2], skiprows = 2)

rfalf = 10.0**np.loadtxt(paths.it1h + 'var/F/fal/tau_ross.dat', usecols = [7], skiprows = 1)

plt.close('all')

fig, ax = plt.subplots(nrows = 2, ncols = 1, figsize = (6.0, 7.0))

fig.tight_layout()

plt.subplots_adjust(hspace = 0.15)

auxplt.figpar(3, 3, 15)

ax[0].plot(hfalq, tfalq, color = 'k', linewidth = 3, label = 'Quiet Sun (FAL99-C)')
ax[0].plot(hfalf, tfalf, color = 'r', label = 'Facula (FAL99-P)')
ax[0].plot(hkurq, tkurq, color = 'c', label = 'Quiet Sun (U99)', linestyle = '--')

ax[0].plot(hkurp, tkurp, color = 'm', label = 'Penumbra (U99)', linestyle = '--')
ax[0].plot(hkuru, tkuru, color = 'b', label = 'Umbra (U99)', linestyle = '--')

#ax[0].tick_params(labelbottom = 'off')

ax[0].set_ylabel('Temperature, [K]')
ax[1].set_ylabel('Temperature, [K]')
ax[1].set_xlabel('Height, [km]')

ax[0].set_xlim(-500, 2500)
ax[0].set_ylim(2500, 10100)

ax[0].xaxis.set_minor_locator(AutoMinorLocator(5))
ax[0].yaxis.set_minor_locator(AutoMinorLocator(4))
ax[1].xaxis.set_minor_locator(AutoMinorLocator(4))
ax[1].yaxis.set_minor_locator(AutoMinorLocator(4))

ax[0].axvline(x = 0.0, linestyle = ':', color = 'k')
ax[1].axvline(x = 0.0, linestyle = ':', color = 'k')

ax[1].plot(hfalq, tfalq, color = 'k', linewidth = 2, label = 'Quiet Sun (FAL99-C)')
ax[1].plot(hfalf, tfalf, color = 'r', linewidth = 2, label = 'Facula (FAL99-P)')
ax[1].plot(hkurq, tkurq, color = 'c', linewidth = 2, label = 'Quiet Sun (U99)',   linestyle = '--')
ax[1].plot(hkurf, tkurf, color = 'g', linewidth = 2, label = 'Facula (U99)', linestyle = '--')

ax[1].set_xlim(-160, 800)
ax[1].set_ylim(3750, 10100)

ax[1].fill_between(np.array([130, 395]), 3750, 10100, facecolor = 'gray', linewidth = 0, alpha = 0.5)

leg0 = ax[0].legend(framealpha = 1, loc = 1, bbox_to_anchor = (0.559, 0.9999), handletextpad = 1, prop = {'size': 10.5})
leg1 = ax[1].legend(framealpha = 1, loc = 1,                                   handletextpad = 1, prop = {'size': 10.5})

for obj in leg0.legendHandles: obj.set_linewidth(3.0)
for obj in leg1.legendHandles: obj.set_linewidth(3.0)

auxplt.savepdf('var/atm_mod_2_od')

sys.exit()

fig, ax = plt.subplots(nrows = 2, ncols = 1, figsize = (12.0, 10.0))

ax[0].set_title('Quiet Sun')

ax[0].plot(hkurq, tkurq, label = 'Kurucz')
ax[0].plot(hfalq, tfalq, label = 'FAL')

ax[1].set_title('Facula')

ax[1].plot(hkurf, tkurf, label = 'Kurucz')
ax[1].plot(hfalf, tfalf, label = 'FAL')

ax[0].set_xlabel('Height, [km]')
ax[1].set_xlabel('Height, [km]')

ax[0].set_ylabel('Temperature, [K]')
ax[1].set_ylabel('Temperature, [K]')

ax[0].set_ylim(3500, 10000)
ax[1].set_ylim(4500, 9000)

ax[0].axvline(x = 0.0, linestyle = '--', color = 'k')
ax[1].axvline(x = 0.0, linestyle = '--', color = 'k')

auxplt.savepdf('atmcomp_falkur')

fig, ax = plt.subplots(nrows = 2, ncols = 1, figsize = (12.0, 10.0))

ax[0].set_title('FAL99 models')

ax[0].plot(hfalq, tfalq, color = 'k', label = 'Quiet Sun')
ax[0].plot(hfalf, tfalf, color = 'r', label = 'Facula')

ax[1].set_title('Kurucz models')

ax[1].plot(hkurq, tkurq, color = 'k', label = 'Quiet Sun')
ax[1].plot(hkurf, tkurf, color = 'r', label = 'Facula')
ax[1].plot(hkuru, tkuru, color = 'b', label = 'Umbra')
ax[1].plot(hkurp, tkurp, color = 'm', label = 'Penumbra')

ax[0].set_xlabel('Height, [km]')
ax[1].set_xlabel('Height, [km]')

ax[0].set_ylabel('Temperature, [K]')
ax[1].set_ylabel('Temperature, [K]')

ax[0].set_ylim(4000, 10000)
ax[1].set_ylim(2500, 10100)

ax[0].axvline(x = 0.0, linestyle = '--', color = 'k')
ax[1].axvline(x = 0.0, linestyle = '--', color = 'k')

leg0 = ax[0].legend(framealpha = 1, loc = 4, handletextpad = 1, prop = {'size': 15.5})
leg1 = ax[1].legend(framealpha = 1, loc = 1, handletextpad = 1, prop = {'size': 15.5})

for obj in leg0.legendHandles: obj.set_linewidth(3.0)
for obj in leg1.legendHandles: obj.set_linewidth(3.0)

auxplt.savepdf('atm_mod_old')

fig, ax = plt.subplots(nrows = 1, ncols = 1, figsize = (6.0, 5.0))

fig.tight_layout()

auxplt.figpar(3, 3, 15)

ax.plot(hfalq, tfalq, color = 'k', linewidth = 3, label = 'Quiet Sun (FAL-C)')
ax.plot(hfalf, tfalf, color = 'r', label = 'Facula (FAL-P)')
ax.plot(hkurq, tkurq, color = 'c', label = 'Quiet Sun (Kurucz)')

ax.plot(hkurp, tkurp, color = 'm', label = 'Penumbra (Kurucz)')
ax.plot(hkuru, tkuru, color = 'b', label = 'Umbra (Kurucz)')

ax.set_xlabel('Height, [km]')
ax.set_ylabel('Temperature, [K]')

ax.set_xlim(-500, 2500)
ax.set_ylim(2500, 10100)

ax.xaxis.set_minor_locator(AutoMinorLocator(5))
ax.yaxis.set_minor_locator(AutoMinorLocator(4))

ax.axvline(x = 0.0, linestyle = '--', color = 'k')

leg = ax.legend(framealpha = 1, loc = 1, bbox_to_anchor = (0.62, 0.99), handletextpad = 1, prop = {'size': 13.0})

for obj in leg.legendHandles: obj.set_linewidth(3.0)

auxplt.savepdf('atm_mod')

fig, ax = plt.subplots(nrows = 2, ncols = 1, figsize = (6.0, 7.0))

fig.tight_layout()

plt.subplots_adjust(hspace = 0.15)

auxplt.figpar(3, 3, 15)

ax[0].plot(rfalq, tfalq[1 : len(tfalq)], color = 'k', linewidth = 3, label = 'Quiet Sun (FAL99-C)')
ax[0].plot(rfalf, tfalf[1 : len(tfalf)], color = 'r', label = 'Facula (FAL99-P)')
ax[0].plot(rkurq, tkurq[1 : len(tkurq)], color = 'c', label = 'Quiet Sun (U99)', linestyle = '--')

ax[0].plot(rkurp, tkurp[1 : len(tkurp)], color = 'm', label = 'Penumbra (U99)', linestyle = '--')
ax[0].plot(rkuru, tkuru[1 : len(tkuru)], color = 'b', label = 'Umbra (U99)', linestyle = '--')

ax[0].set_xscale('log')

#ax[0].tick_params(labelbottom = 'off')

ax[0].set_ylabel('Temperature, [K]')
ax[1].set_ylabel('Temperature, [K]')
ax[1].set_xlabel(r'$\tau_\mathrm{Ross}$')

ax[0].set_xlim(1e2, 1e-10)
ax[0].set_ylim(2500, 10100)

#ax[0].xaxis.set_minor_locator(AutoMinorLocator(5))
#ax[0].yaxis.set_minor_locator(AutoMinorLocator(4))
#ax[1].xaxis.set_minor_locator(AutoMinorLocator(4))
#ax[1].yaxis.set_minor_locator(AutoMinorLocator(4))

ax[0].axvline(x = 1.0, linestyle = '--', color = 'k', linewidth = 0.3)
ax[1].axvline(x = 1.0, linestyle = '--', color = 'k', linewidth = 0.3)

ax[1].plot(rfalq, tfalq[1 : len(tfalq)], color = 'k', linewidth = 2, label = 'Quiet Sun (FAL99-C)')
ax[1].plot(rfalf, tfalf[1 : len(tfalf)], color = 'r', linewidth = 2, label = 'Facula (FAL99-P)')
ax[1].plot(rkurq, tkurq[1 : len(tkurq)], color = 'c', linewidth = 2, label = 'Quiet Sun (U99)',   linestyle = '--')
ax[1].plot(rkurf, tkurf[1 : len(tkurf)], color = 'g', linewidth = 2, label = 'Facula (U99)', linestyle = '--')

ax[1].set_xscale('log')

ax[1].set_xlim(1e2, 1e-10)
ax[1].set_ylim(3750, 10100)

leg0 = ax[0].legend(framealpha = 1, loc = 1, bbox_to_anchor = (0.489, 0.9999), handletextpad = 1, prop = {'size': 8.5})
leg1 = ax[1].legend(framealpha = 1, loc = 1, bbox_to_anchor = (0.559, 0.9999), handletextpad = 1, prop = {'size': 10.0})

for obj in leg0.legendHandles: obj.set_linewidth(3.0)
for obj in leg1.legendHandles: obj.set_linewidth(3.0)

auxplt.savepdf('var/atm_mod_3')
