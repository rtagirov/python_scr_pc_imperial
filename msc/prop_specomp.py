import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

from matplotlib.ticker import AutoMinorLocator
from matplotlib.ticker import MultipleLocator

import importlib
import sys
import os

from tqdm import tqdm

if not '../aux/' in sys.path:

    sys.path.append('../aux/')

import paths;   importlib.reload(paths)
import spec;    importlib.reload(spec)
import nessy;   importlib.reload(nessy)
import auxsys;  importlib.reload(auxsys)
import auxplt;  importlib.reload(auxplt)
import auxfunc; importlib.reload(auxfunc)

#w, Icf = nessy.read_spec(paths.it1f + 'prop/C/full',   wvl1 = 1005, wvl2 = 10000)
#w, Icm = nessy.read_spec(paths.it1f + 'prop/C/nlte11', wvl1 = 1005, wvl2 = 10000)
#w, Icl = nessy.read_spec(paths.it0f + 'prop/C/lte',    wvl1 = 1005, wvl2 = 10000)
#w, Iff = nessy.read_spec(paths.it1f + 'prop/P/full',   wvl1 = 1005, wvl2 = 10000)
#w, Ifm = nessy.read_spec(paths.it1f + 'prop/P/nlte11', wvl1 = 1005, wvl2 = 10000)
#w, Ifl = nessy.read_spec(paths.it0f + 'prop/P/lte',    wvl1 = 1005, wvl2 = 10000)

#w = w / 10.0

#ws, Icfs = spec.mean(w, Icf, 1.0)
#ws, Icms = spec.mean(w, Icm, 1.0)
#ws, Icls = spec.mean(w, Icl, 1.0)
#ws, Iffs = spec.mean(w, Iff, 1.0)
#ws, Ifms = spec.mean(w, Ifm, 1.0)
#ws, Ifls = spec.mean(w, Ifl, 1.0)

#np.savez(paths.npz + 'w', w = ws)

#np.savez(paths.npz + 'Icfs', I = Icfs)
#np.savez(paths.npz + 'Icms', I = Icms)
#np.savez(paths.npz + 'Icls', I = Icls)
#np.savez(paths.npz + 'Iffs', I = Iffs)
#np.savez(paths.npz + 'Ifms', I = Ifms)
#np.savez(paths.npz + 'Ifls', I = Ifls)

w = np.load(paths.npz + 'w.npz')['w']

Icfs = np.load(paths.npz + 'Icfs.npz')['I']
Icms = np.load(paths.npz + 'Icms.npz')['I']
Icls = np.load(paths.npz + 'Icls.npz')['I']
Iffs = np.load(paths.npz + 'Iffs.npz')['I']
Ifms = np.load(paths.npz + 'Ifms.npz')['I']
Ifls = np.load(paths.npz + 'Ifls.npz')['I']

plt.close('all')

auxplt.figpar(fontsize = 20)

fig, ax = plt.subplots(nrows = 1, ncols = 1, figsize = (12.0, 10.0))

fig.tight_layout()

ax.set_xlabel('Wavelength, [nm]', labelpad = 15)
ax.set_ylabel(r'$F_\mathrm{NLTE} / F_\mathrm{LTE}$')

ax.plot(w, Iffs / Ifls, color = 'r', linewidth = 4, label = 'NLTE 30')
ax.plot(w, Ifms / Ifls, color = 'c', linewidth = 2, label = 'NLTE 11')

ax.plot(np.ones(len(Icls)) * 1.0, '--', color = 'k', linewidth = 0.5)

ax.set_xlim(100, 400)

ax.set_ylim(0.0, 1.25)

ax.xaxis.set_minor_locator(ticker.AutoMinorLocator(5))

leg = ax.legend(framealpha = 1, loc = 'best', handletextpad = 1, prop = {'size': 25})

auxplt.savepdf('prop_specomp')
