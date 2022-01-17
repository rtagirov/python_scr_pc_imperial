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

lambda1 = 1005
lambda2 = 10000

#wvl, faldef = nessy.read_spec(paths.it1f + 'runtime/def',               wvl1 = lambda1, wvl2 = lambda2)
#wvl, falacc = nessy.read_spec(paths.it1f + 'runtime/def_fioss_acc',     wvl1 = lambda1, wvl2 = lambda2)
#wvl, kurdef = nessy.read_spec(paths.it1f + 'runtime/def_kur',           wvl1 = lambda1, wvl2 = lambda2)
#wvl, kuracc = nessy.read_spec(paths.it1f + 'runtime/def_kur_fioss_acc', wvl1 = lambda1, wvl2 = lambda2)

#wvl = wvl / 10.0

#wvls, faldef_s = spec.mean_within_delta(wvl, faldef, 1.0)
#wvls, falacc_s = spec.mean_within_delta(wvl, falacc, 1.0)
#wvls, kurdef_s = spec.mean_within_delta(wvl, kurdef, 1.0)
#wvls, kuracc_s = spec.mean_within_delta(wvl, kuracc, 1.0)

outname = 'fioss_acc_dev'

#np.savez(paths.npz + outname, w = wvls, f1 = faldef_s,
#                                        f2 = falacc_s,
#                                        f3 = kurdef_s,
#                                        f4 = kuracc_s)

specs = np.load(paths.npz + outname + '.npz')

wvl =    specs['w']
faldef = specs['f1']
falacc = specs['f2']
kurdef = specs['f3']
kuracc = specs['f4']

plt.close('all')

fig, ax = plt.subplots(nrows = 1, ncols = 1, figsize = (6.0, 3.0))

fig.tight_layout()

ax.plot(wvl, kurdef / kuracc, color = 'r', label = 'Kurucz')
ax.plot(wvl, faldef / falacc, color = 'b', label = 'FAL99-C')

ax.set_xlim(100, 1000)

ax.xaxis.set_major_locator(MultipleLocator(100))
#ax.yaxis.set_minor_locator(AutoMinorLocator(5))
ax.xaxis.set_minor_locator(AutoMinorLocator(4))

ax.set_ylabel('Ratio')
ax.set_xlabel('Wavelength, nm')

leg0 = ax.legend(framealpha = 1, loc = 4, handletextpad = 1, prop = {'size':7.5})

for obj in leg0.legendHandles: obj.set_linewidth(3.0)

auxplt.savepdf('fioss_acc_dev')
