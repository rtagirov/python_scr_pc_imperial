import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as const

from matplotlib.ticker import AutoMinorLocator

import sys

if not '../aux/' in sys.path: sys.path.append('../aux/')

import importlib
import auxplt
import auxsys
import paths

importlib.reload(auxplt)
importlib.reload(auxsys)
importlib.reload(paths)

def read_cfunc(pathtofile, cutoff1, cutoff2):

    f = open(paths.it1f + pathtofile, 'r')

    wvl = []

    cc = []
    cl = []

    for line in f:

        elems = line.split()

        if len(elems) == 1:

            w = float(elems[0]) / 10.0

            if w >= cutoff2: break
#            if w >  cutoff1:  wvl.append(w)
            if w ==  cutoff1:  wvl.append(w)

#        if len(elems) == 3 and w > cutoff1:
        if len(elems) == 3 and w == cutoff1:

            cc.append(float(elems[1]))
            cl.append(float(elems[2]))

        if len(elems) > 3 and w == cutoff1:

            auxsys.abort('Wrong line. Abort.')

    gs = int(len(cc) / len(wvl))

    cfunc_c = np.array([cc[i : i + gs] for i in range(0, len(cc), gs)])
    cfunc_l = np.array([cl[i : i + gs] for i in range(0, len(cl), gs)])

    return wvl, cfunc_c, cfunc_l

#wvl, cq, lq = read_cfunc('var/Q/fal/contr.txt', 12.5, 401.0)
#wvl, cf, lf = read_cfunc('var/F/fal/contr.txt', 12.5, 401.0)

wvl, cq, lq = read_cfunc('var/Q/fal_contr/contr.func', 200.5, 401.0)
wvl, cf, lf = read_cfunc('var/F/fal_contr/contr.func', 200.5, 401.0)

cq[0] = np.flip(cq[0], 0)
cf[0] = np.flip(cf[0], 0)

lq[0] = np.flip(lq[0], 0)
lf[0] = np.flip(lf[0], 0)

hq, Tq, x, x, x = np.loadtxt(paths.it1f + 'var/Q/fal/atm.inp', unpack = True)
hf, Tf, x, x, x = np.loadtxt(paths.it1f + 'var/F/fal/atm.inp', unpack = True)

hq = np.delete(hq, len(hq) - 1)
hf = np.delete(hf, len(hf) - 1)
Tq = np.delete(Tq, len(Tq) - 1)
Tf = np.delete(Tf, len(Tf) - 1)

plt.close('all')

#fig, ax1 = plt.subplots(nrows = 1, ncols = 1, figsize = (const.golden * 5, 5))
fig, ax1 = plt.subplots(nrows = 1, ncols = 2, figsize = (12, 5))

plt.subplots_adjust(wspace = 0.50, hspace = 0.5)

auxplt.figpar(3, 3, 15)

#fig.tight_layout()

ax1[0].plot(hq, cq[0] / max(cq[0]), color = 'k' ,label = 'Quiet Sun')
ax1[0].plot(hf, cf[0] / max(cf[0]), color = 'r', label = 'Facula')

ax1[0].set_xlabel('Height, [km]')
ax1[0].set_ylabel('Normalized contribution function at 200 nm')

ax1[0].set_xlim(0, 1000)
ax1[0].set_ylim(0, 1)

#ax1[0].axvline(x = 270., ymin = 0.225, ymax = 1,  linestyle = ':', color = 'r')
ax1[0].axvline(x = 270., ymin = 0.0, ymax = 1,  linestyle = ':', color = 'r')
ax1[0].axhline(y = 0.22, xmin = 0.27, xmax = 1, linestyle = ':', color = 'r')

#ax1[0].axvline(x = 350., ymin = 0.15, ymax = 1,  linestyle = ':', color = 'k')
ax1[0].axvline(x = 350., ymin = 0.0, ymax = 1,  linestyle = ':', color = 'k')
ax1[0].axhline(y = 0.15, xmin = 0.35, xmax = 1, linestyle = ':', color = 'k')

leg = ax1[0].legend(framealpha = 1, loc = 1, handletextpad = 1, prop = {'size': 15.0})

ax2 = ax1[0].twinx()

ax2.plot(hq, Tq, color = 'k')
ax2.plot(hf, Tf, color = 'r')

ax2.yaxis.set_minor_locator(AutoMinorLocator(10))

ax2.set_xlim(0, 1000)
ax2.set_ylim(4000, 10000)

ax2.set_ylabel('Temperature, [K]')

for obj in leg.legendHandles: obj.set_linewidth(3.0)

ax1[0].set_title(r'NESSY-FAL99, $\mu = 1.0$')

ax1[1].plot(hq, lq[0] / max(lq[0]), color = 'k' ,label = 'Quiet Sun')
ax1[1].plot(hf, lf[0] / max(lf[0]), color = 'r', label = 'Facula')

ax1[1].set_xlabel('Height, [km]')
ax1[1].set_ylabel('Normalized contribution function at 200 nm')

ax1[1].set_xlim(0, 1000)
ax1[1].set_ylim(0, 1)

ax1[1].axvline(x = 390., ymin = 0.0, ymax = 1,  linestyle = ':', color = 'r')
ax1[1].axhline(y = 0.18, xmin = 0.38, xmax = 1, linestyle = ':', color = 'r')

ax1[1].axvline(x = 450., ymin = 0.0, ymax = 1,  linestyle = ':', color = 'k')
ax1[1].axhline(y = 0.115, xmin = 0.45, xmax = 1, linestyle = ':', color = 'k')

leg = ax1[1].legend(framealpha = 1, loc = 1, handletextpad = 1, prop = {'size': 15.0})

ax2 = ax1[1].twinx()

ax2.plot(hq, Tq, color = 'k')
ax2.plot(hf, Tf, color = 'r')

ax2.yaxis.set_minor_locator(AutoMinorLocator(10))

ax2.set_xlim(0, 1000)
ax2.set_ylim(4000, 10000)

ax2.set_ylabel('Temperature, [K]')

for obj in leg.legendHandles: obj.set_linewidth(3.0)

ax1[1].set_title(r'NESSY-FAL99, $\mu = 0.31$')

auxplt.savepdf('var/contrib_func')
