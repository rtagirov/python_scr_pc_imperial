import numpy as np
import matplotlib.pyplot as plt

from matplotlib.ticker import MultipleLocator
from matplotlib.ticker import AutoMinorLocator

import sys

if not '/mnt/SSD/sim/python/src/aux' in sys.path:

    sys.path.append('/mnt/SSD/sim/python/src/aux')

import pandas as pd

import importlib

import auxplt; importlib.reload(auxplt)

def conv(t):

    result = []

    for elem in t:

        (h, m, s) = elem.split(':')

        res = int(h) * 60 + int(m) + int(s) / 60

        result.append(res)

    return np.asarray(result)

runs = ['falc', 'fale', 'falf', 'falp', 'fals', 'vala', 'kurf', 'kurq', 't5400g4', 'fchhtb', 'mur']
#runs = ['vala', 'kurf', 'kurq', 'kurt5400g4', 'fale', 'falf', 'falp', 'fals', 'falc', 'fchhtb', 'mur']

cp = np.zeros((320, 11))
ti = np.zeros((320, 11))

for i, r in enumerate(runs):

    data = np.genfromtxt('./test_fdfh_11/' + r + '/CONV/ALL', dtype = None)

    df = pd.DataFrame(data)

    if i == 0:

        li = df.loc[:, 0].tolist()

        li = [elem.decode('utf-8') for elem in li]

        li.pop(0)

        li = [int(elem) for elem in li]

        li = np.asarray(li)

    time = df.loc[:, 1].tolist()

    time = [elem.decode('utf-8') for elem in time]

    time.pop(0)

    ti[:, i] = conv(time)

    corm = df.loc[:, 2].tolist()

    corm = [elem.decode('utf-8') for elem in corm]

    corm.pop(0)

    corm = [float(elem) for elem in corm]

    cp[:, i] = np.asarray(corm)

plt.close('all')

fig, ax = plt.subplots(nrows = 1, ncols = 1, figsize = (5, 5))

#col = ['gray' for k in range(len(cormax[0, :]) - 2)]

ax.plot(ti[:, 5],  cp[:, 5],  label = 'VAL81-A, (ND = 51)')
ax.plot(ti[:, 7],  cp[:, 7],  label = 'Kurucz (Facula, ND = 56)')
ax.plot(ti[:, 6],  cp[:, 6],  label = 'Kurucz (QS, ND = 63)')
ax.plot(ti[:, 8],  cp[:, 8],  label = 'Kurucz (Teff = 5400, g = 4, ND = 72)')
ax.plot(ti[:, 1],  cp[:, 1],  label = 'FAL99-E, (ND = 79)')
ax.plot(ti[:, 2],  cp[:, 2],  label = 'FAL99-F, (ND = 79)')
ax.plot(ti[:, 3],  cp[:, 3],  label = 'FAL99-P, (ND = 79)')
ax.plot(ti[:, 4],  cp[:, 4],  label = 'FAL99-S, (ND = 79)')
ax.plot(ti[:, 0],  cp[:, 0],  label = 'FAL99-C, (ND = 82)')
ax.plot(ti[:, 10], cp[:, 10], label = 'MURaM, (ND = 86)', color = 'k')
ax.plot(ti[:, 9],  cp[:, 9],  label = 'FCHHT-B, (ND = 91)')

ax.axhline(y = 1e-4, linestyle = '--')
ax.axvline(x = 22.5, linestyle = '-', color = 'k', linewidth = 0.5)

ax.set_yscale('log')

ax.set_xlim(0, 100)
ax.set_ylim(1e-8, 1e+4)

ax.set_xlabel('Time, [min]')
ax.set_ylabel('Convergence')

ax.xaxis.set_major_locator(MultipleLocator(10))
ax.xaxis.set_minor_locator(AutoMinorLocator(10))

fig.suptitle('H, C, O, Mg, Al, Si, S, Ca, Cr, Fe, Ni', y = 0.92)

leg = ax.legend(framealpha = 1, loc = 1, handletextpad = 1, prop = {'size': 7.0})

auxplt.savepdf('conv_times')
