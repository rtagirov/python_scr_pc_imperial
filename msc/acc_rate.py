import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

from matplotlib.ticker import AutoMinorLocator
from matplotlib.ticker import MultipleLocator

import numpy.polynomial.polynomial as poly

import importlib
import sys
import os

from tqdm import tqdm

import paths;     importlib.reload(paths)
import spec;      importlib.reload(spec)
import nessy;     importlib.reload(nessy)
import sysaux;    importlib.reload(sysaux)
import pltaux;    importlib.reload(pltaux)
import auxfunc;   importlib.reload(auxfunc)
import oper_file; importlib.reload(oper_file)

sub = sys.argv[1]

time = []; nlev = []

os.chdir(paths.it1h + 'spec')

for j in range(1, 31):

    if j == 2: continue

    name = 'nlte' + str(j)

    s, f = np.loadtxt(name + '/wall_time.' + sub, unpack = True); t = f - s

    time.append(np.mean(t))

    os.system('grep -irw LEVEL ' + name + '/DATOM_NLTE > lev.out')

    nlev.append(oper_file.num_lines('lev.out'))

    os.system('rm lev.out')

os.chdir(paths.pydir)

c = poly.polyfit(nlev, time, 1)

l = np.arange(1, 104)

pltaux.figpar(fontsize = 20)

fig, ax = plt.subplots(nrows = 1, ncols = 1, figsize = (12.0, 10.0))

fig.tight_layout()

fig.suptitle(sub, y = 1.01)

ax.scatter(nlev, time, color = 'k')

if sub == 'wrcont' or sub == 'como':

    ax.plot(np.arange(1, max(nlev) + 1), np.mean(time) * np.ones(max(nlev)), color = 'g', label = 'Mean', linewidth = 2.0)

else:

    ax.plot(l, poly.polyval(l, c), color = 'r', label = 'Fit')

ax.xaxis.set_minor_locator(ticker.AutoMinorLocator(10))
ax.xaxis.set_major_locator(ticker.MultipleLocator(10))

ax.yaxis.set_minor_locator(ticker.AutoMinorLocator(5))
#ax.yaxis.set_major_locator(ticker.MultipleLocator(1))

ax.set_xlim(5, 105)
#ax.set_ylim(0, 15)

ax.set_xlabel('Number of levels', labelpad = 15)
ax.set_ylabel('Average execution time, [s]')

leg = ax.legend(framealpha = 1, loc = 'best', handletextpad = 1, prop = {'size': 10.5})

for obj in leg.legendHandles: obj.set_linewidth(3.0)

pltaux.savepdf('acc_rate_' + sub)
