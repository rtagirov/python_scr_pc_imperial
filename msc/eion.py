import numpy as np
import matplotlib.pyplot as plt

from matplotlib.ticker import AutoMinorLocator
from matplotlib.ticker import MultipleLocator

import importlib
import sys
import os

import paths;  importlib.reload(paths)
import pltaux; importlib.reload(pltaux)
import sysaux; importlib.reload(sysaux)

def get_elems(filename):

    f = open(filename, 'r')

    elems = []; anums = []

    for line in f:

        if line[0 : 7] == 'ELEMENT':

            parts = line.split()

            elem = parts[2].strip(' () ')

            if parts[3] == ')': anum = int(float(parts[4]))

            if parts[3] != ')': anum = int(float(parts[3]))

            if elem == 'HE':
 
                elem = 'He'

                anum = 2

            elems.append(elem)

            anums.append(anum)

    f.close()

    return elems, anums

#def get_eion(filename):

#    eion = np.zeros((30, 12))

#    for i in range(1, 31):

#        elem_eion = e[np.where(lan == i)]

#        for j in range(0, len(elem_eion)):

#            eion[i - 1, j] = elem_eion[j]

#    return chn, eion

rctoev = 0.00012

prefix = paths.it0h

elems, anums = get_elems(prefix + 'eion_test/DATOM_FULL')

#chn,   eion =  get_eion(prefix + 'eion_test/eion.out')

ln, lan, chn, eion = np.loadtxt(prefix + 'eion_test/eion.out', unpack = True)

chd = np.zeros(len(ln))

for i in range(1, int(np.max(ln))): chd[i] = chn[i] - abs(chn[i - 1])

#print(chn_diff)

#sys.exit()

#for i in range(0, 30):

#    print(i + 1, eion[i, :])

#sys.exit()

plt.close('all')

fig1, ax1 = plt.subplots(nrows = 1, ncols = 1, figsize = (12.0, 10.0))

fig1.tight_layout()

ax1.set_yscale('log')

ax1.set_ylim(0.1, 100)
ax1.set_xlim(0.5, 30.5)

ax1.xaxis.grid(True)

ax1.xaxis.set_major_locator(MultipleLocator(1))

ax1.scatter(lan,                      rctoev * eion,                      color = 'k', label = r'$Z = 0$, exited states')
ax1.scatter(lan[np.where(chd == -1)], rctoev * eion[np.where(chd == -1)], color = 'g', label = r'$Z = 0$, ground state')
ax1.scatter(lan[np.where(chn == 1)],  rctoev * eion[np.where(chn == 1)],  color = 'r', label = r'$Z = 1$')
ax1.scatter(lan[np.where(chn == -1)], rctoev * eion[np.where(chn == -1)], color = 'm', label = r'$Z = -1$')

ax1.set_xlabel('Atomic Number',           labelpad = 10.5)
ax1.set_ylabel('Ionization Energy, [eV]', labelpad = 10.5)

ax2 = ax1.twiny()

ax2.set_xlim(ax1.get_xlim())

ax2.set_xticks(np.arange(1, 31))

ax2.set_xticklabels(elems)

ax2.set_xlabel('Element', labelpad = 10.5)

leg = ax1.legend(framealpha = 1, loc = 4, handletextpad = 1, prop = {'size': 20})

#for obj in leg.legendHandles: obj.set_linewidth(3.0)

#fig1.show()
pltaux.savepdf('eion')
