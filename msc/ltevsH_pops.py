import numpy as np
import matplotlib.pyplot as plt

from matplotlib.ticker import AutoMinorLocator
from matplotlib.ticker import MultipleLocator
from matplotlib.ticker import LogLocator

import importlib
import sys
import os

if not '../aux/' in sys.path: sys.path.append('../aux/')

import paths;   importlib.reload(paths)
import spec;    importlib.reload(spec)
import nessy;   importlib.reload(nessy)
import auxsys;  importlib.reload(auxsys)
import auxplt;  importlib.reload(auxplt)
import phys;    importlib.reload(phys)

lev_atl = ['HMINUS',
           'HI1',
           'HI2',
           'HI3',
           'HI4',
           'HI5',
           'HI6',
           'HI7',
           'HI8',
           'HII',
           'CI1',
           'CI2',
           'CI3',
           'CII1',
           'MgI1',
           'MgII1',
           'AlI1',
           'AlII1',
           'SiI1',
           'SiII1',
           'FeI1',
           'FeI2',
           'FeI3',
           'FeI4',
           'FeI5',
           'FeII1']

ele_atl = []

i = 0

for lev in lev_atl:

    ele_atl.append(lev[0 : 2])

    if ele_atl[i][1] == 'M' or ele_atl[i][1] == 'I': ele_atl[i] = lev[0]

    i += 1

entot  = np.loadtxt(paths.it0h + 'ltevsH/H_C_Mg_Al_Si_Fe/ATM_STR', skiprows = 2, usecols = [3])

xnatom = np.loadtxt(paths.atlruns + 'ltevsH/H_C_Mg_Al_Si_Fe/ltepop/XNATOM')

ndp = len(xnatom)

pop_nes = np.zeros((len(lev_atl), ndp))
pop_atl = np.zeros((len(lev_atl), ndp))

#NESSY populations

#-----------------------------------------------------------

lev_nes, rne, popnum = nessy.read_popnum(paths.it0h + 'ltevsH/lte')

lev_atl_idx = []

for i in range(len(lev_atl)):

    for j in range(len(lev_nes)):

        if lev_nes[j] == lev_atl[i]: lev_atl_idx.append(j)

i = 0

for idx in lev_atl_idx:

    pop_nes[i, :] = popnum[idx, :]

    i += 1

#-----------------------------------------------------------

#ATLAS populations

i = 0

for lev in lev_atl:

    pop_atl[i, :] = np.loadtxt(paths.atlruns + 'ltevsH/H_C_Mg_Al_Si_Fe/ltepop/' + lev)

    i += 1

xne = np.loadtxt(paths.atlruns + 'ltevsH/H_C_Mg_Al_Si_Fe/ltepop/ELECTR')

#-----------------------------------------------------------

plt.close('all')

os.system('mkdir -p ' + paths.figdir + 'ltevsH_pops')

ele_str = ''

for i in range(0, len(lev_atl)):

    if i == 0 or (i != 0 and ele_atl[i] != ele_atl[i - 1]):

        fig, ax = plt.subplots(nrows = 1, ncols = 1, figsize = (12.0, 10.0))

        fig.tight_layout()

        fig.suptitle(ele_atl[i], y = 1.02)

        ax.plot(np.arange(ndp), np.zeros(ndp), 'k--')

        ax.set_xlim(ndp - 1, 0)

        if ele_atl[i] != 'Fe': ax.set_ylim(-20, 6)
#        ax.set_ylim(-20, 120)

        ax.set_xlabel('Depth Index', fontsize = 12.5)
        ax.set_ylabel('(NESSY - ATLAS) * 100 / ATLAS, [%]', fontsize = 12.5)

        ax.yaxis.set_minor_locator(AutoMinorLocator(5))

    ax.plot((pop_nes[i, :] - pop_atl[i, :]) * 100 / pop_atl[i, :], label = lev_atl[i])

    leg = ax.legend(framealpha = 1, loc = 'best', handletextpad = 1, prop = {'size': 14.5})

    for obj in leg.legendHandles: obj.set_linewidth(3.0)

    if i == len(lev_atl) - 1 or (i != len(lev_atl) - 1 and ele_atl[i + 1] != ele_atl[i]):

        auxplt.savepdf('ltevsH_pops/' + ele_atl[i])

        ele_str += ele_atl[i] + '.pdf '

#-----------------------------------------------------------
# Plotting electrons and total particle concentration

fig, ax = plt.subplots(nrows = 1, ncols = 1, figsize = (12.0, 10.0))

fig.tight_layout()

fig.suptitle('ELectrons and heavy particle concentration', y = 1.02)

ax.plot(np.arange(ndp), np.zeros(ndp), 'k--')

ax.set_xlim(ndp - 1, 0)

ax.set_xlabel('Depth Index', fontsize = 12.5)
ax.set_ylabel('(NESSY - ATLAS) * 100 / ATLAS, [%]', fontsize = 12.5)

ax.yaxis.set_minor_locator(AutoMinorLocator(5))

ax.plot((rne - xne) * 100 / xne,         label = 'Electrons')
ax.plot((entot - xnatom) * 100 / xnatom, label = 'Heavy Particle Concentration')

leg = ax.legend(framealpha = 1, loc = 'best', handletextpad = 1, prop = {'size': 14.5})

for obj in leg.legendHandles: obj.set_linewidth(3.0)

auxplt.savepdf('ltevsH_pops/elec')

ele_str += 'elec.pdf '

#-----------------------------------------------------------

cwd = os.getcwd()

os.chdir(paths.figdir + 'ltevsH_pops')

os.system('pdftk ' + ele_str + 'output ' + 'ltepops.pdf')

os.chdir(cwd)
