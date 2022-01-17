import numpy as np
import matplotlib.pyplot as plt

from matplotlib.ticker import AutoMinorLocator
from matplotlib.ticker import MultipleLocator

import importlib
import sys
import os

import paths;     importlib.reload(paths)
import pltaux;    importlib.reload(pltaux)
import sysaux;    importlib.reload(sysaux)
import oper_file; importlib.reload(oper_file)

plt.close('all')

sysaux.clean_dir(paths.figdir + 'comppop', mode = 'noverbose')

name = sys.argv[1]

l, dm, eftd = np.loadtxt(paths.it0h + 'noder/' + name + '/eftd.out', unpack = True)

ndp = oper_file.num_lines(paths.it0h + 'noder/' + name + '/ATM_MOD')

nl  = oper_file.num_lines(paths.it0h + 'noder/' + name + '/eftd.out')

last_eftd = np.loadtxt(paths.it0h + 'noder/' + name + '/eftd.out', usecols = [2], skiprows = nl - ndp)

last_eftd = abs(last_eftd)

last_eftd[np.where(last_eftd == 0)] = np.min(last_eftd[np.nonzero(last_eftd)])

idx = np.arange(1, max(l) + 1)

dir1 = paths.it0h + 'der/'   + name + '/' + paths.lev
dir2 = paths.it0h + 'noder/' + name + '/' + paths.lev

files = os.listdir(dir1)

elec = np.loadtxt(dir1 + 'ELECTR', usecols = [3], skiprows = 2)

names_pdf = ''

for f in files:

#    if f[0 : 2] != 'HE': continue

    pop1 = np.loadtxt(dir1 + f, usecols = [3], skiprows = 2)
    pop2 = np.loadtxt(dir2 + f, usecols = [3], skiprows = 2)

    rel_dev = (pop1  - pop2) * 100.0 / pop2

    if max(abs(rel_dev)) <= 0.0001: continue

    print('Plotting ', f, max(abs(rel_dev)))

    names_pdf = names_pdf + f + '.pdf '

    fig, ax1 = plt.subplots(nrows = 1, ncols = 1, figsize = (12.0, 10.0))

    fig.tight_layout()

    fig.suptitle(f, y = 1.001)

    ax1.plot(idx, rel_dev)

    ax1.set_xlim(91, 1)

#    ax1.set_ylim(-2, 2)

    ax1.set_ylabel('Deviation, [%]')

    ax1.set_xlabel('Depth index')

    ax2 = ax1.twinx()

    ax2.scatter(l, abs(eftd), color = 'k')

    ax2.plot(idx, last_eftd, color = 'r')

    ax2.plot(idx, pop1, color = 'b')
    ax2.plot(idx, elec, color = 'g')

    ax2.set_ylim(1e-30, 1)

    ax2.set_yscale('log')

    ax2.set_ylabel('Derivative/Helium Level Population/Electron population')

    pltaux.savepdf('comppop/' + f)

os.chdir(paths.figdir + 'comppop/')

os.system('pdftk ' + names_pdf + 'output ' + 'joined.pdf')

#os.system('pdftk ' + 'HE* ' + 'output ' + 'He.pdf')

os.chdir(paths.pydir)
