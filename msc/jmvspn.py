import numpy as np

import matplotlib.pyplot as plt

import os

import importlib

import paths;  importlib.reload(paths)
import nessy;  importlib.reload(nessy)

import pltaux; importlib.reload(pltaux)
import sysaux; importlib.reload(sysaux)

plt.close('all')

sysaux.clean_dir(paths.figdir + 'jmvspn', mode = 'noverbose')

dir1 = paths.it0h + 'popcomp/lte/'

dir2 = paths.it0h + 'popcomp/lte_jobmax/'

lev, rne1, popnum1 = nessy.read_popnum(dir1)

fig, ax1 = plt.subplots(nrows = 1, ncols = 1, figsize = (12, 10))

fig.tight_layout

ax1.set_ylim(-0.0001, 0.0001)

j = 0

for l in lev:

    print('Plotting ', l)

    popnum2 = np.loadtxt(dir2 + paths.lev + l, skiprows = 2, usecols = 3)

#    fig, ax1 = plt.subplots(nrows = 1, ncols = 1, figsize = (12, 10))

#    fig.tight_layout

#    fig.suptitle(l, y = 1.001)

    ax1.plot((popnum2 - popnum1[j, :]) * 100 / popnum2)

#    ax1.set_ylim(-0.0001, 0.0001)

#    pltaux.savepdf('jmvspn/' + l)
    j = j + 1

rne2 = np.loadtxt(dir2 + paths.lev + 'ELECTR', skiprows = 2, usecols = 3)

ax1.plot((rne2 - rne1) * 100 / rne2, color = 'k')

pltaux.savepdf('jmvspn')

#os.chdir(paths.figdir + 'jmvspn/')

#os.system('pdftk * output joined.pdf')

#os.chdir(paths.pydir)
