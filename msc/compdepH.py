import numpy as np

import matplotlib.pyplot as plt

from glob import glob

import os
import importlib

import sysaux; importlib.reload(sysaux)
import pltaux; importlib.reload(pltaux)
import paths;  importlib.reload(paths)

plt.close('all')

sysaux.clean_dir(paths.figdir + 'compdepH', mode = 'noverbose')

dir1 = paths.it0h + 'popcomp/H/'
dir2 = paths.it0h + 'popcomp/full/'

height = np.loadtxt(dir1 + 'ATM_MOD', usecols = [0])

files = glob(dir1 + paths.lev + 'HI*'); files.append(dir1 + paths.lev + 'HMINUS')

levs = []

for f in files:

    parts = f.split('/')

    levs.append(parts[len(parts) - 1])

for lev in levs:

    print('Plotting ', lev)

    b1 = np.loadtxt(dir1 + paths.lev + lev, skiprows = 2, usecols = [4])
    b2 = np.loadtxt(dir2 + paths.lev + lev, skiprows = 2, usecols = [4])

    fig, ax = plt.subplots(nrows = 2, ncols = 1, figsize = (12, 10))

    ax[0].set_yscale('log')

    ax[0].plot(height, b1, label = 'H' + '(' + lev + ')')
    ax[0].plot(height, b2, label = 'full' + '(' + lev + ')')

    ax[1].plot(height, b1 / b2)

    ax[1].set_ylabel('H / full')
    ax[1].set_xlabel('Height, [km]')

    ax[0].legend(framealpha = 1, loc = 'best', handletextpad = 1, prop = {'size': 20})

    pltaux.savepdf('compdepH/' + lev)

os.chdir(paths.figdir + 'compdepH/')

os.system('pdftk * output joined.pdf')

os.chdir(paths.pydir)
