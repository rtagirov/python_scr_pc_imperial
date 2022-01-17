import numpy as np

import matplotlib.pyplot as plt

from glob import glob

import os
import importlib

import sysaux; importlib.reload(sysaux)
import pltaux; importlib.reload(pltaux)
import paths;  importlib.reload(paths)

plt.close('all')

#sysaux.clean_dir(paths.figdir + 'compltepopH', mode = 'noverbose')

dir1 = paths.it0h + 'popcomp/H/'
dir2 = paths.it0h + 'popcomp/lte_jobmax/'

height = np.loadtxt(dir1 + 'ATM_MOD', usecols = [0])

files = glob(dir1 + paths.lev + '*')#; files.append(dir1 + paths.lev + 'HMINUS')

levs = []

for f in files:

    parts = f.split('/')

    levs.append(parts[len(parts) - 1])

fig, ax = plt.subplots(nrows = 1, ncols = 1, figsize = (12, 10))

for lev in levs:

    print('Plotting ', lev)

    pop1 = np.loadtxt(dir1 + paths.lev + lev, skiprows = 2, usecols = [2])
    pop2 = np.loadtxt(dir2 + paths.lev + lev, skiprows = 2, usecols = [2])

    ax.plot(height, pop1 / pop2)

rne1 = np.loadtxt(dir1 + paths.lev + 'ELECTR', skiprows = 2, usecols = [2])
rne2 = np.loadtxt(dir2 + paths.lev + 'ELECTR', skiprows = 2, usecols = [2])

ax.plot(height, rne1 / rne2, color = 'k')

ax.set_ylabel('pop1 / pop2')
ax.set_xlabel('Height, [km]')

pltaux.savepdf('compltepopH')
