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

sysaux.clean_dir(paths.figdir + 'comppop', mode = 'noverbose')

name = sys.argv[1]

dir1 = paths.it0h + 'der/'   + name + '/' + paths.lev
dir2 = paths.it0h + 'noder/' + name + '/' + paths.lev

files = os.listdir(dir1)

names_pdf = ''

for f in files:

    print('Plotting ', f)

    names_pdf = names_pdf + f + '.pdf '

    pop1 = np.loadtxt(dir1 + f, usecols = [3], skiprows = 2)

    pop2 = np.loadtxt(dir2 + f, usecols = [3], skiprows = 2)

    fig, ax = plt.subplots(nrows = 1, ncols = 1, figsize = (12.0, 10.0))

    fig.suptitle(f, y = 1.02)

    ax.plot((pop1  - pop2) * 100.0 / pop2)

    ax.set_ylim(-0.5, 0.5)

    pltaux.savepdf('comppop/' + f)

os.chdir(paths.figdir + 'comppop/')

os.system('pdftk ' + names_pdf + 'output ' + 'joined.pdf')

os.chdir(paths.pydir)
