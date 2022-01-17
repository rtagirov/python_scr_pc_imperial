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

#sysaux.clean_dir(paths.figdir + 'comppop', mode = 'noverbose')

run1 = sys.argv[1]
run2 = sys.argv[2]

dir1 = paths.it0h + run1 + '/' + paths.lev
dir2 = paths.it0h + run2 + '/' + paths.lev

fig, ax = plt.subplots(nrows = 1, ncols = 1, figsize = (12.0, 10.0))

files = os.listdir(dir1)

for f in files:

    pop1 = np.loadtxt(dir1 + f, usecols = [2], skiprows = 2)

    pop2 = np.loadtxt(dir2 + f, usecols = [2], skiprows = 2)

    ax.plot((pop1 - pop2) * 100.0 / pop2)

pltaux.savepdf('comppop')

#os.chdir(paths.figdir + 'comppop/')

#os.system('pdftk ' + names_pdf + 'output ' + 'joined.pdf')

#os.chdir(paths.pydir)
