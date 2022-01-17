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

#dir_new = paths.it0h + 'datom_new_master/' + paths.lev
#dir_new = paths.it0h + 'datom_full_wie_master/' + paths.lev
#dir_new = paths.it0h + 'datom_def_ge_master/' + paths.lev
#dir_old = paths.it0h + 'datom_old_master/' + paths.lev

sysaux.clean_dir(paths.figdir + 'comppop', mode = 'noverbose')

inpstr = sys.argv[1]

names = inpstr.split()

names_pdf = ''

for name in names:

    print('Plotting ', name)

    names_pdf = names_pdf + name + '.pdf '

#    dir1 = paths.it0h + 'nltemark_rd/'    + name + '/' + paths.lev
#    dir2 = paths.it0h + 'nltemark_noder/' + name + '/' + paths.lev
    dir1 = paths.it0h + 'der/'   + name + '/' + paths.lev
    dir2 = paths.it0h + 'noder/' + name + '/' + paths.lev

    fig, ax = plt.subplots(nrows = 1, ncols = 1, figsize = (12.0, 10.0))

    fig.suptitle(name)

    files = os.listdir(dir1)

    for f in files:

        pop1 = np.loadtxt(dir1 + f, usecols = [3], skiprows = 2)

        pop2 = np.loadtxt(dir2 + f, usecols = [3], skiprows = 2)

        ax.plot((pop1 - pop2) * 100.0 / pop2)

    pltaux.savepdf('comppop/' + name)

os.chdir(paths.figdir + 'comppop/')

os.system('pdftk ' + names_pdf + 'output ' + 'joined.pdf')

os.chdir(paths.pydir)
