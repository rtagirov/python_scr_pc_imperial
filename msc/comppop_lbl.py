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

plt.close('all')

sysaux.clean_dir(paths.figdir + 'comppop', mode = 'noverbose')

name1 = sys.argv[1]
name2 = sys.argv[2]

odf_it = '0'

if len(sys.argv) > 3: odf_it = sys.argv[3]

if odf_it == '0': prefix = paths.it0h
if odf_it == '1': prefix = paths.it1h
if odf_it == '2': prefix = paths.it2h
if odf_it == '3': prefix = paths.it3h

dir1 = prefix + name1 + '/' + paths.lev
dir2 = prefix + name2 + '/' + paths.lev

files = os.listdir(dir1)

names_pdf = ''

for f in files:

    print('Plotting ', f)

    names_pdf = names_pdf + f + '.pdf '

    fig, ax = plt.subplots(nrows = 1, ncols = 1, figsize = (12.0, 10.0))

    fig.tight_layout()

    fig.suptitle(f, y = 1.001)

    pop1 = np.loadtxt(dir1 + f, usecols = [3], skiprows = 2)
    pop2 = np.loadtxt(dir2 + f, usecols = [3], skiprows = 2)

    ax.plot((pop1  - pop2) * 100.0 / pop2)

    pltaux.savepdf('comppop/' + f)

os.chdir(paths.figdir + 'comppop/')

os.system('pdftk ' + names_pdf + 'output ' + 'joined.pdf')

os.chdir(paths.pydir)
