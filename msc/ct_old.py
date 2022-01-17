import numpy as np
import matplotlib.pyplot as plt

from matplotlib.ticker import AutoMinorLocator
from matplotlib.ticker import MultipleLocator

import importlib
import sys
import os

if not '../aux/' in sys.path: sys.path.append('../aux/')

import paths;  importlib.reload(paths)
import auxplt; importlib.reload(auxplt)
import auxsys; importlib.reload(auxsys)

#time_dir = 'tst_time/'
time_dir = 'prop/'

runs = sys.argv[1]

parts = 'wrcont como etl steal linpop cycle'

col = ['C0', 'C1', 'C2', 'C3', 'C4', 'C5', 'C6']

auxsys.clean_dir(paths.figdir + time_dir, mode = 'noverbose')

pdf_names = ''

for part in parts.split():

    fig, ax = plt.subplots(nrows = 1, ncols = 1, figsize = (12.0, 10.0))

    fig.tight_layout()

    fig.suptitle(part, y = 1.02)

    i = 0

    pdf_names = pdf_names + part + '.pdf' + ' '

    for run in runs.split():

        print(run, ' ', part)

        path = paths.it0h + time_dir + run + '/'

        ctime = np.loadtxt(path + 'cpu_time.'  + part)

#        wts, wtf = np.loadtxt(path + 'wall_time.' + part, unpack = True); wtime = wtf - wts

#        ax.scatter(np.arange(1, len(wtime) + 1), wtime, color = col[i], label = 'NLTE ' + run)
        ax.scatter(np.arange(1, len(ctime) + 1), ctime, color = col[i], label = 'NLTE ' + run)

#        ax.plot(np.ones(len(wtime)) * np.mean(wtime), color = col[i], label = 'Mean')
        ax.plot(np.ones(len(ctime)) * np.mean(ctime), color = col[i], label = 'Mean')

        ax.set_xlabel('Call number')
#        ax.set_ylabel('Execution time, [s]')
        ax.set_ylabel('CPU time, [s]')

        i += 1

    leg = ax.legend(framealpha = 1, loc = 1, handletextpad = 1, prop = {'size': 10.5})

    for obj in leg.legendHandles:

        obj.set_linewidth(3.0)

    auxplt.savepdf(time_dir + part)

d = os.getcwd()

os.chdir(paths.figdir + time_dir)

os.system('pdftk ' + pdf_names + 'output ' + 'joined.pdf')

os.chdir(d)
