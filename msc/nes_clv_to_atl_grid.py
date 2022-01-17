import numpy as np
import matplotlib.pyplot as plt

from matplotlib.ticker import AutoMinorLocator
from matplotlib.ticker import MultipleLocator
from matplotlib.ticker import LogLocator

from tqdm              import tqdm

import importlib
import math
import sys

if not '../aux/' in sys.path: sys.path.append('../aux/')

import paths;   importlib.reload(paths)
import spec;    importlib.reload(spec)
import phys;    importlib.reload(phys)
import nessy;   importlib.reload(nessy)
import auxsys;  importlib.reload(auxsys)
import auxplt;  importlib.reload(auxplt)
import auxfile; importlib.reload(auxfile)
import auxfunc; importlib.reload(auxfunc)

prefix1 = paths.it1f

wa = np.loadtxt(paths.atlruns + 'var_m/Q/spec.out', usecols = [1])

fa = phys.c / (wa * 1.0e-7)

nwn = auxfile.num_lines(prefix1 + 'var_od/Q/fal/CLV')
#nwn = auxfile.num_lines(prefix1 + 'var/Q/fal/CLV')

Ih = np.zeros((5, 11, nwn))

Il = np.zeros((5, 11, len(wa)))

runs = ['Q/fal', 'F/fal', 'Q/kur', 'P/kur', 'U/kur']

for i in tqdm(range(0, len(runs)), ncols = auxfunc.term_width(), desc = 'Reading'):
    
    wn, Ih[i, 0,  :], \
        Ih[i, 1,  :], \
        Ih[i, 2,  :], \
        Ih[i, 3,  :], \
        Ih[i, 4,  :], \
        Ih[i, 5,  :], \
        Ih[i, 6,  :], \
        Ih[i, 7,  :], \
        Ih[i, 8,  :], \
        Ih[i, 9,  :], \
        Ih[i, 10, :]  \
        = np.loadtxt(prefix1 + 'var_od/' + runs[i] + '/CLV', unpack = True)
#        = np.loadtxt(prefix1 + 'var/' + runs[i] + '/CLV', unpack = True)

    wn = wn / 10.0

for i in range(0, len(runs)):

    for j in range(0, 11): Il[i, j, :] = spec.mean_over_grid(Ih[i, j, :], wn, wa, message = 'Averaging ' + runs[i] + ' angle ' + str(j))

itup = ('%4i',)
ftup = ('%9.2f',)
etup = ('%15.7E',)

frmt = itup + ftup + 12 * etup

nmbr = np.arange(1, len(wa) + 1)

for i in tqdm(range(0, len(runs)), ncols = auxfunc.term_width(), desc = 'Saving'):

    np.savetxt(paths.out + runs[i].replace('/', '_'), np.transpose((nmbr, wa, fa, Il[i, 0,  :],\
                                                                                  Il[i, 1,  :],\
                                                                                  Il[i, 2,  :],\
                                                                                  Il[i, 3,  :],\
                                                                                  Il[i, 4,  :],\
                                                                                  Il[i, 5,  :],\
                                                                                  Il[i, 6,  :],\
                                                                                  Il[i, 7,  :],\
                                                                                  Il[i, 8,  :],\
                                                                                  Il[i, 9,  :],\
                                                                                  Il[i, 10, :])), fmt = frmt, delimiter = '  ')
