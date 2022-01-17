import numpy as np

import importlib
import sys

if not '../aux/' in sys.path: sys.path.append('../aux/')

import paths;  importlib.reload(paths)
import spec;   importlib.reload(spec)
import nessy;  importlib.reload(nessy)

wn, fn = nessy.read_spec(paths.it3f + 'NEW_SPEC_ALI_COL_HMI_STI_GAU', wvl1 = 900, wvl2 = 24000)

wn = wn / 10.0

wa = np.loadtxt(paths.inp + 'composite.atl3', usecols = [0])

fnc = spec.mean_within_delta_over_grid(fn, wn, 1.0, wa) / 1e+6

np.savetxt('nessy_spec_2017.dat', np.transpose((wa, fnc)), fmt = ('%10.4f', '%10.4E'), delimiter = '  ')
