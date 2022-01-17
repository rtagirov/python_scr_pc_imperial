import numpy as np

import importlib

import sys

if not '../aux/' in sys.path: sys.path.append('../aux/')

import paths;  importlib.reload(paths)
import nessy;  importlib.reload(nessy)
import auxsys; importlib.reload(auxsys)

if len(sys.argv) < 3: auxsys.abort('Not enough arguments. Abort.')

run = sys.argv[1]

w1 = int(sys.argv[2].split()[0])
w2 = int(sys.argv[2].split()[1])

mode = 'dir'; it = '0'

if len(sys.argv) == 4: it =   sys.argv[3]
if len(sys.argv) == 5: mode = sys.argv[4]

if it == '0': prefix = paths.it0f
if it == '1': prefix = paths.it1f
if it == '2': prefix = paths.it2f
if it == '3': prefix = paths.it3f

wvl, flu = nessy.read_spec(prefix + run, w1, w2, mode)

np.savetxt(prefix + run + '/flux', np.transpose((wvl, flu)), fmt = ('%10.4f', '%11.5E'), delimiter = '  ')
