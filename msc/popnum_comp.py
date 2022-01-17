import numpy as np
import matplotlib.pyplot as plt
import os
import sys
import importlib

if not '../aux/' in sys.path: sys.path.append('../aux/')

import paths;  importlib.reload(paths)
import nessy;  importlib.reload(nessy)
import auxplt; importlib.reload(auxplt)
import auxsys; importlib.reload(auxsys)

plt.close('all')

name1 = sys.argv[1]
name2 = sys.argv[2]

odf_it = '0'

if len(sys.argv) > 3: odf_it = sys.argv[3]

if odf_it == '0': prefix = paths.it0h
if odf_it == '1': prefix = paths.it1h
if odf_it == '2': prefix = paths.it2h
if odf_it == '3': prefix = paths.it3h

lev, rne1, popnum1 = nessy.read_popnum(prefix + name1)
lev, rne2, popnum2 = nessy.read_popnum(prefix + name2)

fig, ax = plt.subplots(nrows = 1, ncols = 1, figsize = (12, 10))

fig.tight_layout()

for i in range(0, 114): ax.plot((popnum2[i, :] - popnum1[i, :]) * 100 / popnum2[i, :])
#for i in range(100, 106): ax.plot((popnum2[i, :] - popnum1[i, :]) * 100 / popnum2[i, :])

ax.plot((rne2 - rne1) * 100 / rne2, color = 'k')

plt.show()
