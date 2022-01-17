import numpy as np
import matplotlib.pyplot as plt

import sys
import importlib

if not '../aux/' in sys.path: sys.path.append('../aux/')

import auxplt; importlib.reload(auxplt)
import paths;  importlib.reload(paths)

wvl = np.loadtxt(paths.atlruns + 'var_m/Q/spec.out', usecols = [1])

delta = np.zeros(len(wvl) - 1)

for i in range(len(wvl) - 1): delta[i] = wvl[i + 1] - wvl[i]

plt.close('all')

fig, ax = plt.subplots(nrows = 1, ncols = 1, figsize = (18.0, 9.0))

auxplt.figpar(3, 3, 23)

ax.scatter(wvl[0 : len(wvl) - 1], delta, marker = 'o', s = 10, color = 'k')

ax.axvline(x = 30.0, linewidth = 0.5, linestyle = '--', color = 'r')
ax.axvline(x = 290., linewidth = 0.5, linestyle = '--', color = 'r')
ax.axvline(x = 1000, linewidth = 0.5, linestyle = '--', color = 'r')
ax.axvline(x = 1600, linewidth = 0.5, linestyle = '--', color = 'r')
ax.axvline(x = 3200, linewidth = 0.5, linestyle = '--', color = 'r')
ax.axvline(x = 6400, linewidth = 0.5, linestyle = '--', color = 'r')

ax.axhline(y = 1,  linewidth = 0.5, linestyle = '--', color = 'r')
ax.axhline(y = 2,  linewidth = 0.5, linestyle = '--', color = 'r')
ax.axhline(y = 5,  linewidth = 0.5, linestyle = '--', color = 'r')
ax.axhline(y = 10, linewidth = 0.5, linestyle = '--', color = 'r')
ax.axhline(y = 20, linewidth = 0.5, linestyle = '--', color = 'r')
ax.axhline(y = 40, linewidth = 0.5, linestyle = '--', color = 'r')

s = 23

ax.text(27.5, 110, '30',   color = 'r', size = s)
ax.text(250., 110, '290',  color = 'r', size = s)
ax.text(825., 110, '1000', color = 'r', size = s)
ax.text(1300, 110, '1600', color = 'r', size = s)
ax.text(2650, 110, '3200', color = 'r', size = s)
ax.text(5400, 110, '6400', color = 'r', size = s)

ax.text(1.05e+4, 0.92, '1',  color = 'r', size = s)
ax.text(1.05e+4, 1.85, '2',  color = 'r', size = s)
ax.text(1.05e+4, 4.60, '5',  color = 'r', size = s)
ax.text(1.05e+4, 9.20, '10', color = 'r', size = s)
ax.text(1.05e+4, 18.5, '20', color = 'r', size = s)
ax.text(1.05e+4, 37.0, '40', color = 'r', size = s)

#ax.set_xlabel(r'$\lambda$, [nm]')
#ax.set_ylabel(r'$\Delta\lambda$, [nm]')
ax.set_xlabel('Wavelength, [nm]')
ax.set_ylabel('Resolution, [nm]')

#plt.suptitle('ATLAS intensity spectrum resolution')

ax.set_xscale('log')
ax.set_yscale('log')

ax.set_xlim(1e+1, 1e+4)
ax.set_ylim(1e-1, 1e+2)

auxplt.savepdf('atlas_res')
