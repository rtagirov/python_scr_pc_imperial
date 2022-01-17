import numpy as np
import matplotlib.pyplot as plt

from matplotlib.ticker import AutoMinorLocator
from matplotlib.ticker import MultipleLocator
from matplotlib.ticker import LogLocator

import importlib
import sys

if not '../aux/' in sys.path: sys.path.append('../aux/')

import paths;  importlib.reload(paths)
import spec;   importlib.reload(spec)
import nessy;  importlib.reload(nessy)
import auxsys; importlib.reload(auxsys)
import auxplt; importlib.reload(auxplt)

dp, to, om_d, di_d, vc_d, svc_d, dvc_d, c1_d, c2_d, c3_d = np.loadtxt(paths.it1f + 'runtime/def_test/linop_cycles.time.default', skiprows = 1, unpack = True)
dp, to, om_r, di_r, vc_r, svc_r, dvc_r, c1_r, c2_r, c3_r = np.loadtxt(paths.it1f + 'runtime/def_test/linop_cycles.time.reduced', skiprows = 1, unpack = True)

pl_d = to - om_d
pl_r = to - om_r

plt.close('all')

fig, ax = plt.subplots(nrows = 3, ncols = 1, figsize = (6.0, 6.0))

fig.tight_layout()

ax[0].plot(dp, pl_d,  color = 'r', linestyle = '-', label = 'lines processed')
ax[0].plot(dp, pl_r,  color = 'r', linestyle = '--')
ax[0].plot(dp, vc_d,  color = 'g', linestyle = '-', label = 'voigt instances')
ax[0].plot(dp, dvc_d, color = 'b', linestyle = '-', label = 'double voigt instances')
ax[0].plot(dp, vc_r,  color = 'g', linestyle = '--')
ax[0].plot(dp, dvc_r, color = 'b', linestyle = '--')

ax[0].set_yscale('log')
ax[0].yaxis.set_minor_locator(LogLocator(base = 10.0, subs = (2, 3, 4, 5, 6, 7, 8, 9)))
#ax[0].yaxis.set_minor_locator(AutoMinorLocator(10))
ax[0].yaxis.set_major_locator(LogLocator(base = 10.0, numticks = 12))

ax[0].axvspan(82, 61, alpha=0.5, color='gray')
#ax[0].axvline(x = 62, color = 'k', linewidth = 0.5)

ax[1].plot(dp, pl_d  / pl_r,  color = 'r', linestyle = '-')
ax[1].plot(dp, vc_d  / vc_r,  color = 'g', linestyle = '-')
ax[1].plot(dp, dvc_d / dvc_r, color = 'b', linestyle = '-')

ax[1].set_yscale('log')

ax[1].axvspan(82, 61, alpha=0.5, color='gray')
#ax[1].axvline(x = 62, color = 'k', linewidth = 0.5)

ax[2].plot(dp, c3_d, color = 'orange', label = 'cycle 3: ' + str(sum(c3_d))[0 : 5] + r' s $\rightarrow$ ' + str(sum(c3_r))[0 : 5] + ' s')
ax[2].plot(dp, c2_d, color = 'black',  label = 'cycle 2: ' + str(sum(c2_d))[0 : 5] + ' s')
ax[2].plot(dp, c3_r, color = 'orange', linestyle = '--')

ax[2].set_yscale('log')

ax[2].axvspan(82, 61, alpha=0.5, color='gray')
#ax[2].axvline(x = 62, color = 'k', linewidth = 0.5)

ax[2].set_xlabel('Depth Index')
ax[0].set_ylabel('Number of instances')
ax[1].set_ylabel('Ratio')
ax[2].set_ylabel('Execution time, s')

ax[0].set_xlim(82, 1)
ax[1].set_xlim(82, 1)
ax[2].set_xlim(82, 1)

ax[0].set_ylim(1e+0, 1e+5)
ax[1].set_ylim(1e+1, 3e+3)

leg0 = ax[0].legend(framealpha = 1, loc = 3, handletextpad = 1, prop = {'size':7.5})
leg2 = ax[2].legend(framealpha = 1, loc = 1, handletextpad = 1, prop = {'size':7.5})

for obj in leg0.legendHandles: obj.set_linewidth(3.0)
for obj in leg2.legendHandles: obj.set_linewidth(3.0)

auxplt.savepdf('linop_cycles_dvc')
