import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as const

import sys

if not '../aux/' in sys.path: sys.path.append('../aux/')

import importlib
import auxplt
import paths

importlib.reload(auxplt)
importlib.reload(paths)

def P(A, mu):

    P = 0.0

    for i, a in enumerate(A):

        P += a * mu**i

    return P

def wvl_selection(wvl, wvl_neckel):

    idxs = []
    wvls = []

    for i, w in enumerate(wvl_neckel):

        idx = np.abs(wvl - w).argmin()

        idxs.append(idx)
        wvls.append(wvl[idx])

    return idxs, wvls

mu = np.array([1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.05])

w = [303.327, 329.897, 401.970, 445.125, 519.930, 669.400, 1046.600]

A = np.array([[0.08011, 0.70695,  0.49910, -0.31080, -0.02177,  0.04642],  # 303.327
              [0.09188, 0.92459,  0.19604, -0.39546,  0.23599, -0.05303],  # 329.897
              [0.12323, 1.08648, -0.43974,  0.45912, -0.32759,  0.09850],  # 401.970
              [0.15248, 1.38517, -1.49615,  1.99886, -1.48155,  0.44119],  # 445.125
              [0.23695, 1.29927, -1.28034,  1.37760, -0.85054,  0.21706],  # 519.930
              [0.34685, 1.37539, -2.04425,  2.70493, -1.94290,  0.55999],  # 669.400
              [0.49870, 1.21429, -2.06976,  2.80703, -2.05247,  0.60221]]) # 1046.600

clv_a = np.loadtxt(paths.atlruns + 'var_m/Q/spec.out')
clv_n = np.loadtxt(paths.it1f + 'var/Q/fal/CLV_UVI')

wvl_a = clv_a[:, 1]
int_a = clv_a[:, 4 : ]

wvl_n = clv_n[:, 0] / 10.0
int_n = clv_n[:, 1 : ]

idx_n, wvl_nessy = wvl_selection(wvl_n, np.array(w))
idx_a, wvl_atlas = wvl_selection(wvl_a, np.array(w))

plt.close('all')

fig, ax = plt.subplots(nrows = 1, ncols = 1, figsize = (const.golden * 5, 5))

col = ['g', 'r', 'c', 'm', 'b', 'y', 'k']

auxplt.figpar(3, 3, 15)

fig.tight_layout()

ss = 10.1

for i in range(len(w) - 1, -1, -1):

    plt.plot(mu, P(A[i, :], mu), label = str(w[i]) + ' nm', color = col[i])

    clv_nessy = int_n[idx_n[i], :] / int_n[idx_n[i], 0]
    clv_atlas = int_a[idx_a[i], :] / int_a[idx_a[i], 0]

    plt.scatter(mu, clv_nessy, s = ss * 2, marker = 'o', facecolors = 'none', edgecolors = col[i])
    plt.scatter(mu, clv_atlas, s = ss, marker = '*', color = col[i])

plt.xlim(1.0, 0.0)
plt.ylim(0.0, 1.0)

plt.xlabel(r'$\mu$')
plt.ylabel(r'$I(\mu) / I_c$')

leg = plt.legend(framealpha = 1, loc = 3, handletextpad = 1, prop = {'size': 17.0})

for obj in leg.legendHandles: obj.set_linewidth(3.0)

auxplt.savepdf('var/neckel_clv_comp')
