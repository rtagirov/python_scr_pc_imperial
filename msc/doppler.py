import numpy as np
import sys
import importlib

import matplotlib.pyplot as plt

if not '../aux/' in sys.path: sys.path.append('../aux/')

import auxplt; importlib.reload(auxplt)

x = np.arange(-100, 101, 1)

y = np.zeros((3, len(x)))
z = np.zeros((3, len(x)))

a = [5, 10, 20]
s = [4, 9, 19]

f1 = 1.0
f2 = 5.0

for i in range(len(a)):

#    y[i, :] = np.exp(-(x / a[i])**2.0) / a[i]
    y[i, :] = f1 * np.exp(-((x - s[i]) / a[i])**2.0) / a[i] + f2 * np.exp(-((x + s[i]) / a[i])**2.0) / a[i]

    z[i, :] = np.sort(y[i, :])

fig, ax = plt.subplots(nrows = 2, ncols = 1, figsize = (6.0, 7.0))

for i in range(len(a)):

    ax[0].plot(x, y[i, :], label = 'a = ' + str(a[i]))
    ax[1].plot(x, z[i, :])

    print(np.mean(z[i, :]))

leg = ax[0].legend(framealpha = 1, loc = 2, handletextpad = 1, prop = {'size': 8.5})

auxplt.savepdf('doppler')
