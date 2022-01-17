import numpy as np

import random

from tqdm import tqdm

st = 'G2V'

mf = 'hydro'
#mf = '300G'

ss = '123000'
#ss = '118000'

atms = '/mnt/HDD/atms/'   + st + '/' + mf + '/' + ss
path = '/mnt/HDD/slices/' + st + '/' + mf + '/' + ss

prefix = st + '_' + mf + '_' + ss + '_rot_'

r1 = random.sample(range(512), 100)
r2 = random.sample(range(512), 100)

nums0 = []
nums1 = []
nums2 = []

n0 = 0

for n1 in tqdm(r1):

    f = open(path + '/' + prefix + str(n1), 'r')

    zs, Ts, ps, ds = np.loadtxt(f, skiprows = 2, unpack = True)

    zs = zs.reshape((512, 348))
    Ts = Ts.reshape((512, 348))
    ps = ps.reshape((512, 348))
    ds = ds.reshape((512, 348))

    for n2 in r2:

        n0 += 1

        z = np.flip(zs[n2, :]) + 240.0

        T = Ts[n2, :]
        p = ps[n2, :]
        d = ds[n2, :]

        np.savetxt(atms + '/atm.' + str(n0), \
        np.column_stack([z, T, p, d]), \
        fmt = ('%6.1f', '%7.1f', '%7.5e', '%7.5e'), delimiter = '  ')

        nums0.append(n0)
        nums1.append(n1)
        nums2.append(n2)

np.savetxt(atms + '/nums.log', \
           np.transpose([nums0, nums1, nums2]), \
           fmt = ('%5i', '%3i', '%3i'), delimiter = '  ')
