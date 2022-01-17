import numpy as np

import importlib

import paths;  importlib.reload(paths)

wvl = np.loadtxt(paths.inp + 'atl_wvl.dat')

infile = open(paths.inp + 'atl_its_86.dat', 'r')

clv = np.zeros((len(wvl), 11))

i = 1; j = 0

for line in infile:

    if i % 2 == 1: clv[j, 0 : 8] = np.float128(line.split())

    if i % 2 == 0:

       clv[j, 8 : 11] = np.float128(line.split())

       j = j + 1

    i = i + 1

infile.close()

np.savetxt(paths.out + 'atl_clv_86.dat', np.transpose((wvl, clv[:, 0],\
                                                            clv[:, 1],\
                                                            clv[:, 2],\
                                                            clv[:, 3],\
                                                            clv[:, 4],\
                                                            clv[:, 5],\
                                                            clv[:, 6],\
                                                            clv[:, 7],\
                                                            clv[:, 8],\
                                                            clv[:, 9],\
                                                            clv[:, 10])), fmt = ('%9.3E'), delimiter = ' ')
