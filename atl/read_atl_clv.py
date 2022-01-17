import sys
import numpy as np

wvl = np.loadtxt('wvl')

fins = []

for arg in sys.argv:

    if arg != sys.argv[0]: fins.append(arg)

for fin in fins:

    fi = open(fin, 'r')

    clv = np.zeros((len(wvl), 11))

    i = 1; j = 0

    for line in fi:

        if i % 2 == 1: clv[j, 0 : 8] = np.float128(line.split())

        if i % 2 == 0:

            clv[j, 8 : 11] = np.float128(line.split())

            j = j + 1

        i = i + 1

    fi.close()

    np.savetxt(fin + '_ord', np.transpose((wvl, clv[:, 0],\
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
