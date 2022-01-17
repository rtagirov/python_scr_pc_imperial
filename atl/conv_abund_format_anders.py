import numpy as np

import matplotlib.pyplot as plt

import importlib
import re
import sys

if not '../aux/' in sys.path: sys.path.append('../aux/')

import paths; importlib.reload(paths)
import phys;  importlib.reload(phys)

from decimal import Decimal

ab_all = np.loadtxt(paths.inp + 'CARDS_DEF_AB_ONLY')

atmean = ab_all[len(ab_all) - 1]

ab_all = np.delete(ab_all, len(ab_all) - 1)

#print(ab_all)

#sys.exit()

ab_he = 1.0 - np.sum(ab_all)

ab_all = np.insert(ab_all, 1, ab_he)

#print(ab_all)

#sys.exit()

fa = open(paths.inp + 'Anders_abundances.txt', 'r')

ab = np.zeros(30)

i = 1

for line in fa:

   if i == 5:

        ab[0] =  float(line.split()[6])
        ab[1] =  float(line.split()[8])

   if i == 6:

        ab[2] =  10.0**float(line.split()[3])
        ab[3] =  10.0**float(line.split()[5])
        ab[4] =  10.0**float(line.split()[7])
        ab[5] =  10.0**float(line.split()[9])
        ab[6] =  10.0**float(line.split()[11])
        ab[7] =  10.0**float(line.split()[13])

   if i == 7:

        ab[8] =  10.0**float(line.split()[3])
        ab[9] =  10.0**float(line.split()[5])
        ab[10] = 10.0**float(line.split()[7])
        ab[11] = 10.0**float(line.split()[9])
        ab[12] = 10.0**float(line.split()[11])
        ab[13] = 10.0**float(line.split()[13])

   if i == 8:

        ab[14] = 10.0**float(line.split()[3])
        ab[15] = 10.0**float(line.split()[5])
        ab[16] = 10.0**float(line.split()[7])
        ab[17] = 10.0**float(line.split()[9])
        ab[18] = 10.0**float(line.split()[11])
        ab[19] = 10.0**float(line.split()[13])

   if i == 9:

        ab[20] = 10.0**float(line.split()[3])
        ab[21] = 10.0**float(line.split()[5])
        ab[22] = 10.0**float(line.split()[7])
        ab[23] = 10.0**float(line.split()[9])
        ab[24] = 10.0**float(line.split()[11])
        ab[25] = 10.0**float(line.split()[13])

   if i == 10:

        ab[26] = 10.0**float(line.split()[3])
        ab[27] = 10.0**float(line.split()[5])
        ab[28] = 10.0**float(line.split()[7])
        ab[29] = 10.0**float(line.split()[9])

   i += 1

fa.close()

atmean = 0.0; i = 0

for value in phys.ptable.values():

    atmean += value['atmass'] * ab[i]

    i += 1

fno = open(paths.inp + 'CARDS_DEF', 'r+')

ablist = []

i = 0

for line in fno:

    i += 1

    if i == 1: continue

    if i > 31: break

    ablist.append(line.split()[0])

ablist.insert(1, 'AHE=')

for i in range(0, 30):

    ablist[i] += '  ' + '{:.2E}'.format(Decimal(str(ab[i])))

    part1 = ablist[i].split('-')[0]
    part2 = ablist[i].split('-')[1]

    if len(part2) == 1:

        part2 = '0' + part2

        ablist[i] = part1 + '-' + part2

ablist[30] += ' ' + str(atmean)[0 : 7]

for i in range(0, 30): print(ablist[i], ' ', ab[i], ' ', ab_all[i])

print(ablist[30])

#plt.close('all')

#plt.scatter(np.arange(len(ab)),     ab,     color = 'k')
#plt.scatter(np.arange(len(ab_all)), ab_all, color = 'r')

#plt.yscale('log')

#plt.show()

sys.exit()

fno.seek(0, 0)

fnn = open(paths.out + 'cards_atl_ab', 'w')

for line in fno:

    flag = 0

    for i in range(0, 31):

        if re.search(ablist[i].split('=')[0] + '=', line):

            flag = 1

            newline = ablist[i]

    if flag == 1: fnn.write(newline + '\n')

    if flag == 0: fnn.write(line)

fno.close()
fnn.close()
