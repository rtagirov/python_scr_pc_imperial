import numpy as np

import importlib
import re
import sys

if not '../aux/' in sys.path: sys.path.append('../aux/')

import paths; importlib.reload(paths)
import phys;  importlib.reload(phys)

from decimal import Decimal

#fa = open(paths.inp + 'atlas9.punched.input', 'r')
#fa = open('./inp/odf_mh_m1.input', 'r')
fa = open('./inp/odf_mh_0.input', 'r')

ab = np.zeros(30)

i = 1

for line in fa:

   if i == 5:

        ab[0] =  float(line.split()[6])
        ab[1] =  float(line.split()[8])

   if i == 6:

        ab[2] =  ab[0] * 10.0**float(line.split()[3])
        ab[3] =  ab[0] * 10.0**float(line.split()[5])
        ab[4] =  ab[0] * 10.0**float(line.split()[7])
        ab[5] =  ab[0] * 10.0**float(line.split()[9])
        ab[6] =  ab[0] * 10.0**float(line.split()[11])
        ab[7] =  ab[0] * 10.0**float(line.split()[13])

   if i == 7:

        ab[8] =  ab[0] * 10.0**float(line.split()[3])
        ab[9] =  ab[0] * 10.0**float(line.split()[5])
        ab[10] = ab[0] * 10.0**float(line.split()[7])
        ab[11] = ab[0] * 10.0**float(line.split()[9])
        ab[12] = ab[0] * 10.0**float(line.split()[11])
        ab[13] = ab[0] * 10.0**float(line.split()[13])

   if i == 8:

        ab[14] = ab[0] * 10.0**float(line.split()[3])
        ab[15] = ab[0] * 10.0**float(line.split()[5])
        ab[16] = ab[0] * 10.0**float(line.split()[7])
        ab[17] = ab[0] * 10.0**float(line.split()[9])
        ab[18] = ab[0] * 10.0**float(line.split()[11])
        ab[19] = ab[0] * 10.0**float(line.split()[13])

   if i == 9:

        ab[20] = ab[0] * 10.0**float(line.split()[3])
        ab[21] = ab[0] * 10.0**float(line.split()[5])
        ab[22] = ab[0] * 10.0**float(line.split()[7])
        ab[23] = ab[0] * 10.0**float(line.split()[9])
        ab[24] = ab[0] * 10.0**float(line.split()[11])
        ab[25] = ab[0] * 10.0**float(line.split()[13])

   if i == 10:

        ab[26] = ab[0] * 10.0**float(line.split()[3])
        ab[27] = ab[0] * 10.0**float(line.split()[5])
        ab[28] = ab[0] * 10.0**float(line.split()[7])
        ab[29] = ab[0] * 10.0**float(line.split()[9])

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

for i in range(0, 30): print(ablist[i], ' ', ab[i])

print(ablist[30])

#sys.exit()

fno.seek(0, 0)

#fnn = open(paths.out + 'cards_atl_ab', 'w')
fnn = open('abund.inp', 'w')

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
