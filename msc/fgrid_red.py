import importlib

import paths; importlib.reload(paths)

import sys

name = sys.argv[1]

f = int(sys.argv[2])

w = float(sys.argv[3])

grid_in = open(paths.inp + name, 'r')

grid_out = open(paths.out + name + '_red_' + str(f), 'w')

i = 1

for line in grid_in:

    wvl = float(line)

    if wvl <= w: grid_out.write(line)

    if wvl > w and i % f == 0:

        grid_out.write(line)

        i = i + 1

    if wvl > w and i % f != 0: i = i + 1

grid_in.close()
grid_out.close()
