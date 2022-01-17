import numpy as np

import sys

qnums = sys.argv[1].split()

energ = sys.argv[2].split()

if (len(qnums) != len(energ)):

    print('Argument lengths do not match. Abort.')

    sys.exit()

g = np.array([])

for j in qnums: g = np.append(g, [2 * float(j) + 1])

E = 0.0

for i in range(len(qnums)): E = E + g[i] * float(energ[i])

E = E / sum(g)

print(sum(g), '      ', E)
