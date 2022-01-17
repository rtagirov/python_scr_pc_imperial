import numpy as np

import os
import sys

if not '/mnt/SSD/sim/python/src/aux/' in sys.path: sys.path.append('/mnt/SSD/sim/python/src/aux/')

import paths

modes = np.arange(13)

for mode in modes:

    command = 'python ' + paths.pydir + 'var/form_temp_weighted_rings.py ' + str(mode)

    os.system(command)
