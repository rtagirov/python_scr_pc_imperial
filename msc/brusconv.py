import numpy as np

import importlib

import nessy_spec; importlib.reload(nessy_spec)

wvl_o, flu_o = nessy_spec.read('/home/rtagirov/Dropbox/Work/spec/old', wvl1 = 905, wvl2 = 40000, mode = 'file')
wvl_n, flu_n = nessy_spec.read('/home/rtagirov/Dropbox/Work/spec/new', wvl1 = 905, wvl2 = 40000, mode = 'file')

np.savetxt('/home/rtagirov/Dropbox/Work/spec/old/SI', np.transpose((wvl_o, flu_o)), fmt = ('%10.4f', '%10.4E'), delimiter = '  ')
np.savetxt('/home/rtagirov/Dropbox/Work/spec/new/SI', np.transpose((wvl_n, flu_n)), fmt = ('%10.4f', '%10.4E'), delimiter = '  ')
