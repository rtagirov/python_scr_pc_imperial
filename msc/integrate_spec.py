import numpy as np
import matplotlib.pyplot as plt

import importlib
import sys

if not '../aux/' in sys.path: sys.path.append('../aux/')

import paths;  importlib.reload(paths)
import spec;   importlib.reload(spec)
import nessy;  importlib.reload(nessy)

w, f1 = nessy.read_spec(paths.it1f + 'ltt/fala',                 1305, 40000)
w, f2 = nessy.read_spec(paths.it1f + 'ltt/falb_correct',         1305, 40000)
w, f3 = nessy.read_spec(paths.it1f + 'ltt/falb_correct_Tplus10', 1305, 40000)
w, f4 = nessy.read_spec(paths.it1f + 'ltt/falc',                 1305, 40000)

w /= 10.0

wc, f1c = spec.mean_within_delta(w, f1, 10.0)
wc, f2c = spec.mean_within_delta(w, f2, 10.0)
wc, f3c = spec.mean_within_delta(w, f3, 10.0)
wc, f4c = spec.mean_within_delta(w, f4, 10.0)

np.savez(paths.npz + 'cb', w =  wc,
                           f1 = f1c,\
                           f2 = f2c,\
                           f3 = f3c,\
                           f4 = f4c)

cb = np.load(paths.npz + 'cb.npz')

w =  cb['w']

f1 = cb['f1']
f2 = cb['f2']
f3 = cb['f3']
f4 = cb['f4']

#wc, f1c = spec.mean_within_delta(w, f1, 10.0)
#wc, f2c = spec.mean_within_delta(w, f2, 10.0)
#wc, f3c = spec.mean_within_delta(w, f3, 10.0)

print(np.trapz(f1, w))
print(np.trapz(f2, w))
print(np.trapz(f3, w))
print(np.trapz(f4, w))

plt.clf()

plt.plot(w, f4 - f1, color = 'k', label = 'C - A')
plt.plot(w, f4 - f2, color = 'r', label = 'C - B')

plt.xlim(130, 1000)

plt.xlabel('Wavelength, [nm]')
plt.ylabel(r'Contrast, [W / m$^2$ / nm]')

leg = plt.legend(framealpha = 1, loc = 1, handletextpad = 1, prop = {'size': 25.0})

plt.show()
