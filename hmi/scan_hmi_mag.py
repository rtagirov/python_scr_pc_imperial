from astropy.io import fits

import numpy               as np
import matplotlib.pyplot   as plt
import matplotlib.gridspec as gridspec

import sys
import importlib

import paths;  importlib.reload(paths)
import pltaux; importlib.reload(pltaux)

date = '2013_11_16'

filename_int = 'mag/' + date + '.fits'

hdu = fits.open(paths.inp + filename_int)[0]

x0 = hdu.header['X0']
y0 = hdu.header['Y0']

rsun = hdu.header['R_SUN']

r = np.sqrt(1 - 0.3**2.0) * rsun

img = hdu.data

sigma = 10.0

x = np.arange(4096)

idx_r = np.where(abs(x - x0) <= r)

uc = y0 + np.sqrt(r**2 - (x[idx_r] - x0)**2)
lc = y0 - np.sqrt(r**2 - (x[idx_r] - x0)**2)

ys = 2164

x1 = x0 - np.sqrt(r**2 - (ys - y0)**2)
x2 = x0 + np.sqrt(r**2 - (ys - y0)**2)

idx_scn = np.where((x >= x1) & (x <= x2))

x_scn = x[idx_scn]

mf_scn = img[ys, idx_scn[0]]

pltaux.figpar(fontsize = 30)

fig1 = plt.figure(figsize = (26, 8.26))

gs = gridspec.GridSpec(nrows=1, ncols=2, left=0.06, bottom=0.00, right=0.92, top=0.99,
                       wspace=0.25, hspace=0.0, width_ratios=[1, 1])

img_ax = plt.subplot(gs[0])
scn_ax = plt.subplot(gs[1])

image = img_ax.imshow(img, cmap = 'gray')

img_ax.plot(x[idx_r], uc, color = 'k')
img_ax.plot(x[idx_r], lc, color = 'k')

img_ax.plot(np.array([x1, x2]), np.array([ys, ys]), color = 'k')

img_ax.set_xlim(0, 4096)
img_ax.set_ylim(0, 4096)

img_ax.set_xlabel('X position, [px]', labelpad = 10)
img_ax.set_ylabel('Y position, [px]', labelpad = 10)

img_ax.set_title('Magnetogram')

img_cbar = fig1.colorbar(image, ax = img_ax, fraction = 0.046, pad = 0.04)

img_cbar.set_label('LOS magnetic field, [Gauss]', labelpad = 10)

scn_ax.plot(x_scn, mf_scn, label = 'LOS magnetic field along the scan')

scn_ax.fill_between(x_scn, - 3 * sigma,  3 * sigma, color = 'gray', label = r'$-30G < B_\mu < 30G$ (noise)', alpha = 0.4)

scn_ax.set_xlim(min(x_scn), max(x_scn))

scn_ax.set_xlabel('X position, [px]', labelpad = 10)
scn_ax.set_ylabel('LOS magnetic field, [Gauss]', labelpad = 10)

scn_ax.yaxis.tick_right()

scn_ax.yaxis.set_label_position('right')

scn_ax.set_title('Scan at Y = 2164')

leg = scn_ax.legend(framealpha = 0, loc = 4, handlelength = 2, handletextpad=0.5, prop={'size': 30})

for obj in leg.legendHandles: obj.set_linewidth(5.0)

fig1.savefig('/home/rtagirov/Dropbox/Work/Y1/fig/mag_scan.png', bbox_inches = 'tight')
