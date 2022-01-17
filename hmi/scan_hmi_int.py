from astropy.io    import fits
from scipy.ndimage import filters

import numpy               as np
import matplotlib.pyplot   as plt
import matplotlib.gridspec as gridspec

import sys
import importlib

import paths; importlib.reload(paths)
import pltaux;  importlib.reload(pltaux)

date = '2013_11_16'
#date = '2016_06_10'

filename_int = 'int/' + date + '.fits'

hdulist = fits.open(paths.inp + filename_int)

hdu = hdulist[0]

x0 =   hdu.header['X0']
y0 =   hdu.header['Y0']
rsun = hdu.header['R_SUN']

clv = np.load(paths.npz + 'clv.npz'); img_ld = clv['img_ld']

r = np.sqrt(1 - 0.3**2.0) * rsun

#img = hdu.data
img = hdu.data / img_ld

#img = filters.gaussian_filter(hdu.data, 10)

x = np.arange(4096)

idx_r = np.where(abs(x - x0) <= r)

uc = y0 + np.sqrt(r**2 - (x[idx_r] - x0)**2)
lc = y0 - np.sqrt(r**2 - (x[idx_r] - x0)**2)

ys = 2164

x1 = x0 - np.sqrt(r**2 - (ys - y0)**2)
x2 = x0 + np.sqrt(r**2 - (ys - y0)**2)

idx_scn = np.where((x >= x1) & (x <= x2))

x_scn = x[idx_scn]

int_scn = img[ys, idx_scn[0]]

pltaux.figpar(fontsize = 30)

fig1 = plt.figure(figsize = (26, 8.26))

gs = gridspec.GridSpec(nrows=1, ncols=2, left=0.06, bottom=0.00, right=0.92, top=0.99,
                       wspace=0.25, hspace=0.0, width_ratios=[1, 1])

img_ax = plt.subplot(gs[0])
scn_ax = plt.subplot(gs[1])

image = img_ax.imshow(img, cmap = 'afmhot')

img_ax.plot(x[idx_r], uc, color = 'k')
img_ax.plot(x[idx_r], lc, color = 'k')

img_ax.plot(np.array([x1, x2]), np.array([ys, ys]), color = 'k')

img_ax.set_xlim(0, 4096)
img_ax.set_ylim(0, 4096)

img_ax.set_xlabel('X position, [px]', labelpad = 10)
img_ax.set_ylabel('Y position, [px]', labelpad = 10)

#img_ax.set_title('Intensity image')
img_ax.set_title('Normalized intensity image')

img_cbar = fig1.colorbar(image, ax = img_ax, fraction = 0.046, pad = 0.04)

#img_cbar.set_label('Intensity, [arbitrary units]', labelpad = 10)
img_cbar.set_label('Normalized intensity', labelpad = 10)

scn_ax.plot(x_scn, int_scn, label = 'Normalized intensity along the scan')

sigma = np.sqrt(np.mean((1.0 - img[np.where(np.logical_not(np.isnan(img)))])**2.0))

scn_ax.plot(x_scn, np.ones(len(x_scn)), color = 'r', linewidth = 2.0, label = 'Mean of the normalized image')

scn_ax.fill_between(x_scn, 1.0 - 3.0 * sigma,  1.0 + 3.0 * sigma, color = 'grey', label = r'$1 - 3\sigma_{I_n} < I_n < 1 + 3\sigma_{I_n}$ (noise)', alpha = 0.4)

#scn_ax.fill_between(x_scn, 1.0 - 10.0 * sigma, 1.0 - 3.0 * sigma, color = 'grey', label = r'$1 - 10\sigma_{I_n} < I_n < 1 - 3\sigma_{I_n}$', alpha = 1.0)

scn_ax.set_xlim(min(x_scn), max(x_scn))
scn_ax.set_ylim(0.0, 1.2)

scn_ax.set_xlabel('X position, [px]', labelpad = 10)
#scn_ax.set_ylabel('Intensity, [arbitrary units]', labelpad = 10)
scn_ax.set_ylabel('Normalized intensity', labelpad = 10)

scn_ax.yaxis.tick_right()

scn_ax.yaxis.set_label_position('right')

scn_ax.set_title('Scan at Y = 2164')

leg = scn_ax.legend(framealpha = 0, loc = 4, handlelength = 2, handletextpad=0.5, prop={'size': 30})

for obj in leg.legendHandles: obj.set_linewidth(5.0)

#fig1.savefig('/home/rtagirov/Dropbox/Work/Y1/fig/int_scan_orig.png', bbox_inches = 'tight')
#fig1.savefig(paths.figdir + 'int_scan_orig.png', bbox_inches = 'tight')
fig1.savefig(paths.figdir + 'int_scan_norm.png', bbox_inches = 'tight')

sys.exit()

f1 = plt.figure(1)

plt.clf()

plt.imshow(img, cmap = 'gray')

plt.xlim(0, 4096)
plt.ylim(0, 4096)

plt.show()

f2 = plt.figure(2)

plt.clf()

plt.imshow(img, cmap = 'gray')

plt.xlim(0, 4096)
plt.ylim(0, 4096)

plt.show()

sys.exit()

x = np.array(range(4096))

plt.xlim(0, 4096)
plt.ylim(0, 4096)

idx_rad = np.where(abs(x - x0) <= r)

uc = y0 + np.sqrt(r**2 - (x[idx_rad] - x0)**2)
lc = y0 - np.sqrt(r**2 - (x[idx_rad] - x0)**2)

plt.plot(x[idx_rad], uc, color = 'k')
plt.plot(x[idx_rad], lc, color = 'k')

ys = 2525

x1 = x0 - np.sqrt(r**2 - (ys - y0)**2)
x2 = x0 + np.sqrt(r**2 - (ys - y0)**2)

plt.plot(np.array([x1, x2]), np.array([ys, ys]), color = 'k')

f3 = plt.figure(2)

plt.clf()

idx_scan = np.where((x >= x1) & (x <= x2))

xs = x[idx_scan]

ints = img[ys, idx_scan[0]]

c = np.polynomial.polynomial.polyfit(xs, ints, 2)

p = c[0] + c[1] * xs + c[2] * xs**2

plt.plot(xs, ints)
plt.plot(xs, p, color = 'r', linewidth = 3.0)

plt.xlim(x1, x2)

plt.show()

f4 = plt.figure(3)

plt.clf()

ints_n = ints / p

mean = sum(ints_n) / len(ints_n)

sigma = np.sqrt(sum((ints_n - mean)**2) / len(ints_n))

plt.fill_between(xs, mean - 3 * sigma, mean + 3 * sigma, color = 'gray', alpha = 0.5)

plt.plot(xs, ints_n)
plt.plot(xs, np.ones(len(xs)), color = 'r', linewidth = 2.0)

plt.xlim(x1, x2)

plt.show()

idx_spot = np.where((ints_n <= mean - 3 * sigma) | (ints_n >= mean + 3 * sigma))

print(xs[idx_spot])
