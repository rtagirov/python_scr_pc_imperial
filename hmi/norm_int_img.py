from astropy.io    import fits
from scipy.ndimage import filters
from tqdm          import tqdm

import numpy                       as np
import matplotlib.pyplot           as plt
import matplotlib.ticker           as ticker
import numpy.polynomial.polynomial as poly

import importlib
import sys

import auxfunc; importlib.reload(auxfunc)
import pltaux;  importlib.reload(pltaux)
import paths;   importlib.reload(paths)

#mode = 'skip'
mode = 'go'

image = 'int/2013_11_16'

filename = paths.inp + image + '.fits'

hdu = fits.open(filename)[0]

img = hdu.data

rad = its = np.array([])

x0 = hdu.header['X0']
y0 = hdu.header['Y0']

rsun = hdu.header['R_SUN']

r = np.sqrt(1 - 0.3**2.0) * rsun

if mode != 'skip':

    y1 = int(np.ceil(y0 - r))
    y2 = int(np.floor(y0 + r))

    for ys in tqdm(range(y1, y2), ncols = auxfunc.term_width(), desc = 'Scanning intensity image'):

        if ys % 10 != 0: continue

        x = np.arange(4096)
 
        x1 = x0 - np.sqrt(r**2 - (ys - y0)**2)
        x2 = x0 + np.sqrt(r**2 - (ys - y0)**2)

        idx_scan = np.where((x >= x1) & (x <= x2) & (x % 10 == 0))
#        idx_scan = np.where((x >= x1) & (x <= x2))

        xs = x[idx_scan]

        rs = np.sqrt((xs - x0)**2.0 + (ys - y0)**2.0)

        ints = img[ys, idx_scan[0]]

        its = np.concatenate((its, ints))

        rad = np.concatenate((rad, rs))

    c = poly.polyfit(rad, its, 6)

#    p = c[0] + c[1] * rad    + \
#               c[2] * rad**2 + \
#               c[3] * rad**3 + \
#               c[4] * rad**4 + \
#               c[5] * rad**5 + \
#               c[6] * rad**6

    p = poly.polyval(rad, c)

    img_ld = np.nan * np.ones((4096, 4096))

    y1 = int(np.ceil(y0 - rsun))
    y2 = int(np.floor(y0 + rsun))

    for ys in tqdm(range(y1, y2), ncols = auxfunc.term_width(), desc = 'Building CLV image'):

        x = np.arange(4096)
 
        x1 = x0 - np.sqrt(rsun**2 - (ys - y0)**2)
        x2 = x0 + np.sqrt(rsun**2 - (ys - y0)**2)

        idx_scan = np.where((x >= x1) & (x <= x2))

        xs = x[idx_scan]

        rs = np.sqrt((xs - x0)**2.0 + (ys - y0)**2.0)

        img_ld[ys, idx_scan[0]] = poly.polyval(rs, c)

#        img_ld[ys, idx_scan[0]] = c[0] + c[1] * rs    + \
#                                         c[2] * rs**2 + \
#                                         c[3] * rs**3 + \
#                                         c[4] * rs**4 + \
#                                         c[5] * rs**5 + \
#                                         c[6] * rs**6

np.savez(paths.npz + 'clv', rad = rad, its = its, p = p, img_ld = img_ld)

clv = np.load(paths.npz + 'clv.npz')

rad =    clv['rad']
its =    clv['its']
p =      clv['p']
img_ld = clv['img_ld']

pltaux.figpar(fontsize = 17)

plt.close('all')

fig, ax = plt.subplots(nrows=1, ncols=2, figsize = (17.0, 7.0))

sct_ax = ax[0]
img_ax = ax[1]

sct_ax.scatter(rad / rsun, its,              label = 'Data')
sct_ax.plot(rad / rsun,    p,   color = 'r', label = 'Fit', linewidth = 4.0)

sct_ax.set_xlim(0, 1)

sct_ax.set_xlabel(r'Impact parameter, $p = r / R_\odot$', labelpad = 10)
sct_ax.set_ylabel(r'Intensity, [arbitrary units]',        labelpad = 10)

sct_ax.xaxis.set_major_locator(ticker.MultipleLocator(0.1))

leg = sct_ax.legend(framealpha = 0, loc = 4, handletextpad=1, prop={'size': 25})

for obj in leg.legendHandles: obj.set_linewidth(5.0)

img = img_ax.imshow(img_ld, cmap = 'afmhot')

img_ax.set_xlabel('X position, [px]', labelpad = 10)
img_ax.set_ylabel('Y position, [px]', labelpad = 10)

img_ax.set_ylim(0, 4096)

cbar = fig.colorbar(img, ax = img_ax, fraction = 0.046, pad = 0.04)

cbar.set_label('Intensity, [arbitrary units]', labelpad = 10)

sct_ax.set_title(r'$6^\mathrm{th}$ order polynomial fit to the data')
img_ax.set_title('CLV image')

#fig.show()

#fig.savefig('/home/rtagirov/Dropbox/Work/Y1/fig/clv_img.png', bbox_inches = 'tight')
fig.savefig(paths.figdir + 'clv_img.png', bbox_inches = 'tight')
