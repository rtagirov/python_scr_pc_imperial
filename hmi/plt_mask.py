from astropy.io    import fits
from tqdm          import tqdm
from scipy.ndimage import filters
from scipy.ndimage import measurements

import numpy                       as np
import matplotlib.pyplot           as plt
import matplotlib.gridspec         as gridspec
import matplotlib.ticker           as ticker
import matplotlib.patches          as patches
import numpy.polynomial.polynomial as poly

import sys
import importlib
import cv2

import paths;   importlib.reload(paths)
import phys;    importlib.reload(phys)
import auxfunc; importlib.reload(auxfunc)
import pltaux;  importlib.reload(pltaux)
import mask;    importlib.reload(mask)

def rad_mu_pxa(x0, y0, rsun):

    rad = np.zeros((4096, 4096))

    for x in tqdm(range(4096), ncols = auxfunc.term_width(), desc = 'Calculating r, mu and pixel areas'):

        for y in range(4096):

            rad[y, x] = (x - x0)**2.0 + (y - y0)**2.0

    rad = np.sqrt(rad)

    mu = np.sqrt(1.0 - (rad / rsun)**2.0)

    pxa_norm = (phys.solar_radius * 1.0e-5 / rsun)**2.0

    pxa = np.zeros((4096, 4096))

    mu_is_nan_idx = np.where(np.isnan(mu))

    mu_is_not_nan_idx = np.where(np.logical_not(np.isnan(mu)))

    pxa[mu_is_nan_idx] = np.nan

    pxa[mu_is_not_nan_idx] = pxa_norm / mu[mu_is_not_nan_idx]

    return rad, mu, pxa

def clv(img, x0, y0, rsun):

    r = np.sqrt(1 - 0.3**2.0) * rsun

    rad = its = np.array([])

    y1 = int(np.ceil(y0 - r))
    y2 = int(np.ceil(y0 + r))

    for ys in tqdm(range(y1, y2), ncols = auxfunc.term_width(), desc = 'Building I(r)'):

        if ys % 2 != 0: continue

        x = np.arange(4096)
 
        x1 = x0 - np.sqrt(r**2 - (ys - y0)**2)
        x2 = x0 + np.sqrt(r**2 - (ys - y0)**2)

        idx_scan = np.where((x >= x1) & (x <= x2) & (x % 2 == 0))

        xs = x[idx_scan]

        rs = np.sqrt((xs - x0)**2.0 + (ys - y0)**2.0)

        ints = img[ys, idx_scan[0]]

        its = np.concatenate((its, ints))

        rad = np.concatenate((rad, rs))

    c = poly.polyfit(rad, its, 6)

    img_ld = np.zeros((4096, 4096))

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

    return img_ld

def get_cnt_coord(cnt):

    x = np.array([])
    y = np.array([])

    for elem in cnt:

        x = np.concatenate((x, [elem[0][0]]))
        y = np.concatenate((y, [elem[0][1]]))

    first_elem = cnt[0]

    x = np.concatenate((x, [first_elem[0][0]]))
    y = np.concatenate((y, [first_elem[0][1]]))

    return x, y

p = np.arange(1, 11) / 10.0

date = '2013_11_16'

year = date[0 : 4]; month = date[5 : 7]; day = date[8 : 10]

date_df = day + '.' + month + '.' + year

filename_int = 'int/' + date + '.fits'
filename_mag = 'mag/' + date + '.fits'

hdulist_int = fits.open(paths.inp + filename_int)
hdulist_mag = fits.open(paths.inp + filename_mag)

hdu_int = hdulist_int[0]
hdu_mag = hdulist_mag[0]

x0 =   hdu_int.header['X0']
y0 =   hdu_int.header['Y0']
rsun = hdu_int.header['R_SUN']

r = np.sqrt(1 - 0.3**2.0) * rsun

I = hdu_int.data
B = hdu_mag.data

#rad, mu, pxa = rad_mu_pxa(x0, y0, rsun)

#Is = filters.gaussian_filter(I, 10)

#In = I / clv(I, x0, y0, rsun)

#Isn = Is / clv(Is, x0, y0, rsun)

#np.savez(paths.npz + date + '_Iradmupxa', rad = rad, mu = mu, pxa = pxa, In = In, Isn = Isn)

Iradmupxa = np.load(paths.npz + date + '_Iradmupxa.npz')

rad = Iradmupxa['rad']
mu =  Iradmupxa['mu']
pxa = Iradmupxa['pxa']
In =  Iradmupxa['In']
Isn = Iradmupxa['Isn']

#spt_mask, umb_mask, pen_mask, sff, uff, pff = mask.intensity(In, r, rad)

#spt_smask, umb_smask, pen_smask, ssff, suff, spff = mask.intensity(Isn, r, rad)

#mag_mask = mask.magnetic_field(B, r, rad)

#fac_mask = mag_mask - spt_mask

#np.savez(paths.npz + date + '_masks', mag_mask  = mag_mask, \
#                                      spt_mask  = spt_mask, \
#                                      umb_mask  = umb_mask, \
#                                      pen_mask  = pen_mask, \
#                                      spt_smask = spt_smask, \
#                                      umb_smask = umb_smask, \
#                                      pen_smask = pen_smask, \
#                                      fac_mask  = fac_mask)

masks = np.load(paths.npz + date + '_masks.npz')

spt_mask  = masks['spt_mask']
umb_mask  = masks['umb_mask']
pen_mask  = masks['pen_mask']
spt_smask = masks['spt_smask']
umb_smask = masks['umb_smask']
pen_smask = masks['pen_smask']
mag_mask  = masks['mag_mask']
fac_mask  = masks['fac_mask']

fff = len(fac_mask[np.where(fac_mask == 1.0)]) / len(fac_mask[np.where(rad <= r)])

spt_shapes = cv2.inRange(spt_mask, 0.9, 1.1)
umb_shapes = cv2.inRange(umb_mask, 0.9, 1.1)

w1, spt_cnts, w2 = cv2.findContours(spt_shapes.copy(), cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
w1, umb_cnts, w2 = cv2.findContours(umb_shapes.copy(), cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)

B = np.zeros((4096, 4096))

B[np.where(mag_mask == 1)] = hdu_mag.data[np.where(mag_mask == 1)]
B[np.where(np.isnan(hdu_mag.data) | (mag_mask == 0))] = np.nan

fac_neg_idx = np.where(np.logical_not(np.isnan(B)) & (fac_mask == 1) & (B < 0.0))
fac_pos_idx = np.where(np.logical_not(np.isnan(B)) & (fac_mask == 1) & (B > 0.0))

umb_neg_idx = np.where(np.logical_not(np.isnan(B)) & (umb_mask == 1) & (umb_mask != -1) & (B < 0.0))
umb_pos_idx = np.where(np.logical_not(np.isnan(B)) & (umb_mask == 1) & (umb_mask != -1) & (B > 0.0))

fac_idx = np.where(np.logical_not(np.isnan(B)) & (fac_mask == 1))

umb_idx = np.where(np.logical_not(np.isnan(B)) & (umb_mask == 1) & (fac_mask != -1))

B_fac_neg_mean = np.mean(B[fac_neg_idx] / mu[fac_neg_idx])
B_fac_pos_mean = np.mean(B[fac_pos_idx] / mu[fac_pos_idx])

B_umb_neg_mean = np.mean(B[umb_neg_idx] / mu[umb_neg_idx])
B_umb_pos_mean = np.mean(B[umb_pos_idx] / mu[umb_pos_idx])

B_fac_mean = np.mean(abs(B[fac_idx]) / mu[fac_idx])
B_umb_mean = np.mean(abs(B[umb_idx]) / mu[umb_idx])

x = np.arange(4096)

idx_ir = np.where(abs(x - x0) <= r)
idx_or = np.where(abs(x - x0) <= rsun)

iuc = y0 + np.sqrt(r**2 -    (x[idx_ir] - x0)**2)
ilc = y0 - np.sqrt(r**2 -    (x[idx_ir] - x0)**2)
ouc = y0 + np.sqrt(rsun**2 - (x[idx_or] - x0)**2)
olc = y0 - np.sqrt(rsun**2 - (x[idx_or] - x0)**2)

spt_lmask, spt_num  = measurements.label(spt_smask)

sn = np.arange(1, spt_num + 1)

spt_diam = np.zeros(spt_num)

for i in sn: spt_diam[i - 1] = 2.0 * np.sqrt(sum(pxa[np.where(spt_lmask == i)]) / np.pi)

pltaux.figpar(fontsize = 17)

fig3, ax = plt.subplots(nrows=1, ncols=2, figsize = (10, 5))

fig3.tight_layout()

plt.subplots_adjust(wspace = 0.20)

img1 = ax[0].imshow(B, cmap = 'gray')

ax[0].plot(x[idx_ir], iuc, '--', color = 'k', linewidth = 0.2)
ax[0].plot(x[idx_ir], ilc, '--', color = 'k', linewidth = 0.2)

ax[0].plot(x[idx_or], ouc, '-', color = 'k', linewidth = 0.6)
ax[0].plot(x[idx_or], olc, '-', color = 'k', linewidth = 0.6)

ax[0].add_patch(patches.Rectangle((700, 1900), 500, 550, fill = False, linewidth = 0.4))

for cnt in spt_cnts:

    xc, yc = get_cnt_coord(cnt)

    ax[0].plot(xc, yc, '-', color = 'r', linewidth = 0.5)

for cnt in umb_cnts:

    xc, yc = get_cnt_coord(cnt)

    ax[0].plot(xc, yc, '-', color = 'b', linewidth = 0.5)

ax[0].set_xlim(0, 4096)
ax[0].set_ylim(0, 4096)

ax[0].set_xlabel('X position, [px]', labelpad = 10)
ax[0].set_ylabel('Y position, [px]', labelpad = 10)

img2 = ax[1].imshow(B, cmap = 'gray')

for cnt in spt_cnts:

    xc, yc = get_cnt_coord(cnt)

    ax[1].plot(xc, yc, '-', color = 'r', linewidth = 0.7)

for cnt in umb_cnts:

    xc, yc = get_cnt_coord(cnt)

    ax[1].plot(xc, yc, '-', color = 'b', linewidth = 0.7)

ax[1].set_xlim(700, 1200)
ax[1].set_ylim(1900, 2450)

ax[1].set_xlabel('X position, [px]', labelpad = 10)
#ax[1].set_ylabel('Y position, [px]', labelpad = 10)

cbar = fig3.colorbar(img2, fraction = 0.050, pad = 0.04)

cbar.set_label('LOS magnetic field, [Gauss]', labelpad = 10)

fig3.savefig('../fig/spt_cnts.png', bbox_inches = 'tight')

#fig3.show()

sys.exit()

fig4, ax = plt.subplots()

ax.imshow(spt_shapes, cmap = 'gray_r')

ax.plot(x[idx_ir], iuc, '--', color = 'k')
ax.plot(x[idx_ir], ilc, '--', color = 'k')

ax.plot(x[idx_or], ouc, '-', color = 'k')
ax.plot(x[idx_or], olc, '-', color = 'k')

ax.set_xlim(0, 4096)
ax.set_ylim(0, 4096)

f3 = plt.figure(figsize = (13, 8.26))

#f3.clf()

#plt.scatter(sn, spt_diam)

#plt.yscale('log')

#plt.imshow(umb_mask + pen_mask, cmap = 'gray_r')
plt.imshow(spt_mask, cmap = 'gray_r')

plt.plot(x[idx_ir], iuc, '--', color = 'k')
plt.plot(x[idx_ir], ilc, '--', color = 'k')

plt.plot(x[idx_or], ouc, '-', color = 'k')
plt.plot(x[idx_or], olc, '-', color = 'k')

plt.ylim(0, 4096)

plt.xlabel('X position, [px]', labelpad = 10)
plt.ylabel('Y position, [px]', labelpad = 10)

#f3.savefig('/home/rtagirov/Dropbox/Work/Y1/fig/spt_mask.png', bbox_inches = 'tight')

#f3.show()

#sys.exit()

plt.imshow(spt_smask, cmap = 'gray_r')

plt.plot(x[idx_ir], iuc, '--', color = 'k')
plt.plot(x[idx_ir], ilc, '--', color = 'k')

plt.plot(x[idx_or], ouc, '-', color = 'k')
plt.plot(x[idx_or], olc, '-', color = 'k')

f2 = plt.figure()

plt.imshow(spt_lmask)

plt.plot(x[idx_ir], iuc, '--', color = 'k')
plt.plot(x[idx_ir], ilc, '--', color = 'k')

plt.plot(x[idx_or], ouc, '-', color = 'k')
plt.plot(x[idx_or], olc, '-', color = 'k')

plt.ylim(0, 4096)

#f2.show()

#sys.exit()

#plt.close('all')

pltaux.figpar(fontsize = 17)

#fig1 = plt.figure(figsize = (15.0, 7.5))

#gs = gridspec.GridSpec(nrows=1, ncols=2, left=0.06, bottom=0.00, right=0.92, top=0.99,
#                       wspace=0.45, hspace=0.0, width_ratios=[1, 1])

#int_ax = plt.subplot(gs[0])
#mag_ax = plt.subplot(gs[1])

fig1, ax = plt.subplots(nrows=1, ncols=2, figsize = (10, 5))

fig1.tight_layout()

plt.subplots_adjust(wspace = 0.25)

#int_img = int_ax.imshow(hdu_int.data, cmap = 'afmhot')
#mag_img = mag_ax.imshow(hdu_mag.data, cmap = 'gray')
zoom1 = ax[0].imshow(In, cmap = 'afmhot')
zoom2 = ax[1].imshow(In, cmap = 'afmhot')

#int_ax.set_ylim(0, 4096)
#mag_ax.set_ylim(0, 4096)

ax[0].set_xlim(750,  2000)
ax[1].set_xlim(2000, 3000)

ax[0].set_ylim(2000, 2500)
ax[1].set_ylim(1300, 1700)

#int_ax.set_xlabel('X position, [px]', labelpad = 10)
#int_ax.set_ylabel('Y position, [px]', labelpad = 10)
#mag_ax.set_xlabel('X position, [px]', labelpad = 10)
#mag_ax.set_ylabel('Y position, [px]', labelpad = 10)

ax[0].set_xlabel('X position, [px]', labelpad = 10)
ax[0].set_ylabel('Y position, [px]', labelpad = 10)
ax[1].set_xlabel('X position, [px]', labelpad = 10)

#mag_ax.xaxis.set_major_locator(ticker.MultipleLocator(250))
#mag_ax.yaxis.set_major_locator(ticker.MultipleLocator(200))

#int_ax.set_title('Intensity image')
#mag_ax.set_title('Magnetogram')

#int_cbar = fig1.colorbar(int_img, ax = int_ax, fraction = 0.045, pad = 0.06)
#mag_cbar = fig1.colorbar(mag_img, ax = mag_ax, fraction = 0.045, pad = 0.06)

#int_cbar.set_label('Normalized intensity', labelpad = 10)
#mag_cbar.set_label('Normalized intensity', labelpad = 10, fontsize = 12.5)

#int_cbar.set_label('Intensity at 6173 \AA, [arbitrary units]', labelpad = 10)
#mag_cbar.set_label('Line of sight magnetic field, [Gauss]',    labelpad = 10)

#cbar = fig1.colorbar(zoom2, ax = ax[1], fraction = 0.019, pad = 0.04)

#cbar.set_label('Normalized intensity', labelpad = 10, fontsize = 13)
#cbar.set_label('Normalized intensity', labelpad = 10)

#fig1.show()

#fig1.savefig('/home/rtagirov/Dropbox/Work/Y1/fig/smth_spot.png', bbox_inches = 'tight')
#fig1.savefig('../fig/sdo_img.png', bbox_inches = 'tight')
fig1.savefig(paths.figdir + 'sdo_int_img_zoom.png', bbox_inches = 'tight')

#sys.exit()

pltaux.figpar(fontsize = 25)

#fig2, ax = plt.subplots(nrows = 1, ncols = 3, figsize = (15.0, 7.5))
fig2, ax = plt.subplots(nrows = 1, ncols = 2, figsize = (20, 10))

fig2.tight_layout()

fig2.subplots_adjust(wspace = 0.20)

ax[0].imshow(mag_mask, cmap = 'gray_r')

ax[0].plot(x[idx_ir], iuc, '--', color = 'k')
ax[0].plot(x[idx_ir], ilc, '--', color = 'k')

ax[0].plot(x[idx_or], ouc, '-', color = 'k')
ax[0].plot(x[idx_or], olc, '-', color = 'k')

ax[0].set_xlim(0, 4096)
ax[0].set_ylim(0, 4096)

ax[1].imshow(fac_mask, cmap = 'gray_r')

ax[1].plot(x[idx_ir], iuc, '--', color = 'k')
ax[1].plot(x[idx_ir], ilc, '--', color = 'k')

ax[1].plot(x[idx_or], ouc, '-', color = 'k')
ax[1].plot(x[idx_or], olc, '-', color = 'k')

ax[1].set_xlim(0, 4096)
ax[1].set_ylim(0, 4096)

ax[0].set_xlabel('X position, [px]', labelpad = 10)
ax[0].set_ylabel('Y position, [px]', labelpad = 10)
ax[1].set_xlabel('X position, [px]', labelpad = 10)
ax[1].set_ylabel('Y position, [px]', labelpad = 10)

ax[0].set_title('Magnetic mask')
ax[1].set_title('Facular mask')

fig2.savefig('/home/rtagirov/Dropbox/Work/Y1/fig/mag_fac_mask.png', bbox_inches = 'tight')

#ax[2].imshow(fac_mask, cmap = 'gray_r')

#ax[2].plot(x[idx_ir], iuc, '--', color = 'k')
#ax[2].plot(x[idx_ir], ilc, '--', color = 'k')

#ax[2].plot(x[idx_or], ouc, '-', color = 'k')
#ax[2].plot(x[idx_or], olc, '-', color = 'k')

#ax[2].set_xlim(0, 4096)
#ax[2].set_ylim(0, 4096)

fig3, ax = plt.subplots(nrows = 1, ncols = 2, figsize = (17, 7))

img = ax[0].imshow(B, cmap = 'gray')
#img = ax[0].imshow(I, cmap = 'gray')
#img = ax[0].imshow(hdu_int.data, cmap = 'gray')

for i in range(len(p)):

    x = np.arange(4096)

    rp = rsun * p[i]

    idx_p = np.where(abs(x - x0) <= rp)

    uc = y0 + np.sqrt(rp**2.0 - (x[idx_p] - x0)**2.0)
    lc = y0 - np.sqrt(rp**2.0 - (x[idx_p] - x0)**2.0)

    if i == len(p) - 1:

        ax[0].plot(x[idx_p], uc, color = 'k', linewidth = 0.6)
        ax[0].plot(x[idx_p], lc, color = 'k', linewidth = 0.6)

    else:

        ax[0].plot(x[idx_p], uc, '--', color = 'k', linewidth = 0.2)
        ax[0].plot(x[idx_p], lc, '--', color = 'k', linewidth = 0.2)

for cnt in spt_cnts:

    xc, yc = get_cnt_coord(cnt)

    ax[0].plot(xc, yc, '-', color = 'r')

for cnt in umb_cnts:

    xc, yc = get_cnt_coord(cnt)

    ax[0].plot(xc, yc, '-', color = 'b')

ax[0].set_xlim(0, 4096)
ax[0].set_ylim(0, 4096)

cbar = fig3.colorbar(img, ax = ax[0], fraction = 0.046, pad = 0.04)

idx_fac = np.where(fac_mask == 1)
idx_umb = np.where((umb_mask == 1) & (fac_mask != -1))

ax[1].scatter(rad[idx_fac] / rsun, B[idx_fac] / mu[idx_fac], color = 'gray', s = 12)
ax[1].scatter(rad[idx_umb] / rsun, B[idx_umb] / mu[idx_umb], color = 'r',    s = 0.5, alpha = 0.2)

ax[1].plot(np.array([p[0], p[8]]), np.zeros(2), '--', color = 'w', linewidth = 1)

ax[1].plot(np.array([p[0], p[8]]), B_fac_pos_mean * np.ones(2), color = 'k')
ax[1].plot(np.array([p[0], p[8]]), B_fac_neg_mean * np.ones(2), color = 'k')

ax[1].plot(np.array([p[0], p[8]]), B_umb_pos_mean * np.ones(2), color = 'r')
ax[1].plot(np.array([p[0], p[8]]), B_umb_neg_mean * np.ones(2), color = 'r')

ax[1].yaxis.tick_right()

ax[1].xaxis.grid(True)

fig3.savefig('../fig/spt_fac_mf_scat.png')

fig4, ax = plt.subplots()

ax.imshow(spt_shapes, cmap = 'gray_r')

ax.plot(x[idx_ir], iuc, '--', color = 'k')
ax.plot(x[idx_ir], ilc, '--', color = 'k')

ax.plot(x[idx_or], ouc, '-', color = 'k')
ax.plot(x[idx_or], olc, '-', color = 'k')

ax.set_xlim(0, 4096)
ax.set_ylim(0, 4096)
