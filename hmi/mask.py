import numpy as np

def intensity(I, r, rad):

#    sigma = np.std(I[np.where(np.logical_not(np.isnan(I)))])

    sigma = np.sqrt(np.mean((1.0 - I[np.where(np.logical_not(np.isnan(I)))])**2.0))

    I[np.where(np.isnan(I))] = 0.0

    Iu = 0.53

    mask_s = np.zeros((4096, 4096))
    mask_u = np.zeros((4096, 4096))
    mask_p = np.zeros((4096, 4096))

    mask_s[np.where((rad <= r) & (I <= 1.0 - 3.0 * sigma) & (I != 0.0))] = 1.0

    mask_u[np.where((rad <= r) & (I <= Iu) & (I != 0.0))] = 1.0

    mask_p[np.where((rad <= r) & (I > Iu) & (I <= 1.0 - 3.0 * sigma) & (I != 0.0))] = 0.5

    sff = len(mask_s[np.where(mask_s == 1.0)]) / len(mask_s[np.where(rad <= r)])
    uff = len(mask_u[np.where(mask_u == 1.0)]) / len(mask_u[np.where(rad <= r)])
    pff = len(mask_p[np.where(mask_p == 0.5)]) / len(mask_p[np.where(rad <= r)])

    return mask_s, mask_u, mask_p, sff, uff, pff

def magnetic_field(B, r, rad):

    B[np.where(np.isnan(B))] = 0.0

    sigma = 10.0

    mask = np.zeros((4096, 4096))

    mask[np.where((rad <= r) & ((B <= -3.0 * sigma) | (B >= 3.0 * sigma)) & (B != 0.0))] = 1.0

    return mask
