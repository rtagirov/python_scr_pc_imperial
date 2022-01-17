import oper_file as of
import numpy     as np

import spec; reload(spec)

wvl = np.loadtxt('../../idl/input/whi_uv.dat', skiprows = 700, usecols = [0])
irr = np.loadtxt('../../idl/input/whi_uv.dat', skiprows = 700, usecols = [3])

idx_euv = np.where((wvl >= 30.0) & (wvl < 116.0))

idx_sol = np.where((wvl >= 116.0) & (wvl < 310.0))

irr_euv_l = irr[idx_euv] * (1.0 - 15.0 / 100.0)
irr_euv_u = irr[idx_euv] * (1.0 + 15.0 / 100.0)

irr_sol_l = irr[idx_sol] * (1.0 - 7.0 / 100.0)
irr_sol_u = irr[idx_sol] * (1.0 + 7.0 / 100.0)

irr_l = np.concatenate((irr_euv_l, irr_sol_l))
irr_u = np.concatenate((irr_euv_u, irr_sol_u))

#wvls, irrs =   spec.running_mean(wvl, irr,   1.0)
#wvls, irrs_l = spec.running_mean(wvl, irr_l, 1.0)
#wvls, irrs_u = spec.running_mean(wvl, irr_u, 1.0)

#of.write4('../../idl/output/whi_uv_unc.dat', wvls, irrs_l, irrs, irrs_u)
of.write4('../../idl/output/whi_uv_unc.dat', wvl, irr_l, irr, irr_u)
