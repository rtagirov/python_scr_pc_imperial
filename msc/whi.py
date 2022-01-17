import oper_file as of
import numpy     as np
import math      as m

import spec; reload(spec)

wvl_whi = np.loadtxt('../../idl/input/whi.dat', skiprows = 1700, usecols = [0])
irr_whi = np.loadtxt('../../idl/input/whi.dat', skiprows = 1700, usecols = [3])

idx_sol = np.where(wvl_whi < 310.0)

wvls_sol, irrs_sol = spec.running_mean(wvl_whi[idx_sol], irr_whi[idx_sol], 1.0)

idx_sim  = np.where(wvl_whi >= 310.0)

irr_sol_l = irrs_sol * (1.0 - 3.0 / 100.0)
irr_sol_u = irrs_sol * (1.0 + 3.0 / 100.0)

irr_sim_l = irr_whi[idx_sim] * (1.0 - 1.5 / 100.0)
irr_sim_u = irr_whi[idx_sim] * (1.0 + 1.5 / 100.0)

irr_l = np.concatenate((irr_sol_l, irr_sim_l))
irr_u = np.concatenate((irr_sol_u, irr_sim_u))

wvl = np.concatenate((wvls_sol, wvl_whi[idx_sim]))

irr = np.concatenate((irrs_sol, irr_whi[idx_sim]))

of.write4('../../idl/output/whi_unc_min.dat', wvl, irr_l, irr, irr_u)

irr_sol_l = irrs_sol * (1.0 - 7.0 / 100.0)
irr_sol_u = irrs_sol * (1.0 + 7.0 / 100.0)

irr_sim_l = irr_whi[idx_sim] * (1.0 - 3.5 / 100.0)
irr_sim_u = irr_whi[idx_sim] * (1.0 + 3.5 / 100.0)

irr_l = np.concatenate((irr_sol_l, irr_sim_l))
irr_u = np.concatenate((irr_sol_u, irr_sim_u))

of.write4('../../idl/output/whi_unc_max.dat', wvl, irr_l, irr, irr_u)
