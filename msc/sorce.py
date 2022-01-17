import oper_file as of
import numpy     as np

wvl_sorce = np.loadtxt('../../idl/input/sorce.dat', skiprows = 50, usecols = [0])
irr_sorce = np.loadtxt('../../idl/input/sorce.dat', skiprows = 50, usecols = [1])

idx_sol = np.where(wvl_sorce < 320.0)

idx_sim  = np.where(wvl_sorce >= 320.0)

irr_sol_l = irr_sorce[idx_sol] * (1.0 - 3.0 / 100.0)
irr_sol_u = irr_sorce[idx_sol] * (1.0 + 3.0 / 100.0)

irr_sim_l = irr_sorce[idx_sim] * (1.0 - 1.5 / 100.0)
irr_sim_u = irr_sorce[idx_sim] * (1.0 + 1.5 / 100.0)

irr_l = np.concatenate((irr_sol_l, irr_sim_l))
irr_u = np.concatenate((irr_sol_u, irr_sim_u))

of.write4('../../idl/output/sorce_unc_min.dat', wvl_sorce, irr_l, irr_sorce, irr_u)

irr_sol_l = irr_sorce[idx_sol] * (1.0 - 7.0 / 100.0)
irr_sol_u = irr_sorce[idx_sol] * (1.0 + 7.0 / 100.0)

irr_sim_l = irr_sorce[idx_sim] * (1.0 - 3.5 / 100.0)
irr_sim_u = irr_sorce[idx_sim] * (1.0 + 3.5 / 100.0)

irr_l = np.concatenate((irr_sol_l, irr_sim_l))
irr_u = np.concatenate((irr_sol_u, irr_sim_u))

of.write4('../../idl/output/sorce_unc_max.dat', wvl_sorce, irr_l, irr_sorce, irr_u)
