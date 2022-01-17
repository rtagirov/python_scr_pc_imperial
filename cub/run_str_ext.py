import str_ext_new

from tqdm import tqdm

reload(str_ext_new)

#str_ext_new.Slice_rot_3D('/mnt/HDD/cubes', '/mnt/HDD/extr_atms', 'G2V', 'hydro', '123000', 10.0, 24, 0)
#str_ext_new.Slice_rot_3D('/mnt/HDD/cubes', '/mnt/HDD/extr_atms', 'G2V', '300G',  '118000', 10.0, 24, 0)

#str_ext_new.ext_cube_write('/mnt/HDD/extr_atms', '/mnt/HDD/cubes', 'G2V', 'hydro', '123000')

#str_ext_new.ext_cube('/mnt/HDD/cubes', '/mnt/HDD/extr_atms', 'G2V', 'hydro', '123000', 10.0, 0, 0)

for i in tqdm(range(0, 25)):
#for i in tqdm(range(0, 1)):
#for i in range(10, 11):
#for i in range(1, 2):
#for i in range(0, 1):
#for i in tqdm(range(20, 21)):

#    str_ext_new.ext_cube('/mnt/HDD/cubes', '/mnt/HDD/extr_atms', 'G2V', 'hydro', '123000', 10.0, i, 0)
    str_ext_new.ext_cube('/mnt/HDD/cubes', '/mnt/HDD/extr_atms', 'G2V', '300G', '118000', 10.0, i, 0)
