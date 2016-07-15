import os
regressionParameters = {os.environ["PETSC_TO_POINTBEM"]+'/src/testSrfOnSurfacePoints':
                            [# Sphere tests [sphere_0-sphere_10]
                             {'num': 'sphere_0',  'numProcs': 1, 'args': '-srf_base ../geometry/sphere_R6_vdens0 -num_charges 10 -nmax 25 -srf_num 0125 -density 0.125 -ksp_rtol 1.0e-10 -snes_linesearch_type basic -snes_type ksponly'},
                             {'num': 'sphere_1',  'numProcs': 1, 'args': '-srf_base ../geometry/sphere_R6_vdens0 -num_charges 10 -nmax 25 -srf_num 025 -density 0.25 -ksp_rtol 1.0e-10 -snes_linesearch_type basic -snes_type ksponly'},
                             {'num': 'sphere_2',  'numProcs': 1, 'args': '-srf_base ../geometry/sphere_R6_vdens0 -num_charges 10 -nmax 25 -srf_num 05 -density 0.5 -ksp_rtol 1.0e-10 -snes_linesearch_type basic -snes_type ksponly'},
                             {'num': 'sphere_3',  'numProcs': 1, 'args': '-srf_base ../geometry/sphere_R6_vdens  -num_charges 10 -nmax 25 -srf_num 1 -density 1.0 -ksp_rtol 1.0e-10 -snes_linesearch_type basic -snes_type ksponly'},
                             {'num': 'sphere_4',  'numProcs': 1, 'args': '-srf_base ../geometry/sphere_R6_vdens  -num_charges 10 -nmax 25 -srf_num 2 -density 2.0 -ksp_rtol 1.0e-10 -snes_linesearch_type basic -snes_type ksponly'},
                             {'num': 'sphere_5',  'numProcs': 1, 'args': '-srf_base ../geometry/sphere_R6_vdens  -num_charges 10 -nmax 25 -srf_num 3 -density 3.0 -ksp_rtol 1.0e-10 -snes_linesearch_type basic -snes_type ksponly'},
                             {'num': 'sphere_6',  'numProcs': 1, 'args': '-srf_base ../geometry/sphere_R6_vdens  -num_charges 10 -nmax 25 -srf_num 4 -density 4.0 -ksp_rtol 1.0e-10 -snes_linesearch_type basic -snes_type ksponly'},
                             {'num': 'sphere_7',  'numProcs': 1, 'args': '-srf_base ../geometry/sphere_R6_vdens  -num_charges 10 -nmax 25 -srf_num 5 -density 5.0 -ksp_rtol 1.0e-10 -snes_linesearch_type basic -snes_type ksponly'},
                             {'num': 'sphere_8',  'numProcs': 1, 'args': '-srf_base ../geometry/sphere_R6_vdens  -num_charges 10 -nmax 25 -srf_num 6 -density 6.0 -ksp_rtol 1.0e-10 -snes_linesearch_type basic -snes_type ksponly'},
                             {'num': 'sphere_9',  'numProcs': 1, 'args': '-srf_base ../geometry/sphere_R6_vdens  -num_charges 10 -nmax 25 -srf_num 8 -density 8.0 -ksp_rtol 1.0e-10 -snes_linesearch_type basic -snes_type ksponly'},
                             {'num': 'sphere_10', 'numProcs': 1, 'args': '-srf_base ../geometry/sphere_R6_vdens  -num_charges 10 -nmax 25 -srf_num 10 -density 10.0 -ksp_rtol 1.0e-10 -snes_linesearch_type basic -snes_type ksponly'},
                             # Arginine tests [arg_0-arg_13]
                             {'num': 'arg_0',  'numProcs': 1, 'args': '-is_sphere 0 -pdb_filename ../geometry/arg.pdb -crg_filename ../geometry/jr1.crg -srf_base ../geometry/arg_scaledcharmm_ -srf_num 1'},
                             {'num': 'arg_1',  'numProcs': 1, 'args': '-is_sphere 0 -pdb_filename ../geometry/arg.pdb -crg_filename ../geometry/jr1.crg -srf_base ../geometry/arg_scaledcharmm_ -srf_num 2'},
                             {'num': 'arg_2',  'numProcs': 1, 'args': '-is_sphere 0 -pdb_filename ../geometry/arg.pdb -crg_filename ../geometry/jr1.crg -srf_base ../geometry/arg_scaledcharmm_ -srf_num 3'},
                             {'num': 'arg_3',  'numProcs': 1, 'args': '-is_sphere 0 -pdb_filename ../geometry/arg.pdb -crg_filename ../geometry/jr1.crg -srf_base ../geometry/arg_scaledcharmm_ -srf_num 4'},
                             {'num': 'arg_4',  'numProcs': 1, 'args': '-is_sphere 0 -pdb_filename ../geometry/arg.pdb -crg_filename ../geometry/jr1.crg -srf_base ../geometry/arg_scaledcharmm_ -srf_num 5'},
                             {'num': 'arg_5',  'numProcs': 1, 'args': '-is_sphere 0 -pdb_filename ../geometry/arg.pdb -crg_filename ../geometry/jr1.crg -srf_base ../geometry/arg_scaledcharmm_ -srf_num 6'},
                             {'num': 'arg_6',  'numProcs': 1, 'args': '-is_sphere 0 -pdb_filename ../geometry/arg.pdb -crg_filename ../geometry/jr1.crg -srf_base ../geometry/arg_scaledcharmm_ -srf_num 7'},
                             {'num': 'arg_7',  'numProcs': 1, 'args': '-is_sphere 0 -pdb_filename ../geometry/arg.pdb -crg_filename ../geometry/jr1.crg -srf_base ../geometry/arg_scaledcharmm_ -srf_num 8'},
                             {'num': 'arg_8',  'numProcs': 1, 'args': '-is_sphere 0 -pdb_filename ../geometry/arg.pdb -crg_filename ../geometry/jr1.crg -srf_base ../geometry/arg_scaledcharmm_ -srf_num 12'},
                             {'num': 'arg_9',  'numProcs': 1, 'args': '-is_sphere 0 -pdb_filename ../geometry/arg.pdb -crg_filename ../geometry/jr1.crg -srf_base ../geometry/arg_scaledcharmm_ -srf_num 16'},
                             {'num': 'arg_10', 'numProcs': 1, 'args': '-is_sphere 0 -pdb_filename ../geometry/arg.pdb -crg_filename ../geometry/jr1.crg -srf_base ../geometry/arg_scaledcharmm_ -srf_num 20'},
                             {'num': 'arg_11', 'numProcs': 1, 'args': '-is_sphere 0 -pdb_filename ../geometry/arg.pdb -crg_filename ../geometry/jr1.crg -srf_base ../geometry/arg_scaledcharmm_ -srf_num 24'},
                             {'num': 'arg_12', 'numProcs': 1, 'args': '-is_sphere 0 -pdb_filename ../geometry/arg.pdb -crg_filename ../geometry/jr1.crg -srf_base ../geometry/arg_scaledcharmm_ -srf_num 28'},
                             {'num': 'arg_13', 'numProcs': 1, 'args': '-is_sphere 0 -pdb_filename ../geometry/arg.pdb -crg_filename ../geometry/jr1.crg -srf_base ../geometry/arg_scaledcharmm_ -srf_num 32'},
                             # Aspartic acid tests [asp_0-asp_13]
                             {'num': 'asp_0',  'numProcs': 1, 'args': '-is_sphere 0 -pdb_filename ../geometry/asp.pdb -crg_filename ../geometry/jd1.crg -srf_base ../geometry/asp_scaledcharmm_ -srf_num 1'},
                             {'num': 'asp_1',  'numProcs': 1, 'args': '-is_sphere 0 -pdb_filename ../geometry/asp.pdb -crg_filename ../geometry/jd1.crg -srf_base ../geometry/asp_scaledcharmm_ -srf_num 2'},
                             {'num': 'asp_2',  'numProcs': 1, 'args': '-is_sphere 0 -pdb_filename ../geometry/asp.pdb -crg_filename ../geometry/jd1.crg -srf_base ../geometry/asp_scaledcharmm_ -srf_num 3'},
                             {'num': 'asp_3',  'numProcs': 1, 'args': '-is_sphere 0 -pdb_filename ../geometry/asp.pdb -crg_filename ../geometry/jd1.crg -srf_base ../geometry/asp_scaledcharmm_ -srf_num 4'},
                             {'num': 'asp_4',  'numProcs': 1, 'args': '-is_sphere 0 -pdb_filename ../geometry/asp.pdb -crg_filename ../geometry/jd1.crg -srf_base ../geometry/asp_scaledcharmm_ -srf_num 5'},
                             {'num': 'asp_5',  'numProcs': 1, 'args': '-is_sphere 0 -pdb_filename ../geometry/asp.pdb -crg_filename ../geometry/jd1.crg -srf_base ../geometry/asp_scaledcharmm_ -srf_num 6'},
                             {'num': 'asp_6',  'numProcs': 1, 'args': '-is_sphere 0 -pdb_filename ../geometry/asp.pdb -crg_filename ../geometry/jd1.crg -srf_base ../geometry/asp_scaledcharmm_ -srf_num 7'},
                             {'num': 'asp_7',  'numProcs': 1, 'args': '-is_sphere 0 -pdb_filename ../geometry/asp.pdb -crg_filename ../geometry/jd1.crg -srf_base ../geometry/asp_scaledcharmm_ -srf_num 8'},
                             {'num': 'asp_8',  'numProcs': 1, 'args': '-is_sphere 0 -pdb_filename ../geometry/asp.pdb -crg_filename ../geometry/jd1.crg -srf_base ../geometry/asp_scaledcharmm_ -srf_num 12'},
                             {'num': 'asp_9',  'numProcs': 1, 'args': '-is_sphere 0 -pdb_filename ../geometry/asp.pdb -crg_filename ../geometry/jd1.crg -srf_base ../geometry/asp_scaledcharmm_ -srf_num 16'},
                             {'num': 'asp_10', 'numProcs': 1, 'args': '-is_sphere 0 -pdb_filename ../geometry/asp.pdb -crg_filename ../geometry/jd1.crg -srf_base ../geometry/asp_scaledcharmm_ -srf_num 20'},
                             {'num': 'asp_11', 'numProcs': 1, 'args': '-is_sphere 0 -pdb_filename ../geometry/asp.pdb -crg_filename ../geometry/jd1.crg -srf_base ../geometry/asp_scaledcharmm_ -srf_num 24'},
                             {'num': 'asp_12', 'numProcs': 1, 'args': '-is_sphere 0 -pdb_filename ../geometry/asp.pdb -crg_filename ../geometry/jd1.crg -srf_base ../geometry/asp_scaledcharmm_ -srf_num 28'},
                             {'num': 'asp_13', 'numProcs': 1, 'args': '-is_sphere 0 -pdb_filename ../geometry/asp.pdb -crg_filename ../geometry/jd1.crg -srf_base ../geometry/asp_scaledcharmm_ -srf_num 32'},
                             ],
                        }
