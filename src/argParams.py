import os
petscToPointbem = os.path.relpath(os.path.abspath('..'), os.environ['PETSC_DIR'])
print petscToPointbem
regressionParameters = {os.environ.get("PETSC_TO_POINTBEM", petscToPointbem)+'/src/testSrfOnSurfacePoints':
                        [# Arginine tests [arg_0-arg_13]
                            {'num': 'arg_0',  'numProcs': 1, 'args': '-is_sphere 0 -pdb_filename ../geometry/arg.pdb -crg_filename ../geometry/jr1.crg -srf_base ../geometry/arg_scaledcharmm_ -srf_num 1 -ksp_rtol 1.0e-10'},
                            {'num': 'arg_1',  'numProcs': 1, 'args': '-is_sphere 0 -pdb_filename ../geometry/arg.pdb -crg_filename ../geometry/jr1.crg -srf_base ../geometry/arg_scaledcharmm_ -srf_num 2 -ksp_rtol 1.0e-10'}, 
                            {'num': 'arg_2',  'numProcs': 1, 'args': '-is_sphere 0 -pdb_filename ../geometry/arg.pdb -crg_filename ../geometry/jr1.crg -srf_base ../geometry/arg_scaledcharmm_ -srf_num 3 -ksp_rtol 1.0e-10'},
                            {'num': 'arg_3',  'numProcs': 1, 'args': '-is_sphere 0 -pdb_filename ../geometry/arg.pdb -crg_filename ../geometry/jr1.crg -srf_base ../geometry/arg_scaledcharmm_ -srf_num 4 -ksp_rtol 1.0e-10'},
                            {'num': 'arg_4',  'numProcs': 1, 'args': '-is_sphere 0 -pdb_filename ../geometry/arg.pdb -crg_filename ../geometry/jr1.crg -srf_base ../geometry/arg_scaledcharmm_ -srf_num 5 -ksp_rtol 1.0e-10'},
                            {'num': 'arg_5',  'numProcs': 1, 'args': '-is_sphere 0 -pdb_filename ../geometry/arg.pdb -crg_filename ../geometry/jr1.crg -srf_base ../geometry/arg_scaledcharmm_ -srf_num 6 -ksp_rtol 1.0e-10'},
                            {'num': 'arg_6',  'numProcs': 1, 'args': '-is_sphere 0 -pdb_filename ../geometry/arg.pdb -crg_filename ../geometry/jr1.crg -srf_base ../geometry/arg_scaledcharmm_ -srf_num 7 -ksp_rtol 1.0e-10'},
                            {'num': 'arg_7',  'numProcs': 1, 'args': '-is_sphere 0 -pdb_filename ../geometry/arg.pdb -crg_filename ../geometry/jr1.crg -srf_base ../geometry/arg_scaledcharmm_ -srf_num 8 -ksp_rtol 1.0e-10'},
                            {'num': 'arg_8',  'numProcs': 1, 'args': '-is_sphere 0 -pdb_filename ../geometry/arg.pdb -crg_filename ../geometry/jr1.crg -srf_base ../geometry/arg_scaledcharmm_ -srf_num 12 -ksp_rtol 1.0e-10'},
                            {'num': 'arg_9',  'numProcs': 1, 'args': '-is_sphere 0 -pdb_filename ../geometry/arg.pdb -crg_filename ../geometry/jr1.crg -srf_base ../geometry/arg_scaledcharmm_ -srf_num 16 -ksp_rtol 1.0e-10'},
                            {'num': 'arg_10', 'numProcs': 1, 'args': '-is_sphere 0 -pdb_filename ../geometry/arg.pdb -crg_filename ../geometry/jr1.crg -srf_base ../geometry/arg_scaledcharmm_ -srf_num 20 -ksp_rtol 1.0e-10'},
                            {'num': 'arg_11', 'numProcs': 1, 'args': '-is_sphere 0 -pdb_filename ../geometry/arg.pdb -crg_filename ../geometry/jr1.crg -srf_base ../geometry/arg_scaledcharmm_ -srf_num 24 -ksp_rtol 1.0e-10'},
                            {'num': 'arg_12', 'numProcs': 1, 'args': '-is_sphere 0 -pdb_filename ../geometry/arg.pdb -crg_filename ../geometry/jr1.crg -srf_base ../geometry/arg_scaledcharmm_ -srf_num 28 -ksp_rtol 1.0e-10'},
                            {'num': 'arg_13', 'numProcs': 1, 'args': '-is_sphere 0 -pdb_filename ../geometry/arg.pdb -crg_filename ../geometry/jr1.crg -srf_base ../geometry/arg_scaledcharmm_ -srf_num 32 -ksp_rtol 1.0e-10'},
                        ],
}
