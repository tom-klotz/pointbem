import os
petscToPointbem = os.path.relpath(os.path.abspath('..'), os.environ['PETSC_DIR'])
print petscToPointbem
regressionParameters = {os.environ.get("PETSC_TO_POINTBEM", petscToPointbem)+'/src/testSrfOnSurfacePoints':
                        [# Aspartic acid tests [asp_0-asp_13]
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
