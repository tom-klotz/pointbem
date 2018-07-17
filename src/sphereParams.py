import os
petscToPointbem = os.path.relpath(os.path.abspath('..'), os.environ['PETSC_DIR'])
print petscToPointbem
regressionParameters = {os.environ.get("PETSC_TO_POINTBEM", petscToPointbem)+'/src/testSrfOnSurfacePoints':
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
                        ],
}
                        
