CFLAGS = -I.

testSrfOnSurfacePoints: testSrfOnSurfacePoints.o surface.o sphere.o molecule.o BEM.o
	${CLINKER} -o $@ $^ ${PETSC_LIB}
	${DSYMUTIL} $@

calcSolvationEnergy: calcSolvationEnergy.o BEM.o surface.o molecule.o sphere.o
	${CLINKER} -o $@ $^ ${PETSC_LIB}
	${DSYMUTIL} $@

calcEigs: calcEigs.o
	${CLINKER} -o $@ $^ ${PETSC_LIB} ${SLEPC_EPS_LIB}
	${DSYMUTIL} $@

ex7: ex7.o chkopts
	-${CLINKER} -o ex7 ex7.o ${SLEPC_EPS_LIB}
	${RM} ex7.o

spheremat:
	./calcSolvationEnergy -srf_base ../geometry/sphere_R6_vdens0 -num_charges 10 -nmax 25 -srf_num 0125 -density 0.125 -ksp_rtol 1.0e-10 -snes_linesearch_type basic -snes_type ksponly -potential_petscspace_degree 1 -alpha 0.5 -beta -60 -gamma -0.5 -k_view output/Ksphere0125.mat
	./calcSolvationEnergy -srf_base ../geometry/sphere_R6_vdens0 -num_charges 10 -nmax 25 -srf_num 025 -density 0.25 -ksp_rtol 1.0e-10 -snes_linesearch_type basic -snes_type ksponly -potential_petscspace_degree 1 -alpha 0.5 -beta -60 -gamma -0.5 -k_view output/Ksphere025.mat
	./calcSolvationEnergy -srf_base ../geometry/sphere_R6_vdens0 -num_charges 10 -nmax 25 -srf_num 05 -density 0.5 -ksp_rtol 1.0e-10 -snes_linesearch_type basic -snes_type ksponly -potential_petscspace_degree 1 -alpha 0.5 -beta -60 -gamma -0.5 -k_view output/Ksphere05.mat
	./calcSolvationEnergy -srf_base ../geometry/sphere_R6_vdens -num_charges 10 -nmax 25 -srf_num 1 -density 1.0 -ksp_rtol 1.0e-10 -snes_linesearch_type basic -snes_type ksponly -potential_petscspace_degree 1 -alpha 0.5 -beta -60 -gamma -0.5 -k_view output/Ksphere1.mat
	./calcSolvationEnergy -srf_base ../geometry/sphere_R6_vdens -num_charges 10 -nmax 25 -srf_num 2 -density 2.0 -ksp_rtol 1.0e-10 -snes_linesearch_type basic -snes_type ksponly -potential_petscspace_degree 1 -alpha 0.5 -beta -60 -gamma -0.5 -k_view output/Ksphere2.mat
	./calcSolvationEnergy -srf_base ../geometry/sphere_R6_vdens -num_charges 10 -nmax 25 -srf_num 3 -density 3.0 -ksp_rtol 1.0e-10 -snes_linesearch_type basic -snes_type ksponly -potential_petscspace_degree 1 -alpha 0.5 -beta -60 -gamma -0.5 -k_view output/Ksphere3.mat
	./calcSolvationEnergy -srf_base ../geometry/sphere_R6_vdens -num_charges 10 -nmax 25 -srf_num 4 -density 4.0 -ksp_rtol 1.0e-10 -snes_linesearch_type basic -snes_type ksponly -potential_petscspace_degree 1 -alpha 0.5 -beta -60 -gamma -0.5 -k_view output/Ksphere4.mat
	./calcSolvationEnergy -srf_base ../geometry/sphere_R6_vdens -num_charges 10 -nmax 25 -srf_num 5 -density 5.0 -ksp_rtol 1.0e-10 -snes_linesearch_type basic -snes_type ksponly -potential_petscspace_degree 1 -alpha 0.5 -beta -60 -gamma -0.5 -k_view output/Ksphere5.mat
	./calcSolvationEnergy -srf_base ../geometry/sphere_R6_vdens -num_charges 10 -nmax 25 -srf_num 6 -density 6.0 -ksp_rtol 1.0e-10 -snes_linesearch_type basic -snes_type ksponly -potential_petscspace_degree 1 -alpha 0.5 -beta -60 -gamma -0.5 -k_view output/Ksphere6.mat
	./calcSolvationEnergy -srf_base ../geometry/sphere_R6_vdens -num_charges 10 -nmax 25 -srf_num 8 -density 8.0 -ksp_rtol 1.0e-10 -snes_linesearch_type basic -snes_type ksponly -potential_petscspace_degree 1 -alpha 0.5 -beta -60 -gamma -0.5 -k_view output/Ksphere8.mat
	./calcSolvationEnergy -srf_base ../geometry/sphere_R6_vdens -num_charges 10 -nmax 25 -srf_num 10 -density 10.0 -ksp_rtol 1.0e-10 -snes_linesearch_type basic -snes_type ksponly -potential_petscspace_degree 1 -alpha 0.5 -beta -60 -gamma -0.5 -k_view output/Ksphere10.mat

argmat:
	./calcSolvationEnergy -is_sphere 0 -pdb_filename ../geometry/arg.pdb -crg_filename ../geometry/jr1.crg -srf_base ../geometry/arg_scaledcharmm_ -srf_num 1 -ksp_rtol 1.0e-10 -potential_petscspace_degree 1 -k_view output/Karg1.mat
	./calcSolvationEnergy -is_sphere 0 -pdb_filename ../geometry/arg.pdb -crg_filename ../geometry/jr1.crg -srf_base ../geometry/arg_scaledcharmm_ -srf_num 2 -ksp_rtol 1.0e-10 -potential_petscspace_degree 1 -k_view output/Karg2.mat
	./calcSolvationEnergy -is_sphere 0 -pdb_filename ../geometry/arg.pdb -crg_filename ../geometry/jr1.crg -srf_base ../geometry/arg_scaledcharmm_ -srf_num 3 -ksp_rtol 1.0e-10 -potential_petscspace_degree 1 -k_view output/Karg3.mat
	./calcSolvationEnergy -is_sphere 0 -pdb_filename ../geometry/arg.pdb -crg_filename ../geometry/jr1.crg -srf_base ../geometry/arg_scaledcharmm_ -srf_num 4 -ksp_rtol 1.0e-10 -potential_petscspace_degree 1 -k_view output/Karg4.mat
	./calcSolvationEnergy -is_sphere 0 -pdb_filename ../geometry/arg.pdb -crg_filename ../geometry/jr1.crg -srf_base ../geometry/arg_scaledcharmm_ -srf_num 5 -ksp_rtol 1.0e-10 -potential_petscspace_degree 1 -k_view output/Karg5.mat
	./calcSolvationEnergy -is_sphere 0 -pdb_filename ../geometry/arg.pdb -crg_filename ../geometry/jr1.crg -srf_base ../geometry/arg_scaledcharmm_ -srf_num 6 -ksp_rtol 1.0e-10 -potential_petscspace_degree 1 -k_view output/Karg6.mat
	./calcSolvationEnergy -is_sphere 0 -pdb_filename ../geometry/arg.pdb -crg_filename ../geometry/jr1.crg -srf_base ../geometry/arg_scaledcharmm_ -srf_num 7 -ksp_rtol 1.0e-10 -potential_petscspace_degree 1 -k_view output/Karg7.mat
	./calcSolvationEnergy -is_sphere 0 -pdb_filename ../geometry/arg.pdb -crg_filename ../geometry/jr1.crg -srf_base ../geometry/arg_scaledcharmm_ -srf_num 8 -ksp_rtol 1.0e-10 -potential_petscspace_degree 1 -k_view output/Karg8.mat
	./calcSolvationEnergy -is_sphere 0 -pdb_filename ../geometry/arg.pdb -crg_filename ../geometry/jr1.crg -srf_base ../geometry/arg_scaledcharmm_ -srf_num 12 -ksp_rtol 1.0e-10 -potential_petscspace_degree 1 -k_view output/Karg12.mat
	./calcSolvationEnergy -is_sphere 0 -pdb_filename ../geometry/arg.pdb -crg_filename ../geometry/jr1.crg -srf_base ../geometry/arg_scaledcharmm_ -srf_num 16 -ksp_rtol 1.0e-10 -potential_petscspace_degree 1 -k_view output/Karg16.mat
	./calcSolvationEnergy -is_sphere 0 -pdb_filename ../geometry/arg.pdb -crg_filename ../geometry/jr1.crg -srf_base ../geometry/arg_scaledcharmm_ -srf_num 20 -ksp_rtol 1.0e-10 -potential_petscspace_degree 1 -k_view output/Karg20.mat
	./calcSolvationEnergy -is_sphere 0 -pdb_filename ../geometry/arg.pdb -crg_filename ../geometry/jr1.crg -srf_base ../geometry/arg_scaledcharmm_ -srf_num 24 -ksp_rtol 1.0e-10 -potential_petscspace_degree 1 -k_view output/Karg24.mat
	./calcSolvationEnergy -is_sphere 0 -pdb_filename ../geometry/arg.pdb -crg_filename ../geometry/jr1.crg -srf_base ../geometry/arg_scaledcharmm_ -srf_num 28 -ksp_rtol 1.0e-10 -potential_petscspace_degree 1 -k_view output/Karg28.mat
	./calcSolvationEnergy -is_sphere 0 -pdb_filename ../geometry/arg.pdb -crg_filename ../geometry/jr1.crg -srf_base ../geometry/arg_scaledcharmm_ -srf_num 32 -ksp_rtol 1.0e-10 -potential_petscspace_degree 1 -k_view output/Karg32.mat

aspmat:
	./calcSolvationEnergy -is_sphere 0 -pdb_filename ../geometry/asp.pdb -crg_filename ../geometry/jd1.crg -srf_base ../geometry/asp_scaledcharmm_ -srf_num 1 -potential_petscspace_degree 1 -k_view output/Kasp1.mat
	./calcSolvationEnergy -is_sphere 0 -pdb_filename ../geometry/asp.pdb -crg_filename ../geometry/jd1.crg -srf_base ../geometry/asp_scaledcharmm_ -srf_num 2 -potential_petscspace_degree 1 -k_view output/Kasp2.mat
	./calcSolvationEnergy -is_sphere 0 -pdb_filename ../geometry/asp.pdb -crg_filename ../geometry/jd1.crg -srf_base ../geometry/asp_scaledcharmm_ -srf_num 3 -potential_petscspace_degree 1 -k_view output/Kasp3.mat
	./calcSolvationEnergy -is_sphere 0 -pdb_filename ../geometry/asp.pdb -crg_filename ../geometry/jd1.crg -srf_base ../geometry/asp_scaledcharmm_ -srf_num 4 -potential_petscspace_degree 1 -k_view output/Kasp4.mat
	./calcSolvationEnergy -is_sphere 0 -pdb_filename ../geometry/asp.pdb -crg_filename ../geometry/jd1.crg -srf_base ../geometry/asp_scaledcharmm_ -srf_num 5 -potential_petscspace_degree 1 -k_view output/Kasp5.mat
	./calcSolvationEnergy -is_sphere 0 -pdb_filename ../geometry/asp.pdb -crg_filename ../geometry/jd1.crg -srf_base ../geometry/asp_scaledcharmm_ -srf_num 6 -potential_petscspace_degree 1 -k_view output/Kasp6.mat
	./calcSolvationEnergy -is_sphere 0 -pdb_filename ../geometry/asp.pdb -crg_filename ../geometry/jd1.crg -srf_base ../geometry/asp_scaledcharmm_ -srf_num 7 -potential_petscspace_degree 1 -k_view output/Kasp7.mat
	./calcSolvationEnergy -is_sphere 0 -pdb_filename ../geometry/asp.pdb -crg_filename ../geometry/jd1.crg -srf_base ../geometry/asp_scaledcharmm_ -srf_num 8 -potential_petscspace_degree 1 -k_view output/Kasp8.mat
	./calcSolvationEnergy -is_sphere 0 -pdb_filename ../geometry/asp.pdb -crg_filename ../geometry/jd1.crg -srf_base ../geometry/asp_scaledcharmm_ -srf_num 12 -potential_petscspace_degree 1 -k_view output/Kasp12.mat
	./calcSolvationEnergy -is_sphere 0 -pdb_filename ../geometry/asp.pdb -crg_filename ../geometry/jd1.crg -srf_base ../geometry/asp_scaledcharmm_ -srf_num 16 -potential_petscspace_degree 1 -k_view output/Kasp16.mat
	./calcSolvationEnergy -is_sphere 0 -pdb_filename ../geometry/asp.pdb -crg_filename ../geometry/jd1.crg -srf_base ../geometry/asp_scaledcharmm_ -srf_num 20 -potential_petscspace_degree 1 -k_view output/Kasp20.mat
	./calcSolvationEnergy -is_sphere 0 -pdb_filename ../geometry/asp.pdb -crg_filename ../geometry/jd1.crg -srf_base ../geometry/asp_scaledcharmm_ -srf_num 24 -potential_petscspace_degree 1 -k_view output/Kasp24.mat
	./calcSolvationEnergy -is_sphere 0 -pdb_filename ../geometry/asp.pdb -crg_filename ../geometry/jd1.crg -srf_base ../geometry/asp_scaledcharmm_ -srf_num 28 -potential_petscspace_degree 1 -k_view output/Kasp28.mat
	./calcSolvationEnergy -is_sphere 0 -pdb_filename ../geometry/asp.pdb -crg_filename ../geometry/jd1.crg -srf_base ../geometry/asp_scaledcharmm_ -srf_num 32 -potential_petscspace_degree 1 -k_view output/Kasp32.mat
regression:
	PYTHONPATH=${PYTHONPATH}:.. ${PETSC_DIR}/config/builder2.py check [testSrfOnSurfacePoints.c,surface.c,sphere.c,molecule.c] --regParams=testParams ${ARGS}

test_sphere:
	PYTHONPATH=${PYTHONPATH}:.. ${PETSC_DIR}/config/builder2.py check [testSrfOnSurfacePoints.c,surface.c,sphere.c,molecule.c] --regParams=sphereParams ${ARGS}

test_arg:
	PYTHONPATH=${PYTHONPATH}:.. ${PETSC_DIR}/config/builder2.py check [testSrfOnSurfacePoints.c,surface.c,sphere.c,molecule.c] --regParams=argParams ${ARGS}

test_asp:
	PYTHONPATH=${PYTHONPATH}:.. ${PETSC_DIR}/config/builder2.py check [testSrfOnSurfacePoints.c,surface.c,sphere.c,molecule.c] --regParams=aspParams ${ARGS}

asym_sphere:
	PYTHONPATH=${PYTHONPATH}:.. ${PETSC_DIR}/config/builder2.py check [testSrfOnSurfacePoints.c,surface.c,sphere.c,molecule.c] --regParams=asymParams ${ARGS}

figures_sphere:
	./meshConvergence.py --type=sphere --output=meshConvSphere
	./workPrecision.py --type=sphere --output=workPrecSphere

figures_arg:
	./meshConvergence.py --type=arg --output=meshConvARG
	./workPrecision.py --type=arg --output=workPrecARG

figures_asp:
	./meshConvergence.py --type=asp --output=meshConvASP
	./workPrecision.py --type=asp --output=workPrecASP

figures: figures_sphere figures_arg figures_asp

figures_mv: figures
	mv meshConvSphere_RelErr.png /PETSc3/papers/imece15-asme/figures/
	mv workPrecSphere.png /PETSc3/papers/imece15-asme/figures/
	mv meshConvARG_RelErr.png /PETSc3/papers/imece15-asme/figures/
	mv workPrecARG.png /PETSc3/papers/imece15-asme/figures/
	mv meshConvASP_RelErr.png /PETSc3/papers/imece15-asme/figures/
	mv meshConvASP_SolvEng.png /PETSc3/papers/imece15-asme/figures/
	mv workPrecASP.png /PETSc3/papers/imece15-asme/figures/

#  -y asp.pat -r asp.rad
msp:
	-../../msp/bin/msroll -m ../../jay-asymmetry/saltresidues/asp/asp.pdb -p 0.1 -q asp.pqms -j asp.pnt
	-../../msp/bin/msdraw -i asp.script -r asp.sgi

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules
include ${SLEPC_DIR}/lib/slepc/conf/slepc_common
