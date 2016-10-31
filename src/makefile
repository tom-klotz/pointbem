CFLAGS = -I. -std=c99

testSrfOnSurfacePoints: testSrfOnSurfacePoints.o surface.o
	${CLINKER} -o $@ $^ ${PETSC_LIB} ${SLEPC_SYS_LIB}
	${DSYMUTIL} $@

regression:
	${PETSC_DIR}/config/builder2.py check [testSrfOnSurfacePoints.c,surface.c] --regParams=testParams ${ARGS}

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
