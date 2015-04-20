CFLAGS = -I.

testSrfOnSurfacePoints: testSrfOnSurfacePoints.o surface.o
	${CLINKER} -o $@ $^ ${PETSC_LIB}
	${DSYMUTIL} $@

regression:
	${PETSC_DIR}/config/builder2.py check [testSrfOnSurfacePoints.c,surface.c] --regParams=testParams ${ARGS}

figures:
	./meshConvergence.py ./meshConvSphere.png
	./workPrecision.py ./workPrecSphere.png

msp:
	../../msp/bin/msroll -m ../../jay-asymmetry/saltresidues/asp/asp.pdb -q asp.pqms -j asp.pnt
	../../msp/bin/msdraw -i asp.script -r asp.sgi

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules
