CFLAGS = -I.

testSrfOnSurfacePoints: testSrfOnSurfacePoints.o surface.o
	${CLINKER} -o $@ $^ ${PETSC_LIB}
	${DSYMUTIL} $@

regression:
	${PETSC_DIR}/config/builder2.py check [testSrfOnSurfacePoints.c,surface.c] --regParams=testParams ${ARGS}

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules
