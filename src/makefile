CFLAGS = -I.

testSrfOnSurfacePoints: testSrfOnSurfacePoints.o
	${CLINKER} -o $@ $< ${PETSC_LIB}
	${DSYMUTIL} $@

include ${PETSC_DIR}/lib/petsc-conf/variables
include ${PETSC_DIR}/lib/petsc-conf/rules
