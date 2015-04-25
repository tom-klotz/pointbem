MESHMAKER = /PETSc3/biology/fftsvd/meshmaker
CFLAGS    = -I.

meshmaker: QuadratureRule.o Vector.o Vector3D.o Matrix.o VertFace.o Panel.o FlatPanel.o FlatIntegration.o Integration.o GST.o TOR.o GreensFunction.o calcpc_GST.o calcpc_TOR.o meshutils.o meshmaker.o
	${CLINKER} -o $@ $^ ${PETSC_LIB}

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules
