MSMS_PATH = /PETSc3/biology/msms/msms
MESHMAKER = /PETSc3/biology/fftsvd/meshmaker
CFLAGS = -I.

testSrfOnSurfacePoints: testSrfOnSurfacePoints.o surface.o
	${CLINKER} -o $@ $^ ${PETSC_LIB}
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

#  -y asp.pat -r asp.rad
msp:
	-../../msp/bin/msroll -m ../../jay-asymmetry/saltresidues/asp/asp.pdb -p 0.1 -q asp.pqms -j asp.pnt
	-../../msp/bin/msdraw -i asp.script -r asp.sgi

srf_sphere:
	cd ../../jay-pointbem/geometry; MSMS_PATH=${MSMS_PATH} ${MESHMAKER} sphere_R6.pdb radii.siz sphere_R6_vdens05.srf 1.4 2.0 0.5 1 1 0 .
	cd ../../jay-pointbem/geometry; MSMS_PATH=${MSMS_PATH} ${MESHMAKER} sphere_R6.pdb radii.siz sphere_R6_vdens025.srf 1.4 2.0 0.25 1 1 0 .
	cd ../../jay-pointbem/geometry; MSMS_PATH=${MSMS_PATH} ${MESHMAKER} sphere_R6.pdb radii.siz sphere_R6_vdens0125.srf 1.4 2.0 0.125 1 1 0 .

# https://en.wikipedia.org/wiki/Amino_acid#Table_of_standard_amino_acid_abbreviations_and_properties
srf_asp:
	cd ../../jay-asymmetry/saltresidues/asp; MSMS_PATH=${MSMS_PATH} ${MESHMAKER} asp.pdb ../scaled_charmm22.siz asp_scaledcharmm_12.srf 1.4 2.0 12 1 1 0 .
	cd ../../jay-asymmetry/saltresidues/asp; MSMS_PATH=${MSMS_PATH} ${MESHMAKER} asp.pdb ../scaled_charmm22.siz asp_scaledcharmm_16.srf 1.4 2.0 16 1 1 0 .
	cd ../../jay-asymmetry/saltresidues/asp; MSMS_PATH=${MSMS_PATH} ${MESHMAKER} asp.pdb ../scaled_charmm22.siz asp_scaledcharmm_20.srf 1.4 2.0 20 1 1 0 .
	cd ../../jay-asymmetry/saltresidues/asp; MSMS_PATH=${MSMS_PATH} ${MESHMAKER} asp.pdb ../scaled_charmm22.siz asp_scaledcharmm_24.srf 1.4 2.0 24 1 1 0 .
	cd ../../jay-asymmetry/saltresidues/asp; MSMS_PATH=${MSMS_PATH} ${MESHMAKER} asp.pdb ../scaled_charmm22.siz asp_scaledcharmm_28.srf 1.4 2.0 28 1 1 0 .
	cd ../../jay-asymmetry/saltresidues/asp; MSMS_PATH=${MSMS_PATH} ${MESHMAKER} asp.pdb ../scaled_charmm22.siz asp_scaledcharmm_32.srf 1.4 2.0 32 1 1 0 .

srf_arg:
	cd ../../jay-asymmetry/saltresidues/arg; MSMS_PATH=${MSMS_PATH} ${MESHMAKER} arg.pdb ../scaled_charmm22.siz arg_scaledcharmm_12.srf 1.4 2.0 12 1 1 0 .
	cd ../../jay-asymmetry/saltresidues/arg; MSMS_PATH=${MSMS_PATH} ${MESHMAKER} arg.pdb ../scaled_charmm22.siz arg_scaledcharmm_16.srf 1.4 2.0 16 1 1 0 .
	cd ../../jay-asymmetry/saltresidues/arg; MSMS_PATH=${MSMS_PATH} ${MESHMAKER} arg.pdb ../scaled_charmm22.siz arg_scaledcharmm_20.srf 1.4 2.0 20 1 1 0 .
	cd ../../jay-asymmetry/saltresidues/arg; MSMS_PATH=${MSMS_PATH} ${MESHMAKER} arg.pdb ../scaled_charmm22.siz arg_scaledcharmm_24.srf 1.4 2.0 24 1 1 0 .
	cd ../../jay-asymmetry/saltresidues/arg; MSMS_PATH=${MSMS_PATH} ${MESHMAKER} arg.pdb ../scaled_charmm22.siz arg_scaledcharmm_28.srf 1.4 2.0 28 1 1 0 .
	cd ../../jay-asymmetry/saltresidues/arg; MSMS_PATH=${MSMS_PATH} ${MESHMAKER} arg.pdb ../scaled_charmm22.siz arg_scaledcharmm_32.srf 1.4 2.0 32 1 1 0 .

srf: srf_sphere srf_arg srf_asp

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules
