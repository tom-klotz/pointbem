MSMS_PATH = /PETSc3/biology/msms/msms

meshmaker:
	-@mkdir bin
	cd src/meshmaker; ${MAKE}
	cp src/meshmaker/meshmaker bin

srf_setup:
	-@mkdir geometry
	rm -f *.srf geometry/surf*

srf_sphere:
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker sphere_R6.pdb radii.siz sphere_R6_vdens0125.srf 1.4 2.0 0.125 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker sphere_R6.pdb radii.siz sphere_R6_vdens025.srf 1.4 2.0 0.25 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker sphere_R6.pdb radii.siz sphere_R6_vdens05.srf 1.4 2.0 0.5 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker sphere_R6.pdb radii.siz sphere_R6_vdens1.srf 1.4 2.0 1 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker sphere_R6.pdb radii.siz sphere_R6_vdens2.srf 1.4 2.0 2 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker sphere_R6.pdb radii.siz sphere_R6_vdens3.srf 1.4 2.0 3 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker sphere_R6.pdb radii.siz sphere_R6_vdens4.srf 1.4 2.0 4 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker sphere_R6.pdb radii.siz sphere_R6_vdens5.srf 1.4 2.0 5 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker sphere_R6.pdb radii.siz sphere_R6_vdens6.srf 1.4 2.0 6 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker sphere_R6.pdb radii.siz sphere_R6_vdens7.srf 1.4 2.0 7 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker sphere_R6.pdb radii.siz sphere_R6_vdens8.srf 1.4 2.0 8 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker sphere_R6.pdb radii.siz sphere_R6_vdens10.srf 1.4 2.0 10 1 1 0 .

# https://en.wikipedia.org/wiki/Amino_acid#Table_of_standard_amino_acid_abbreviations_and_properties
srf_asp:
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker asp.pdb scaled_charmm22.siz asp_scaledcharmm_1.srf 1.4 2.0 1 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker asp.pdb scaled_charmm22.siz asp_scaledcharmm_2.srf 1.4 2.0 2 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker asp.pdb scaled_charmm22.siz asp_scaledcharmm_3.srf 1.4 2.0 3 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker asp.pdb scaled_charmm22.siz asp_scaledcharmm_4.srf 1.4 2.0 4 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker asp.pdb scaled_charmm22.siz asp_scaledcharmm_5.srf 1.4 2.0 5 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker asp.pdb scaled_charmm22.siz asp_scaledcharmm_6.srf 1.4 2.0 6 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker asp.pdb scaled_charmm22.siz asp_scaledcharmm_7.srf 1.4 2.0 7 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker asp.pdb scaled_charmm22.siz asp_scaledcharmm_8.srf 1.4 2.0 8 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker asp.pdb scaled_charmm22.siz asp_scaledcharmm_12.srf 1.4 2.0 12 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker asp.pdb scaled_charmm22.siz asp_scaledcharmm_16.srf 1.4 2.0 16 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker asp.pdb scaled_charmm22.siz asp_scaledcharmm_20.srf 1.4 2.0 20 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker asp.pdb scaled_charmm22.siz asp_scaledcharmm_24.srf 1.4 2.0 24 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker asp.pdb scaled_charmm22.siz asp_scaledcharmm_28.srf 1.4 2.0 28 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker asp.pdb scaled_charmm22.siz asp_scaledcharmm_32.srf 1.4 2.0 32 1 1 0 .

srf_arg:
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker arg.pdb scaled_charmm22.siz arg_scaledcharmm_1.srf 1.4 2.0 1 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker arg.pdb scaled_charmm22.siz arg_scaledcharmm_2.srf 1.4 2.0 2 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker arg.pdb scaled_charmm22.siz arg_scaledcharmm_3.srf 1.4 2.0 3 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker arg.pdb scaled_charmm22.siz arg_scaledcharmm_4.srf 1.4 2.0 4 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker arg.pdb scaled_charmm22.siz arg_scaledcharmm_5.srf 1.4 2.0 5 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker arg.pdb scaled_charmm22.siz arg_scaledcharmm_6.srf 1.4 2.0 6 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker arg.pdb scaled_charmm22.siz arg_scaledcharmm_7.srf 1.4 2.0 7 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker arg.pdb scaled_charmm22.siz arg_scaledcharmm_8.srf 1.4 2.0 8 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker arg.pdb scaled_charmm22.siz arg_scaledcharmm_12.srf 1.4 2.0 12 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker arg.pdb scaled_charmm22.siz arg_scaledcharmm_16.srf 1.4 2.0 16 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker arg.pdb scaled_charmm22.siz arg_scaledcharmm_20.srf 1.4 2.0 20 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker arg.pdb scaled_charmm22.siz arg_scaledcharmm_24.srf 1.4 2.0 24 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker arg.pdb scaled_charmm22.siz arg_scaledcharmm_28.srf 1.4 2.0 28 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker arg.pdb scaled_charmm22.siz arg_scaledcharmm_32.srf 1.4 2.0 32 1 1 0 .

srf: meshmaker srf_setup srf_sphere srf_arg srf_asp
