meshmaker:
	-@mkdir bin
	cd src/meshmaker; ${MAKE}
	cp src/meshmaker/meshmaker bin

srf_setup:
	-@mkdir geometry
	rm -f *.srf geometry/surf*

srf_sphere:
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker sphere_R6.pdb radii.siz charmm.crg sphere_R6_vdens0125.srf sphere.crg 1.4 2.0 0.125 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker sphere_R6.pdb radii.siz charmm.crg sphere_R6_vdens025.srf sphere.crg 1.4 2.0 0.25 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker sphere_R6.pdb radii.siz charmm.crg sphere_R6_vdens05.srf sphere.crg 1.4 2.0 0.5 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker sphere_R6.pdb radii.siz charmm.crg sphere_R6_vdens1.srf sphere.crg 1.4 2.0 1 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker sphere_R6.pdb radii.siz charmm.crg sphere_R6_vdens2.srf sphere.crg 1.4 2.0 2 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker sphere_R6.pdb radii.siz charmm.crg sphere_R6_vdens3.srf sphere.crg 1.4 2.0 3 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker sphere_R6.pdb radii.siz charmm.crg sphere_R6_vdens4.srf sphere.crg 1.4 2.0 4 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker sphere_R6.pdb radii.siz charmm.crg sphere_R6_vdens5.srf sphere.crg 1.4 2.0 5 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker sphere_R6.pdb radii.siz charmm.crg sphere_R6_vdens6.srf sphere.crg 1.4 2.0 6 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker sphere_R6.pdb radii.siz charmm.crg sphere_R6_vdens7.srf sphere.crg 1.4 2.0 7 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker sphere_R6.pdb radii.siz charmm.crg sphere_R6_vdens8.srf sphere.crg 1.4 2.0 8 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker sphere_R6.pdb radii.siz charmm.crg sphere_R6_vdens10.srf sphere.crg 1.4 2.0 10 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker sphere_R6.pdb radii.siz charmm.crg sphere_R6_vdens20.srf sphere.crg 1.4 2.0 20 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker sphere_R6.pdb radii.siz charmm.crg sphere_R6_vdens30.srf sphere.crg 1.4 2.0 30 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker sphere_R6.pdb radii.siz charmm.crg sphere_R6_vdens40.srf sphere.crg 1.4 2.0 40 1 1 0 .

srf_sphere_noprobe:
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker sphere_R6.pdb radii.siz charmm.crg sphere_R6_vdens_noprobe_0125.srf sphere.crg 0 2.0 0.125 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker sphere_R6.pdb radii.siz charmm.crg sphere_R6_vdens_noprobe_025.srf sphere.crg 0 2.0 0.25 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker sphere_R6.pdb radii.siz charmm.crg sphere_R6_vdens_noprobe_05.srf sphere.crg 0 2.0 0.5 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker sphere_R6.pdb radii.siz charmm.crg sphere_R6_vdens_noprobe_1.srf sphere.crg 0 2.0 1 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker sphere_R6.pdb radii.siz charmm.crg sphere_R6_vdens_noprobe_2.srf sphere.crg 0 2.0 2 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker sphere_R6.pdb radii.siz charmm.crg sphere_R6_vdens_noprobe_3.srf sphere.crg 0 2.0 3 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker sphere_R6.pdb radii.siz charmm.crg sphere_R6_vdens_noprobe_4.srf sphere.crg 0 2.0 4 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker sphere_R6.pdb radii.siz charmm.crg sphere_R6_vdens_noprobe_5.srf sphere.crg 0 2.0 5 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker sphere_R6.pdb radii.siz charmm.crg sphere_R6_vdens_noprobe_6.srf sphere.crg 0 2.0 6 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker sphere_R6.pdb radii.siz charmm.crg sphere_R6_vdens_noprobe_7.srf sphere.crg 0 2.0 7 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker sphere_R6.pdb radii.siz charmm.crg sphere_R6_vdens_noprobe_8.srf sphere.crg 0 2.0 8 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker sphere_R6.pdb radii.siz charmm.crg sphere_R6_vdens_noprobe_10.srf sphere.crg 0 2.0 10 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker sphere_R6.pdb radii.siz charmm.crg sphere_R6_vdens_noprobe_20.srf sphere.crg 0 2.0 20 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker sphere_R6.pdb radii.siz charmm.crg sphere_R6_vdens_noprobe_30.srf sphere.crg 0 2.0 30 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker sphere_R6.pdb radii.siz charmm.crg sphere_R6_vdens_noprobe_40.srf sphere.crg 0 2.0 40 1 1 0 .

# https://en.wikipedia.org/wiki/Amino_acid#Table_of_standard_amino_acid_abbreviations_and_properties
srf_asp:
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker asp.pdb scaled_charmm22.siz charmm.crg asp_scaledcharmm_1.srf asp.crg 1.4 2.0 1 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker asp.pdb scaled_charmm22.siz charmm.crg asp_scaledcharmm_2.srf asp.crg 1.4 2.0 2 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker asp.pdb scaled_charmm22.siz charmm.crg asp_scaledcharmm_3.srf asp.crg 1.4 2.0 3 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker asp.pdb scaled_charmm22.siz charmm.crg asp_scaledcharmm_4.srf asp.crg 1.4 2.0 4 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker asp.pdb scaled_charmm22.siz charmm.crg asp_scaledcharmm_5.srf asp.crg 1.4 2.0 5 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker asp.pdb scaled_charmm22.siz charmm.crg asp_scaledcharmm_6.srf asp.crg 1.4 2.0 6 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker asp.pdb scaled_charmm22.siz charmm.crg asp_scaledcharmm_7.srf asp.crg 1.4 2.0 7 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker asp.pdb scaled_charmm22.siz charmm.crg asp_scaledcharmm_8.srf asp.crg 1.4 2.0 8 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker asp.pdb scaled_charmm22.siz charmm.crg asp_scaledcharmm_12.srf asp.crg 1.4 2.0 12 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker asp.pdb scaled_charmm22.siz charmm.crg asp_scaledcharmm_16.srf asp.crg 1.4 2.0 16 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker asp.pdb scaled_charmm22.siz charmm.crg asp_scaledcharmm_20.srf asp.crg 1.4 2.0 20 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker asp.pdb scaled_charmm22.siz charmm.crg asp_scaledcharmm_24.srf asp.crg 1.4 2.0 24 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker asp.pdb scaled_charmm22.siz charmm.crg asp_scaledcharmm_28.srf asp.crg 1.4 2.0 28 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker asp.pdb scaled_charmm22.siz charmm.crg asp_scaledcharmm_32.srf asp.crg 1.4 2.0 32 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker asp.pdb scaled_charmm22.siz charmm.crg asp_scaledcharmm_50.srf asp.crg 1.4 2.0 50 1 1 0 .

srf_asp_noprobe:
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker asp.pdb scaled_charmm22.siz charmm.crg asp_scaledcharmm_noprobe_1.srf asp.crg 0 0.0 1 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker asp.pdb scaled_charmm22.siz charmm.crg asp_scaledcharmm_noprobe_2.srf asp.crg 0 0.0 2 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker asp.pdb scaled_charmm22.siz charmm.crg asp_scaledcharmm_noprobe_3.srf asp.crg 0 0.0 3 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker asp.pdb scaled_charmm22.siz charmm.crg asp_scaledcharmm_noprobe_4.srf asp.crg 0 0.0 4 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker asp.pdb scaled_charmm22.siz charmm.crg asp_scaledcharmm_noprobe_5.srf asp.crg 0 0.0 5 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker asp.pdb scaled_charmm22.siz charmm.crg asp_scaledcharmm_noprobe_6.srf asp.crg 0 0.0 6 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker asp.pdb scaled_charmm22.siz charmm.crg asp_scaledcharmm_noprobe_7.srf asp.crg 0 0.0 7 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker asp.pdb scaled_charmm22.siz charmm.crg asp_scaledcharmm_noprobe_8.srf asp.crg 0 0.0 8 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker asp.pdb scaled_charmm22.siz charmm.crg asp_scaledcharmm_noprobe_12.srf asp.crg 0 0.0 12 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker asp.pdb scaled_charmm22.siz charmm.crg asp_scaledcharmm_noprobe_16.srf asp.crg 0 0.0 16 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker asp.pdb scaled_charmm22.siz charmm.crg asp_scaledcharmm_noprobe_20.srf asp.crg 0 0.0 20 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker asp.pdb scaled_charmm22.siz charmm.crg asp_scaledcharmm_noprobe_24.srf asp.crg 0 0.0 24 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker asp.pdb scaled_charmm22.siz charmm.crg asp_scaledcharmm_noprobe_28.srf asp.crg 0 0.0 28 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker asp.pdb scaled_charmm22.siz charmm.crg asp_scaledcharmm_noprobe_32.srf asp.crg 0 0.0 32 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker asp.pdb scaled_charmm22.siz charmm.crg asp_scaledcharmm_noprobe_50.srf asp.crg 0 0.0 50 1 1 0 .
srf_arg:
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker arg.pdb scaled_charmm22.siz charmm.crg arg_scaledcharmm_1.srf arg.crg 1.4 2.0 1 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker arg.pdb scaled_charmm22.siz charmm.crg arg_scaledcharmm_2.srf arg.crg 1.4 2.0 2 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker arg.pdb scaled_charmm22.siz charmm.crg arg_scaledcharmm_3.srf arg.crg 1.4 2.0 3 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker arg.pdb scaled_charmm22.siz charmm.crg arg_scaledcharmm_4.srf arg.crg 1.4 2.0 4 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker arg.pdb scaled_charmm22.siz charmm.crg arg_scaledcharmm_5.srf arg.crg 1.4 2.0 5 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker arg.pdb scaled_charmm22.siz charmm.crg arg_scaledcharmm_6.srf arg.crg 1.4 2.0 6 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker arg.pdb scaled_charmm22.siz charmm.crg arg_scaledcharmm_7.srf arg.crg 1.4 2.0 7 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker arg.pdb scaled_charmm22.siz charmm.crg arg_scaledcharmm_8.srf arg.crg 1.4 2.0 8 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker arg.pdb scaled_charmm22.siz charmm.crg arg_scaledcharmm_12.srf arg.crg 1.4 2.0 12 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker arg.pdb scaled_charmm22.siz charmm.crg arg_scaledcharmm_16.srf arg.crg 1.4 2.0 16 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker arg.pdb scaled_charmm22.siz charmm.crg arg_scaledcharmm_20.srf arg.crg 1.4 2.0 20 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker arg.pdb scaled_charmm22.siz charmm.crg arg_scaledcharmm_24.srf arg.crg 1.4 2.0 24 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker arg.pdb scaled_charmm22.siz charmm.crg arg_scaledcharmm_28.srf arg.crg 1.4 2.0 28 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker arg.pdb scaled_charmm22.siz charmm.crg arg_scaledcharmm_32.srf arg.crg 1.4 2.0 32 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker arg.pdb scaled_charmm22.siz charmm.crg arg_scaledcharmm_50.srf arg.crg 1.4 2.0 50 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker arg.pdb scaled_charmm22.siz charmm.crg arg_scaledcharmm_100.srf arg.crg 1.4 2.0 100 1 1 0 .

srf_arg_noprobe:
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker arg.pdb scaled_charmm22.siz charmm.crg arg_scaledcharmm_noprobe_1.srf arg.crg 0 2.0 1 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker arg.pdb scaled_charmm22.siz charmm.crg arg_scaledcharmm_noprobe_2.srf arg.crg 0 2.0 2 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker arg.pdb scaled_charmm22.siz charmm.crg arg_scaledcharmm_noprobe_3.srf arg.crg 0 2.0 3 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker arg.pdb scaled_charmm22.siz charmm.crg arg_scaledcharmm_noprobe_4.srf arg.crg 0 2.0 4 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker arg.pdb scaled_charmm22.siz charmm.crg arg_scaledcharmm_noprobe_5.srf arg.crg 0 2.0 5 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker arg.pdb scaled_charmm22.siz charmm.crg arg_scaledcharmm_noprobe_6.srf arg.crg 0 2.0 6 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker arg.pdb scaled_charmm22.siz charmm.crg arg_scaledcharmm_noprobe_7.srf arg.crg 0 2.0 7 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker arg.pdb scaled_charmm22.siz charmm.crg arg_scaledcharmm_noprobe_8.srf arg.crg 0 2.0 8 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker arg.pdb scaled_charmm22.siz charmm.crg arg_scaledcharmm_noprobe_12.srf arg.crg 0 2.0 12 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker arg.pdb scaled_charmm22.siz charmm.crg arg_scaledcharmm_noprobe_16.srf arg.crg 0 2.0 16 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker arg.pdb scaled_charmm22.siz charmm.crg arg_scaledcharmm_noprobe_20.srf arg.crg 0 2.0 20 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker arg.pdb scaled_charmm22.siz charmm.crg arg_scaledcharmm_noprobe_24.srf arg.crg 0 2.0 24 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker arg.pdb scaled_charmm22.siz charmm.crg arg_scaledcharmm_noprobe_28.srf arg.crg 0 2.0 28 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker arg.pdb scaled_charmm22.siz charmm.crg arg_scaledcharmm_noprobe_32.srf arg.crg 0 2.0 32 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker arg.pdb scaled_charmm22.siz charmm.crg arg_scaledcharmm_noprobe_50.srf arg.crg 0 2.0 50 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker arg.pdb scaled_charmm22.siz charmm.crg arg_scaledcharmm_noprobe_100.srf arg.crg 0 2.0 100 1 1 0 .

srf: meshmaker srf_setup srf_sphere srf_arg srf_asp srf_sphere_noprobe srf_arg_noprobe srf_asp_noprobe
