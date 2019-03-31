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
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker asp.pdb scaled_charmm22.siz charmm.crg asp_scaledcharmm_1.srf asp.crg 0 2.0 .5 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker asp.pdb scaled_charmm22.siz charmm.crg asp_scaledcharmm_2.srf asp.crg 0 2.0 .6 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker asp.pdb scaled_charmm22.siz charmm.crg asp_scaledcharmm_3.srf asp.crg 0 2.0 .7 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker asp.pdb scaled_charmm22.siz charmm.crg asp_scaledcharmm_4.srf asp.crg 0 2.0 .8 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker asp.pdb scaled_charmm22.siz charmm.crg asp_scaledcharmm_5.srf asp.crg 0 2.0 .9 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker asp.pdb scaled_charmm22.siz charmm.crg asp_scaledcharmm_6.srf asp.crg 0 2.0 1.0 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker asp.pdb scaled_charmm22.siz charmm.crg asp_scaledcharmm_7.srf asp.crg 0 2.0 1.2 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker asp.pdb scaled_charmm22.siz charmm.crg asp_scaledcharmm_8.srf asp.crg 0 2.0 1.4 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker asp.pdb scaled_charmm22.siz charmm.crg asp_scaledcharmm_9.srf asp.crg 0 2.0 1.6 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker asp.pdb scaled_charmm22.siz charmm.crg asp_scaledcharmm_10.srf asp.crg 0 2.0 1.8 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker asp.pdb scaled_charmm22.siz charmm.crg asp_scaledcharmm_11.srf asp.crg 0 2.0 2.0 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker asp.pdb scaled_charmm22.siz charmm.crg asp_scaledcharmm_12.srf asp.crg 0 2.0 2.25 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker asp.pdb scaled_charmm22.siz charmm.crg asp_scaledcharmm_13.srf asp.crg 0 2.0 2.5 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker asp.pdb scaled_charmm22.siz charmm.crg asp_scaledcharmm_14.srf asp.crg 0 2.0 2.75 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker asp.pdb scaled_charmm22.siz charmm.crg asp_scaledcharmm_15.srf asp.crg 0 2.0 3 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker asp.pdb scaled_charmm22.siz charmm.crg asp_scaledcharmm_16.srf asp.crg 0 2.0 3.5 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker asp.pdb scaled_charmm22.siz charmm.crg asp_scaledcharmm_17.srf asp.crg 0 2.0 4 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker asp.pdb scaled_charmm22.siz charmm.crg asp_scaledcharmm_18.srf asp.crg 0 2.0 4.5 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker asp.pdb scaled_charmm22.siz charmm.crg asp_scaledcharmm_19.srf asp.crg 0 2.0 5 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker asp.pdb scaled_charmm22.siz charmm.crg asp_scaledcharmm_20.srf asp.crg 0 2.0 6 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker asp.pdb scaled_charmm22.siz charmm.crg asp_scaledcharmm_21.srf asp.crg 0 2.0 7 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker asp.pdb scaled_charmm22.siz charmm.crg asp_scaledcharmm_22.srf asp.crg 0 2.0 8 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker asp.pdb scaled_charmm22.siz charmm.crg asp_scaledcharmm_23.srf asp.crg 0 2.0 9 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker asp.pdb scaled_charmm22.siz charmm.crg asp_scaledcharmm_24.srf asp.crg 0 2.0 10 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker asp.pdb scaled_charmm22.siz charmm.crg asp_scaledcharmm_25.srf asp.crg 0 2.0 12 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker asp.pdb scaled_charmm22.siz charmm.crg asp_scaledcharmm_26.srf asp.crg 0 2.0 14 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker asp.pdb scaled_charmm22.siz charmm.crg asp_scaledcharmm_27.srf asp.crg 0 2.0 16 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker asp.pdb scaled_charmm22.siz charmm.crg asp_scaledcharmm_28.srf asp.crg 0 2.0 18 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker asp.pdb scaled_charmm22.siz charmm.crg asp_scaledcharmm_29.srf asp.crg 0 2.0 20 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker asp.pdb scaled_charmm22.siz charmm.crg asp_scaledcharmm_30.srf asp.crg 0 2.0 22 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker asp.pdb scaled_charmm22.siz charmm.crg asp_scaledcharmm_31.srf asp.crg 0 2.0 25 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker asp.pdb scaled_charmm22.siz charmm.crg asp_scaledcharmm_32.srf asp.crg 0 2.0 30 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker asp.pdb scaled_charmm22.siz charmm.crg asp_scaledcharmm_33.srf asp.crg 0 2.0 35 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker asp.pdb scaled_charmm22.siz charmm.crg asp_scaledcharmm_34.srf asp.crg 0 2.0 50 1 1 0 .



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
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker arg.pdb scaled_charmm22.siz charmm.crg arg_scaledcharmm_1.srf arg.crg 0 2.0 .35 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker arg.pdb scaled_charmm22.siz charmm.crg arg_scaledcharmm_2.srf arg.crg 0 2.0 .4 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker arg.pdb scaled_charmm22.siz charmm.crg arg_scaledcharmm_3.srf arg.crg 0 2.0 .44 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker arg.pdb scaled_charmm22.siz charmm.crg arg_scaledcharmm_4.srf arg.crg 0 2.0 .5 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker arg.pdb scaled_charmm22.siz charmm.crg arg_scaledcharmm_5.srf arg.crg 0 2.0 .55 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker arg.pdb scaled_charmm22.siz charmm.crg arg_scaledcharmm_6.srf arg.crg 0 2.0 .6 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker arg.pdb scaled_charmm22.siz charmm.crg arg_scaledcharmm_7.srf arg.crg 0 2.0 .7 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker arg.pdb scaled_charmm22.siz charmm.crg arg_scaledcharmm_8.srf arg.crg 0 2.0 .8 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker arg.pdb scaled_charmm22.siz charmm.crg arg_scaledcharmm_9.srf arg.crg 0 2.0 .9 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker arg.pdb scaled_charmm22.siz charmm.crg arg_scaledcharmm_10.srf arg.crg 0 2.0 1 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker arg.pdb scaled_charmm22.siz charmm.crg arg_scaledcharmm_11.srf arg.crg 0 2.0 1.25 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker arg.pdb scaled_charmm22.siz charmm.crg arg_scaledcharmm_12.srf arg.crg 0 2.0 1.5 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker arg.pdb scaled_charmm22.siz charmm.crg arg_scaledcharmm_13.srf arg.crg 0 2.0 1.75 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker arg.pdb scaled_charmm22.siz charmm.crg arg_scaledcharmm_14.srf arg.crg 0 2.0 2 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker arg.pdb scaled_charmm22.siz charmm.crg arg_scaledcharmm_15.srf arg.crg 0 2.0 2.33 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker arg.pdb scaled_charmm22.siz charmm.crg arg_scaledcharmm_16.srf arg.crg 0 2.0 2.66 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker arg.pdb scaled_charmm22.siz charmm.crg arg_scaledcharmm_17.srf arg.crg 0 2.0 3 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker arg.pdb scaled_charmm22.siz charmm.crg arg_scaledcharmm_18.srf arg.crg 0 2.0 3.33 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker arg.pdb scaled_charmm22.siz charmm.crg arg_scaledcharmm_19.srf arg.crg 0 2.0 3.66 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker arg.pdb scaled_charmm22.siz charmm.crg arg_scaledcharmm_20.srf arg.crg 0 2.0 4 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker arg.pdb scaled_charmm22.siz charmm.crg arg_scaledcharmm_21.srf arg.crg 0 2.0 5 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker arg.pdb scaled_charmm22.siz charmm.crg arg_scaledcharmm_22.srf arg.crg 0 2.0 6 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker arg.pdb scaled_charmm22.siz charmm.crg arg_scaledcharmm_23.srf arg.crg 0 2.0 7 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker arg.pdb scaled_charmm22.siz charmm.crg arg_scaledcharmm_24.srf arg.crg 0 2.0 8 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker arg.pdb scaled_charmm22.siz charmm.crg arg_scaledcharmm_25.srf arg.crg 0 2.0 10 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker arg.pdb scaled_charmm22.siz charmm.crg arg_scaledcharmm_26.srf arg.crg 0 2.0 12 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker arg.pdb scaled_charmm22.siz charmm.crg arg_scaledcharmm_27.srf arg.crg 0 2.0 14 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker arg.pdb scaled_charmm22.siz charmm.crg arg_scaledcharmm_28.srf arg.crg 0 2.0 16 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker arg.pdb scaled_charmm22.siz charmm.crg arg_scaledcharmm_29.srf arg.crg 0 2.0 20 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker arg.pdb scaled_charmm22.siz charmm.crg arg_scaledcharmm_30.srf arg.crg 0 2.0 24 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker arg.pdb scaled_charmm22.siz charmm.crg arg_scaledcharmm_31.srf arg.crg 0 2.0 28 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker arg.pdb scaled_charmm22.siz charmm.crg arg_scaledcharmm_32.srf arg.crg 0 2.0 32 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker arg.pdb scaled_charmm22.siz charmm.crg arg_scaledcharmm_33.srf arg.crg 0 2.0 50 1 1 0 .


srf_argpqr:
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker arg.pqr lolwhocares.siz dontneedthis.crg argpqr_scaledcharmm_1.srf argpqr.crg 1.4 2.0 1 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker arg.pqr lolwhocares.siz dontneedthis.crg argpqr_scaledcharmm_2.srf argpqr.crg 1.4 2.0 2 1 1 0 .

srf_6dz0:
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker 6dz0.pdb charmm.siz charmm.crg 6dz0_scaledcharmm_1.srf 6dz0.crg 0 2.0 .5 1 1 0 .

srf_2rh3:
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker 2rh3.pdb charmm.siz charmm.crg 2rh3_scaledcharmm_1.srf 2rh3.crg 0 2.0 .6 1 1 0 .

srf_3nir:
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker 3nir.pqr lolwhocares.siz dontneedthis.crg 3nir_scaledcharmm_1.srf 3nir.crg 0 2.0 .5 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker 3nir.pqr lolwhocares.siz dontneedthis.crg 3nir_scaledcharmm_15.srf 3nir.crg 0 2.0 .75 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker 3nir.pqr lolwhocares.siz dontneedthis.crg 3nir_scaledcharmm_2.srf 3nir.crg 0 2.0 1 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker 3nir.pqr lolwhocares.siz dontneedthis.crg 3nir_scaledcharmm_25.srf 3nir.crg 0 2.0 1.5 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker 3nir.pqr lolwhocares.siz dontneedthis.crg 3nir_scaledcharmm_3.srf 3nir.crg 0 2.0 2 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker 3nir.pqr lolwhocares.siz dontneedthis.crg 3nir_scaledcharmm_35.srf 3nir.crg 0 2.0 2.5 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker 3nir.pqr lolwhocares.siz dontneedthis.crg 3nir_scaledcharmm_4.srf 3nir.crg 0 2.0 3 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker 3nir.pqr lolwhocares.siz dontneedthis.crg 3nir_scaledcharmm_5.srf 3nir.crg 0 2.0 4 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker 3nir.pqr lolwhocares.siz dontneedthis.crg 3nir_scaledcharmm_6.srf 3nir.crg 0 2.0 5 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker 3nir.pqr lolwhocares.siz dontneedthis.crg 3nir_scaledcharmm_7.srf 3nir.crg 0 2.0 6 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker 3nir.pqr lolwhocares.siz dontneedthis.crg 3nir_scaledcharmm_10.srf 3nir.crg 0 2.0 10 1 1 0 .


srf_1hco:
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker 1hco.pqr charmm.siz charmm.crg 1hco_scaledcharmm_1.srf 1hco.crg 1.4 2.0 .5 1 1 0 .	
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker 1hco.pqr charmm.siz charmm.crg 1hco_scaledcharmm_2.srf 1hco.crg 1.4 2.0 1 1 1 0 .
	cd geometry; MSMS_PATH=${MSMS_PATH} ../bin/meshmaker 1hco.pqr charmm.siz charmm.crg 1hco_scaledcharmm_3.srf 1hco.crg 1.4 2.0 2 1 1 0 .


srf: meshmaker srf_setup srf_sphere srf_arg srf_asp srf_sphere_noprobe srf_arg_noprobe srf_asp_noprobe
