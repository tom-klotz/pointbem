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

1hco: 1hco_point_picard 1hco_point_anderson 1hco_panel_picard 1hco_panel_anderson


spheremat_point:
	./calcSolvationEnergy -srf_base ../geometry/sphere_R6_vdens0 -num_charges 10 -nmax 25 -srf_num 0125 -density 0.125 -ksp_rtol 1.0e-10 -snes_linesearch_type basic -snes_type ksponly -potential_petscspace_degree 1 -k_view output/K_sphere_point_0125.mat
	./calcSolvationEnergy -srf_base ../geometry/sphere_R6_vdens0 -num_charges 10 -nmax 25 -srf_num 025 -density 0.25 -ksp_rtol 1.0e-10 -snes_linesearch_type basic -snes_type ksponly -potential_petscspace_degree 1 -k_view output/K_sphere_point_025.mat
	./calcSolvationEnergy -srf_base ../geometry/sphere_R6_vdens0 -num_charges 10 -nmax 25 -srf_num 05 -density 0.5 -ksp_rtol 1.0e-10 -snes_linesearch_type basic -snes_type ksponly -potential_petscspace_degree 1 -k_view output/K_sphere_point_05.mat
	./calcSolvationEnergy -srf_base ../geometry/sphere_R6_vdens -num_charges 10 -nmax 25 -srf_num 1 -density 1.0 -ksp_rtol 1.0e-10 -snes_linesearch_type basic -snes_type ksponly -potential_petscspace_degree 1 -k_view output/K_sphere_point_1.mat
	./calcSolvationEnergy -srf_base ../geometry/sphere_R6_vdens -num_charges 10 -nmax 25 -srf_num 2 -density 2.0 -ksp_rtol 1.0e-10 -snes_linesearch_type basic -snes_type ksponly -potential_petscspace_degree 1 -k_view output/K_sphere_point_2.mat
	./calcSolvationEnergy -srf_base ../geometry/sphere_R6_vdens -num_charges 10 -nmax 25 -srf_num 3 -density 3.0 -ksp_rtol 1.0e-10 -snes_linesearch_type basic -snes_type ksponly -potential_petscspace_degree 1 -k_view output/K_sphere_point_3.mat
	./calcSolvationEnergy -srf_base ../geometry/sphere_R6_vdens -num_charges 10 -nmax 25 -srf_num 4 -density 4.0 -ksp_rtol 1.0e-10 -snes_linesearch_type basic -snes_type ksponly -potential_petscspace_degree 1 -k_view output/K_sphere_point_4.mat
	./calcSolvationEnergy -srf_base ../geometry/sphere_R6_vdens -num_charges 10 -nmax 25 -srf_num 5 -density 5.0 -ksp_rtol 1.0e-10 -snes_linesearch_type basic -snes_type ksponly -potential_petscspace_degree 1 -k_view output/K_sphere_point_5.mat
	./calcSolvationEnergy -srf_base ../geometry/sphere_R6_vdens -num_charges 10 -nmax 25 -srf_num 6 -density 6.0 -ksp_rtol 1.0e-10 -snes_linesearch_type basic -snes_type ksponly -potential_petscspace_degree 1 -k_view output/K_sphere_point_6.mat
	./calcSolvationEnergy -srf_base ../geometry/sphere_R6_vdens -num_charges 10 -nmax 25 -srf_num 8 -density 8.0 -ksp_rtol 1.0e-10 -snes_linesearch_type basic -snes_type ksponly -potential_petscspace_degree 1 -k_view output/K_sphere_point_8.mat
#	./calcSolvationEnergy -srf_base ../geometry/sphere_R6_vdens -num_charges 10 -nmax 25 -srf_num 10 -density 10.0 -ksp_rtol 1.0e-10 -snes_linesearch_type basic -snes_type ksponly -potential_petscspace_degree 1 -k_view output/K_sphere_point_10.mat
#	./calcSolvationEnergy -srf_base ../geometry/sphere_R6_vdens -num_charges 10 -nmax 25 -srf_num 20 -density 20.0 -ksp_rtol 1.0e-10 -snes_linesearch_type basic -snes_type ksponly -potential_petscspace_degree 1 -k_view output/K_sphere_point_20.mat
#	./calcSolvationEnergy -srf_base ../geometry/sphere_R6_vdens -num_charges 10 -nmax 25 -srf_num 30 -density 30.0 -ksp_rtol 1.0e-10 -snes_linesearch_type basic -snes_type ksponly -potential_petscspace_degree 1 -k_view output/K_sphere_point_30.mat
#	./calcSolvationEnergy -srf_base ../geometry/sphere_R6_vdens -num_charges 10 -nmax 25 -srf_num 40 -density 40.0 -ksp_rtol 1.0e-10 -snes_linesearch_type basic -snes_type ksponly -potential_petscspace_degree 1 -k_view output/K_sphere_point_40.mat

spheremat_panel:
	./calcSolvationEnergy -srf_base ../geometry/sphere_R6_vdens0 -num_charges 10 -nmax 25 -srf_num 0125 -density 0.125 -ksp_rtol 1.0e-10 -snes_linesearch_type basic -snes_type ksponly -potential_petscspace_degree 1 -k_view output/K_sphere_panel_0125.mat -use_panels 1 -force_nonlinear 1
	./calcSolvationEnergy -srf_base ../geometry/sphere_R6_vdens0 -num_charges 10 -nmax 25 -srf_num 025 -density 0.25 -ksp_rtol 1.0e-10 -snes_linesearch_type basic -snes_type ksponly -potential_petscspace_degree 1 -k_view output/K_sphere_panel_025.mat -use_panels 1 -force_nonlinear 1
	./calcSolvationEnergy -srf_base ../geometry/sphere_R6_vdens0 -num_charges 10 -nmax 25 -srf_num 05 -density 0.5 -ksp_rtol 1.0e-10 -snes_linesearch_type basic -snes_type ksponly -potential_petscspace_degree 1 -k_view output/K_sphere_panel_05.mat -use_panels 1
	./calcSolvationEnergy -srf_base ../geometry/sphere_R6_vdens -num_charges 10 -nmax 25 -srf_num 1 -density 1.0 -ksp_rtol 1.0e-10 -snes_linesearch_type basic -snes_type ksponly -potential_petscspace_degree 1 -k_view output/K_sphere_panel_1.mat -use_panels 1
	./calcSolvationEnergy -srf_base ../geometry/sphere_R6_vdens -num_charges 10 -nmax 25 -srf_num 2 -density 2.0 -ksp_rtol 1.0e-10 -snes_linesearch_type basic -snes_type ksponly -potential_petscspace_degree 1 -k_view output/K_sphere_panel_2.mat -use_panels 1
	./calcSolvationEnergy -srf_base ../geometry/sphere_R6_vdens -num_charges 10 -nmax 25 -srf_num 3 -density 3.0 -ksp_rtol 1.0e-10 -snes_linesearch_type basic -snes_type ksponly -potential_petscspace_degree 1 -k_view output/K_sphere_panel_3.mat -use_panels 1
	./calcSolvationEnergy -srf_base ../geometry/sphere_R6_vdens -num_charges 10 -nmax 25 -srf_num 4 -density 4.0 -ksp_rtol 1.0e-10 -snes_linesearch_type basic -snes_type ksponly -potential_petscspace_degree 1 -k_view output/K_sphere_panel_4.mat -use_panels 1
	./calcSolvationEnergy -srf_base ../geometry/sphere_R6_vdens -num_charges 10 -nmax 25 -srf_num 5 -density 5.0 -ksp_rtol 1.0e-10 -snes_linesearch_type basic -snes_type ksponly -potential_petscspace_degree 1 -k_view output/K_sphere_panel_5.mat -use_panels 1
	./calcSolvationEnergy -srf_base ../geometry/sphere_R6_vdens -num_charges 10 -nmax 25 -srf_num 6 -density 6.0 -ksp_rtol 1.0e-10 -snes_linesearch_type basic -snes_type ksponly -potential_petscspace_degree 1 -k_view output/K_sphere_panel_6.mat -use_panels 1
	./calcSolvationEnergy -srf_base ../geometry/sphere_R6_vdens -num_charges 10 -nmax 25 -srf_num 8 -density 8.0 -ksp_rtol 1.0e-10 -snes_linesearch_type basic -snes_type ksponly -potential_petscspace_degree 1 -k_view output/K_sphere_panel_8.mat -use_panels 1
#	./calcSolvationEnergy -srf_base ../geometry/sphere_R6_vdens -num_charges 10 -nmax 25 -srf_num 10 -density 10.0 -ksp_rtol 1.0e-10 -snes_linesearch_type basic -snes_type ksponly -potential_petscspace_degree 1 -k_view output/K_sphere_panel_10.mat -use_panels 1
#	./calcSolvationEnergy -srf_base ../geometry/sphere_R6_vdens -num_charges 10 -nmax 25 -srf_num 20 -density 20.0 -ksp_rtol 1.0e-10 -snes_linesearch_type basic -snes_type ksponly -potential_petscspace_degree 1 -k_view output/K_sphere_panel_20.mat -use_panels 1
#	./calcSolvationEnergy -srf_base ../geometry/sphere_R6_vdens -num_charges 10 -nmax 25 -srf_num 30 -density 30.0 -ksp_rtol 1.0e-10 -snes_linesearch_type basic -snes_type ksponly -potential_petscspace_degree 1 -k_view output/K_sphere_panel_30.mat -use_panels 1
#	./calcSolvationEnergy -srf_base ../geometry/sphere_R6_vdens -num_charges 10 -nmax 25 -srf_num 40 -density 40.0 -ksp_rtol 1.0e-10 -snes_linesearch_type basic -snes_type ksponly -potential_petscspace_degree 1 -k_view output/K_sphere_panel_40.mat -use_panels 1

spheremat_noprobe_point:
	./calcSolvationEnergy -srf_base ../geometry/sphere_R6_vdens_noprobe_0 -num_charges 10 -nmax 25 -srf_num 0125 -density 0.125 -ksp_rtol 1.0e-10 -snes_linesearch_type basic -snes_type ksponly -potential_petscspace_degree 1 -k_view output/K_sphere_noprobe_point_0125.mat
	./calcSolvationEnergy -srf_base ../geometry/sphere_R6_vdens_noprobe_0 -num_charges 10 -nmax 25 -srf_num 025 -density 0.25 -ksp_rtol 1.0e-10 -snes_linesearch_type basic -snes_type ksponly -potential_petscspace_degree 1 -k_view output/K_sphere_noprobe_point_025.mat
	./calcSolvationEnergy -srf_base ../geometry/sphere_R6_vdens_noprobe_0 -num_charges 10 -nmax 25 -srf_num 05 -density 0.5 -ksp_rtol 1.0e-10 -snes_linesearch_type basic -snes_type ksponly -potential_petscspace_degree 1 -k_view output/K_sphere_noprobe_point_05.mat
	./calcSolvationEnergy -srf_base ../geometry/sphere_R6_vdens_noprobe_ -num_charges 10 -nmax 25 -srf_num 1 -density 1.0 -ksp_rtol 1.0e-10 -snes_linesearch_type basic -snes_type ksponly -potential_petscspace_degree 1 -k_view output/K_sphere_noprobe_point_1.mat
	./calcSolvationEnergy -srf_base ../geometry/sphere_R6_vdens_noprobe_ -num_charges 10 -nmax 25 -srf_num 2 -density 2.0 -ksp_rtol 1.0e-10 -snes_linesearch_type basic -snes_type ksponly -potential_petscspace_degree 1 -k_view output/K_sphere_noprobe_point_2.mat
	./calcSolvationEnergy -srf_base ../geometry/sphere_R6_vdens_noprobe_ -num_charges 10 -nmax 25 -srf_num 3 -density 3.0 -ksp_rtol 1.0e-10 -snes_linesearch_type basic -snes_type ksponly -potential_petscspace_degree 1 -k_view output/K_sphere_noprobe_point_3.mat
	./calcSolvationEnergy -srf_base ../geometry/sphere_R6_vdens_noprobe_ -num_charges 10 -nmax 25 -srf_num 4 -density 4.0 -ksp_rtol 1.0e-10 -snes_linesearch_type basic -snes_type ksponly -potential_petscspace_degree 1 -k_view output/K_sphere_noprobe_point_4.mat
	./calcSolvationEnergy -srf_base ../geometry/sphere_R6_vdens_noprobe_ -num_charges 10 -nmax 25 -srf_num 5 -density 5.0 -ksp_rtol 1.0e-10 -snes_linesearch_type basic -snes_type ksponly -potential_petscspace_degree 1 -k_view output/K_sphere_noprobe_point_5.mat
	./calcSolvationEnergy -srf_base ../geometry/sphere_R6_vdens_noprobe_ -num_charges 10 -nmax 25 -srf_num 6 -density 6.0 -ksp_rtol 1.0e-10 -snes_linesearch_type basic -snes_type ksponly -potential_petscspace_degree 1 -k_view output/K_sphere_noprobe_point_6.mat
	./calcSolvationEnergy -srf_base ../geometry/sphere_R6_vdens_noprobe_ -num_charges 10 -nmax 25 -srf_num 8 -density 8.0 -ksp_rtol 1.0e-10 -snes_linesearch_type basic -snes_type ksponly -potential_petscspace_degree 1 -k_view output/K_sphere_noprobe_point_8.mat
#	./calcSolvationEnergy -srf_base ../geometry/sphere_R6_vdens_noprobe_ -num_charges 10 -nmax 25 -srf_num 10 -density 10.0 -ksp_rtol 1.0e-10 -snes_linesearch_type basic -snes_type ksponly -potential_petscspace_degree 1 -k_view output/K_sphere_noprobe_point_10.mat
#	./calcSolvationEnergy -srf_base ../geometry/sphere_R6_vdens_noprobe_ -num_charges 10 -nmax 25 -srf_num 20 -density 20.0 -ksp_rtol 1.0e-10 -snes_linesearch_type basic -snes_type ksponly -potential_petscspace_degree 1 -k_view output/K_sphere_noprobe_point_20.mat
#	./calcSolvationEnergy -srf_base ../geometry/sphere_R6_vdens_noprobe_ -num_charges 10 -nmax 25 -srf_num 30 -density 30.0 -ksp_rtol 1.0e-10 -snes_linesearch_type basic -snes_type ksponly -potential_petscspace_degree 1 -k_view output/K_sphere_noprobe_point_30.mat
#	./calcSolvationEnergy -srf_base ../geometry/sphere_R6_vdens_noprobe_ -num_charges 10 -nmax 25 -srf_num 40 -density 40.0 -ksp_rtol 1.0e-10 -snes_linesearch_type basic -snes_type ksponly -potential_petscspace_degree 1 -k_view output/K_sphere_noprobe_point_40.mat

spheremat_noprobe_panel:
	./calcSolvationEnergy -srf_base ../geometry/sphere_R6_vdens_noprobe_0 -num_charges 10 -nmax 25 -srf_num 0125 -density 0.125 -ksp_rtol 1.0e-10 -snes_linesearch_type basic -snes_type ksponly -potential_petscspace_degree 1 -k_view output/K_sphere_noprobe_panel_0125.mat -use_panels 1
	./calcSolvationEnergy -srf_base ../geometry/sphere_R6_vdens_noprobe_0 -num_charges 10 -nmax 25 -srf_num 025 -density 0.25 -ksp_rtol 1.0e-10 -snes_linesearch_type basic -snes_type ksponly -potential_petscspace_degree 1 -k_view output/K_sphere_noprobe_panel_025.mat -use_panels 1
	./calcSolvationEnergy -srf_base ../geometry/sphere_R6_vdens_noprobe_0 -num_charges 10 -nmax 25 -srf_num 05 -density 0.5 -ksp_rtol 1.0e-10 -snes_linesearch_type basic -snes_type ksponly -potential_petscspace_degree 1 -k_view output/K_sphere_noprobe_panel_05.mat -use_panels 1
	./calcSolvationEnergy -srf_base ../geometry/sphere_R6_vdens_noprobe_ -num_charges 10 -nmax 25 -srf_num 1 -density 1.0 -ksp_rtol 1.0e-10 -snes_linesearch_type basic -snes_type ksponly -potential_petscspace_degree 1 -k_view output/K_sphere_noprobe_panel_1.mat -use_panels 1
	./calcSolvationEnergy -srf_base ../geometry/sphere_R6_vdens_noprobe_ -num_charges 10 -nmax 25 -srf_num 2 -density 2.0 -ksp_rtol 1.0e-10 -snes_linesearch_type basic -snes_type ksponly -potential_petscspace_degree 1 -k_view output/K_sphere_noprobe_panel_2.mat -use_panels 1
	./calcSolvationEnergy -srf_base ../geometry/sphere_R6_vdens_noprobe_ -num_charges 10 -nmax 25 -srf_num 3 -density 3.0 -ksp_rtol 1.0e-10 -snes_linesearch_type basic -snes_type ksponly -potential_petscspace_degree 1 -k_view output/K_sphere_noprobe_panel_3.mat -use_panels 1
	./calcSolvationEnergy -srf_base ../geometry/sphere_R6_vdens_noprobe_ -num_charges 10 -nmax 25 -srf_num 4 -density 4.0 -ksp_rtol 1.0e-10 -snes_linesearch_type basic -snes_type ksponly -potential_petscspace_degree 1 -k_view output/K_sphere_noprobe_panel_4.mat -use_panels 1
	./calcSolvationEnergy -srf_base ../geometry/sphere_R6_vdens_noprobe_ -num_charges 10 -nmax 25 -srf_num 5 -density 5.0 -ksp_rtol 1.0e-10 -snes_linesearch_type basic -snes_type ksponly -potential_petscspace_degree 1 -k_view output/K_sphere_noprobe_panel_5.mat -use_panels 1
	./calcSolvationEnergy -srf_base ../geometry/sphere_R6_vdens_noprobe_ -num_charges 10 -nmax 25 -srf_num 6 -density 6.0 -ksp_rtol 1.0e-10 -snes_linesearch_type basic -snes_type ksponly -potential_petscspace_degree 1 -k_view output/K_sphere_noprobe_panel_6.mat -use_panels 1
	./calcSolvationEnergy -srf_base ../geometry/sphere_R6_vdens_noprobe_ -num_charges 10 -nmax 25 -srf_num 8 -density 8.0 -ksp_rtol 1.0e-10 -snes_linesearch_type basic -snes_type ksponly -potential_petscspace_degree 1 -k_view output/K_sphere_noprobe_panel_8.mat -use_panels 1
#	./calcSolvationEnergy -srf_base ../geometry/sphere_R6_vdens_noprobe_ -num_charges 10 -nmax 25 -srf_num 10 -density 10.0 -ksp_rtol 1.0e-10 -snes_linesearch_type basic -snes_type ksponly -potential_petscspace_degree 1 -k_view output/K_sphere_noprobe_panel_10.mat -use_panels 1
#	./calcSolvationEnergy -srf_base ../geometry/sphere_R6_vdens_noprobe_ -num_charges 10 -nmax 25 -srf_num 20 -density 20.0 -ksp_rtol 1.0e-10 -snes_linesearch_type basic -snes_type ksponly -potential_petscspace_degree 1 -k_view output/K_sphere_noprobe_panel_20.mat -use_panels 1
#	./calcSolvationEnergy -srf_base ../geometry/sphere_R6_vdens_noprobe_ -num_charges 10 -nmax 25 -srf_num 30 -density 30.0 -ksp_rtol 1.0e-10 -snes_linesearch_type basic -snes_type ksponly -potential_petscspace_degree 1 -k_view output/K_sphere_noprobe_panel_30.mat -use_panels 1
#	./calcSolvationEnergy -srf_base ../geometry/sphere_R6_vdens_noprobe_ -num_charges 10 -nmax 25 -srf_num 40 -density 40.0 -ksp_rtol 1.0e-10 -snes_linesearch_type basic -snes_type ksponly -potential_petscspace_degree 1 -k_view output/K_sphere_noprobe_panel_40.mat -use_panels 1

argpqr_point:
	./calcSolvationEnergy -is_sphere 0 -pqr_filename ../geometry/arg.pqr -srf_base ../geometry/argpqr_scaledcharmm_ -srf_num 1 -ksp_rtol 1.0e-10 -potential_petscspace_degree 1 -k_view output/K_argpqr_point_1.mat -dm_view hdf5:sol.h5 -charge_view hdf5:sol.h5::append -force_nonlinear 1 -alpha 0.5 -beta -60 -gamma -0.5

6dz0_point:
	./calcSolvationEnergy -is_sphere 0 -pqr_filename ../geometry/6dz0.pqr -srf_base ../geometry/6dz0_scaledcharmm_ -srf_num 1 -ksp_rtol 1.0e-16 -potential_petscspace_degree 1 -k_view output/K_6dz0_point_1.mat -dm_view hdf5:sol6dz01.h5 -charge_view hdf5:sol6dz01.h5::append -force_nonlinear 1 -alpha 0.5 -beta -60 -gamma -0.5 -info_out output/6dz0_point_1_data.txt -flops_out output/6dz0_point_1_flops.txt -ksp_gmres_restart 200

3nir_point:
	./calcSolvationEnergy -is_sphere 0 -pqr_filename ../geometry/3nir.pqr -srf_base ../geometry/3nir_scaledcharmm_ -srf_num 1 -ksp_rtol 1.0e-16 -potential_petscspace_degree 1 -k_view output/K_3nir_point_1.mat -dm_view hdf5:sol3nir1.h5 -charge_view hdf5:sol3nir1.h5::append -force_nonlinear 1 -alpha 0.5 -beta -60 -gamma -0.5 -info_out output/3nir_point_1_data.txt -flops_out output/3nir_point_1_flops.txt -ksp_gmres_restart 200
	./calcSolvationEnergy -is_sphere 0 -pqr_filename ../geometry/3nir.pqr -srf_base ../geometry/3nir_scaledcharmm_ -srf_num 15 -ksp_rtol 1.0e-16 -potential_petscspace_degree 1 -k_view output/K_3nir_point_15.mat -dm_view hdf5:sol3nir15.h5 -charge_view hdf5:sol3nir15.h5::append -force_nonlinear 1 -alpha 0.5 -beta -60 -gamma -0.5 -info_out output/3nir_point_15_data.txt -flops_out output/3nir_point_15_flops.txt -ksp_gmres_restart 200
	./calcSolvationEnergy -is_sphere 0 -pqr_filename ../geometry/3nir.pqr -srf_base ../geometry/3nir_scaledcharmm_ -srf_num 2 -ksp_rtol 1.0e-16 -potential_petscspace_degree 1 -k_view output/K_3nir_point_2.mat -dm_view hdf5:sol3nir2.h5 -charge_view hdf5:sol3nir2.h5::append -force_nonlinear 1 -alpha 0.5 -beta -60 -gamma -0.5 -info_out output/3nir_point_2_data.txt -flops_out output/3nir_point_2_flops.txt -ksp_gmres_restart 200
	./calcSolvationEnergy -is_sphere 0 -pqr_filename ../geometry/3nir.pqr -srf_base ../geometry/3nir_scaledcharmm_ -srf_num 25 -ksp_rtol 1.0e-16 -potential_petscspace_degree 1 -k_view output/K_3nir_point_25.mat -dm_view hdf5:sol3nir25.h5 -charge_view hdf5:sol3nir25.h5::append -force_nonlinear 1 -alpha 0.5 -beta -60 -gamma -0.5 -info_out output/3nir_point_25_data.txt -flops_out output/3nir_point_25_flops.txt -ksp_gmres_restart 200
	./calcSolvationEnergy -is_sphere 0 -pqr_filename ../geometry/3nir.pqr -srf_base ../geometry/3nir_scaledcharmm_ -srf_num 3 -ksp_rtol 1.0e-16 -potential_petscspace_degree 1 -k_view output/K_3nir_point_3.mat -dm_view hdf5:sol3nir3.h5 -charge_view hdf5:sol3nir3.h5::append -force_nonlinear 1 -alpha 0.5 -beta -60 -gamma -0.5 -info_out output/3nir_point_3_data.txt -flops_out output/3nir_point_3_flops.txt -ksp_gmres_restart 200
	./calcSolvationEnergy -is_sphere 0 -pqr_filename ../geometry/3nir.pqr -srf_base ../geometry/3nir_scaledcharmm_ -srf_num 35 -ksp_rtol 1.0e-16 -potential_petscspace_degree 1 -k_view output/K_3nir_point_35.mat -dm_view hdf5:sol3nir35.h5 -charge_view hdf5:sol3nir35.h5::append -force_nonlinear 1 -alpha 0.5 -beta -60 -gamma -0.5 -info_out output/3nir_point_35_data.txt -flops_out output/3nir_point_35_flops.txt -ksp_gmres_restart 200
	./calcSolvationEnergy -is_sphere 0 -pqr_filename ../geometry/3nir.pqr -srf_base ../geometry/3nir_scaledcharmm_ -srf_num 4 -ksp_rtol 1.0e-16 -potential_petscspace_degree 1 -k_view output/K_3nir_point_4.mat -dm_view hdf5:sol3nir4.h5 -charge_view hdf5:sol3nir4.h5::append -force_nonlinear 1 -alpha 0.5 -beta -60 -gamma -0.5 -info_out output/3nir_point_4_data.txt -flops_out output/3nir_point_4_flops.txt -ksp_gmres_restart 200
	./calcSolvationEnergy -is_sphere 0 -pqr_filename ../geometry/3nir.pqr -srf_base ../geometry/3nir_scaledcharmm_ -srf_num 5 -ksp_rtol 1.0e-16 -potential_petscspace_degree 1 -k_view output/K_3nir_point_5.mat -dm_view hdf5:sol3nir5.h5 -charge_view hdf5:sol3nir5.h5::append -force_nonlinear 1 -alpha 0.5 -beta -60 -gamma -0.5 -info_out output/3nir_point_5_data.txt -flops_out output/3nir_point_5_flops.txt -ksp_gmres_restart 200
	./calcSolvationEnergy -is_sphere 0 -pqr_filename ../geometry/3nir.pqr -srf_base ../geometry/3nir_scaledcharmm_ -srf_num 6 -ksp_rtol 1.0e-16 -potential_petscspace_degree 1 -k_view output/K_3nir_point_6.mat -dm_view hdf5:sol3nir6.h5 -charge_view hdf5:sol3nir6.h5::append -force_nonlinear 1 -alpha 0.5 -beta -60 -gamma -0.5 -info_out output/3nir_point_6_data.txt -flops_out output/3nir_point_6_flops.txt -ksp_gmres_restart 200
	./calcSolvationEnergy -is_sphere 0 -pqr_filename ../geometry/3nir.pqr -srf_base ../geometry/3nir_scaledcharmm_ -srf_num 7 -ksp_rtol 1.0e-16 -potential_petscspace_degree 1 -k_view output/K_3nir_point_7.mat -dm_view hdf5:sol3nir7.h5 -charge_view hdf5:sol3nir7.h5::append -force_nonlinear 1 -alpha 0.5 -beta -60 -gamma -0.5 -info_out output/3nir_point_7_data.txt -flops_out output/3nir_point_7_flops.txt -ksp_gmres_restart 200

3nir_point_anderson:
	./calcSolvationEnergy -is_sphere 0 -pqr_filename ../geometry/3nir.pqr -srf_base ../geometry/3nir_scaledcharmm_ -srf_num 1 -ksp_rtol 1.0e-16 -potential_petscspace_degree 1 -force_nonlinear 1 -alpha 0.5 -beta -60 -gamma -0.5 -info_out output/3nir_point_anderson_1_data.txt -flops_out output/3nir_point_anderson_1_flops.txt -ksp_gmres_restart 200 -anderson 1
	./calcSolvationEnergy -is_sphere 0 -pqr_filename ../geometry/3nir.pqr -srf_base ../geometry/3nir_scaledcharmm_ -srf_num 2 -ksp_rtol 1.0e-16 -potential_petscspace_degree 1 -force_nonlinear 1 -alpha 0.5 -beta -60 -gamma -0.5 -info_out output/3nir_point_anderson_2_data.txt -flops_out output/3nir_point_anderson_2_flops.txt -ksp_gmres_restart 200 -anderson 1
	./calcSolvationEnergy -is_sphere 0 -pqr_filename ../geometry/3nir.pqr -srf_base ../geometry/3nir_scaledcharmm_ -srf_num 3 -ksp_rtol 1.0e-16 -potential_petscspace_degree 1 -force_nonlinear 1 -alpha 0.5 -beta -60 -gamma -0.5 -info_out output/3nir_point_anderson_3_data.txt -flops_out output/3nir_point_anderson_3_flops.txt -ksp_gmres_restart 200 -anderson 1
	./calcSolvationEnergy -is_sphere 0 -pqr_filename ../geometry/3nir.pqr -srf_base ../geometry/3nir_scaledcharmm_ -srf_num 4 -ksp_rtol 1.0e-16 -potential_petscspace_degree 1 -force_nonlinear 1 -alpha 0.5 -beta -60 -gamma -0.5 -info_out output/3nir_point_anderson_4_data.txt -flops_out output/3nir_point_anderson_4_flops.txt -ksp_gmres_restart 200 -anderson 1
	./calcSolvationEnergy -is_sphere 0 -pqr_filename ../geometry/3nir.pqr -srf_base ../geometry/3nir_scaledcharmm_ -srf_num 5 -ksp_rtol 1.0e-16 -potential_petscspace_degree 1 -force_nonlinear 1 -alpha 0.5 -beta -60 -gamma -0.5 -info_out output/3nir_point_anderson_5_data.txt -flops_out output/3nir_point_anderson_5_flops.txt -ksp_gmres_restart 200 -anderson 1


3nir_bigboy:
	./calcSolvationEnergy -is_sphere 0 -pqr_filename ../geometry/3nir.pqr -srf_base ../geometry/3nir_scaledcharmm_ -srf_num 10 -ksp_rtol 1.0e-16 -potential_petscspace_degree 1 -k_view output/K_3nir_point_10.mat -dm_view hdf5:sol3nir10.h5 -charge_view hdf5:sol3nir10.h5::append -force_nonlinear 1 -alpha 0.5 -beta -60 -gamma -0.5 -info_out output/3nir_point_10_data.txt -flops_out output/3nir_point_10_flops.txt -ksp_gmres_restart 200

3nir_panel:
	./calcSolvationEnergy -is_sphere 0 -pqr_filename ../geometry/3nir.pqr -srf_base ../geometry/3nir_scaledcharmm_ -srf_num 1 -ksp_rtol 1.0e-10 -potential_petscspace_degree 1 -k_view output/K_3nir_panel_1.mat -force_nonlinear 1 -alpha 0.5 -beta -60 -gamma -0.5 -use_panels 1 -info_out output/3nir_panel_1_data.txt -flops_out output/3nir_panel_1_flops.txt
	./calcSolvationEnergy -is_sphere 0 -pqr_filename ../geometry/3nir.pqr -srf_base ../geometry/3nir_scaledcharmm_ -srf_num 15 -ksp_rtol 1.0e-10 -potential_petscspace_degree 1 -k_view output/K_3nir_panel_15.mat -force_nonlinear 1 -alpha 0.5 -beta -60 -gamma -0.5 -use_panels 1 -info_out output/3nir_panel_15_data.txt -flops_out output/3nir_panel_15_flops.txt
	./calcSolvationEnergy -is_sphere 0 -pqr_filename ../geometry/3nir.pqr -srf_base ../geometry/3nir_scaledcharmm_ -srf_num 2 -ksp_rtol 1.0e-10 -potential_petscspace_degree 1 -k_view output/K_3nir_panel_2.mat -force_nonlinear 1 -alpha 0.5 -beta -60 -gamma -0.5 -use_panels 1 -info_out output/3nir_panel_2_data.txt -flops_out output/3nir_panel_2_flops.txt
	./calcSolvationEnergy -is_sphere 0 -pqr_filename ../geometry/3nir.pqr -srf_base ../geometry/3nir_scaledcharmm_ -srf_num 25 -ksp_rtol 1.0e-10 -potential_petscspace_degree 1 -k_view output/K_3nir_panel_25.mat -force_nonlinear 1 -alpha 0.5 -beta -60 -gamma -0.5 -use_panels 1 -info_out output/3nir_panel_25_data.txt -flops_out output/3nir_panel_25_flops.txt
	./calcSolvationEnergy -is_sphere 0 -pqr_filename ../geometry/3nir.pqr -srf_base ../geometry/3nir_scaledcharmm_ -srf_num 3 -ksp_rtol 1.0e-10 -potential_petscspace_degree 1 -k_view output/K_3nir_panel_3.mat -force_nonlinear 1 -alpha 0.5 -beta -60 -gamma -0.5 -use_panels 1 -info_out output/3nir_panel_3_data.txt -flops_out output/3nir_panel_3_flops.txt
	./calcSolvationEnergy -is_sphere 0 -pqr_filename ../geometry/3nir.pqr -srf_base ../geometry/3nir_scaledcharmm_ -srf_num 35 -ksp_rtol 1.0e-10 -potential_petscspace_degree 1 -k_view output/K_3nir_panel_35.mat -force_nonlinear 1 -alpha 0.5 -beta -60 -gamma -0.5 -use_panels 1 -info_out output/3nir_panel_35_data.txt -flops_out output/3nir_panel_35_flops.txt
	./calcSolvationEnergy -is_sphere 0 -pqr_filename ../geometry/3nir.pqr -srf_base ../geometry/3nir_scaledcharmm_ -srf_num 4 -ksp_rtol 1.0e-10 -potential_petscspace_degree 1 -k_view output/K_3nir_panel_4.mat -force_nonlinear 1 -alpha 0.5 -beta -60 -gamma -0.5 -use_panels 1 -info_out output/3nir_panel_4_data.txt -flops_out output/3nir_panel_4_flops.txt
	./calcSolvationEnergy -is_sphere 0 -pqr_filename ../geometry/3nir.pqr -srf_base ../geometry/3nir_scaledcharmm_ -srf_num 5 -ksp_rtol 1.0e-10 -potential_petscspace_degree 1 -k_view output/K_3nir_panel_5.mat -force_nonlinear 1 -alpha 0.5 -beta -60 -gamma -0.5 -use_panels 1 -info_out output/3nir_panel_5_data.txt -flops_out output/3nir_panel_5_flops.txt



1hco_point_nonlinear:
	./calcSolvationEnergy -is_sphere 0 -pdb_filename ../geometry/1hco.pqr -srf_base ../geometry/1hco_scaledcharmm_ -srf_num 1 -ksp_rtol 1.0e-10 -potential_petscspace_degree 1 -k_view output/K_1hco_point_1.mat -dm_view hdf5:sol1hco1.h5 -charge_view hdf5:sol1hco1.h5::append -force_nonlinear 1 -alpha 0.5 -beta -60 -gamma -0.5


1hco_panel_picard:
	./calcSolvationEnergy -is_sphere 0 -pdb_filename ../geometry/1hco.pqr -srf_base ../geometry/1hco_scaledcharmm_ -srf_num 1 -ksp_rtol 1.0e-10 -potential_petscspace_degree 1 -k_view output/K_1hco_panel_1.mat -dm_view hdf5:sol1hco1.h5 -charge_view hdf5:sol1hco1.h5::append -force_nonlinear 1 -alpha 0.5 -beta -60 -gamma -0.5 -use_panels 1 -flops_out output/1hco_panel_1_picard_flops.txt -info_out output/1hco_panel_1_picard_data.txt
	./calcSolvationEnergy -is_sphere 0 -pdb_filename ../geometry/1hco.pqr -srf_base ../geometry/1hco_scaledcharmm_ -srf_num 2 -ksp_rtol 1.0e-10 -potential_petscspace_degree 1 -k_view output/K_1hco_panel_2.mat -dm_view hdf5:sol1hco2.h5 -charge_view hdf5:sol1hco2.h5::append -force_nonlinear 1 -alpha 0.5 -beta -60 -gamma -0.5 -use_panels 1 -flops_out output/1hco_panel_2_picard_flops.txt -info_out output/1hco_panel_2_picard_data.txt
	./calcSolvationEnergy -is_sphere 0 -pdb_filename ../geometry/1hco.pqr -srf_base ../geometry/1hco_scaledcharmm_ -srf_num 3 -ksp_rtol 1.0e-10 -potential_petscspace_degree 1 -k_view output/K_1hco_panel_3.mat -dm_view hdf5:sol1hco3.h5 -charge_view hdf5:sol1hco3.h5::append -force_nonlinear 1 -alpha 0.5 -beta -60 -gamma -0.5 -use_panels 1 -flops_out output/1hco_panel_2_picard_flops.txt -info_out output/1hco_panel_1_picard_data.txt

1hco_panel_anderson:
	./calcSolvationEnergy -is_sphere 0 -pdb_filename ../geometry/1hco.pqr -srf_base ../geometry/1hco_scaledcharmm_ -srf_num 1 -ksp_rtol 1.0e-10 -potential_petscspace_degree 1 -k_view output/K_1hco_panel_1.mat -dm_view hdf5:sol1hco1.h5 -charge_view hdf5:sol1hco1.h5::append -force_nonlinear 1 -alpha 0.5 -beta -60 -gamma -0.5 -use_panels 1 -flops_out output/1hco_panel_1_anderson_flops.txt -info_out output/1hco_panel_1_anderson_data.txt -anderson 1
	./calcSolvationEnergy -is_sphere 0 -pdb_filename ../geometry/1hco.pqr -srf_base ../geometry/1hco_scaledcharmm_ -srf_num 2 -ksp_rtol 1.0e-10 -potential_petscspace_degree 1 -k_view output/K_1hco_panel_2.mat -dm_view hdf5:sol1hco2.h5 -charge_view hdf5:sol1hco2.h5::append -force_nonlinear 1 -alpha 0.5 -beta -60 -gamma -0.5 -use_panels 1 -flops_out output/1hco_panel_2_anderson_flops.txt -info_out output/1hco_panel_2_anderson_data.txt -anderson 1
	./calcSolvationEnergy -is_sphere 0 -pdb_filename ../geometry/1hco.pqr -srf_base ../geometry/1hco_scaledcharmm_ -srf_num 3 -ksp_rtol 1.0e-10 -potential_petscspace_degree 1 -k_view output/K_1hco_panel_3.mat -dm_view hdf5:sol1hco3.h5 -charge_view hdf5:sol1hco3.h5::append -force_nonlinear 1 -alpha 0.5 -beta -60 -gamma -0.5 -use_panels 1 -flops_out output/1hco_panel_2_anderson_flops.txt -info_out output/1hco_panel_1_anderson_data.txt -anderson 1


1hco_point_picard:
	./calcSolvationEnergy -is_sphere 0 -pdb_filename ../geometry/1hco.pqr -srf_base ../geometry/1hco_scaledcharmm_ -srf_num 1 -ksp_rtol 1.0e-10 -potential_petscspace_degree 1 -k_view output/K_1hco_point_1.mat -dm_view hdf5:sol1hco1.h5 -charge_view hdf5:sol1hco1.h5::append -force_nonlinear 1 -alpha 0.5 -beta -60 -gamma -0.5 -flops_out output/1hco_point_1_picard_flops.txt -info_out output/1hco_point_1_picard_data.txt
	./calcSolvationEnergy -is_sphere 0 -pdb_filename ../geometry/1hco.pqr -srf_base ../geometry/1hco_scaledcharmm_ -srf_num 2 -ksp_rtol 1.0e-10 -potential_petscspace_degree 1 -k_view output/K_1hco_point_2.mat -dm_view hdf5:sol1hco2.h5 -charge_view hdf5:sol1hco2.h5::append -force_nonlinear 1 -alpha 0.5 -beta -60 -gamma -0.5 -flops_out output/1hco_point_2_picard_flops.txt -info_out output/1hco_point_2_picard_data.txt
	./calcSolvationEnergy -is_sphere 0 -pdb_filename ../geometry/1hco.pqr -srf_base ../geometry/1hco_scaledcharmm_ -srf_num 3 -ksp_rtol 1.0e-10 -potential_petscspace_degree 1 -k_view output/K_1hco_point_3.mat -dm_view hdf5:sol1hco3.h5 -charge_view hdf5:sol1hco3.h5::append -force_nonlinear 1 -alpha 0.5 -beta -60 -gamma -0.5 -flops_out output/1hco_point_3_picard_flops.txt -info_out output/1hco_point_3_picard_data.txt

1hco_point_anderson:
	./calcSolvationEnergy -is_sphere 0 -pdb_filename ../geometry/1hco.pqr -srf_base ../geometry/1hco_scaledcharmm_ -srf_num 1 -ksp_rtol 1.0e-10 -potential_petscspace_degree 1 -k_view output/K_1hco_point_1.mat -dm_view hdf5:sol1hco1.h5 -charge_view hdf5:sol1hco1.h5::append -force_nonlinear 1 -alpha 0.5 -beta -60 -gamma -0.5 -flops_out output/1hco_point_1_anderson_flops.txt -info_out output/1hco_point_1_anderson_data.txt -anderson 1
	./calcSolvationEnergy -is_sphere 0 -pdb_filename ../geometry/1hco.pqr -srf_base ../geometry/1hco_scaledcharmm_ -srf_num 2 -ksp_rtol 1.0e-10 -potential_petscspace_degree 1 -k_view output/K_1hco_point_2.mat -dm_view hdf5:sol1hco2.h5 -charge_view hdf5:sol1hco2.h5::append -force_nonlinear 1 -alpha 0.5 -beta -60 -gamma -0.5 -flops_out output/1hco_point_2_anderson_flops.txt -info_out output/1hco_point_2_anderson_data.txt -anderson 1
	./calcSolvationEnergy -is_sphere 0 -pdb_filename ../geometry/1hco.pqr -srf_base ../geometry/1hco_scaledcharmm_ -srf_num 3 -ksp_rtol 1.0e-10 -potential_petscspace_degree 1 -k_view output/K_1hco_point_3.mat -dm_view hdf5:sol1hco3.h5 -charge_view hdf5:sol1hco3.h5::append -force_nonlinear 1 -alpha 0.5 -beta -60 -gamma -0.5 -flops_out output/1hco_point_3_anderson_flops.txt -info_out output/1hco_point_3_anderson_data.txt -anderson 1

arg_point_solution:
	./calcSolvationEnergy -is_sphere 0 -pdb_filename ../geometry/arg.pdb -crg_filename ../geometry/jr1.crg -srf_base ../geometry/arg_scaledcharmm_ -srf_num 1 -ksp_rtol 1.0e-10 -potential_petscspace_degree 1 -k_view output/K_arg_point_1.mat -dm_view hdf5:sol.h5 -charge_view hdf5:sol.h5::append

arg_point_anderson:
	./calcSolvationEnergy -is_sphere 0 -pdb_filename ../geometry/arg.pdb -crg_filename ../geometry/jr1.crg -srf_base ../geometry/arg_scaledcharmm_ -srf_num 1 -ksp_rtol 1.0e-10 -potential_petscspace_degree 1 -k_view output/K_arg_point_1.mat -dm_view hdf5:output/arg_1_point_anderson.h5 -charge_view hdf5:output/arg_1_point_anderson.h5::append -force_nonlinear 1 -alpha .5 -beta -60 -gamma -0.5 -flops_out output/arg_point_1_anderson_flops.txt -anderson 1 -info_out output/arg_point_1_anderson_data.txt
	./calcSolvationEnergy -is_sphere 0 -pdb_filename ../geometry/arg.pdb -crg_filename ../geometry/jr1.crg -srf_base ../geometry/arg_scaledcharmm_ -srf_num 2 -ksp_rtol 1.0e-10 -potential_petscspace_degree 1 -k_view output/K_arg_point_2.mat -force_nonlinear 1 -alpha .5 -beta -60 -gamma -0.5 -flops_out output/arg_point_2_anderson_flops.txt -anderson 1 -info_out output/arg_point_2_anderson_data.txt -dm_view hdf5:output/arg_2_point_anderson.h5 -charge_view hdf5:output/arg_2_point_anderson.h5::append
	./calcSolvationEnergy -is_sphere 0 -pdb_filename ../geometry/arg.pdb -crg_filename ../geometry/jr1.crg -srf_base ../geometry/arg_scaledcharmm_ -srf_num 3 -ksp_rtol 1.0e-10 -potential_petscspace_degree 1 -k_view output/K_arg_point_3.mat -force_nonlinear 1 -alpha .5 -beta -60 -gamma -0.5 -flops_out output/arg_point_3_anderson_flops.txt -anderson 1 -info_out output/arg_point_3_anderson_data.txt
	./calcSolvationEnergy -is_sphere 0 -pdb_filename ../geometry/arg.pdb -crg_filename ../geometry/jr1.crg -srf_base ../geometry/arg_scaledcharmm_ -srf_num 4 -ksp_rtol 1.0e-10 -potential_petscspace_degree 1 -k_view output/K_arg_point_4.mat -force_nonlinear 1 -alpha .5 -beta -60 -gamma -0.5 -flops_out output/arg_point_4_anderson_flops.txt -anderson 1 -info_out output/arg_point_4_anderson_data.txt

arg_panel_anderson:
	./calcSolvationEnergy -is_sphere 0 -pdb_filename ../geometry/arg.pdb -crg_filename ../geometry/jr1.crg -srf_base ../geometry/arg_scaledcharmm_ -srf_num 1 -ksp_rtol 1.0e-10 -potential_petscspace_degree 1 -k_view output/K_arg_panel_1.mat -dm_view hdf5:output/arg_1_panel_anderson.h5 -charge_view hdf5:output/arg_1_panel_anderson.h5::append -force_nonlinear 1 -alpha .5 -beta -60 -gamma -0.5 -flops_out output/arg_panel_1_anderson_flops.txt -anderson 1 -info_out output/arg_panel_1_anderson_data.txt -use_panels 1
	./calcSolvationEnergy -is_sphere 0 -pdb_filename ../geometry/arg.pdb -crg_filename ../geometry/jr1.crg -srf_base ../geometry/arg_scaledcharmm_ -srf_num 2 -ksp_rtol 1.0e-10 -potential_petscspace_degree 1 -k_view output/K_arg_panel_2.mat -force_nonlinear 1 -alpha .5 -beta -60 -gamma -0.5 -flops_out output/arg_panel_2_anderson_flops.txt -anderson 1 -info_out output/arg_panel_2_anderson_data.txt -dm_view hdf5:output/arg_2_panel_anderson.h5 -charge_view hdf5:output/arg_2_panel_anderson.h5::append -use_panels 1
	./calcSolvationEnergy -is_sphere 0 -pdb_filename ../geometry/arg.pdb -crg_filename ../geometry/jr1.crg -srf_base ../geometry/arg_scaledcharmm_ -srf_num 3 -ksp_rtol 1.0e-10 -potential_petscspace_degree 1 -k_view output/K_arg_panel_3.mat -force_nonlinear 1 -alpha .5 -beta -60 -gamma -0.5 -flops_out output/arg_panel_3_anderson_flops.txt -anderson 1 -info_out output/arg_panel_3_anderson_data.txt -use_panels 1
	./calcSolvationEnergy -is_sphere 0 -pdb_filename ../geometry/arg.pdb -crg_filename ../geometry/jr1.crg -srf_base ../geometry/arg_scaledcharmm_ -srf_num 4 -ksp_rtol 1.0e-10 -potential_petscspace_degree 1 -k_view output/K_arg_panel_4.mat -force_nonlinear 1 -alpha .5 -beta -60 -gamma -0.5 -flops_out output/arg_panel_4_anderson_flops.txt -anderson 1 -info_out output/arg_panel_4_anderson_data.txt -use_panels 1




argmat_point:
	./calcSolvationEnergy -is_sphere 0 -pdb_filename ../geometry/arg.pdb -crg_filename ../geometry/jr1.crg -srf_base ../geometry/arg_scaledcharmm_ -srf_num 1 -ksp_rtol 1.0e-10 -potential_petscspace_degree 1 -k_view output/K_arg_point_1.mat -force_nonlinear 1 -alpha .5 -beta -60 -gamma -0.5 -flops_out output/arg_point_1_picard_flops.txt -info_out output/arg_point_1_picard_data.txt
	./calcSolvationEnergy -is_sphere 0 -pdb_filename ../geometry/arg.pdb -crg_filename ../geometry/jr1.crg -srf_base ../geometry/arg_scaledcharmm_ -srf_num 2 -ksp_rtol 1.0e-10 -potential_petscspace_degree 1 -k_view output/K_arg_point_2.mat -force_nonlinear 1 -alpha .5 -beta -60 -gamma -0.5 -flops_out output/arg_point_2_picard_flops.txt -info_out output/arg_point_2_picard_data.txt
	./calcSolvationEnergy -is_sphere 0 -pdb_filename ../geometry/arg.pdb -crg_filename ../geometry/jr1.crg -srf_base ../geometry/arg_scaledcharmm_ -srf_num 3 -ksp_rtol 1.0e-10 -potential_petscspace_degree 1 -k_view output/K_arg_point_3.mat -force_nonlinear 1 -alpha .5 -beta -60 -gamma -0.5 -flops_out output/arg_point_3_picard_flops.txt -info_out output/arg_point_3_picard_data.txt
	./calcSolvationEnergy -is_sphere 0 -pdb_filename ../geometry/arg.pdb -crg_filename ../geometry/jr1.crg -srf_base ../geometry/arg_scaledcharmm_ -srf_num 4 -ksp_rtol 1.0e-10 -potential_petscspace_degree 1 -k_view output/K_arg_point_4.mat -force_nonlinear 1 -alpha .5 -beta -60 -gamma -0.5 -flops_out output/arg_point_4_picard_flops.txt
	./calcSolvationEnergy -is_sphere 0 -pdb_filename ../geometry/arg.pdb -crg_filename ../geometry/jr1.crg -srf_base ../geometry/arg_scaledcharmm_ -srf_num 5 -ksp_rtol 1.0e-10 -potential_petscspace_degree 1 -k_view output/K_arg_point_5.mat -force_nonlinear 1 -alpha .5 -beta -60 -gamma -0.5
	./calcSolvationEnergy -is_sphere 0 -pdb_filename ../geometry/arg.pdb -crg_filename ../geometry/jr1.crg -srf_base ../geometry/arg_scaledcharmm_ -srf_num 6 -ksp_rtol 1.0e-10 -potential_petscspace_degree 1 -k_view output/K_arg_point_6.mat -force_nonlinear 1 -alpha .5 -beta -60 -gamma -0.5
	./calcSolvationEnergy -is_sphere 0 -pdb_filename ../geometry/arg.pdb -crg_filename ../geometry/jr1.crg -srf_base ../geometry/arg_scaledcharmm_ -srf_num 7 -ksp_rtol 1.0e-10 -potential_petscspace_degree 1 -k_view output/K_arg_point_7.mat -force_nonlinear 1 -alpha .5 -beta -60 -gamma -0.5
	./calcSolvationEnergy -is_sphere 0 -pdb_filename ../geometry/arg.pdb -crg_filename ../geometry/jr1.crg -srf_base ../geometry/arg_scaledcharmm_ -srf_num 8 -ksp_rtol 1.0e-10 -potential_petscspace_degree 1 -k_view output/K_arg_point_8.mat -force_nonlinear 1 -alpha .5 -beta -60 -gamma -0.5
	./calcSolvationEnergy -is_sphere 0 -pdb_filename ../geometry/arg.pdb -crg_filename ../geometry/jr1.crg -srf_base ../geometry/arg_scaledcharmm_ -srf_num 12 -ksp_rtol 1.0e-10 -potential_petscspace_degree 1 -k_view output/K_arg_point_12.mat -force_nonlinear 1 -alpha .5 -beta -60 -gamma -0.5
	./calcSolvationEnergy -is_sphere 0 -pdb_filename ../geometry/arg.pdb -crg_filename ../geometry/jr1.crg -srf_base ../geometry/arg_scaledcharmm_ -srf_num 16 -ksp_rtol 1.0e-10 -potential_petscspace_degree 1 -k_view output/K_arg_point_16.mat -force_nonlinear 1 -alpha .5 -beta -60 -gamma -0.5
#	./calcSolvationEnergy -is_sphere 0 -pdb_filename ../geometry/arg.pdb -crg_filename ../geometry/jr1.crg -srf_base ../geometry/arg_scaledcharmm_ -srf_num 20 -ksp_rtol 1.0e-10 -potential_petscspace_degree 1 -k_view output/K_arg_point_20.mat
#	./calcSolvationEnergy -is_sphere 0 -pdb_filename ../geometry/arg.pdb -crg_filename ../geometry/jr1.crg -srf_base ../geometry/arg_scaledcharmm_ -srf_num 24 -ksp_rtol 1.0e-10 -potential_petscspace_degree 1 -k_view output/K_arg_point_24.mat
#	./calcSolvationEnergy -is_sphere 0 -pdb_filename ../geometry/arg.pdb -crg_filename ../geometry/jr1.crg -srf_base ../geometry/arg_scaledcharmm_ -srf_num 28 -ksp_rtol 1.0e-10 -potential_petscspace_degree 1 -k_view output/K_arg_point_28.mat
#	./calcSolvationEnergy -is_sphere 0 -pdb_filename ../geometry/arg.pdb -crg_filename ../geometry/jr1.crg -srf_base ../geometry/arg_scaledcharmm_ -srf_num 32 -ksp_rtol 1.0e-10 -potential_petscspace_degree 1 -k_view output/K_arg_point_32.mat
#	./calcSolvationEnergy -is_sphere 0 -pdb_filename ../geometry/arg.pdb -crg_filename ../geometry/jr1.crg -srf_base ../geometry/arg_scaledcharmm_ -srf_num 50 -ksp_rtol 1.0e-10 -potential_petscspace_degree 1 -k_view output/K_arg_point_50.mat

argmat_panel:
	./calcSolvationEnergy -is_sphere 0 -pdb_filename ../geometry/arg.pdb -crg_filename ../geometry/jr1.crg -srf_base ../geometry/arg_scaledcharmm_ -srf_num 1 -ksp_rtol 1.0e-10 -potential_petscspace_degree 1 -k_view output/K_arg_panel_1.mat -use_panels 1 -dm_view hdf5:solarg1.h5 -charge_view hdf5:solarg1.h5::append -force_nonlinear 1 -alpha .5 -beta -60 -gamma -0.5
	./calcSolvationEnergy -is_sphere 0 -pdb_filename ../geometry/arg.pdb -crg_filename ../geometry/jr1.crg -srf_base ../geometry/arg_scaledcharmm_ -srf_num 2 -ksp_rtol 1.0e-10 -potential_petscspace_degree 1 -k_view output/K_arg_panel_2.mat -use_panels 1 -dm_view hdf5:solarg2.h5 -charge_view hdf5:solarg2.h5::append -force_nonlinear 1 -alpha .5 -beta -60 -gamma -0.5
	./calcSolvationEnergy -is_sphere 0 -pdb_filename ../geometry/arg.pdb -crg_filename ../geometry/jr1.crg -srf_base ../geometry/arg_scaledcharmm_ -srf_num 3 -ksp_rtol 1.0e-10 -potential_petscspace_degree 1 -k_view output/K_arg_panel_3.mat -use_panels 1 -dm_view hdf5:solarg3.h5 -charge_view hdf5:solarg3.h5::append -force_nonlinear 1 -alpha .5 -beta -60 -gamma -0.5
	./calcSolvationEnergy -is_sphere 0 -pdb_filename ../geometry/arg.pdb -crg_filename ../geometry/jr1.crg -srf_base ../geometry/arg_scaledcharmm_ -srf_num 4 -ksp_rtol 1.0e-10 -potential_petscspace_degree 1 -k_view output/K_arg_panel_4.mat -use_panels 1 -dm_view hdf5:solarg4.h5 -charge_view hdf5:solarg4.h5::append -force_nonlinear 1 -alpha .5 -beta -60 -gamma -0.5
	./calcSolvationEnergy -is_sphere 0 -pdb_filename ../geometry/arg.pdb -crg_filename ../geometry/jr1.crg -srf_base ../geometry/arg_scaledcharmm_ -srf_num 5 -ksp_rtol 1.0e-10 -potential_petscspace_degree 1 -k_view output/K_arg_panel_5.mat -use_panels 1 -dm_view hdf5:solarg5.h5 -charge_view hdf5:solarg5.h5::append -force_nonlinear 1 -alpha .5 -beta -60 -gamma -0.5
	./calcSolvationEnergy -is_sphere 0 -pdb_filename ../geometry/arg.pdb -crg_filename ../geometry/jr1.crg -srf_base ../geometry/arg_scaledcharmm_ -srf_num 6 -ksp_rtol 1.0e-10 -potential_petscspace_degree 1 -k_view output/K_arg_panel_6.mat -use_panels 1 -dm_view hdf5:solarg6.h5 -charge_view hdf5:solarg6.h5::append -force_nonlinear 1 -alpha .5 -beta -60 -gamma -0.5
	./calcSolvationEnergy -is_sphere 0 -pdb_filename ../geometry/arg.pdb -crg_filename ../geometry/jr1.crg -srf_base ../geometry/arg_scaledcharmm_ -srf_num 7 -ksp_rtol 1.0e-10 -potential_petscspace_degree 1 -k_view output/K_arg_panel_7.mat -use_panels 1 -dm_view hdf5:solarg7.h5 -charge_view hdf5:solarg7.h5::append -force_nonlinear 1 -alpha .5 -beta -60 -gamma -0.5
	./calcSolvationEnergy -is_sphere 0 -pdb_filename ../geometry/arg.pdb -crg_filename ../geometry/jr1.crg -srf_base ../geometry/arg_scaledcharmm_ -srf_num 8 -ksp_rtol 1.0e-10 -potential_petscspace_degree 1 -k_view output/K_arg_panel_8.mat -use_panels 1 -dm_view hdf5:solarg8.h5 -charge_view hdf5:solarg8.h5::append -force_nonlinear 1 -alpha .5 -beta -60 -gamma -0.5
	./calcSolvationEnergy -is_sphere 0 -pdb_filename ../geometry/arg.pdb -crg_filename ../geometry/jr1.crg -srf_base ../geometry/arg_scaledcharmm_ -srf_num 12 -ksp_rtol 1.0e-10 -potential_petscspace_degree 1 -k_view output/K_arg_panel_12.mat -use_panels 1
	./calcSolvationEnergy -is_sphere 0 -pdb_filename ../geometry/arg.pdb -crg_filename ../geometry/jr1.crg -srf_base ../geometry/arg_scaledcharmm_ -srf_num 16 -ksp_rtol 1.0e-10 -potential_petscspace_degree 1 -k_view output/K_arg_panel_16.mat -use_panels 1
#	./calcSolvationEnergy -is_sphere 0 -pdb_filename ../geometry/arg.pdb -crg_filename ../geometry/jr1.crg -srf_base ../geometry/arg_scaledcharmm_ -srf_num 20 -ksp_rtol 1.0e-10 -potential_petscspace_degree 1 -k_view output/K_arg_panel_20.mat -use_panels 1
#	./calcSolvationEnergy -is_sphere 0 -pdb_filename ../geometry/arg.pdb -crg_filename ../geometry/jr1.crg -srf_base ../geometry/arg_scaledcharmm_ -srf_num 24 -ksp_rtol 1.0e-10 -potential_petscspace_degree 1 -k_view output/K_arg_panel_24.mat -use_panels 1
#	./calcSolvationEnergy -is_sphere 0 -pdb_filename ../geometry/arg.pdb -crg_filename ../geometry/jr1.crg -srf_base ../geometry/arg_scaledcharmm_ -srf_num 28 -ksp_rtol 1.0e-10 -potential_petscspace_degree 1 -k_view output/K_arg_panel_28.mat -use_panels 1
#	./calcSolvationEnergy -is_sphere 0 -pdb_filename ../geometry/arg.pdb -crg_filename ../geometry/jr1.crg -srf_base ../geometry/arg_scaledcharmm_ -srf_num 32 -ksp_rtol 1.0e-10 -potential_petscspace_degree 1 -k_view output/K_arg_panel_32.mat -use_panels 1
#	./calcSolvationEnergy -is_sphere 0 -pdb_filename ../geometry/arg.pdb -crg_filename ../geometry/jr1.crg -srf_base ../geometry/arg_scaledcharmm_ -srf_num 50 -ksp_rtol 1.0e-10 -potential_petscspace_degree 1 -k_view output/K_arg_panel_50.mat -use_panels 1

argmat_noprobe_point:
	./calcSolvationEnergy -is_sphere 0 -pdb_filename ../geometry/arg.pdb -crg_filename ../geometry/jr1.crg -srf_base ../geometry/arg_scaledcharmm_noprobe_ -srf_num 1 -ksp_rtol 1.0e-10 -potential_petscspace_degree 1 -k_view output/K_arg_noprobe_point_1.mat
	./calcSolvationEnergy -is_sphere 0 -pdb_filename ../geometry/arg.pdb -crg_filename ../geometry/jr1.crg -srf_base ../geometry/arg_scaledcharmm_noprobe_ -srf_num 2 -ksp_rtol 1.0e-10 -potential_petscspace_degree 1 -k_view output/K_arg_noprobe_point_2.mat
	./calcSolvationEnergy -is_sphere 0 -pdb_filename ../geometry/arg.pdb -crg_filename ../geometry/jr1.crg -srf_base ../geometry/arg_scaledcharmm_noprobe_ -srf_num 3 -ksp_rtol 1.0e-10 -potential_petscspace_degree 1 -k_view output/K_arg_noprobe_point_3.mat
	./calcSolvationEnergy -is_sphere 0 -pdb_filename ../geometry/arg.pdb -crg_filename ../geometry/jr1.crg -srf_base ../geometry/arg_scaledcharmm_noprobe_ -srf_num 4 -ksp_rtol 1.0e-10 -potential_petscspace_degree 1 -k_view output/K_arg_noprobe_point_4.mat
	./calcSolvationEnergy -is_sphere 0 -pdb_filename ../geometry/arg.pdb -crg_filename ../geometry/jr1.crg -srf_base ../geometry/arg_scaledcharmm_noprobe_ -srf_num 5 -ksp_rtol 1.0e-10 -potential_petscspace_degree 1 -k_view output/K_arg_noprobe_point_5.mat
	./calcSolvationEnergy -is_sphere 0 -pdb_filename ../geometry/arg.pdb -crg_filename ../geometry/jr1.crg -srf_base ../geometry/arg_scaledcharmm_noprobe_ -srf_num 6 -ksp_rtol 1.0e-10 -potential_petscspace_degree 1 -k_view output/K_arg_noprobe_point_6.mat
	./calcSolvationEnergy -is_sphere 0 -pdb_filename ../geometry/arg.pdb -crg_filename ../geometry/jr1.crg -srf_base ../geometry/arg_scaledcharmm_noprobe_ -srf_num 7 -ksp_rtol 1.0e-10 -potential_petscspace_degree 1 -k_view output/K_arg_noprobe_point_7.mat
	./calcSolvationEnergy -is_sphere 0 -pdb_filename ../geometry/arg.pdb -crg_filename ../geometry/jr1.crg -srf_base ../geometry/arg_scaledcharmm_noprobe_ -srf_num 8 -ksp_rtol 1.0e-10 -potential_petscspace_degree 1 -k_view output/K_arg_noprobe_point_8.mat
	./calcSolvationEnergy -is_sphere 0 -pdb_filename ../geometry/arg.pdb -crg_filename ../geometry/jr1.crg -srf_base ../geometry/arg_scaledcharmm_noprobe_ -srf_num 12 -ksp_rtol 1.0e-10 -potential_petscspace_degree 1 -k_view output/K_arg_noprobe_point_12.mat
	./calcSolvationEnergy -is_sphere 0 -pdb_filename ../geometry/arg.pdb -crg_filename ../geometry/jr1.crg -srf_base ../geometry/arg_scaledcharmm_noprobe_ -srf_num 16 -ksp_rtol 1.0e-10 -potential_petscspace_degree 1 -k_view output/K_arg_noprobe_point_16.mat
#	./calcSolvationEnergy -is_sphere 0 -pdb_filename ../geometry/arg.pdb -crg_filename ../geometry/jr1.crg -srf_base ../geometry/arg_scaledcharmm_noprobe_ -srf_num 20 -ksp_rtol 1.0e-10 -potential_petscspace_degree 1 -k_view output/K_arg_noprobe_point_20.mat
#	./calcSolvationEnergy -is_sphere 0 -pdb_filename ../geometry/arg.pdb -crg_filename ../geometry/jr1.crg -srf_base ../geometry/arg_scaledcharmm_noprobe_ -srf_num 24 -ksp_rtol 1.0e-10 -potential_petscspace_degree 1 -k_view output/K_arg_noprobe_point_24.mat
#	./calcSolvationEnergy -is_sphere 0 -pdb_filename ../geometry/arg.pdb -crg_filename ../geometry/jr1.crg -srf_base ../geometry/arg_scaledcharmm_noprobe_ -srf_num 28 -ksp_rtol 1.0e-10 -potential_petscspace_degree 1 -k_view output/K_arg_noprobe_point_28.mat
#	./calcSolvationEnergy -is_sphere 0 -pdb_filename ../geometry/arg.pdb -crg_filename ../geometry/jr1.crg -srf_base ../geometry/arg_scaledcharmm_noprobe_ -srf_num 32 -ksp_rtol 1.0e-10 -potential_petscspace_degree 1 -k_view output/K_arg_noprobe_point_32.mat
#	./calcSolvationEnergy -is_sphere 0 -pdb_filename ../geometry/arg.pdb -crg_filename ../geometry/jr1.crg -srf_base ../geometry/arg_scaledcharmm_noprobe_ -srf_num 50 -ksp_rtol 1.0e-10 -potential_petscspace_degree 1 -k_view output/K_arg_noprobe_point_50.mat

argmat_noprobe_panel:
	./calcSolvationEnergy -is_sphere 0 -pdb_filename ../geometry/arg.pdb -crg_filename ../geometry/jr1.crg -srf_base ../geometry/arg_scaledcharmm_noprobe_ -srf_num 1 -ksp_rtol 1.0e-10 -potential_petscspace_degree 1 -k_view output/K_arg_noprobe_panel_1.mat -use_panels 1
	./calcSolvationEnergy -is_sphere 0 -pdb_filename ../geometry/arg.pdb -crg_filename ../geometry/jr1.crg -srf_base ../geometry/arg_scaledcharmm_noprobe_ -srf_num 2 -ksp_rtol 1.0e-10 -potential_petscspace_degree 1 -k_view output/K_arg_noprobe_panel_2.mat -use_panels 1
	./calcSolvationEnergy -is_sphere 0 -pdb_filename ../geometry/arg.pdb -crg_filename ../geometry/jr1.crg -srf_base ../geometry/arg_scaledcharmm_noprobe_ -srf_num 3 -ksp_rtol 1.0e-10 -potential_petscspace_degree 1 -k_view output/K_arg_noprobe_panel_3.mat -use_panels 1
	./calcSolvationEnergy -is_sphere 0 -pdb_filename ../geometry/arg.pdb -crg_filename ../geometry/jr1.crg -srf_base ../geometry/arg_scaledcharmm_noprobe_ -srf_num 4 -ksp_rtol 1.0e-10 -potential_petscspace_degree 1 -k_view output/K_arg_noprobe_panel_4.mat -use_panels 1
	./calcSolvationEnergy -is_sphere 0 -pdb_filename ../geometry/arg.pdb -crg_filename ../geometry/jr1.crg -srf_base ../geometry/arg_scaledcharmm_noprobe_ -srf_num 5 -ksp_rtol 1.0e-10 -potential_petscspace_degree 1 -k_view output/K_arg_noprobe_panel_5.mat -use_panels 1
	./calcSolvationEnergy -is_sphere 0 -pdb_filename ../geometry/arg.pdb -crg_filename ../geometry/jr1.crg -srf_base ../geometry/arg_scaledcharmm_noprobe_ -srf_num 6 -ksp_rtol 1.0e-10 -potential_petscspace_degree 1 -k_view output/K_arg_noprobe_panel_6.mat -use_panels 1
	./calcSolvationEnergy -is_sphere 0 -pdb_filename ../geometry/arg.pdb -crg_filename ../geometry/jr1.crg -srf_base ../geometry/arg_scaledcharmm_noprobe_ -srf_num 7 -ksp_rtol 1.0e-10 -potential_petscspace_degree 1 -k_view output/K_arg_noprobe_panel_7.mat -use_panels 1
	./calcSolvationEnergy -is_sphere 0 -pdb_filename ../geometry/arg.pdb -crg_filename ../geometry/jr1.crg -srf_base ../geometry/arg_scaledcharmm_noprobe_ -srf_num 8 -ksp_rtol 1.0e-10 -potential_petscspace_degree 1 -k_view output/K_arg_noprobe_panel_8.mat -use_panels 1
	./calcSolvationEnergy -is_sphere 0 -pdb_filename ../geometry/arg.pdb -crg_filename ../geometry/jr1.crg -srf_base ../geometry/arg_scaledcharmm_noprobe_ -srf_num 12 -ksp_rtol 1.0e-10 -potential_petscspace_degree 1 -k_view output/K_arg_noprobe_panel_12.mat -use_panels 1
	./calcSolvationEnergy -is_sphere 0 -pdb_filename ../geometry/arg.pdb -crg_filename ../geometry/jr1.crg -srf_base ../geometry/arg_scaledcharmm_noprobe_ -srf_num 16 -ksp_rtol 1.0e-10 -potential_petscspace_degree 1 -k_view output/K_arg_noprobe_panel_16.mat -use_panels 1
#	./calcSolvationEnergy -is_sphere 0 -pdb_filename ../geometry/arg.pdb -crg_filename ../geometry/jr1.crg -srf_base ../geometry/arg_scaledcharmm_noprobe_ -srf_num 20 -ksp_rtol 1.0e-10 -potential_petscspace_degree 1 -k_view output/K_arg_noprobe_panel_20.mat -use_panels 1
#	./calcSolvationEnergy -is_sphere 0 -pdb_filename ../geometry/arg.pdb -crg_filename ../geometry/jr1.crg -srf_base ../geometry/arg_scaledcharmm_noprobe_ -srf_num 24 -ksp_rtol 1.0e-10 -potential_petscspace_degree 1 -k_view output/K_arg_noprobe_panel_24.mat -use_panels 1
#	./calcSolvationEnergy -is_sphere 0 -pdb_filename ../geometry/arg.pdb -crg_filename ../geometry/jr1.crg -srf_base ../geometry/arg_scaledcharmm_noprobe_ -srf_num 28 -ksp_rtol 1.0e-10 -potential_petscspace_degree 1 -k_view output/K_arg_noprobe_panel_28.mat -use_panels 1
#	./calcSolvationEnergy -is_sphere 0 -pdb_filename ../geometry/arg.pdb -crg_filename ../geometry/jr1.crg -srf_base ../geometry/arg_scaledcharmm_noprobe_ -srf_num 32 -ksp_rtol 1.0e-10 -potential_petscspace_degree 1 -k_view output/K_arg_noprobe_panel_32.mat -use_panels 1
#	./calcSolvationEnergy -is_sphere 0 -pdb_filename ../geometry/arg.pdb -crg_filename ../geometry/jr1.crg -srf_base ../geometry/arg_scaledcharmm_noprobe_ -srf_num 50 -ksp_rtol 1.0e-10 -potential_petscspace_degree 1 -k_view output/K_arg_noprobe_panel_50.mat -use_panels 1


aspmat_point:
	./calcSolvationEnergy -is_sphere 0 -pdb_filename ../geometry/asp.pdb -crg_filename ../geometry/jd1.crg -srf_base ../geometry/asp_scaledcharmm_ -srf_num 1 -potential_petscspace_degree 1 -k_view output/K_asp_point_1.mat
	./calcSolvationEnergy -is_sphere 0 -pdb_filename ../geometry/asp.pdb -crg_filename ../geometry/jd1.crg -srf_base ../geometry/asp_scaledcharmm_ -srf_num 2 -potential_petscspace_degree 1 -k_view output/K_asp_point_2.mat
	./calcSolvationEnergy -is_sphere 0 -pdb_filename ../geometry/asp.pdb -crg_filename ../geometry/jd1.crg -srf_base ../geometry/asp_scaledcharmm_ -srf_num 3 -potential_petscspace_degree 1 -k_view output/K_asp_point_3.mat
	./calcSolvationEnergy -is_sphere 0 -pdb_filename ../geometry/asp.pdb -crg_filename ../geometry/jd1.crg -srf_base ../geometry/asp_scaledcharmm_ -srf_num 4 -potential_petscspace_degree 1 -k_view output/K_asp_point_4.mat
	./calcSolvationEnergy -is_sphere 0 -pdb_filename ../geometry/asp.pdb -crg_filename ../geometry/jd1.crg -srf_base ../geometry/asp_scaledcharmm_ -srf_num 5 -potential_petscspace_degree 1 -k_view output/K_asp_point_5.mat
	./calcSolvationEnergy -is_sphere 0 -pdb_filename ../geometry/asp.pdb -crg_filename ../geometry/jd1.crg -srf_base ../geometry/asp_scaledcharmm_ -srf_num 6 -potential_petscspace_degree 1 -k_view output/K_asp_point_6.mat
	./calcSolvationEnergy -is_sphere 0 -pdb_filename ../geometry/asp.pdb -crg_filename ../geometry/jd1.crg -srf_base ../geometry/asp_scaledcharmm_ -srf_num 7 -potential_petscspace_degree 1 -k_view output/K_asp_point_7.mat
	./calcSolvationEnergy -is_sphere 0 -pdb_filename ../geometry/asp.pdb -crg_filename ../geometry/jd1.crg -srf_base ../geometry/asp_scaledcharmm_ -srf_num 8 -potential_petscspace_degree 1 -k_view output/K_asp_point_8.mat
	./calcSolvationEnergy -is_sphere 0 -pdb_filename ../geometry/asp.pdb -crg_filename ../geometry/jd1.crg -srf_base ../geometry/asp_scaledcharmm_ -srf_num 12 -potential_petscspace_degree 1 -k_view output/K_asp_point_12.mat
	./calcSolvationEnergy -is_sphere 0 -pdb_filename ../geometry/asp.pdb -crg_filename ../geometry/jd1.crg -srf_base ../geometry/asp_scaledcharmm_ -srf_num 16 -potential_petscspace_degree 1 -k_view output/K_asp_point_16.mat
#	./calcSolvationEnergy -is_sphere 0 -pdb_filename ../geometry/asp.pdb -crg_filename ../geometry/jd1.crg -srf_base ../geometry/asp_scaledcharmm_ -srf_num 20 -potential_petscspace_degree 1 -k_view output/K_asp_point_20.mat
#	./calcSolvationEnergy -is_sphere 0 -pdb_filename ../geometry/asp.pdb -crg_filename ../geometry/jd1.crg -srf_base ../geometry/asp_scaledcharmm_ -srf_num 24 -potential_petscspace_degree 1 -k_view output/K_asp_point_24.mat
#	./calcSolvationEnergy -is_sphere 0 -pdb_filename ../geometry/asp.pdb -crg_filename ../geometry/jd1.crg -srf_base ../geometry/asp_scaledcharmm_ -srf_num 28 -potential_petscspace_degree 1 -k_view output/K_asp_point_28.mat
#	./calcSolvationEnergy -is_sphere 0 -pdb_filename ../geometry/asp.pdb -crg_filename ../geometry/jd1.crg -srf_base ../geometry/asp_scaledcharmm_ -srf_num 32 -potential_petscspace_degree 1 -k_view output/K_asp_point_32.mat
#	./calcSolvationEnergy -is_sphere 0 -pdb_filename ../geometry/asp.pdb -crg_filename ../geometry/jd1.crg -srf_base ../geometry/asp_scaledcharmm_ -srf_num 50 -potential_petscspace_degree 1 -k_view output/K_asp_point_50.mat

aspmat_panel:
	./calcSolvationEnergy -is_sphere 0 -pdb_filename ../geometry/asp.pdb -crg_filename ../geometry/jd1.crg -srf_base ../geometry/asp_scaledcharmm_ -srf_num 1 -potential_petscspace_degree 1 -k_view output/K_asp_panel_1.mat -use_panels 1
	./calcSolvationEnergy -is_sphere 0 -pdb_filename ../geometry/asp.pdb -crg_filename ../geometry/jd1.crg -srf_base ../geometry/asp_scaledcharmm_ -srf_num 2 -potential_petscspace_degree 1 -k_view output/K_asp_panel_2.mat -use_panels 1
	./calcSolvationEnergy -is_sphere 0 -pdb_filename ../geometry/asp.pdb -crg_filename ../geometry/jd1.crg -srf_base ../geometry/asp_scaledcharmm_ -srf_num 3 -potential_petscspace_degree 1 -k_view output/K_asp_panel_3.mat -use_panels 1
	./calcSolvationEnergy -is_sphere 0 -pdb_filename ../geometry/asp.pdb -crg_filename ../geometry/jd1.crg -srf_base ../geometry/asp_scaledcharmm_ -srf_num 4 -potential_petscspace_degree 1 -k_view output/K_asp_panel_4.mat -use_panels 1
	./calcSolvationEnergy -is_sphere 0 -pdb_filename ../geometry/asp.pdb -crg_filename ../geometry/jd1.crg -srf_base ../geometry/asp_scaledcharmm_ -srf_num 5 -potential_petscspace_degree 1 -k_view output/K_asp_panel_5.mat -use_panels 1
	./calcSolvationEnergy -is_sphere 0 -pdb_filename ../geometry/asp.pdb -crg_filename ../geometry/jd1.crg -srf_base ../geometry/asp_scaledcharmm_ -srf_num 6 -potential_petscspace_degree 1 -k_view output/K_asp_panel_6.mat -use_panels 1
	./calcSolvationEnergy -is_sphere 0 -pdb_filename ../geometry/asp.pdb -crg_filename ../geometry/jd1.crg -srf_base ../geometry/asp_scaledcharmm_ -srf_num 7 -potential_petscspace_degree 1 -k_view output/K_asp_panel_7.mat -use_panels 1
	./calcSolvationEnergy -is_sphere 0 -pdb_filename ../geometry/asp.pdb -crg_filename ../geometry/jd1.crg -srf_base ../geometry/asp_scaledcharmm_ -srf_num 8 -potential_petscspace_degree 1 -k_view output/K_asp_panel_8.mat -use_panels 1
	./calcSolvationEnergy -is_sphere 0 -pdb_filename ../geometry/asp.pdb -crg_filename ../geometry/jd1.crg -srf_base ../geometry/asp_scaledcharmm_ -srf_num 12 -potential_petscspace_degree 1 -k_view output/K_asp_panel_12.mat -use_panels 1
	./calcSolvationEnergy -is_sphere 0 -pdb_filename ../geometry/asp.pdb -crg_filename ../geometry/jd1.crg -srf_base ../geometry/asp_scaledcharmm_ -srf_num 16 -potential_petscspace_degree 1 -k_view output/K_asp_panel_16.mat -use_panels 1
#	./calcSolvationEnergy -is_sphere 0 -pdb_filename ../geometry/asp.pdb -crg_filename ../geometry/jd1.crg -srf_base ../geometry/asp_scaledcharmm_ -srf_num 20 -potential_petscspace_degree 1 -k_view output/K_asp_panel_20.mat -use_panels 1
#	./calcSolvationEnergy -is_sphere 0 -pdb_filename ../geometry/asp.pdb -crg_filename ../geometry/jd1.crg -srf_base ../geometry/asp_scaledcharmm_ -srf_num 24 -potential_petscspace_degree 1 -k_view output/K_asp_panel_24.mat -use_panels 1
#	./calcSolvationEnergy -is_sphere 0 -pdb_filename ../geometry/asp.pdb -crg_filename ../geometry/jd1.crg -srf_base ../geometry/asp_scaledcharmm_ -srf_num 28 -potential_petscspace_degree 1 -k_view output/K_asp_panel_28.mat -use_panels 1
#	./calcSolvationEnergy -is_sphere 0 -pdb_filename ../geometry/asp.pdb -crg_filename ../geometry/jd1.crg -srf_base ../geometry/asp_scaledcharmm_ -srf_num 32 -potential_petscspace_degree 1 -k_view output/K_asp_panel_32.mat -use_panels 1
#	./calcSolvationEnergy -is_sphere 0 -pdb_filename ../geometry/asp.pdb -crg_filename ../geometry/jd1.crg -srf_base ../geometry/asp_scaledcharmm_ -srf_num 50 -potential_petscspace_degree 1 -k_view output/K_asp_panel_50.mat -use_panels 1

aspmat_noprobe_point:
	./calcSolvationEnergy -is_sphere 0 -pdb_filename ../geometry/asp.pdb -crg_filename ../geometry/jd1.crg -srf_base ../geometry/asp_scaledcharmm_noprobe_ -srf_num 1 -potential_petscspace_degree 1 -k_view output/K_asp_noprobe_point_1.mat
	./calcSolvationEnergy -is_sphere 0 -pdb_filename ../geometry/asp.pdb -crg_filename ../geometry/jd1.crg -srf_base ../geometry/asp_scaledcharmm_noprobe_ -srf_num 2 -potential_petscspace_degree 1 -k_view output/K_asp_noprobe_point_2.mat
	./calcSolvationEnergy -is_sphere 0 -pdb_filename ../geometry/asp.pdb -crg_filename ../geometry/jd1.crg -srf_base ../geometry/asp_scaledcharmm_noprobe_ -srf_num 3 -potential_petscspace_degree 1 -k_view output/K_asp_noprobe_point_3.mat
	./calcSolvationEnergy -is_sphere 0 -pdb_filename ../geometry/asp.pdb -crg_filename ../geometry/jd1.crg -srf_base ../geometry/asp_scaledcharmm_noprobe_ -srf_num 4 -potential_petscspace_degree 1 -k_view output/K_asp_noprobe_point_4.mat
	./calcSolvationEnergy -is_sphere 0 -pdb_filename ../geometry/asp.pdb -crg_filename ../geometry/jd1.crg -srf_base ../geometry/asp_scaledcharmm_noprobe_ -srf_num 5 -potential_petscspace_degree 1 -k_view output/K_asp_noprobe_point_5.mat
	./calcSolvationEnergy -is_sphere 0 -pdb_filename ../geometry/asp.pdb -crg_filename ../geometry/jd1.crg -srf_base ../geometry/asp_scaledcharmm_noprobe_ -srf_num 6 -potential_petscspace_degree 1 -k_view output/K_asp_noprobe_point_6.mat
	./calcSolvationEnergy -is_sphere 0 -pdb_filename ../geometry/asp.pdb -crg_filename ../geometry/jd1.crg -srf_base ../geometry/asp_scaledcharmm_noprobe_ -srf_num 7 -potential_petscspace_degree 1 -k_view output/K_asp_noprobe_point_7.mat
	./calcSolvationEnergy -is_sphere 0 -pdb_filename ../geometry/asp.pdb -crg_filename ../geometry/jd1.crg -srf_base ../geometry/asp_scaledcharmm_noprobe_ -srf_num 8 -potential_petscspace_degree 1 -k_view output/K_asp_noprobe_point_8.mat
	./calcSolvationEnergy -is_sphere 0 -pdb_filename ../geometry/asp.pdb -crg_filename ../geometry/jd1.crg -srf_base ../geometry/asp_scaledcharmm_noprobe_ -srf_num 12 -potential_petscspace_degree 1 -k_view output/K_asp_noprobe_point_12.mat
	./calcSolvationEnergy -is_sphere 0 -pdb_filename ../geometry/asp.pdb -crg_filename ../geometry/jd1.crg -srf_base ../geometry/asp_scaledcharmm_noprobe_ -srf_num 16 -potential_petscspace_degree 1 -k_view output/K_asp_noprobe_point_16.mat
#	./calcSolvationEnergy -is_sphere 0 -pdb_filename ../geometry/asp.pdb -crg_filename ../geometry/jd1.crg -srf_base ../geometry/asp_scaledcharmm_noprobe_ -srf_num 20 -potential_petscspace_degree 1 -k_view output/K_asp_noprobe_point_20.mat
#	./calcSolvationEnergy -is_sphere 0 -pdb_filename ../geometry/asp.pdb -crg_filename ../geometry/jd1.crg -srf_base ../geometry/asp_scaledcharmm_noprobe_ -srf_num 24 -potential_petscspace_degree 1 -k_view output/K_asp_noprobe_point_24.mat
#	./calcSolvationEnergy -is_sphere 0 -pdb_filename ../geometry/asp.pdb -crg_filename ../geometry/jd1.crg -srf_base ../geometry/asp_scaledcharmm_noprobe_ -srf_num 28 -potential_petscspace_degree 1 -k_view output/K_asp_noprobe_point_28.mat
#	./calcSolvationEnergy -is_sphere 0 -pdb_filename ../geometry/asp.pdb -crg_filename ../geometry/jd1.crg -srf_base ../geometry/asp_scaledcharmm_noprobe_ -srf_num 32 -potential_petscspace_degree 1 -k_view output/K_asp_noprobe_point_32.mat
#	./calcSolvationEnergy -is_sphere 0 -pdb_filename ../geometry/asp.pdb -crg_filename ../geometry/jd1.crg -srf_base ../geometry/asp_scaledcharmm_noprobe_ -srf_num 50 -potential_petscspace_degree 1 -k_view output/K_asp_noprobe_point_50.mat

aspmat_noprobe_panel:
	./calcSolvationEnergy -is_sphere 0 -pdb_filename ../geometry/asp.pdb -crg_filename ../geometry/jd1.crg -srf_base ../geometry/asp_scaledcharmm_noprobe_ -srf_num 1 -potential_petscspace_degree 1 -k_view output/K_asp_noprobe_panel_1.mat -use_panels 1
	./calcSolvationEnergy -is_sphere 0 -pdb_filename ../geometry/asp.pdb -crg_filename ../geometry/jd1.crg -srf_base ../geometry/asp_scaledcharmm_noprobe_ -srf_num 2 -potential_petscspace_degree 1 -k_view output/K_asp_noprobe_panel_2.mat -use_panels 1
	./calcSolvationEnergy -is_sphere 0 -pdb_filename ../geometry/asp.pdb -crg_filename ../geometry/jd1.crg -srf_base ../geometry/asp_scaledcharmm_noprobe_ -srf_num 3 -potential_petscspace_degree 1 -k_view output/K_asp_noprobe_panel_3.mat -use_panels 1
	./calcSolvationEnergy -is_sphere 0 -pdb_filename ../geometry/asp.pdb -crg_filename ../geometry/jd1.crg -srf_base ../geometry/asp_scaledcharmm_noprobe_ -srf_num 4 -potential_petscspace_degree 1 -k_view output/K_asp_noprobe_panel_4.mat -use_panels 1
	./calcSolvationEnergy -is_sphere 0 -pdb_filename ../geometry/asp.pdb -crg_filename ../geometry/jd1.crg -srf_base ../geometry/asp_scaledcharmm_noprobe_ -srf_num 5 -potential_petscspace_degree 1 -k_view output/K_asp_noprobe_panel_5.mat -use_panels 1
	./calcSolvationEnergy -is_sphere 0 -pdb_filename ../geometry/asp.pdb -crg_filename ../geometry/jd1.crg -srf_base ../geometry/asp_scaledcharmm_noprobe_ -srf_num 6 -potential_petscspace_degree 1 -k_view output/K_asp_noprobe_panel_6.mat -use_panels 1
	./calcSolvationEnergy -is_sphere 0 -pdb_filename ../geometry/asp.pdb -crg_filename ../geometry/jd1.crg -srf_base ../geometry/asp_scaledcharmm_noprobe_ -srf_num 7 -potential_petscspace_degree 1 -k_view output/K_asp_noprobe_panel_7.mat -use_panels 1
	./calcSolvationEnergy -is_sphere 0 -pdb_filename ../geometry/asp.pdb -crg_filename ../geometry/jd1.crg -srf_base ../geometry/asp_scaledcharmm_noprobe_ -srf_num 8 -potential_petscspace_degree 1 -k_view output/K_asp_noprobe_panel_8.mat -use_panels 1
	./calcSolvationEnergy -is_sphere 0 -pdb_filename ../geometry/asp.pdb -crg_filename ../geometry/jd1.crg -srf_base ../geometry/asp_scaledcharmm_noprobe_ -srf_num 12 -potential_petscspace_degree 1 -k_view output/K_asp_noprobe_panel_12.mat -use_panels 1
	./calcSolvationEnergy -is_sphere 0 -pdb_filename ../geometry/asp.pdb -crg_filename ../geometry/jd1.crg -srf_base ../geometry/asp_scaledcharmm_noprobe_ -srf_num 16 -potential_petscspace_degree 1 -k_view output/K_asp_noprobe_panel_16.mat -use_panels 1
#	./calcSolvationEnergy -is_sphere 0 -pdb_filename ../geometry/asp.pdb -crg_filename ../geometry/jd1.crg -srf_base ../geometry/asp_scaledcharmm_noprobe_ -srf_num 20 -potential_petscspace_degree 1 -k_view output/K_asp_noprobe_panel_20.mat -use_panels 1
#	./calcSolvationEnergy -is_sphere 0 -pdb_filename ../geometry/asp.pdb -crg_filename ../geometry/jd1.crg -srf_base ../geometry/asp_scaledcharmm_noprobe_ -srf_num 24 -potential_petscspace_degree 1 -k_view output/K_asp_noprobe_panel_24.mat -use_panels 1
#	./calcSolvationEnergy -is_sphere 0 -pdb_filename ../geometry/asp.pdb -crg_filename ../geometry/jd1.crg -srf_base ../geometry/asp_scaledcharmm_noprobe_ -srf_num 28 -potential_petscspace_degree 1 -k_view output/K_asp_noprobe_panel_28.mat -use_panels 1
#	./calcSolvationEnergy -is_sphere 0 -pdb_filename ../geometry/asp.pdb -crg_filename ../geometry/jd1.crg -srf_base ../geometry/asp_scaledcharmm_noprobe_ -srf_num 32 -potential_petscspace_degree 1 -k_view output/K_asp_noprobe_panel_32.mat -use_panels 1
#	./calcSolvationEnergy -is_sphere 0 -pdb_filename ../geometry/asp.pdb -crg_filename ../geometry/jd1.crg -srf_base ../geometry/asp_scaledcharmm_noprobe_ -srf_num 50 -potential_petscspace_degree 1 -k_view output/K_asp_noprobe_panel_50.mat -use_panels 1


matsph: spheremat_point spheremat_panel spheremat_noprobe_point spheremat_noprobe_panel

matarg: argmat_point argmat_panel argmat_noprobe_point argmat_noprobe_panel

matasp: aspmat_point aspmat_panel aspmat_noprobe_point aspmat_noprobe_panel

mats: matsph matarg matasp



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
