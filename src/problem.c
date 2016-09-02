#include <petsc.h>
#include <petscmat.h>
#include <petsc/private/dmpleximpl.h>
#include "constants.h"
#include "surface.h"
#include "ellipsoid.h"
#include "problem.h"

/*

DistToOrigin - Calculates distance to origin from surface of ellipsoid in spherical coordinates

*/
#undef __FUNCT__
#define __FUNCT__ "DistToOrigin"
PetscErrorCode DistToOrigin(Ellipsoid *ell, PetscReal theta, PetscReal phi, PetscReal *dist)
{

  PetscFunctionBegin;
  
  PetscReal cos2theta = PetscCosReal(theta)*PetscCosReal(theta);
  PetscReal sin2theta = PetscSinReal(theta)*PetscSinReal(theta);
  PetscReal cos2phi   = PetscCosReal(phi)*PetscCosReal(phi);
  PetscReal sin2phi   = PetscSinReal(phi)*PetscSinReal(phi);
 
  *dist = ell->a*ell->a*sin2theta*cos2phi +
    ell->b*ell->b*sin2theta*sin2phi + 
    ell->c*ell->c*cos2theta;
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "CalcInteractionHessian"
PetscErrorCode CalcInteractionHessian(Tao tao, Vec x, Mat H, Mat Hpre, InteractionContext *ctx)
{
  const PetscReal del = .001;
  const PetscReal del2 = 2*del;
  Vec x1, x2;
  Vec df1, df2;
  PetscReal *vals;
  PetscInt ndim;
  //PetscBool assembled;
  PetscErrorCode ierr;
  PetscFunctionBegin;

  //initialize hessian entries
  //ierr = MatAssembled(H, &assembled);
  //if(assembled) { MatZeroEntries(H); }

  MatZeroEntries(H);
  
  //get dimension of x
  ierr = VecGetSize(x, &ndim); CHKERRQ(ierr);

  //copy x to x1 and x2
  ierr = VecDuplicate(x, &x1); CHKERRQ(ierr);
  ierr = VecDuplicate(x, &x2); CHKERRQ(ierr);
  ierr = VecCopy(x, x1); CHKERRQ(ierr);
  ierr = VecCopy(x, x2); CHKERRQ(ierr);  

  //initialize gradient vectors
  ierr = VecDuplicate(x, &df1); CHKERRQ(ierr);
  ierr = VecDuplicate(x, &df2); CHKERRQ(ierr);

  //loop over columns
  for(PetscInt j=0; j<ndim; ++j) {
    //add/subtract del to row row j
    ierr = VecSetValue(x1, j, -del, ADD_VALUES); CHKERRQ(ierr);
    ierr = VecSetValue(x2, j,  del, ADD_VALUES); CHKERRQ(ierr);
    ierr = VecAssemblyBegin(x1); CHKERRQ(ierr); ierr = VecAssemblyEnd(x1); CHKERRQ(ierr);
    ierr = VecAssemblyBegin(x2); CHKERRQ(ierr); ierr = VecAssemblyEnd(x2); CHKERRQ(ierr);
    //calculate gradient for x1, x2
    ierr = CalcInteractionObjectiveGradient(tao, x1, NULL, df1, ctx); CHKERRQ(ierr);
    ierr = CalcInteractionObjectiveGradient(tao, x2, NULL, df2, ctx); CHKERRQ(ierr);
    ierr = VecAXPY(x2, -1.0, x1); CHKERRQ(ierr);
    ierr = VecScale(x2, 1./del2); CHKERRQ(ierr);
    
    ierr = VecGetArray(x2, &vals); CHKERRQ(ierr);
    //loop over rows and enter into matrix
    for(PetscInt i=0; i<ndim; ++i) {
      ierr = MatSetValue(H, i, j, vals[i], INSERT_VALUES); CHKERRQ(ierr);
      ierr = MatSetValue(Hpre, i, j, vals[i], INSERT_VALUES); CHKERRQ(ierr);
    }

  }
  ierr = MatAssemblyBegin(H, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyBegin(Hpre, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(H, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(Hpre, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  
  ierr = VecDestroy(&x1); CHKERRQ(ierr);
  ierr = VecDestroy(&x2); CHKERRQ(ierr);
  ierr = VecDestroy(&df1); CHKERRQ(ierr);
  ierr = VecDestroy(&df2); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}


/*
Calculates the function and the gradient at position x
if g is NULL then only the function is evaluated

 */
#undef __FUNCT__
#define __FUNCT__ "CalcInteractionObjectiveGradient"
PetscErrorCode CalcInteractionObjectiveGradient(Tao tao, Vec x, PetscReal *f, Vec g, InteractionContext *ctx)
{
  const PetscReal del = .001;
  PetscReal evals[2*9];
  PetscReal storeVal;
  PetscReal inter;
  PetscReal der;
  PetscErrorCode ierr;
  Vec xLocal;
  PetscFunctionBegin;

  //calculate function value
  if(f) {
    ierr = EllipsoidInteractionInterface(tao, x, f, ctx); CHKERRQ(ierr);
  }
  
  if(g) {

    
    ierr = VecDuplicate(x, &xLocal); CHKERRQ(ierr);
    ierr = VecCopy(x, xLocal); CHKERRQ(ierr);
    
    PetscReal stencil[2] = {-del, 2*del}; //stencil
    
    for(PetscInt i=0; i<9; ++i) { //over each ellipsoid parameter
      ierr = VecGetValues(xLocal, 1, &i, &storeVal); CHKERRQ(ierr);
      for(PetscInt j=0; j<2; ++j) { //finite difference evaluations
	
	ierr = VecSetValues(xLocal, 1, &i, stencil+j, ADD_VALUES); CHKERRQ(ierr);
	ierr = VecAssemblyBegin(xLocal); CHKERRQ(ierr);
	ierr = VecAssemblyEnd(xLocal); CHKERRQ(ierr);
	ierr = EllipsoidInteractionInterface(tao, xLocal, &inter, ctx); CHKERRQ(ierr);
	evals[2*i+j] = inter;
	
      }
      //restore x values
      ierr = VecSetValues(xLocal, 1, &i, &storeVal, INSERT_VALUES); CHKERRQ(ierr);
      ierr = VecAssemblyBegin(xLocal); CHKERRQ(ierr);
      ierr = VecAssemblyEnd(xLocal); CHKERRQ(ierr);
      
      der = (evals[i*2+1] - evals[i*2])/(2.0*del);
      ierr = VecSetValue(g, i, der, INSERT_VALUES); CHKERRQ(ierr);
    }
    
    ierr = VecAssemblyBegin(g); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(g); CHKERRQ(ierr);
  }

  PetscFunctionReturn(0);
}

/*
  x[0] = a
  x[1] = b
  x[2] = c
  x[3] = translation[0]
  x[4] = translation[1]
  x[5] = translation[2]
  x[6] = rotation[0]
  x[7] = rotation[1]
  x[8] = rotation[2]
  
 */
#undef __FUNCT__
#define __FUNCT__ "EllipsoidInteractionInterface"
PetscErrorCode EllipsoidInteractionInterface(Tao tao, Vec x, PetscReal *val, InteractionContext *ctx)
{
  Ellipsoid locEll;
  PetscErrorCode ierr;
  PetscFunctionBegin;
  PetscInt i = 0;
  ierr = VecGetValues(x, 1, &i, &(locEll.a)); CHKERRQ(ierr); i++;
  ierr = VecGetValues(x, 1, &i, &(locEll.b)); CHKERRQ(ierr); i++;
  ierr = VecGetValues(x, 1, &i, &(locEll.c)); CHKERRQ(ierr); i++;
  ierr = VecGetValues(x, 1, &i, &(locEll.origin[0])); CHKERRQ(ierr); i++;
  ierr = VecGetValues(x, 1, &i, &(locEll.origin[1])); CHKERRQ(ierr); i++;
  ierr = VecGetValues(x, 1, &i, &(locEll.origin[2])); CHKERRQ(ierr); i++;
  ierr = VecGetValues(x, 1, &i, &(locEll.rotation[0])); CHKERRQ(ierr); i++;
  ierr = VecGetValues(x, 1, &i, &(locEll.rotation[1])); CHKERRQ(ierr); i++;
  ierr = VecGetValues(x, 1, &i, &(locEll.rotation[2])); CHKERRQ(ierr); i++;

  ierr = CalcEllipsoidInteraction(&locEll, ctx->pqr, val); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "CalcEllipsoidInteraction"
PetscErrorCode CalcEllipsoidInteraction(Ellipsoid *ell, PQRData *pqr, PetscReal *val)
{
  Mat         P, Ptemp;
  Mat         rotX, rotY, rotZ;
  Vec         xyzCopy;
  PetscScalar *data;
  PetscInt    size, npts;
  PetscViewer out;
  PetscErrorCode ierr;

  PetscFunctionBegin;

  //form matrix with xyz data in columns
  ierr = VecGetSize(pqr->xyz, &size); CHKERRQ(ierr);
  npts = size/3;
  ierr = VecDuplicate(pqr->xyz, &xyzCopy); CHKERRQ(ierr);
  ierr = VecCopy(pqr->xyz, xyzCopy); CHKERRQ(ierr);
  ierr = VecGetArray(xyzCopy, &data); CHKERRQ(ierr);
  ierr = MatCreateSeqDense(PETSC_COMM_SELF, 3, npts, data, &P); CHKERRQ(ierr);

  
  //output P matrix before any transformations
  ierr = PetscViewerASCIIOpen(PETSC_COMM_SELF, "out1.txt", &out); CHKERRQ(ierr);
  ierr = MatView(P, out); CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&out);

  PetscReal t1, t2, t3;
  t1 = ell->rotation[0];
  t2 = ell->rotation[1];
  t3 = ell->rotation[2];
  PetscScalar rotXvals[9] = {1,       0,        0,
			     0,  PetscCosReal(t1), PetscSinReal(t1),
			     0, -PetscSinReal(t1), PetscCosReal(t1) };
  PetscScalar rotYvals[9] = {PetscCosReal(t2),  0, -PetscSinReal(t2),
			     0,        1,        0,
			     PetscSinReal(t2),  0,  PetscCosReal(t2) };
  PetscScalar rotZvals[9] = {PetscCosReal(t3) , PetscSinReal(t3), 0,
			     -PetscSinReal(t3), PetscCosReal(t3), 0,
			     0,        0,       1 };


  //create rotation matrices
  ierr = MatCreateSeqDense(PETSC_COMM_SELF, 3, 3, &rotXvals[0], &rotX); CHKERRQ(ierr);
  ierr = MatCreateSeqDense(PETSC_COMM_SELF, 3, 3, &rotYvals[0], &rotY); CHKERRQ(ierr);
  ierr = MatCreateSeqDense(PETSC_COMM_SELF, 3, 3, &rotZvals[0], &rotZ); CHKERRQ(ierr);
  

  //translate all points
  PetscInt rows[3] = {0, 1, 2};
  ierr = MatAssemblyBegin(P, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  for(PetscInt i=0; i<npts; ++i)
    ierr = MatSetValues(P, 3, rows, 1, &i, &ell->origin[0], ADD_VALUES); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(P, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

  //rotate all points
  ierr = MatDuplicate(P, MAT_COPY_VALUES, &Ptemp); CHKERRQ(ierr);
  ierr = MatMatMult(rotX, Ptemp, MAT_REUSE_MATRIX, PETSC_DEFAULT, &P); CHKERRQ(ierr);     //x-rotation first
  ierr = MatMatMult(rotY, P    , MAT_REUSE_MATRIX, PETSC_DEFAULT, &Ptemp); CHKERRQ(ierr); //y-rotation second
  ierr = MatMatMult(rotZ, Ptemp, MAT_REUSE_MATRIX, PETSC_DEFAULT, &P); CHKERRQ(ierr);     //z-rotation last
  
  ierr = PetscViewerASCIIOpen(PETSC_COMM_SELF, "out2.txt", &out); CHKERRQ(ierr);
  ierr = MatView(P, out); CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&out); CHKERRQ(ierr);
  

  /////
  //INTERACTION ENERGY - using variables similar to Ford, Wang publication
  ////
  PetscReal rj;     //radii parameters from Allinger
  PetscReal *ri;
  PetscReal epsj; //energy parameters from Allinger
  PetscReal *epsi;
  PetscReal theta, phi;
  PetscReal ellX, ellY, ellZ;
  PetscReal solventX, solventY, solventZ;
  PetscReal *soluteX, *soluteY, *soluteZ;
  PetscReal s;
  PetscReal ellA, ellB, ellC;
  PetscReal ellA4, ellB4, ellC4;
  PetscReal Aconst = 2.9e5;
  PetscReal Bconst = 12.50;
  PetscReal Cconst = 2.25;
  rj = 1.82;    //coefficients for oxygen molecule, ignoring hydrogen in water
  epsj = 0.059; //

  ellA = ell->a; ellB = ell->b; ellC = ell->c;
  ellA4 = PetscPowReal(ellA, 4);
  ellB4 = PetscPowReal(ellB, 4);
  ellC4 = PetscPowReal(ellC, 4);
  //PetscInt nT = (b-a)/(nT-1);
  PetscInt nT = 10; // FOR NOW nP must be 2*nT due to how delS is calculated
  PetscInt nP = 2*nT; //
  PetscReal delT = PETSC_PI/(nT);
  PetscReal delP = 2*PETSC_PI/(nP);

  //Matrix containing distance from ellipsoid surface to origin for each m,n
  Mat R;
  ierr = MatCreateSeqDense(PETSC_COMM_SELF, nT, nP, NULL, &R); CHKERRQ(ierr);

  ierr = PetscMalloc1(npts, &soluteX); CHKERRQ(ierr);
  ierr = PetscMalloc1(npts, &soluteY); CHKERRQ(ierr);
  ierr = PetscMalloc1(npts, &soluteZ); CHKERRQ(ierr);
  ierr = PetscMalloc1(npts, &epsi); CHKERRQ(ierr);
  ierr = PetscMalloc1(npts, &ri); CHKERRQ(ierr);

  PetscReal V = 0;
  PetscReal surfA = 0;
  //loop over theta angles
  for(PetscInt m=0; m<nT; ++m) {
    theta = delT*m;
    //loop over phi angles
    for(PetscInt n=0; n<nP; ++n) {
      phi = delP*n;
      //calculate point on ellipsoid for theta, phi
      ellX = ellA*sin(theta)*cos(phi);
      ellY = ellB*sin(theta)*sin(phi);
      ellZ = ellC*cos(theta);
      s = rj/((ellX*ellX/ellA4) + (ellY*ellY/ellB4) + (ellZ*ellZ/ellC4));
      solventX = ellX*(1 + s/(ellA*ellA));
      solventY = ellY*(1 + s/(ellB*ellB));
      solventZ = ellZ*(1 + s/(ellC*ellC));
      PetscReal rmn;
      ierr = DistToOrigin(ell, theta, phi, &rmn); CHKERRQ(ierr);
      //store distance to origin in R matrix
      ierr = MatSetValue(R, m, n, rmn, INSERT_VALUES); CHKERRQ(ierr);
      ierr = MatAssemblyBegin(R, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
      ierr = MatAssemblyEnd(R, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
     
      PetscReal rm1n, rmn1;  //previous radius values for calculating delS
      PetscInt m1 = m-1;
      PetscInt n1 = n-1;
      if(m==0 || n==0) {
	ierr = DistToOrigin(ell, delT*(m1), delP*n, &rm1n); CHKERRQ(ierr);
	ierr = DistToOrigin(ell, delT*m, delP*n1, &rmn1); CHKERRQ(ierr);
      }
      else {
	ierr = MatGetValues(R, 1, &m1, 1, &n, &rm1n); CHKERRQ(ierr);
	ierr = MatGetValues(R, 1, &m, 1, &n1, &rmn1); CHKERRQ(ierr);
      }

      //calculate delS
      PetscReal delS = PetscSqrtReal(rmn*rmn + rm1n*rm1n - 2*rmn*rm1n*cos(delT))
	* PetscSqrtReal(rmn*rmn + rmn1*rmn1 - 2*rmn*rmn1*cos(delT));

      //loop over solute molecules
      for(PetscInt i=0; i<npts; ++i) {

	if(m==0 && n==0) {
	  PetscInt val = 0;
	  ierr = MatGetValues(P, 1, &val, 1, &i, soluteX+i); CHKERRQ(ierr); val++;
	  ierr = MatGetValues(P, 1, &val, 1, &i, soluteY+i); CHKERRQ(ierr); val++;
	  ierr = MatGetValues(P, 1, &val, 1, &i, soluteZ+i); CHKERRQ(ierr);
	  
	  ierr = VecGetValues(pqr->MM3rad, 1, &i, ri+i); CHKERRQ(ierr);
	  ierr = VecGetValues(pqr->MM3eps, 1, &i, epsi+i); CHKERRQ(ierr);

	}
	//get solute radius and energy parameters

	PetscReal eij = PetscSqrtReal(epsi[i]*epsj);
	PetscReal rij = PetscSqrtReal((solventX - soluteX[i])*(solventX - soluteX[i]) + (solventY - soluteY[i])*(solventY - soluteY[i]) + (solventZ - soluteZ[i])*(solventZ - soluteZ[i]));
	PetscReal Pij = (ri[i] + rj)/rij;
	PetscReal Pij6 = Pij*Pij*Pij*Pij*Pij*Pij;
	PetscReal Vij = eij*(Aconst*PetscExpReal(-Bconst/Pij) - Cconst*Pij6);
	V += Vij*delS;
	surfA += delS;

	
	
      } 
    } 
      
  }
  *val = V/surfA;
  ierr = VecRestoreArray(xyzCopy, &data); CHKERRQ(ierr);
  ierr = MatDestroy(&P); CHKERRQ(ierr);
  ierr = MatDestroy(&R); CHKERRQ(ierr);
  ierr = MatDestroy(&Ptemp); CHKERRQ(ierr);
  ierr = MatDestroy(&rotX); CHKERRQ(ierr);
  ierr = MatDestroy(&rotY); CHKERRQ(ierr);
  ierr = MatDestroy(&rotZ); CHKERRQ(ierr);

  ierr = VecDestroy(&xyzCopy); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ProcessOptions"
PetscErrorCode ProcessOptions(MPI_Comm comm, SolvationContext *ctx)
{
  PetscErrorCode ierr;

  PetscFunctionBeginUser;
  ctx->epsIn      = 4;
  ctx->epsOut     = 80;
  ctx->srfNum     = 1;
  ctx->isSphere   = PETSC_TRUE;
  ctx->R          = 6.0;
  ctx->origin[0]  = 0.0;
  ctx->origin[1]  = 0.0;
  ctx->origin[2]  = 0.0;
  ctx->numCharges = 100;
  ctx->h          = 1.0;
  ctx->Nmax       = 100;
  ctx->density    = 1.0;


  ierr = PetscStrcpy(ctx->basename, "../geometry/sphere_R6_vdens");CHKERRQ(ierr);
  ierr = PetscOptionsBegin(comm, "", "Solvation Problem Options", "BIBEE");CHKERRQ(ierr);
    ierr = PetscOptionsReal("-epsilon_solute", "The dielectric coefficient of the solute", "testSrfOnSurfacePoints", ctx->epsIn, &ctx->epsIn, NULL);CHKERRQ(ierr);
    ierr = PetscOptionsReal("-epsilon_solvent", "The dielectric coefficient of the solvent", "testSrfOnSurfacePoints", ctx->epsOut, &ctx->epsOut, NULL);CHKERRQ(ierr);
    ierr = PetscOptionsString("-pdb_filename", "The filename for the .pdb file", "testSrfOnSurfacePoints", ctx->pdbFile, ctx->pdbFile, sizeof(ctx->pdbFile), NULL);CHKERRQ(ierr);
    ierr = PetscOptionsString("-pdb_filename", "The filename for the .pdb file", "testSrfOnSurfacePoints", ctx->pdbFile, ctx->pdbFile, sizeof(ctx->pdbFile), NULL); CHKERRQ(ierr);
    ierr = PetscOptionsString("-crg_filename", "The filename for the .crg file", "testSrfOnSurfacePoints", ctx->crgFile, ctx->crgFile, sizeof(ctx->crgFile), NULL); CHKERRQ(ierr);
    ierr = PetscOptionsString("-MM3_filename", "The filename for the .mm3 file", "testSrfOnSurfacePoints", ctx->MM3File, ctx->MM3File, sizeof(ctx->MM3File), NULL); CHKERRQ(ierr);
    ierr = PetscOptionsBool("-is_sphere", "Use a spherical test case", "testSrfOnSurfacePoints", ctx->isSphere, &ctx->isSphere, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsBool("-do_ellipsoid", "Form approximate ellipsoid, do computation on it", "testSrfOnSurfacePoints", ctx->do_ellipsoid, &ctx->do_ellipsoid, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsInt("-num_charges", "The number of atomic charges in the solute", "testSrfOnSurfacePoints", ctx->numCharges, &ctx->numCharges, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsString("-srf_base", "The basename for the .srf file", "testSrfOnSurfacePoints", ctx->basename, ctx->basename, sizeof(ctx->basename), NULL); CHKERRQ(ierr);
    ierr = PetscOptionsInt("-srf_num", "The resolution number of the mesh", "testSrfOnSurfacePoints", ctx->srfNum, &ctx->srfNum, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsInt("-nmax", "The order of the multipole expansion", "testSrfOnSurfacePoints", ctx->Nmax, &ctx->Nmax, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsReal("-density", "The density of points for BEM", "testSrfOnSurfacePoints", ctx->density, &ctx->density, NULL); CHKERRQ(ierr);
  ierr = PetscOptionsEnd();

  ierr = PetscSNPrintf(ctx->srfFile, PETSC_MAX_PATH_LEN-1, "%s%d.srf", ctx->basename, (int) ctx->srfNum);CHKERRQ(ierr);
  ierr = PetscSNPrintf(ctx->pntFile, PETSC_MAX_PATH_LEN-1, "%s%d.pnt", ctx->basename, (int) ctx->srfNum);CHKERRQ(ierr);

  ierr = PetscLogEventRegister("IntegratePanel",   DM_CLASSID, &IntegratePanel_Event);CHKERRQ(ierr);
  ierr = PetscLogEventRegister("CalcSurfToSurf",   DM_CLASSID, &CalcStoS_Event);CHKERRQ(ierr);
  ierr = PetscLogEventRegister("CalcSurfToCharge", DM_CLASSID, &CalcStoQ_Event);CHKERRQ(ierr);
  ierr = PetscLogEventRegister("CalcLMatrix",      DM_CLASSID, &CalcL_Event);CHKERRQ(ierr);
  ierr = PetscLogEventRegister("CalcReactPot",     DM_CLASSID, &CalcR_Event);CHKERRQ(ierr);
  ierr = PetscLogEventRegister("CalcSolvEnergy",   DM_CLASSID, &CalcE_Event);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "PQRCreateFromPDB"
PetscErrorCode PQRCreateFromPDB(MPI_Comm comm, const char pdbFile[], const char crgFile[], PetscBool do_ellipsoid, const char MM3File[], PQRData *pqr)
{
  PetscViewer    viewerPDB, viewerCRG, viewerMM3;
  PetscScalar   *q, *x, *MM3rad, *MM3eps;
  PetscReal     *charges, *coords;
  PetscInt       n = 0, i, d;
  PetscMPIInt    rank;
  //PetscBool      do_ellipsoid;
  PetscErrorCode ierr;

  PetscFunctionBeginUser;

  //do_ellipsoid checks whether the MM3 data is provided
  //ierr = PetscOptionsGetBool(NULL, NULL, "-do_ellipsoid", NULL, &do_ellipsoid);

  
  ierr = MPI_Comm_rank(comm, &rank);CHKERRQ(ierr);
  ierr = PetscViewerCreate(comm, &viewerPDB);CHKERRQ(ierr);
  ierr = PetscViewerSetType(viewerPDB, PETSCVIEWERASCII);CHKERRQ(ierr);
  ierr = PetscViewerFileSetMode(viewerPDB, FILE_MODE_READ);CHKERRQ(ierr);
  ierr = PetscViewerFileSetName(viewerPDB, pdbFile);CHKERRQ(ierr);
  if (!rank) {
    char     buf[128];
    PetscInt line = 0, maxSize = 1024, cnt = 1;

    ierr = PetscMalloc1(maxSize, &charges); CHKERRQ(ierr);
    ierr = PetscMalloc1(maxSize*3, &coords); CHKERRQ(ierr);
    while (cnt) {
      PetscInt c = 0;

      /* Read line */
      do {ierr = PetscViewerRead(viewerPDB, &buf[c++], 1, &cnt, PETSC_CHAR);CHKERRQ(ierr);}
      while (buf[c-1] != '\n' && buf[c-1] != '\0' && cnt);
      /* Parse line */
      if (c > 6 &&
          ((buf[0] == 'A' && buf[1] == 'T' && buf[2] == 'O' && buf[3] == 'M') ||
           (buf[0] == 'H' && buf[1] == 'E' && buf[2] == 'T' && buf[3] == 'A' && buf[4] == 'T' && buf[5] == 'M'))) {
        double tmp;

        if (n >= maxSize) {
          /* Reallocate and copy */
        }
        buf[66] = '\0';
        ierr = sscanf(&buf[60], "%lg", &tmp); if (ierr != 1) SETERRQ2(comm, PETSC_ERR_ARG_WRONG, "Could not read charge for line %d of PDB file %s", line, pdbFile);
        charges[n] = tmp;
	buf[54] = '\0';
        ierr = sscanf(&buf[46], "%lg", &tmp); if (ierr != 1) SETERRQ2(comm, PETSC_ERR_ARG_WRONG, "Could not read z coordinate for line %d of PDB file %s", line, pdbFile);
        coords[n*3+2] = tmp;
        buf[46] = '\0';
        ierr = sscanf(&buf[38], "%lg", &tmp); if (ierr != 1) SETERRQ2(comm, PETSC_ERR_ARG_WRONG, "Could not read y coordinate for line %d of PDB file %s", line, pdbFile);
        coords[n*3+1] = tmp;
        buf[38] = '\0';
        ierr = sscanf(&buf[31], "%lg", &tmp); if (ierr != 1) SETERRQ2(comm, PETSC_ERR_ARG_WRONG, "Could not read x coordinate for line %d of PDB file %s", line, pdbFile);
        coords[n*3+0] = tmp;
        /* Residue id [23-27] */
        /* Segment id [21-22] */
        /* Residue name [17-20] */
        /* Atom name [12-15] */
        ++n;
      }
      ++line;
    }
  }
  ierr = PetscViewerDestroy(&viewerPDB);CHKERRQ(ierr);

  //create q and xyz vectors
  ierr = VecCreate(comm, &pqr->q);CHKERRQ(ierr);
  ierr = VecSetSizes(pqr->q, n, PETSC_DETERMINE);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) pqr->q, "Atomic Charges");CHKERRQ(ierr);
  ierr = VecSetFromOptions(pqr->q);CHKERRQ(ierr);
  ierr = VecCreate(comm, &pqr->xyz);CHKERRQ(ierr);
  ierr = VecSetSizes(pqr->xyz, n*3, PETSC_DETERMINE);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) pqr->xyz, "Atomic XYZ");CHKERRQ(ierr);
  ierr = VecSetBlockSize(pqr->xyz, 3);CHKERRQ(ierr);
  ierr = VecSetFromOptions(pqr->xyz);CHKERRQ(ierr);
  //create MM3 parameter vectors
  ierr = VecCreate(comm, &pqr->MM3rad); CHKERRQ(ierr);
  ierr = VecCreate(comm, &pqr->MM3eps); CHKERRQ(ierr);
  ierr = VecSetSizes(pqr->MM3rad, n, PETSC_DETERMINE); CHKERRQ(ierr);
  ierr = VecSetSizes(pqr->MM3eps, n, PETSC_DETERMINE); CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) pqr->MM3rad, "MM3 Radius"); CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) pqr->MM3eps, "MM3 Energy Parameter"); CHKERRQ(ierr);
  ierr = VecSetFromOptions(pqr->MM3rad); CHKERRQ(ierr);
  ierr = VecSetFromOptions(pqr->MM3eps); CHKERRQ(ierr);


  //set local pointers to point towards pqr arrays
  ierr = VecGetArray(pqr->q, &q);CHKERRQ(ierr);
  ierr = VecGetArray(pqr->xyz, &x);CHKERRQ(ierr);

  for (i = 0; i < n; ++i) {
    q[i] = charges[i];
    for (d = 0; d < 3; ++d) x[i*3+d] = coords[i*3+d];
  }
  ierr = VecRestoreArray(pqr->q, &q);CHKERRQ(ierr);
  ierr = VecRestoreArray(pqr->xyz, &x);CHKERRQ(ierr);
  ierr = PetscFree2(charges, coords);CHKERRQ(ierr);

  if (crgFile) {
    ierr = PetscViewerCreate(comm, &viewerCRG);CHKERRQ(ierr);
    ierr = PetscViewerSetType(viewerCRG, PETSCVIEWERASCII);CHKERRQ(ierr);
    ierr = PetscViewerFileSetMode(viewerCRG, FILE_MODE_READ);CHKERRQ(ierr);
    ierr = PetscViewerFileSetName(viewerCRG, crgFile);CHKERRQ(ierr);
    if (!rank) {
      char     buf[128];
      PetscInt cnt = 0;

      /* The CRG file is required to have the same nubmer of atoms in the same order as the PDB */
      ierr = VecGetArray(pqr->q, &q); CHKERRQ(ierr);
      for (i = -1; i <= n; ++i) {
        double    tmp;
        PetscInt  c = 0;

        /* Read line */
        do {ierr = PetscViewerRead(viewerCRG, &buf[c++], 1, &cnt, PETSC_CHAR); CHKERRQ(ierr);}
        while (buf[c-1] != '\n' && buf[c-1] != '\0' && cnt);
        if (!cnt) break;
        if (i < 0) continue;
        buf[22] = '\0';
        ierr = sscanf(&buf[14], "%lg", &tmp); if (ierr != 1) SETERRQ2(comm, PETSC_ERR_ARG_WRONG, "Could not read charge for line %d of CRG file %", i+1, crgFile);
        q[i] = tmp;
	
        /* Segment id [13] */
        /* Residue number [9-12] */
        /* Residue name [6-8] */
        /* Atom name [0-5] */
      }
      ierr = VecRestoreArray(pqr->q, &q); CHKERRQ(ierr);
    }
    ierr = PetscViewerDestroy(&viewerCRG);CHKERRQ(ierr);
  }
  
  if(do_ellipsoid) {
    ierr = PetscViewerCreate(comm, &viewerMM3); CHKERRQ(ierr);
    ierr = PetscViewerSetType(viewerMM3, PETSCVIEWERASCII); CHKERRQ(ierr);
    ierr = PetscViewerFileSetMode(viewerMM3, FILE_MODE_READ); CHKERRQ(ierr);
    ierr = PetscViewerFileSetName(viewerMM3, MM3File); CHKERRQ(ierr);
    if(!rank) {
      char     buf[128];
      PetscInt cnt = 0;
      
      /* The MM3 file is required to have the same nubmer of atoms in the same order as the PDB */
      ierr = VecGetArray(pqr->MM3eps, &MM3eps); CHKERRQ(ierr);
      ierr = VecGetArray(pqr->MM3rad, &MM3rad); CHKERRQ(ierr);
      for (i = 0; i <= n; ++i) {
        double    tmp;
        PetscInt  c = 0;

        /* Read line */
        do {ierr = PetscViewerRead(viewerMM3, &buf[c++], 1, &cnt, PETSC_CHAR); CHKERRQ(ierr);}
        while (buf[c-1] != '\n' && buf[c-1] != '\0' && cnt);
        if (!cnt) break;
        //if (i < 0) continue;
        buf[22] = '\0';
        ierr = sscanf(&buf[5] , "%lg", &tmp); if (ierr != 1) SETERRQ2(comm, PETSC_ERR_ARG_WRONG, "Could not read charge for line %d of MM3 file %", i+1, MM3File);
	MM3rad[i] = tmp;
        ierr = sscanf(&buf[13], "%lg", &tmp); if (ierr != 1) SETERRQ2(comm, PETSC_ERR_ARG_WRONG, "Could not read charge for line %d of MM3 file %", i+1, MM3File);
	MM3eps[i] = tmp;

        //q[i] = tmp;

        /* Segment id [13] */
        /* Residue number [9-12] */
        /* Residue name [6-8] */
        /* Atom name [0-5] */
      }
      ierr = VecRestoreArray(pqr->MM3rad, &MM3rad); CHKERRQ(ierr);
      ierr = VecRestoreArray(pqr->MM3eps, &MM3eps); CHKERRQ(ierr);
    }
    ierr = PetscViewerDestroy(&viewerMM3); CHKERRQ(ierr);
  }
  
  ierr = VecDuplicate(pqr->q, &pqr->R);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) pqr->R, "Atomic radii");CHKERRQ(ierr);
  ierr = VecSet(pqr->R, 0.0);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "PQRViewFromOptions"
PetscErrorCode PQRViewFromOptions(PQRData *pqr)
{
  PetscErrorCode ierr;

  PetscFunctionBeginUser;
  ierr = VecViewFromOptions(pqr->xyz, NULL, "-pqr_vec_view");CHKERRQ(ierr);
  ierr = VecViewFromOptions(pqr->q,   NULL, "-pqr_vec_view");CHKERRQ(ierr);
  ierr = VecViewFromOptions(pqr->R,   NULL, "-pqr_vec_view");CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "PQRDestroy"
PetscErrorCode PQRDestroy(PQRData *pqr)
{
  PetscErrorCode ierr;

  PetscFunctionBeginUser;
  ierr = VecDestroy(&pqr->xyz);CHKERRQ(ierr);
  ierr = VecDestroy(&pqr->q);CHKERRQ(ierr);
  ierr = VecDestroy(&pqr->R);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

