#include <petsc.h>
#include <petsc/private/dmpleximpl.h>
#include "constants.h"
#include "surface.h"
#include "ellipsoid.h"
#include "problem.h"


#undef __FUNCT__
#define __FUNCT__ "CalcEllipsoidInteraction"
PetscErrorCode CalcEllipsoidInteraction(Ellipsoid *ell, PQRData *pqr, PetscReal *val)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;

  ierr = 1;
  *val = 1.0;



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
PetscErrorCode PQRCreateFromPDB(MPI_Comm comm, const char pdbFile[], const char crgFile[], const char MM3File[], PQRData *pqr)
{
  PetscViewer    viewerPDB, viewerCRG, viewerMM3;
  PetscScalar   *q, *x, *MM3rad, *MM3eps;
  PetscReal     *charges, *coords;
  PetscInt       n = 0, i, d;
  PetscMPIInt    rank;
  PetscBool      do_ellipsoid;
  PetscErrorCode ierr;

  PetscFunctionBeginUser;

  //do_ellipsoid checks whether the MM3 data is provided
  ierr = PetscOptionsGetBool(NULL, NULL, "-MM3File", NULL, &do_ellipsoid);

  
  ierr = MPI_Comm_rank(comm, &rank);CHKERRQ(ierr);
  ierr = PetscViewerCreate(comm, &viewerPDB);CHKERRQ(ierr);
  ierr = PetscViewerSetType(viewerPDB, PETSCVIEWERASCII);CHKERRQ(ierr);
  ierr = PetscViewerFileSetMode(viewerPDB, FILE_MODE_READ);CHKERRQ(ierr);
  ierr = PetscViewerFileSetName(viewerPDB, pdbFile);CHKERRQ(ierr);
  if (!rank) {
    char     buf[128];
    PetscInt line = 0, maxSize = 1024, cnt = 1;

    ierr = PetscMalloc2(maxSize, &charges, maxSize*3, &coords);CHKERRQ(ierr);
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

