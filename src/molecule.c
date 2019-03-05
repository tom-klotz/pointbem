#include <petsc.h>
#include "molecule.h"



#undef __FUNCT__
#define __FUNCT__ "PQRCreateFromPDB"
PetscErrorCode PQRCreateFromPDB(MPI_Comm comm, const char pdbFile[], const char crgFile[], PQRData *pqr)
{
  PetscViewer    viewerPDB, viewerCRG;
  PetscScalar   *q, *x;
  PetscReal     *charges, *coords;
  PetscInt       n = 0, i, d;
  PetscMPIInt    rank;
  PetscErrorCode ierr;

  PetscFunctionBeginUser;
  ierr = MPI_Comm_rank(comm, &rank);CHKERRQ(ierr);
  ierr = PetscViewerCreate(comm, &viewerPDB);CHKERRQ(ierr);
  ierr = PetscViewerSetType(viewerPDB, PETSCVIEWERASCII);CHKERRQ(ierr);
  ierr = PetscViewerFileSetMode(viewerPDB, FILE_MODE_READ);CHKERRQ(ierr);
  ierr = PetscViewerFileSetName(viewerPDB, pdbFile);CHKERRQ(ierr);
  if (!rank) {
    char     buf[128];
    PetscInt line = 0, maxSize = 10000, cnt = 1;

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

  ierr = VecCreate(comm, &pqr->q);CHKERRQ(ierr);
  ierr = VecSetSizes(pqr->q, n, PETSC_DETERMINE);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) pqr->q, "Atomic Charges");CHKERRQ(ierr);
  ierr = VecSetFromOptions(pqr->q);CHKERRQ(ierr);
  ierr = VecCreate(comm, &pqr->xyz);CHKERRQ(ierr);
  ierr = VecSetSizes(pqr->xyz, n*3, PETSC_DETERMINE);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) pqr->xyz, "Atomic XYZ");CHKERRQ(ierr);
  ierr = VecSetBlockSize(pqr->xyz, 3);CHKERRQ(ierr);
  ierr = VecSetFromOptions(pqr->xyz);CHKERRQ(ierr);

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
      PetscInt cnt = 1;

      /* The CRG file is required to have the same nubmer of atoms in the same order as the PDB */
      ierr = VecGetArray(pqr->q, &q);CHKERRQ(ierr);
      for (i = -1; i < n; ++i) {
        double    tmp;
        PetscInt  c = 0;

        /* Read line */
        do {ierr = PetscViewerRead(viewerCRG, &buf[c++], 1, &cnt, PETSC_CHAR);CHKERRQ(ierr);}
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
      ierr = VecRestoreArray(pqr->q, &q);CHKERRQ(ierr);
    }
    ierr = PetscViewerDestroy(&viewerCRG);CHKERRQ(ierr);
  }

  ierr = VecDuplicate(pqr->q, &pqr->R);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) pqr->R, "Atomic radii");CHKERRQ(ierr);
  ierr = VecSet(pqr->R, 0.0);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "PQRCreateFromPQR"
PetscErrorCode PQRCreateFromPQR(MPI_Comm comm, const char pqrFile[], PQRData *pqr)
{
  PetscViewer    viewerPDB, viewerCRG;
  PetscScalar   *q, *x, *r;
  PetscReal     *charges, *coords, *radius;
  PetscInt       n = 0, i, d;
  PetscMPIInt    rank;
  PetscErrorCode ierr;

  PetscFunctionBeginUser;
  ierr = MPI_Comm_rank(comm, &rank);CHKERRQ(ierr);
  ierr = PetscViewerCreate(comm, &viewerPDB);CHKERRQ(ierr);
  ierr = PetscViewerSetType(viewerPDB, PETSCVIEWERASCII);CHKERRQ(ierr);
  ierr = PetscViewerFileSetMode(viewerPDB, FILE_MODE_READ);CHKERRQ(ierr);
  ierr = PetscViewerFileSetName(viewerPDB, pqrFile);CHKERRQ(ierr);
  if (!rank) {
    char     buf[128];
    PetscInt line = 0, maxSize = 10000, cnt = 1;

    ierr = PetscMalloc2(maxSize, &charges, maxSize*3, &coords);CHKERRQ(ierr);
    ierr = PetscMalloc(maxSize, &radius);CHKERRQ(ierr);
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
	buf[69] = '\0';
	ierr = sscanf(&buf[62], "%lg", &tmp); if(ierr != 1) SETERRQ2(comm, PETSC_ERR_ARG_WRONG, "Could not read radius for line %d of PQR file %s", line, pqrFile);
	radius[n] = tmp;
        buf[62] = '\0';
        ierr = sscanf(&buf[54], "%lg", &tmp); if (ierr != 1) SETERRQ2(comm, PETSC_ERR_ARG_WRONG, "Could not read charge for line %d of PQR file %s", line, pqrFile);
        charges[n] = tmp;
        buf[54] = '\0';
        ierr = sscanf(&buf[46], "%lg", &tmp); if (ierr != 1) SETERRQ2(comm, PETSC_ERR_ARG_WRONG, "Could not read z coordinate for line %d of PQR file %s", line, pqrFile);
        coords[n*3+2] = tmp;
        buf[46] = '\0';
        ierr = sscanf(&buf[38], "%lg", &tmp); if (ierr != 1) SETERRQ2(comm, PETSC_ERR_ARG_WRONG, "Could not read y coordinate for line %d of PQR file %s", line, pqrFile);
        coords[n*3+1] = tmp;
        buf[38] = '\0';
        ierr = sscanf(&buf[31], "%lg", &tmp); if (ierr != 1) SETERRQ2(comm, PETSC_ERR_ARG_WRONG, "Could not read x coordinate for line %d of PQR file %s", line, pqrFile);
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

  ierr = VecCreate(comm, &pqr->q);CHKERRQ(ierr);
  ierr = VecSetSizes(pqr->q, n, PETSC_DETERMINE);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) pqr->q, "Atomic Charges");CHKERRQ(ierr);
  ierr = VecSetFromOptions(pqr->q);CHKERRQ(ierr);
  ierr = VecCreate(comm, &pqr->xyz);CHKERRQ(ierr);
  ierr = VecSetSizes(pqr->xyz, n*3, PETSC_DETERMINE);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) pqr->xyz, "Atomic XYZ");CHKERRQ(ierr);
  ierr = VecSetBlockSize(pqr->xyz, 3);CHKERRQ(ierr);
  ierr = VecSetFromOptions(pqr->xyz);CHKERRQ(ierr);
  
  ierr = VecDuplicate(pqr->q, &pqr->R);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) pqr->R, "Atomic radii");CHKERRQ(ierr);
  ierr = VecSet(pqr->R, 0.0);CHKERRQ(ierr);
  ierr = VecGetArray(pqr->q, &q);CHKERRQ(ierr);
  ierr = VecGetArray(pqr->xyz, &x);CHKERRQ(ierr);
  ierr = VecGetArray(pqr->R, &r);CHKERRQ(ierr);
  for (i = 0; i < n; ++i) {
    r[i] = radius[i];

    q[i] = charges[i];
    printf("CHARGE: %5.5f\n", charges[i]);
    for (d = 0; d < 3; ++d) x[i*3+d] = coords[i*3+d];
  }
  ierr = VecRestoreArray(pqr->q, &q);CHKERRQ(ierr);
  ierr = VecRestoreArray(pqr->xyz, &x);CHKERRQ(ierr);
  ierr = VecRestoreArray(pqr->R, &r);CHKERRQ(ierr);
  ierr = PetscFree2(charges, coords);CHKERRQ(ierr);
  ierr = PetscFree(radius);CHKERRQ(ierr);




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
