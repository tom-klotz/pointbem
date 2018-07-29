#if !defined(__MOLECULE_H)
#define __MOLECULE_H

typedef struct {
  Vec q;   /* Charge values */
  Vec xyz; /* Charge coordinates, always 3D */
  Vec R;   /* Charge radii */
} PQRData;

PetscErrorCode PQRCreateFromPDB(MPI_Comm comm, const char pdbFile[], const char crgFile[], PQRData *pqr);
PetscErrorCode PQRViewFromOptions(PQRData *pqr);
PetscErrorCode PQRDestroy(PQRData *pqr);

#endif
