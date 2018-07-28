#if !defined(__SPHERE_H)
#define __SPHERE_H

PetscErrorCode convertToSpherical(const PetscScalar xyz[], PetscScalar r[]);
PetscErrorCode legendre2(PetscInt l, PetscInt nz, PetscScalar z, PetscScalar leg[]);
PetscErrorCode computeEnm(PetscReal b, PetscReal epsIn, PQRData *pqr, Vec qVec, PetscInt Nmax, Vec Enm);
PetscErrorCode computeBnm(PetscReal b, PetscReal epsIn, PetscReal epsOut, PetscInt Nmax, Vec Enm, Vec Bnm);
PetscErrorCode legendre(PetscInt l, PetscInt m, PetscScalar x, PetscScalar *leg, PetscScalar *err);
PetscErrorCode computePotentialSpherical(PQRData *pqr, PetscInt Nmax, Vec Bnm, Vec phi);

#endif
