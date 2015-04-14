#if !defined(__SURFACE_H)
#define __SURFACE_H

PETSC_EXTERN PetscErrorCode DMPlexCreateBardhanFromFile(MPI_Comm, const char[], PetscBool, Vec *, DM *);
PETSC_EXTERN PetscErrorCode DMPlexCreateBardhan(MPI_Comm, PetscViewer, PetscViewer, PetscBool, Vec *, DM *);

PETSC_EXTERN PetscErrorCode loadSrfIntoSurfacePoints(MPI_Comm, const char[], Vec *, Vec *, Vec *, PetscReal *, DM *);
PETSC_EXTERN PetscErrorCode makeSphereSurface(MPI_Comm, PetscReal[], PetscReal, PetscInt, Vec *, Vec *, Vec *, DM *);
#endif
