#if !defined(__SURFACE_H)
#define __SURFACE_H

/* TODO Should this just move to the FVM geometry structs? */
typedef struct {
  DM  dm;      /* DM with surface topology and coordinates */
  Vec weights; /* Areas associated with surface pieces */
  Vec normals; /* Normals associated with surface pieces */
} PetscSurface;

PETSC_EXTERN PetscErrorCode DMPlexCreateBardhanFromFile(MPI_Comm, const char[], PetscBool, Vec *, DM *);
PETSC_EXTERN PetscErrorCode DMPlexCreateBardhan(MPI_Comm, PetscViewer, PetscViewer, PetscBool, Vec *, DM *);

PETSC_EXTERN PetscErrorCode loadSrfIntoSurfacePoints(MPI_Comm, const char[], Vec *, Vec *, Vec *, PetscReal *, DM *);
PETSC_EXTERN PetscErrorCode makeSphereSurface(MPI_Comm, PetscReal[], PetscReal, PetscInt, Vec *, Vec *, Vec *, DM *);

PETSC_EXTERN PetscErrorCode PetscSurfaceCreateMSP(MPI_Comm, const char[], PetscSurface *);
PETSC_EXTERN PetscErrorCode PetscSurfaceDestroy(PetscSurface *);
#endif
