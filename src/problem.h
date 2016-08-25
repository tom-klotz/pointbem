#if !defined(__PROBLEM_H)
#define __PROBLEM_H

typedef enum {BEM_POINT, BEM_PANEL, BEM_POINT_MF, BEM_PANEL_MF} BEMType;

/* Performance characterization */
PetscLogEvent CalcE_Event, CalcL_Event, CalcR_Event, CalcStoQ_Event, CalcStoS_Event, IntegratePanel_Event;

typedef struct {
  Vec q;   /* Charge values */
  Vec xyz; /* Charge coordinates, always 3D */
  Vec R;   /* Charge radii */
} PQRData;

typedef struct {
  /* Physical parameters */
  PetscReal epsIn;      /* solute dielectric coefficient */
  PetscReal epsOut;     /* solvent dielectric coefficient */
  char      pdbFile[PETSC_MAX_PATH_LEN]; /* Chemists are crazy and have never heard of normalized data */
  char      crgFile[PETSC_MAX_PATH_LEN];
  /* Surface file */
  PetscInt  srfNum;     /* Resolution of mesh file */
  char      basename[PETSC_MAX_PATH_LEN];
  char      srfFile[PETSC_MAX_PATH_LEN];
  char      pntFile[PETSC_MAX_PATH_LEN];
  /* Point BEM parameters */
  PetscReal density;    /* Density of points on surface */
  /* Sphere setup */
  PetscBool isSphere;   /* Indicates we are running the sphere test */
  PetscReal R;          /* Sphere radius */
  PetscReal origin[3];  /* Sphere center */
  PetscInt  numCharges; /* Number of atomic charges in the solute */
  PetscReal h;          /* Charge spacing */
  /* Ellipsoid setup */
  PetscBool doEllipsoid;  /* Indicates whethr we want to form approximate ellipsoid */
  PetscReal Ea, Eb, Ec;   /* Lengths of semi-axes of ellipsoid */
  PetscReal Eorigin[3];   /* Center of ellipsoid */
  PetscReal Erotation[3]; /* Axes Rotations */
  /* Analytical parameters */
  PetscInt  Nmax;       /* Order of the multipole expansion */
} SolvationContext;

PETSC_EXTERN PetscErrorCode ProcessOptions(MPI_Comm, SolvationContext*);
PETSC_EXTERN PetscErrorCode PQRCreateFromPDB(MPI_Comm, const char[], const char[], PQRData*);
PETSC_EXTERN PetscErrorCode PQRViewFromOptions(PQRData*);
PETSC_EXTERN PetscErrorCode PQRDestroy(PQRData*);
#endif
