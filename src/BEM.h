#if !defined(__BEM_H)
#define __BEM_H

typedef enum {BEM_POINT, BEM_PANEL, BEM_POINT_MF, BEM_PANEL_MF} BEMType;

typedef struct {
  PetscReal alpha;
  PetscReal beta;
  PetscReal gamma;
} HContext;

typedef struct {
  PQRData   *pqr;
  PetscReal epsIn;
  PetscReal epsOut;
  Mat*      B;
  Mat*      K;
  Vec*      Bq;
  Vec*      w;
  HContext* hctx;
} NonlinearContext;

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
  /* Analytical parameters */
  PetscInt  Nmax;       /* Order of the multipole expansion */
  /* SLIC model parameters */
  PetscReal alpha, beta, gamma;
} SolvationContext;


PetscErrorCode doAnalytical(PetscReal b, PetscReal epsIn, PetscReal epsOut, PQRData *pqr, PetscInt Nmax, Mat *L);
PetscErrorCode IntegratePanel(PetscInt numCorners, const PetscReal npanel[], const PetscReal point[], const PetscReal normal[], PetscScalar *fss, PetscScalar *fds, PetscScalar *fess, PetscScalar *feds);
PetscErrorCode makeSurfaceToSurfacePanelOperators_Laplace(DM dm, Vec w, Vec n, Mat *V, Mat *K);
PetscErrorCode makeSurfaceToChargePanelOperators(DM dm, Vec w, Vec n, PQRData *pqr, Mat *potential, Mat *field, Mat *singleLayer, Mat *doubleLayer);
PetscErrorCode makeSurfaceToChargePointOperators(Vec coordinates, Vec w, Vec n, PQRData *pqr, Mat *potential, Mat *field, Mat *singleLayer, Mat *doubleLayer);
PetscErrorCode makeSurfaceToSurfacePointOperators_Laplace(Vec coordinates, Vec w, Vec n, Mat *V, Mat *K);
PetscErrorCode makeBEMPcmQualMatrices(DM dm, BEMType bem, PetscReal epsIn, PetscReal epsOut, PQRData *pqr, Vec coordinates, Vec w, Vec n, Mat *L);
PetscErrorCode ComputeBEMResidual(SNES snes, Vec x, Vec r, void *ctx);
PetscErrorCode ComputeBEMJacobian(SNES snes, Vec x, Mat J, Mat P, void *ctx);
PetscErrorCode nonlinearH(Vec E, HContext *ctx, Vec *hEn);
PetscErrorCode ASCBq(Vec sigma, Vec *Bq, NonlinearContext *ctx);
PetscErrorCode FormASCNonlinearMatrix(Vec sigma, Mat *A, NonlinearContext *ctx);
PetscErrorCode NonlinearPicard(PetscErrorCode (*lhs)(Vec, Mat*, void*), PetscErrorCode (*rhs)(Vec, Vec*, void*), Vec guess, void *ctx, Vec *sol);
PetscErrorCode makeBEMPcmQualReactionPotentialNonlinear(DM dm, BEMType bem, HContext params, PetscReal epsIn, PetscReal epsOut, PQRData *pqr, Vec coordinates, Vec w, Vec n, Vec react);
PetscErrorCode makeBEMPcmQualReactionPotential(DM dm, BEMType bem, PetscReal epsIn, PetscReal epsOut, PQRData *pqr, Vec coordinates, Vec w, Vec n, Vec react);
PetscErrorCode CalculateAnalyticSolvationEnergy(PetscReal epsIn, PetscReal epsOut, PQRData *pqr, PetscReal R, PetscInt Nmax, Vec react, PetscReal *E);
PetscErrorCode CalculateBEMSolvationEnergy(DM dm, const char prefix[], BEMType bem, HContext params, PetscReal epsIn, PetscReal epsOut, PQRData *pqr, Vec w, Vec n, Vec react, PetscReal *E);
PetscErrorCode ProcessOptions(MPI_Comm comm, SolvationContext *ctx);
PetscErrorCode workprectests(int argc, char **argv);
#endif
