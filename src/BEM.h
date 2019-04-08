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
  Mat*      A; //left hand side matrix for nonlinear problem
  Mat*      Elec;
  Vec*      Bq;
  Vec*      w;
  Vec      En;
  Vec      hEn;
  HContext* hctx;
} NonlinearContext;

typedef struct {
  /* Physical parameters */
  PetscReal epsIn;      /* solute dielectric coefficient */
  PetscReal epsOut;     /* solvent dielectric coefficient */
  char      pdbFile[PETSC_MAX_PATH_LEN]; /* Chemists are crazy and have never heard of normalized data */
  char      crgFile[PETSC_MAX_PATH_LEN];
  char      pqrFile[PETSC_MAX_PATH_LEN];
  /* Surface file */
  PetscInt  srfNum;     /* Resolution of mesh file */
  char      basename[PETSC_MAX_PATH_LEN];
  char      srfFile[PETSC_MAX_PATH_LEN];
  char      pntFile[PETSC_MAX_PATH_LEN];
  char      kFile[PETSC_MAX_PATH_LEN];
  char      flopsFile[PETSC_MAX_PATH_LEN];
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
  /* User flags to specify solver */
  PetscBool forceNonlinear;
  PetscBool usePanels;
} SolvationContext;


PetscErrorCode outputK(SolvationContext *ctx, Mat K);
PetscErrorCode doAnalytical(PetscReal b, PetscReal epsIn, PetscReal epsOut, PQRData *pqr, PetscInt Nmax, Mat *L);
PetscErrorCode IntegratePanel(PetscInt numCorners, const PetscReal npanel[], const PetscReal point[], const PetscReal normal[], PetscScalar *fss, PetscScalar *fds, PetscScalar *fess, PetscScalar *feds);
PetscErrorCode makeSurfaceToSurfacePanelOperators_Laplace(DM dm, Vec w, Vec n, Mat *V, Mat *K);
PetscErrorCode makeSurfaceToChargePanelOperators(DM dm, Vec w, Vec n, PQRData *pqr, Mat *potential, Mat *field, Mat *singleLayer, Mat *doubleLayer);
PetscErrorCode makeSurfaceToChargePointOperators(Vec coordinates, Vec w, Vec n, PQRData *pqr, Mat *potential, Mat *field, Mat *singleLayer, Mat *doubleLayer);
PetscErrorCode makeSurfaceToSurfacePointOperators_Laplace(Vec coordinates, Vec w, Vec n, Mat *V, Mat *K, Mat *Elec);
PetscErrorCode makeBEMPcmQualMatrices(DM dm, BEMType bem, PetscReal epsIn, PetscReal epsOut, PQRData *pqr, Vec coordinates, Vec w, Vec n, Mat *L);
PetscErrorCode ComputeBEMResidual(SNES snes, Vec x, Vec r, void *ctx);
PetscErrorCode ComputeBEMJacobian(SNES snes, Vec x, Mat J, Mat P, void *ctx);
PetscErrorCode nonlinearH(Vec E, HContext *ctx, Vec *hEn);
PetscErrorCode ASCBq(Vec sigma, Vec *Bq, NonlinearContext *ctx);
PetscErrorCode FormASCNonlinearMatrix(Vec sigma, Mat *A, NonlinearContext *ctx);
PetscErrorCode FastRHS(Vec sigma, Vec *out, NonlinearContext *ctx);
PetscErrorCode NonlinearPicard(PetscErrorCode (*lhs)(Vec, Mat*, void*), PetscErrorCode (*rhs)(Vec, Vec*, void*), Vec guess, Vec weights, void *ctx, Vec *sol, PetscReal *estError);
PetscErrorCode NonlinearAnderson(PetscErrorCode (*lhs)(Vec, Mat*, void*), PetscErrorCode (*rhs)(Vec, Vec*, void*), Vec guess, Vec weights, void *ctx, Vec *sol, PetscReal *estError);
PetscErrorCode FastPicard(PetscErrorCode (*rhs)(Vec, Vec*, void*), Vec guess, Vec weights, void *ctx, Vec *sol);
PetscErrorCode makeBEMPcmQualReactionPotentialNonlinear(DM dm, BEMType bem, HContext params, SolvationContext *ctx, PetscReal epsIn, PetscReal epsOut, PQRData *pqr, Vec coordinates, Vec w, Vec n, Vec react, PetscReal *estError);
PetscErrorCode makeBEMPcmQualReactionPotential(DM dm, BEMType bem, SolvationContext *ctx, PetscReal epsIn, PetscReal epsOut, PQRData *pqr, Vec coordinates, Vec w, Vec n, Vec react);
PetscErrorCode CalculateAnalyticSolvationEnergy(PetscReal epsIn, PetscReal epsOut, PQRData *pqr, PetscReal R, PetscInt Nmax, Vec react, PetscReal *E);
PetscErrorCode CalculateBEMSolvationEnergy(DM dm, SolvationContext *ctx, const char prefix[], BEMType bem, HContext params, PetscReal epsIn, PetscReal epsOut, PQRData *pqr, Vec w, Vec n, Vec react, PetscReal *E, PetscReal *estError);
PetscErrorCode ProcessOptions(MPI_Comm comm, SolvationContext *ctx);
PetscErrorCode workprectests(int argc, char **argv);
PetscErrorCode CalcASCResidual(SNES snes, Vec x, Vec resid, NonlinearContext *ctx);
PetscErrorCode ComputeBEMNonlinearJacobian(SNES snes, Vec x, Mat J, Mat P, NonlinearContext *ctx);
PetscErrorCode MatFormViewOptions(Vec En, Vec hEn, NonlinearContext *ctx);
#endif
