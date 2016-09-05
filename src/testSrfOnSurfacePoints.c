#include <petsc.h>
#include <petsc/private/dmpleximpl.h>
#include "constants.h"
#include "surface.h"

typedef enum {BEM_POINT, BEM_PANEL, BEM_POINT_MF, BEM_PANEL_MF} BEMType;

/* Performance characterization */
PetscLogEvent CalcE_Event, CalcL_Event, CalcR_Event, CalcStoQ_Event, CalcStoS_Event, IntegratePanel_Event;

typedef struct {
  Vec q;   /* Charge values */
  Vec xyz; /* Charge coordinates, always 3D */
  Vec R;   /* Charge radii */
} PQRData;

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
} SolvationContext;

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
    ierr = PetscOptionsString("-crg_filename", "The filename for the .crg file", "testSrfOnSurfacePoints", ctx->crgFile, ctx->crgFile, sizeof(ctx->crgFile), NULL);CHKERRQ(ierr);
    ierr = PetscOptionsBool("-is_sphere", "Use a spherical test case", "testSrfOnSurfacePoints", ctx->isSphere, &ctx->isSphere, NULL);CHKERRQ(ierr);
    ierr = PetscOptionsInt("-num_charges", "The number of atomic charges in the solute", "testSrfOnSurfacePoints", ctx->numCharges, &ctx->numCharges, NULL);CHKERRQ(ierr);
    ierr = PetscOptionsString("-srf_base", "The basename for the .srf file", "testSrfOnSurfacePoints", ctx->basename, ctx->basename, sizeof(ctx->basename), NULL);CHKERRQ(ierr);
    ierr = PetscOptionsInt("-srf_num", "The resolution number of the mesh", "testSrfOnSurfacePoints", ctx->srfNum, &ctx->srfNum, NULL);CHKERRQ(ierr);
    ierr = PetscOptionsInt("-nmax", "The order of the multipole expansion", "testSrfOnSurfacePoints", ctx->Nmax, &ctx->Nmax, NULL);CHKERRQ(ierr);
    ierr = PetscOptionsReal("-density", "The density of points for BEM", "testSrfOnSurfacePoints", ctx->density, &ctx->density, NULL);CHKERRQ(ierr);
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

#undef __FUNCT__
#define __FUNCT__ "factorial"
PETSC_STATIC_INLINE PetscErrorCode factorial(PetscInt n, PetscReal *fact)
{
  PetscReal f = 1.0;
  PetscInt  i;

  PetscFunctionBeginUser;
  for (i = 2; i <= n; ++i) f *= i;
  *fact = f;
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "convertToSpherical"
PETSC_STATIC_INLINE PetscErrorCode convertToSpherical(const PetscScalar xyz[], PetscScalar r[])
{
  const PetscReal x = PetscRealPart(xyz[0]);
  const PetscReal y = PetscRealPart(xyz[1]);
  const PetscReal z = PetscRealPart(xyz[2]);

  PetscFunctionBeginUser;
  r[0] = PetscSqrtReal(x*x + y*y + z*z);
  r[1] = atan2(y, x);
  r[2] = 0.0;
  if (PetscAbsReal(r[0]) > 0.0) r[2] = acos(z/r[0]);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "legendre"
/*@
  legendre - Calculate the Associated Legendre polynomial P^l_m(x).

  Input Parameters:
+ l - The order
. m - The suborder
. x - The argument

  Output Parameter:
. leg - An array of the associated Legendre polynomials evaluated at the points x

  Note: Speed is obtained by direct calculation of polynomial coefficients rather than recursion.
  Polynomial coefficients can increase in magnitude very quickly with polynomial degree,
  leading to decreased accuracy (estimated by err). If you need higher degrees for the polynomials,
  use recursion-based algorithms.

  Level: intermediate

.seealso: 
@*/
PetscErrorCode legendre(PetscInt l, PetscInt m, PetscScalar x, PetscScalar *leg, PetscScalar *err)
{
  /* The error estimate is based on worst case scenario and the significant digits, and thus
     based on the largest polynomial coefficient and machine error, "eps" */
  PetscInt       maxcf;    /* largest polynomial coefficient */
  PetscReal      cfnm = 1; /* proportionality constant for m < 0 polynomials compared to m > 0 */
  PetscReal      cl, x2, p, f1, f2, f3;
  PetscInt       px, j;
  PetscErrorCode ierr;

  PetscFunctionBeginUser;
  /* The polynomials are not defined for |x| > 1 */
  if (PetscAbsScalar(x) > 1) SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE, "Invalid input |x %g| > 1", PetscRealPart(x));
  /* Could also define this to be 0 */
  if (PetscAbsInt(m) > PetscAbsInt(l)) SETERRQ2(PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE, "Invalid input m %d > l %d", m, l);
  if (l < 0) l = -(l+1);
  if (m < 0) {
    PetscReal num, den;

    m    = -m;
    ierr = factorial(l-m, &num);CHKERRQ(ierr);
    ierr = factorial(l+m, &den);CHKERRQ(ierr);
    cfnm = PetscPowInt(-1, m)*num/den;
  }
  /* Calculate coef of maximum degree in x from the explicit analytical formula */
  ierr  = factorial(2*l, &f1);CHKERRQ(ierr);
  ierr  = factorial(l,   &f2);CHKERRQ(ierr);
  ierr  = factorial(l-m, &f3);CHKERRQ(ierr);
  cl    = PetscPowInt(-1, m) * cfnm * f1/((1 << l)*f2*f3);
  maxcf = PetscAbsInt(cl);
  px    = l-m;
  /* Power of x changes from one term to the next by 2. Also needed for sqrt(1-x^2). */
  x2    = x*x; /* TODO make pointwise square */
  /* Calculate efficiently P_l^m (x)/sqrt(1-x^2)^(m/2) - that is, only the polynomial part.
     At least one coefficient is guaranteed to exist - there is no null Legendre polynomial. */
  p     = cl; /* TODO make an array of cl */
  for (j = l-1; j >= 0; --j) {
    /* Check the exponent of x for current coefficient, px. If it is 0 or 1, just exit the loop */
    if (px < 2) break;
    /* If current exponent is >=2, there is a "next" coefficient; multiply p by x2 and add it. Calculate the current coefficient */
    cl = -(j+j+2-l-m)*(j+j+1-l-m)/(2*(j+j+1)*(l-j))*cl;
    
    if (maxcf < PetscAbsReal(cl)) maxcf = PetscAbsReal(cl);
    /* ...and add to the polynomial */
    p = p*x2 + cl; /* TODO make this pointwise multiply */
    /* Decrease the exponent of x - this is the exponent of x corresponding to the newly added coefficient */
    px -= 2;
  }
  /* Estimate the error */
  if (err) *err = maxcf*PETSC_MACHINE_EPSILON;

  /* Now we're done adding coefficients. However, if the exponent of x
     corresponding to the last added coefficient is 1 (polynomial is odd),
     multiply the polynomial by x */
  if (px == 1) p = p*x;

  /* All that's left is to multiply the whole thing with sqrt(1-x^2)^(m/2). No further calculations are needed if m = 0. */
  if (m == 0) {*leg = p; PetscFunctionReturn(0);}

  x2 = 1-x2;
  /* First, multiply by the integer part of m/2 */
  for (j = 1; j < PetscFloorReal(m/2.0); ++j) p = p*x2; /* TODO make this pointwise multiply */
  /* If m is odd, there is an additional factor sqrt(1-x^2) */
  if (m % 2) p = p*PetscSqrtReal(x2); /* TODO make this pointwise multiply */
  *leg = p;
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "legendre2"
/*
  This is code to compute P^m_n(z) = (-1)^m (1 - z^2)^{m/2} \frac{d^m P_n(z)}{dz^m}

  leg is an arry of length nz*(l+1)

  Note: http://www.accefyn.org.co/revista/Vol_37/145/541-544.pdf
*/
PetscErrorCode legendre2(PetscInt l, PetscInt nz, PetscScalar z, PetscScalar leg[])
{
  PetscReal sqz2   = PetscSqrtReal(1.0 - PetscSqr(PetscRealPart(z)));
  PetscReal hsqz2  = 0.5*sqz2;
  PetscReal ihsqz2 = PetscRealPart(z)/hsqz2;
  PetscReal fac    = 1.0;
  PetscInt  pre    = l % 2 ? -1 : 1;
  PetscInt  m;

  PetscFunctionBeginUser;
  for (m = 2; m <= l; ++m) fac *= m;
  if (!l) {
    leg[0] = 1.0;
  } else if (l == 1) {
    leg[0] = -hsqz2;
    leg[1] = z;
    leg[2] = sqz2;
  } else {
    leg[0] = (1.0 - 2.0*PetscAbsReal(l - 2.0*PetscFloorReal(l/2.0)))*PetscPowReal(hsqz2, l)/fac;
    leg[1] = -leg[0]*l*ihsqz2;
    for (m = 1; m < 2*l; ++m) leg[m+1] = (m - l)*ihsqz2*leg[m] - (2*l - m + 1)*m*leg[m-1];
  }
  for (m = 0; m <= 2*l; ++m, pre = -pre) leg[m] *= pre;
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "makeSphereChargeDistribution"
/*
  Select a set of point charges from a grid with spacing dx which are inside sphere of radius R, and delta away from the surface
*/
PetscErrorCode makeSphereChargeDistribution(PetscReal R, PetscInt numCharges, PetscReal dx, PetscReal delta, PQRData *data)
{
  PetscRandom     rand;
  const PetscReal maxChargeValue = 0.85;
  PetscInt        numPoints      = 0, *select, c;
  PetscReal       x, y, z;
  PetscErrorCode  ierr;

  PetscFunctionBeginUser;
  {
    PetscReal vals[8];
    PetscInt  nmax = 8, i;
    PetscBool flg;

    ierr = PetscOptionsGetRealArray(NULL, NULL, "-test", vals, &nmax, &flg);CHKERRQ(ierr);
    if (flg) {
      numCharges = nmax/4;
      ierr = VecCreate(PETSC_COMM_WORLD, &data->q);CHKERRQ(ierr);
      ierr = VecSetSizes(data->q, numCharges, PETSC_DETERMINE);CHKERRQ(ierr);
      ierr = PetscObjectSetName((PetscObject) data->q, "Atomic Charges");CHKERRQ(ierr);
      ierr = VecSetFromOptions(data->q);CHKERRQ(ierr);
      ierr = VecCreate(PETSC_COMM_WORLD, &data->xyz);CHKERRQ(ierr);
      ierr = VecSetSizes(data->xyz, numCharges*3, PETSC_DETERMINE);CHKERRQ(ierr);
      ierr = PetscObjectSetName((PetscObject) data->xyz, "Atomic XYZ");CHKERRQ(ierr);
      ierr = VecSetBlockSize(data->xyz, 3);CHKERRQ(ierr);
      ierr = VecSetFromOptions(data->xyz);CHKERRQ(ierr);
      ierr = VecDuplicate(data->q, &data->R);CHKERRQ(ierr);
      ierr = PetscObjectSetName((PetscObject) data->R, "Atomic radii");CHKERRQ(ierr);
      ierr = VecSet(data->R, 0.0);CHKERRQ(ierr);
      for (i = 0; i < numCharges; ++i) {
        ierr = VecSetValues(data->q, 1, &i, &vals[i*4], INSERT_VALUES);CHKERRQ(ierr);
        ierr = VecSetValuesBlocked(data->xyz, 1, &i, &vals[i*4+1], INSERT_VALUES);CHKERRQ(ierr);
      }
      PetscFunctionReturn(0);
    }
  }
  if (delta < 0.0) delta = dx;

  /* Form a grid of points [-R, R]^3 with spacing dx */
  for (z = -R; z < R; z += dx) {
    for (y = -R; y < R; y += dx) {
      for (x = -R; x < R; x += dx) {
        const PetscReal dist = sqrt(x*x + y*y + z*z);

        if (dist < R - delta) ++numPoints;
      }
    }
  }

  ierr = VecCreate(PETSC_COMM_WORLD, &data->q);CHKERRQ(ierr);
  ierr = VecSetSizes(data->q, numCharges, PETSC_DETERMINE);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) data->q, "Atomic Charges");CHKERRQ(ierr);
  ierr = VecSetFromOptions(data->q);CHKERRQ(ierr);
  ierr = VecCreate(PETSC_COMM_WORLD, &data->xyz);CHKERRQ(ierr);
  ierr = VecSetSizes(data->xyz, numCharges*3, PETSC_DETERMINE);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) data->xyz, "Atomic XYZ");CHKERRQ(ierr);
  ierr = VecSetBlockSize(data->xyz, 3);CHKERRQ(ierr);
  ierr = VecSetFromOptions(data->xyz);CHKERRQ(ierr);

  ierr = PetscCalloc1(numPoints, &select);CHKERRQ(ierr);
  if ((numCharges >= 0) && (numCharges < numPoints)) {
    ierr = PetscRandomCreate(PETSC_COMM_WORLD, &rand);CHKERRQ(ierr);
    ierr = PetscRandomSetFromOptions(rand);CHKERRQ(ierr);
    ierr = PetscRandomSetInterval(rand, 0, numPoints);CHKERRQ(ierr);
    for (c = 0; c < numCharges; ++c) {
      ierr = PetscRandomGetValueReal(rand, &x);CHKERRQ(ierr);
      if (select[(PetscInt) PetscFloorReal(x)]) --c;
      select[(PetscInt) PetscFloorReal(x)] = 1;
    }
    ierr = PetscRandomDestroy(&rand);CHKERRQ(ierr);
  }
  
  ierr = PetscRandomCreate(PETSC_COMM_WORLD, &rand);CHKERRQ(ierr);
  ierr = PetscRandomSetFromOptions(rand);CHKERRQ(ierr);
  ierr = PetscRandomSetInterval(rand, -maxChargeValue, maxChargeValue);CHKERRQ(ierr);
  numPoints = 0; c = 0;
  for (z = -R; z < R; z += dx) {
    for (y = -R; y < R; y += dx) {
      for (x = -R; x < R; x += dx) {
        const PetscReal dist   = sqrt(x*x + y*y + z*z);
        PetscReal       pos[3] = {x, y, z}, q;

        if (dist < R - delta) {
          if (select[numPoints]) {
            ierr = PetscRandomGetValueReal(rand, &q);CHKERRQ(ierr);
            ierr = VecSetValues(data->q, 1, &c, &q, INSERT_VALUES);CHKERRQ(ierr);
            ierr = VecSetValuesBlocked(data->xyz, 1, &c, pos, INSERT_VALUES);CHKERRQ(ierr);
            ++c;
          }
          ++numPoints;
        }
      }
    }
  }
  ierr = PetscRandomDestroy(&rand);CHKERRQ(ierr);
  ierr = PetscFree(select);CHKERRQ(ierr);

  ierr = VecDuplicate(data->q, &data->R);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) data->R, "Atomic radii");CHKERRQ(ierr);
  ierr = VecSet(data->R, 0.0);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "computeEnm"
/*@
  computeEnm - Compute the multipole coefficients for the protein charges

  Input Parameters:
+ b - the sphere radius, in Angstroms
. epsIn - the dielectric constant inside the protein
. pqrData - the PQRData context
. qVec - The charge vector to use instead of the pqrData vector
- Nmax - the maximum multipole order to use

  Output Parameters:
. Enm - The vector of multipole coefficients (packed real part, imaginary part)

  Level: beginner

.seealso: computeBnm(), doAnalytical()
@*/
PetscErrorCode computeEnm(PetscReal b, PetscReal epsIn, PQRData *pqr, Vec qVec, PetscInt Nmax, Vec Enm)
{
  PetscScalar   *xyz, *q;
  PetscReal     *P;
  PetscInt       Nq, n, m, k, idx;
  PetscErrorCode ierr;

  PetscFunctionBeginUser;
  ierr = PetscMalloc1(2*Nmax+1,&P);CHKERRQ(ierr);
  ierr = VecGetLocalSize(qVec, &Nq);CHKERRQ(ierr);
  ierr = VecSet(Enm, 0.0);CHKERRQ(ierr);
  ierr = VecGetArray(pqr->xyz, &xyz);CHKERRQ(ierr);
  ierr = VecGetArray(qVec, &q);CHKERRQ(ierr);
  for (k = 0; k < Nq; ++k) {
    PetscScalar r[3], val;

    ierr = convertToSpherical(&xyz[k*3], r);CHKERRQ(ierr);
    for (n = 0, idx = 0; n <= Nmax; ++n) {
      ierr = legendre2(n, 1, cos(r[2]), P);CHKERRQ(ierr);
      for (m = -n ; m <= n; ++m, ++idx) {
        const PetscReal Pnm  = P[PetscAbsInt(m)+n];
        PetscReal       ff; /* (n - |m|)! / (n + |m|)! */
        PetscReal       num, den;

        //ierr = PetscPrintf(PETSC_COMM_SELF, "Charge %d P(%d, |%d|) %g\n", k, n, m, Pnm);CHKERRQ(ierr);
        ierr = factorial(n - PetscAbsInt(m), &num);CHKERRQ(ierr);
        ierr = factorial(n + PetscAbsInt(m), &den);CHKERRQ(ierr);
        ff   = num/den;
		val  = ff * q[k] * PetscPowScalar(r[0], n) * Pnm;
        ierr = VecSetValue(Enm, idx*2+0,  val*cos(m*r[1]), ADD_VALUES);CHKERRQ(ierr);
        ierr = VecSetValue(Enm, idx*2+1, -val*sin(m*r[1]), ADD_VALUES);CHKERRQ(ierr);
      }
    }
  }
  ierr = VecRestoreArray(pqr->xyz, &xyz);CHKERRQ(ierr);
  ierr = VecRestoreArray(qVec, &q);CHKERRQ(ierr);
  ierr = PetscFree(P);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "computeBnm"
/*@
  computeBnm - Compute the multipole coefficients for the sphere interior (reaction potential)

  Input Parameters:
+ b - the sphere radius, in Angstroms
. epsIn - the dielectric constant inside the protein
. epsOut - the dielectric constant outside the protein
. Nmax - the maximum multipole order to use
- Enm - The vector of multipole coefficients for the charge field

  Output Parameters:
. Bnm - The vector of multipole coefficients (packed real part, imaginary part)

  Level: beginner

.seealso: computeEnm(), doAnalytical()
@*/
PetscErrorCode computeBnm(PetscReal b, PetscReal epsIn, PetscReal epsOut, PetscInt Nmax, Vec Enm, Vec Bnm)
{
#if 0
  PetscReal      epsHat = 2.0*(epsIn - epsOut)/(epsIn + epsOut);
#endif
  PetscScalar   *bnm;
  PetscInt       n, m, idx;
  PetscErrorCode ierr;

  PetscFunctionBeginUser;
  ierr = VecCopy(Enm, Bnm);CHKERRQ(ierr);
  ierr = VecGetArray(Bnm, &bnm);CHKERRQ(ierr);
  for (n = 0, idx = 0; n <= Nmax; ++n) {
    const PetscReal val = ((epsIn - epsOut)*(n+1))/(epsIn * (n*epsIn + (n+1)*epsOut)) * (1.0/PetscPowReal(b, (2*n+1)));

    for (m = -n; m <= n; ++m, ++idx) {
      bnm[idx*2+0] *= val;
      bnm[idx*2+1] *= val;
#if 0
     /* I forget what these coefficients are for */
	 Vlambda = b/(1+2*n);
	 Klambda = -1/(2*(1+2*n));
	 Snm(iIndex,jIndex) = Bnm(iIndex,jIndex) / Vlambda;
	 Snm2(iIndex,jIndex) = epsHat/(1 + epsHat*Klambda) * (n+1)/b^(2*n+2)* Enm(iIndex,jIndex);
#endif
    }
  }
  ierr = VecRestoreArray(Bnm, &bnm);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "computePotential"
/*@
  computePotential - Compute the reaction potential at the charge locations

  Input Parameters:
+ b - the sphere radius, in Angstroms
. pqrData - the PQRData context
. epsIn - the dielectric constant inside the protein
. epsOut - the dielectric constant outside the protein
. Nmax - the maximum multipole order to use
- Enm - The vector of multipole coefficients for the charge field

  Output Parameters:
. phi - The reaction potential values at the charge locations

  Level: beginner

.seealso: computeEnm(), computeBnm(), doAnalytical()
@*/
PetscErrorCode computePotential(PQRData *pqr, PetscInt Nmax, Vec Bnm, Vec phi)
{
  PetscScalar   *xyz, *bnm, *p;
  PetscReal     *P;
  PetscInt       Nq, n, m, k, idx = 0;
  PetscErrorCode ierr;

  PetscFunctionBeginUser;
  ierr = VecGetLocalSize(pqr->q, &Nq);CHKERRQ(ierr);
  ierr = PetscMalloc1(2*Nmax+1,&P);CHKERRQ(ierr);
  ierr = VecSet(phi, 0.0);CHKERRQ(ierr);
  ierr = VecGetArray(pqr->xyz, &xyz);CHKERRQ(ierr);
  ierr = VecGetArray(Bnm, &bnm);CHKERRQ(ierr);
  ierr = VecGetArray(phi, &p);CHKERRQ(ierr);
  for (k = 0; k < Nq; ++k) {
    PetscScalar r[3], val;

    ierr = convertToSpherical(&xyz[k*3], r);CHKERRQ(ierr);
    for (n = 0, idx = 0; n <= Nmax; ++n) {
      ierr = legendre2(n, 1, cos(r[2]), P);CHKERRQ(ierr);
      for (m = -n ; m <= n; ++m, ++idx) {
        const PetscReal Pnm  = P[PetscAbsInt(m)+n];

        p[k] += PetscPowReal(r[0], n) * Pnm * (bnm[idx*2+0] * cos(m*r[1]) - bnm[idx*2+1] * sin(m*r[1]));
        val = bnm[idx*2+1] * cos(m*r[1]) + bnm[idx*2+0] * sin(m*r[1]);
        if (val > 1e-2) SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_PLIB, "Imaginary part of potential is nonzero %g", (double) val);
      }
    }
  }
  ierr = VecRestoreArray(pqr->xyz, &xyz);CHKERRQ(ierr);
  ierr = VecRestoreArray(Bnm, &bnm);CHKERRQ(ierr);
  ierr = VecRestoreArray(phi, &p);CHKERRQ(ierr);
  ierr = PetscFree(P);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "doAnalytical"
/*@
  doAnalytical - Compute the analytical solvation matrix L

  Input Parameters:
+ b - the sphere radius, in Angstroms
. epsIn - the dielectric constant inside the protein
. epsOut - the dielectric constant outside the protein
. pqrData - the PQRData context
- Nmax - the maximum multipole order to use

  Output Parameters:
+ L    - the actual solvation matrix (Hessian)
- Lbib - the BIBEE/CFA solvation matrix (Hessian)

  Level: beginner

  Note: In order to get kcal/mol energies, you need to multiply by 332.112 outside this function

.seealso: computeEnm(), computeBnm(), computePotential()
@*/
PetscErrorCode doAnalytical(PetscReal b, PetscReal epsIn, PetscReal epsOut, PQRData *pqr, PetscInt Nmax, Mat *L)
{
  Vec            Enm, Bnm, tmpq;
  PetscScalar   *a;
  PetscInt       Nq, q;
  PetscErrorCode ierr;

  PetscFunctionBeginUser;
  ierr = VecGetLocalSize(pqr->q, &Nq);CHKERRQ(ierr);
  ierr = VecDuplicate(pqr->q, &tmpq);CHKERRQ(ierr);
  ierr = MatCreateSeqDense(PETSC_COMM_SELF, Nq, Nq, NULL, L);CHKERRQ(ierr);
  ierr = PetscObjectSetOptionsPrefix((PetscObject) *L, "lref_");CHKERRQ(ierr);
  ierr = MatDenseGetArray(*L, &a);CHKERRQ(ierr);

  /* \sum^{N_{max}}_{l=0} 2l + 1 = 2 * (Nmax)(Nmax+1)/2 + Nmax+1 = (Nmax+1)^2 */
  ierr = VecCreateSeq(PETSC_COMM_SELF, PetscSqr(Nmax+1) * 2, &Enm);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) Enm, "Enm Coefficients");CHKERRQ(ierr);
  ierr = VecDuplicate(Enm, &Bnm);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) Bnm, "Bnm Coefficients");CHKERRQ(ierr);
  for (q = 0; q < Nq; ++q) {
    Vec phi;

    ierr = VecSet(tmpq, 0.0);CHKERRQ(ierr);
    ierr = VecSetValue(tmpq, q, 1.0, INSERT_VALUES);CHKERRQ(ierr);

    /* The vector phi should be L(:,q) */
    ierr = VecCreateSeqWithArray(PETSC_COMM_SELF, 1, Nq, &a[Nq*q], &phi);CHKERRQ(ierr);
    ierr = PetscObjectSetName((PetscObject) phi, "Reaction Potential");CHKERRQ(ierr);
    ierr = computeEnm(b, epsIn, pqr, tmpq, Nmax, Enm);CHKERRQ(ierr);
    ierr = computeBnm(b, epsIn, epsOut, Nmax, Enm, Bnm);CHKERRQ(ierr);
    ierr = computePotential(pqr, Nmax, Bnm, phi);CHKERRQ(ierr);
    ierr = VecDestroy(&phi);CHKERRQ(ierr);
  }
  ierr = MatDenseRestoreArray(*L, &a);CHKERRQ(ierr);
  ierr = MatAssemblyBegin(*L, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(*L, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = VecDestroy(&Enm);CHKERRQ(ierr);
  ierr = VecDestroy(&Bnm);CHKERRQ(ierr);
  ierr = VecDestroy(&tmpq);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "IntegratePanel"
/*@
  IntegratePanel - Returns potential at evaluation point due to unit monopole and unit dipole uniformly distributed on a panel.

  Input Parameters:
. panel - The vertex coordinates for this panel in the panel coordinate system
. point - The evaluation point in the panel coordinate system
. normal - [Optional] The evaluation direction in the panel coordinate system

  Output Parameters:
. fss - the potential at evalpnt due to a panel monopole
. fds - the potential at evalpnt due to a panel normal dipole distribution
. fess - the derivative of the monopole potential at evalpnt along direction
. feds - the derivative of the dipole potential at evalpnt along direction

  Note: This is called calcp() in Matlab and FFTSVD. All calculations take place in the panel coordinate system,
  in which the face lies in the x-y plane.
@*/
PetscErrorCode IntegratePanel(PetscInt numCorners, const PetscReal npanel[], const PetscReal point[], const PetscReal normal[], PetscScalar *fss, PetscScalar *fds, PetscScalar *fess, PetscScalar *feds)
{
  PetscScalar    fs, fsx, fsy,      fes;
  PetscScalar    fd, fdx, fdy, fdz, fed;
  PetscReal      zn  = point[2], znabs = PetscAbsReal(zn);
  PetscReal      elen[4]; /* the length of each edge in the panel */
  PetscReal      ct[4];   /* cos(th) where th is the angle that panel edge i makes with the x-axis */
  PetscReal      st[4];   /* sin(th) where th is the angle that panel edge i makes with the x-axis */
  PetscReal      fe[4];   /* x-z plane square distance from evalpnt to each vertex */
  PetscReal      r[4];    /* distance from evalpnt to each vertex */
  PetscReal      xmxv[4]; /* x distance from evalpnt to each vertex */
  PetscReal      ymyv[4]; /* y distance from evalpnt to each vertex */
  PetscReal      xri[4];  /* cos(th) where th is the angle to panel made by vector from each vertex to evalpnt */
  PetscReal      yri[4];  /* sin(th) where th is the angle to panel made by vector from each vertex to evalpnt */
  PetscBool      isNormal = PETSC_FALSE; /* The evaluation point lies along the line passing through a vertex oriented along the normal */
  PetscInt       c;
  PetscErrorCode ierr;

  PetscFunctionBeginUser;
  ierr = PetscLogEventBegin(IntegratePanel_Event, 0, 0, 0, 0);CHKERRQ(ierr);
  for (c = 0; c < numCorners; ++c) {
    /* Jay used a left-handed coordinate system, so iterate backwards */
    const PetscInt curr = (numCorners - c)%numCorners;
    const PetscInt next = (numCorners*2 - c - 1)%numCorners;
    PetscReal      dx[3];
    PetscInt       d;

    elen[c] = 0.0;
    for (d = 0; d < 3; ++d) elen[c] += PetscSqr(npanel[next*3+d] - npanel[curr*3+d]);
    elen[c] = PetscSqrtReal(elen[c]);
    /* My coordinate system seems rotated compared to Jay's */
    ct[c] = (npanel[next*3+0] - npanel[curr*3+0])/elen[c];
    st[c] = (npanel[next*3+1] - npanel[curr*3+1])/elen[c];

    for (d = 0; d < 3; ++d) dx[d] = point[d] - npanel[curr*3+d];
    xmxv[c] = dx[0]; 
    ymyv[c] = dx[1]; 
    fe[c]   = PetscSqr(dx[0]) + PetscSqr(dx[2]);
    r[c]    = PetscSqrtReal(PetscSqr(dx[1]) + fe[c]);
    if (r[c] < 1.005*znabs) isNormal = PETSC_TRUE;
    if (normal) {
      xri[c] = xmxv[c]/r[c];
      yri[c] = ymyv[c]/r[c];
    }
  }

  if (feds && *((PetscInt *) feds) == 0) {
    PetscInt d;
    for (d = 0; d < 3; ++d) {ierr = PetscPrintf(PETSC_COMM_SELF, "elen[%d] %g\n", d, elen[d]);CHKERRQ(ierr);}
    for (d = 0; d < 3; ++d) {ierr = PetscPrintf(PETSC_COMM_SELF, "ct[%d] %g\n", d, ct[d]);CHKERRQ(ierr);}
    for (d = 0; d < 3; ++d) {ierr = PetscPrintf(PETSC_COMM_SELF, "st[%d] %g\n", d, st[d]);CHKERRQ(ierr);}
    for (d = 0; d < 3; ++d) {ierr = PetscPrintf(PETSC_COMM_SELF, "xmxv[%d] %g\n", d, xmxv[d]);CHKERRQ(ierr);}
    for (d = 0; d < 3; ++d) {ierr = PetscPrintf(PETSC_COMM_SELF, "ymyv[%d] %g\n", d, ymyv[d]);CHKERRQ(ierr);}
    for (d = 0; d < 3; ++d) {ierr = PetscPrintf(PETSC_COMM_SELF, "fe[%d] %g\n", d, fe[d]);CHKERRQ(ierr);}
    for (d = 0; d < 3; ++d) {ierr = PetscPrintf(PETSC_COMM_SELF, "r[%d] %g\n", d, r[d]);CHKERRQ(ierr);}
  }

  /* The potential and dipole contributions are made by summing up a contribution from each edge */
  fs = 0;
  fd = 0;
  if (normal) {
    fsx = 0; fsy = 0;
    fdx = 0; fdy = 0; fdz = 0;
  }


  for (c = 0; c < numCorners; ++c) {
    const PetscInt next = (c+1)%numCorners;
    PetscReal      v, arg;

    /* v is the projection of the eval-i edge on the perpend to the side-i:  
       Exploits the fact that corner points in panel coordinates. */
    v = xmxv[c]*st[c] - ymyv[c]*ct[c];

    /* arg == zero if eval on next-i edge, but then v = 0. */
    arg = (r[c]+r[next] - elen[c])/(r[c]+r[next] + elen[c]);
    if (arg > 0.0) {
      PetscReal fln;

      fln = -PetscLogReal(arg);
      fs  = fs + v * fln;
      if (normal) {
        PetscReal fac;

        fac = (r[c]+r[next]-elen[c]) * (r[c]+r[next]+elen[c]);
        fac = v*(elen[c]+ elen[c])/fac;
        fsx = fsx + (fln*st[c] - fac*(xri[c] + xri[next]));
        fsy = fsy - (fln*ct[c] + fac*(yri[c] + yri[next]));
        fdz = fdz - (fac*( 1.0/r[c] + 1.0/r[next]));
      }
    }
    PetscReal s1, c1, s2, c2, s12, c12, val;

    if (!isNormal) {
      /* eval not near a vertex normal, use Hess-Smith */
      s1 = v*r[c];
      c1 = znabs*(xmxv[c]*ct[c]+ymyv[c]*st[c]);
      s2 = v*r[next];
      c2 = znabs*(xmxv[next]*ct[c]+ymyv[next]*st[c]);
    } else {
      /* eval near a vertex normal, use Newman */
      s1 = (fe[c]*st[c])-(xmxv[c]*ymyv[c]*ct[c]);
      c1 = znabs*r[c]*ct[c];
      s2 = (fe[next]*st[c])-(xmxv[next]*ymyv[next]*ct[c]);
      c2 = znabs*r[next]*ct[c];
    }
  
    s12 = (s1*c2)-(s2*c1);
    c12 = (c1*c2)+(s1*s2);
    val = atan2(s12, c12);
    fd  = fd+val;
    if (normal) {
      PetscReal fac, u1, u2, rr, fh1, fh2;

      u1 = xmxv[c]*ct[c] + ymyv[c]*st[c];
      u2 = xmxv[next]*ct[c]+ymyv[next]*st[c];
      if (isNormal) {
        rr  = r[c]*r[c];
        fh1 = xmxv[c]*ymyv[c];
        fh2 = xmxv[next]*ymyv[next];
        fac = c1/((c1*c1+s1*s1)*rr );
        if (zn < 0.0) fac = -fac;
        fdx = fdx + ((rr*v+fh1*u1)*fac);
        fdy = fdy - (fe[c]*u1*fac);
        rr  = r[next]*r[next];
        fac = c2/((c2*c2+s2*s2)*rr);
        if (zn < 0.0) fac = -fac;
        fdx = fdx - ((rr*v+fh2*u2)*fac);
        fdy = fdy + fe[next]*u2*fac;
      } else {
        fac = zn/(c1*c1+s1*s1);
        fdx = fdx + (u1*v*xri[c]+r[c]*ymyv[c])*fac;
        fdy = fdy + (u1*v*yri[c]-r[c]*xmxv[c])*fac;
        fac = zn/(c2*c2+s2*s2);
        fdx = fdx - ((u2*v*xri[next]+r[next]*ymyv[next])*fac);
        fdy = fdy - ((u2*v*yri[next]-r[next]*xmxv[next])*fac);
      }
    }
  }

  /* I do not understand this line, and it is screwing up */
  // if (fd < 0.0) fd = fd + 2*PETSC_PI;
  if (fd < -1.0e-7) fd = fd + 2*PETSC_PI;
  if (zn < 0.0) fd = -fd;
 
  fs = fs - zn*fd;

  if (normal) {
    fsx = fsx - zn*fdx;
    fsy = fsy - zn*fdy;
    fes = normal[0]*fsx + normal[1]*fsy - normal[2]*fd;
    fed = normal[0]*fdx + normal[1]*fdy + normal[2]*fdz;
    ierr = PetscLogFlops((2 + 61) * numCorners + 14);CHKERRQ(ierr);
  }
  ierr = PetscLogFlops((24 + 29) * numCorners + 2);CHKERRQ(ierr);

  /* No area normalization */
  *fss = fs;
  *fds = fd;
  if (normal) *fess = fes;
  if (normal) *feds = fed;
  ierr = PetscLogEventEnd(IntegratePanel_Event, 0, 0, 0, 0);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "makeSurfaceToSurfacePanelOperators_Laplace"
/*@
  makeSurfaceToSurfacePanelOperators_Laplace - Make A matrix mapping the surface to itself

  Input Parameters:
+ coordinates - The vertex coordinates
. w - The vertex weights
- n - The vertex normals

  Output Parameters:
+ V - The single layer operator
- K - The double layer operator

  Level: beginner

.seealso: doAnalytical()
@*/
PetscErrorCode makeSurfaceToSurfacePanelOperators_Laplace(DM dm, Vec w, Vec n, Mat *V, Mat *K)
{
  const PetscReal fac = 1.0/4.0/PETSC_PI;
  Vec             coordinates;
  PetscSection    coordSection;
  PetscInt        Np;
  PetscInt        i, j;
  PetscErrorCode  ierr;

  PetscFunctionBeginUser;
  ierr = PetscLogEventBegin(CalcStoS_Event, 0, 0, 0, 0);CHKERRQ(ierr);
  ierr = DMGetCoordinateSection(dm, &coordSection);CHKERRQ(ierr);
  ierr = DMGetCoordinatesLocal(dm, &coordinates);CHKERRQ(ierr);
  ierr = DMPlexGetHeightStratum(dm, 0, NULL, &Np);CHKERRQ(ierr);
  if (V) {ierr = MatCreateSeqDense(PETSC_COMM_SELF, Np, Np, NULL, V);CHKERRQ(ierr);}
  if (K) {ierr = MatCreateSeqDense(PETSC_COMM_SELF, Np, Np, NULL, K);CHKERRQ(ierr);}
  for (i = 0; i < Np; ++i) {
    PetscScalar *coords = NULL;
    PetscReal    panel[12], R[9], v0[3];
    PetscInt     numCorners, coordSize, d, e;

    ierr = DMPlexGetConeSize(dm, i, &numCorners);CHKERRQ(ierr);
    ierr = DMPlexVecGetClosure(dm, coordSection, coordinates, i, &coordSize, &coords);CHKERRQ(ierr);
    for (d = 0; d < 3; ++d) v0[d] = coords[d];
    ierr = DMPlexComputeProjection3Dto2D(coordSize, coords, R);CHKERRQ(ierr); /* 28 + 36 + 27 = 91 flops */
    for (d = 0; d < numCorners; ++d) {
      panel[d*3+0] = PetscRealPart(coords[d*2+0]);
      panel[d*3+1] = PetscRealPart(coords[d*2+1]);
      panel[d*3+2] = 0.0;
    }
    ierr = DMPlexVecRestoreClosure(dm, coordSection, coordinates, i, &coordSize, &coords);CHKERRQ(ierr);
    for (j = 0; j < Np; ++j) {
      PetscScalar *tcoords = NULL;
      PetscReal    centroid[3], cloc[3];
      PetscScalar  fss, fds;

      ierr = DMPlexVecGetClosure(dm, coordSection, coordinates, j, NULL, &tcoords);CHKERRQ(ierr);
      for (d = 0; d < 3; ++d) {
        centroid[d] = 0.0;
        for (e = 0; e < numCorners; ++e) centroid[d] += tcoords[e*3+d];
        centroid[d] /= numCorners;
      }
      /* Rotate centroid into panel coordinate system */
      for (d = 0; d < 3; ++d) {
        cloc[d] = 0.0;
        for (e = 0; e < 3; ++e) {
          cloc[d] += R[e*3+d] * (centroid[e] - v0[e]);
        }
      }
      ierr = DMPlexVecRestoreClosure(dm, coordSection, coordinates, j, NULL, &tcoords);CHKERRQ(ierr);
      /* 'panel' is the coordinates of the panel vertices in the panel coordinate system */
      /*  TODO pass normals if we want fess for Kp */
      ierr = IntegratePanel(numCorners, panel, cloc, NULL, &fss, &fds, NULL, NULL);CHKERRQ(ierr);

      if (V) {ierr = MatSetValue(*V, j, i, fss*fac, INSERT_VALUES);CHKERRQ(ierr);}
      if (K) {ierr = MatSetValue(*K, j, i, fds*fac, INSERT_VALUES);CHKERRQ(ierr);}
      /* if (Kp) {ierr = MatSetValue(*singleLayer, j, i, fess/4/PETSC_PI, INSERT_VALUES);CHKERRQ(ierr);} */
    }
  }
  ierr = PetscLogFlops(37 * Np*Np + 91 * Np + 2);CHKERRQ(ierr);
  if (V) {ierr = PetscLogFlops(Np*Np);CHKERRQ(ierr);}
  if (K) {ierr = PetscLogFlops(Np*Np);CHKERRQ(ierr);}
  if (V) {ierr = MatAssemblyBegin(*V, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);ierr = MatAssemblyEnd(*V, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);}
  if (K) {ierr = MatAssemblyBegin(*K, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);ierr = MatAssemblyEnd(*K, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);}
  ierr = PetscLogEventEnd(CalcStoS_Event, 0, 0, 0, 0);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "makeSurfaceToChargePanelOperators"
/*@
  makeSurfaceToChargePanelOperators - Make B and C matrices mapping point charges to the surface points

  Input Parameters:
+ coordinates - The vertex coordinates
. w - The vertex weights
. n - The vertex normals
- pqrData - the PQRData context

  Output Parameters:
+ potential - 
. field - 
. singleLayer -
- doubleLayer - 

  Level: beginner

.seealso: doAnalytical()
@*/
PetscErrorCode makeSurfaceToChargePanelOperators(DM dm, Vec w, Vec n, PQRData *pqr, Mat *potential, Mat *field, Mat *singleLayer, Mat *doubleLayer)
{
  const PetscReal    fac = 1.0/4.0/PETSC_PI;
  Vec                coordinates;
  PetscSection       coordSection;
  const PetscScalar *xyz;
  PetscInt           Nq, Np;
  PetscInt           i, j;
  PetscErrorCode     ierr;

  PetscFunctionBeginUser;
  ierr = PetscLogEventBegin(CalcStoQ_Event, 0, 0, 0, 0);CHKERRQ(ierr);
  ierr = DMGetCoordinateSection(dm, &coordSection);CHKERRQ(ierr);
  ierr = DMGetCoordinatesLocal(dm, &coordinates);CHKERRQ(ierr);
  ierr = VecGetLocalSize(pqr->q, &Nq);CHKERRQ(ierr);
  ierr = DMPlexGetHeightStratum(dm, 0, NULL, &Np);CHKERRQ(ierr);
  if (potential || field) SETERRQ(PetscObjectComm((PetscObject) dm), PETSC_ERR_SUP, "Do not currently make the potential or field operators");
  if (potential)   {ierr = MatCreateSeqDense(PETSC_COMM_SELF, Np, Nq, NULL, potential);CHKERRQ(ierr);}
  if (field)       {ierr = MatCreateSeqDense(PETSC_COMM_SELF, Np, Nq, NULL, field);CHKERRQ(ierr);}
  if (singleLayer) {ierr = MatCreateSeqDense(PETSC_COMM_SELF, Nq, Np, NULL, singleLayer);CHKERRQ(ierr);}
  if (doubleLayer) {ierr = MatCreateSeqDense(PETSC_COMM_SELF, Nq, Np, NULL, doubleLayer);CHKERRQ(ierr);}
  ierr = VecGetArrayRead(pqr->xyz, &xyz);CHKERRQ(ierr);
  for (i = 0; i < Np; ++i) {
    PetscScalar *coords = NULL;
    PetscReal    panel[12], R[9], v0[3];
    PetscInt     numCorners, coordSize, d, e;

    ierr = DMPlexGetConeSize(dm, i, &numCorners);CHKERRQ(ierr);
    ierr = DMPlexVecGetClosure(dm, coordSection, coordinates, i, &coordSize, &coords);CHKERRQ(ierr);
    for (d = 0; d < 3; ++d) v0[d] = coords[d];
    ierr = DMPlexComputeProjection3Dto2D(coordSize, coords, R);CHKERRQ(ierr); /* 28 + 36 + 27 = 91 flops */
    for (d = 0; d < numCorners; ++d) {
      panel[d*3+0] = PetscRealPart(coords[d*2+0]);
      panel[d*3+1] = PetscRealPart(coords[d*2+1]);
      panel[d*3+2] = 0.0;
    }
    ierr = DMPlexVecRestoreClosure(dm, coordSection, coordinates, i, &coordSize, &coords);CHKERRQ(ierr);
    for (j = 0; j < Nq; ++j) {
      PetscReal   qloc[3];
      PetscScalar fss, fds;

      /* Rotate charge location into panel coordinate system */
      for (d = 0; d < 3; ++d) {
        qloc[d] = 0.0;
        for (e = 0; e < 3; ++e) {
          qloc[d] += R[e*3+d] * (xyz[j*3+e] - v0[e]);
        }
      }
      /* 'panel' is the coordinates of the panel vertices in the panel coordinate system */
      ierr = IntegratePanel(numCorners, panel, qloc, NULL, &fss, &fds, NULL, NULL);CHKERRQ(ierr);

#if 0
      if (!i) {
        for (d = 0; d < numCorners; ++d) {
          ierr = PetscPrintf(PETSC_COMM_SELF, "v%d (%g, %g, %g)\n", d, panel[d*3+0], panel[d*3+1], panel[d*3+2]);CHKERRQ(ierr);
        }
        ierr = PetscPrintf(PETSC_COMM_SELF, "q (%g, %g, %g)\n", d, qloc[0], qloc[1], qloc[2]);CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_SELF, "fss %g fds %g\n", fss, fds);CHKERRQ(ierr);
      }
#endif

      if (potential)   {ierr = MatSetValue(*potential,   i, j, 0.0,     INSERT_VALUES);CHKERRQ(ierr);}
      if (field)       {ierr = MatSetValue(*field,       i, j, 0.0,     INSERT_VALUES);CHKERRQ(ierr);}
      if (singleLayer) {ierr = MatSetValue(*singleLayer, j, i, fss*fac, INSERT_VALUES);CHKERRQ(ierr);}
      if (doubleLayer) {ierr = MatSetValue(*doubleLayer, j, i, fds*fac, INSERT_VALUES);CHKERRQ(ierr);}
    }
  }
  ierr = PetscLogFlops(27 * Np*Nq + 91 * Np + 2);CHKERRQ(ierr);
  if (singleLayer) {ierr = PetscLogFlops(Np*Nq);CHKERRQ(ierr);}
  if (doubleLayer) {ierr = PetscLogFlops(Np*Nq);CHKERRQ(ierr);}
  ierr = VecRestoreArrayRead(pqr->xyz, &xyz);CHKERRQ(ierr);
  if (potential)   {ierr = MatAssemblyBegin(*potential,   MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);ierr = MatAssemblyEnd(*potential,   MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);}
  if (field)       {ierr = MatAssemblyBegin(*field,       MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);ierr = MatAssemblyEnd(*field,       MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);}
  if (singleLayer) {ierr = MatAssemblyBegin(*singleLayer, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);ierr = MatAssemblyEnd(*singleLayer, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);}
  if (doubleLayer) {ierr = MatAssemblyBegin(*doubleLayer, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);ierr = MatAssemblyEnd(*doubleLayer, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);}
  ierr = PetscLogEventEnd(CalcStoQ_Event, 0, 0, 0, 0);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "makeSurfaceToChargePointOperators"
/*@
  makeSurfaceToChargePointOperators - Make B and C matrices mapping point charges to the surface points

  Input Parameters:
+ coordinates - The vertex coordinates
. w - The vertex weights
. n - The vertex normals
- pqrData - the PQRData context

  Output Parameters:
+ potential - 
. field - 
. singleLayer -
- doubleLayer - 

  Level: beginner

.seealso: doAnalytical()
@*/
PetscErrorCode makeSurfaceToChargePointOperators(Vec coordinates, Vec w, Vec n, PQRData *pqr, Mat *potential, Mat *field, Mat *singleLayer, Mat *doubleLayer)
{
  const PetscReal    fac = 1.0/4.0/PETSC_PI;
  const PetscScalar *coords, *xyz, *weights, *normals;
  PetscInt           Nq, Np;
  PetscInt           i, j, d;
  PetscErrorCode     ierr;

  PetscFunctionBeginUser;
  ierr = PetscLogEventBegin(CalcStoQ_Event, 0, 0, 0, 0);CHKERRQ(ierr);
  ierr = VecGetLocalSize(pqr->q, &Nq);CHKERRQ(ierr);
  ierr = VecGetLocalSize(w, &Np);CHKERRQ(ierr);
  if (potential)   {ierr = MatCreateSeqDense(PETSC_COMM_SELF, Np, Nq, NULL, potential);CHKERRQ(ierr);}
  if (field)       {ierr = MatCreateSeqDense(PETSC_COMM_SELF, Np, Nq, NULL, field);CHKERRQ(ierr);}
  if (singleLayer) {ierr = MatCreateSeqDense(PETSC_COMM_SELF, Nq, Np, NULL, singleLayer);CHKERRQ(ierr);}
  if (doubleLayer) {ierr = MatCreateSeqDense(PETSC_COMM_SELF, Nq, Np, NULL, doubleLayer);CHKERRQ(ierr);}
  ierr = VecGetArrayRead(coordinates, &coords);CHKERRQ(ierr);
  ierr = VecGetArrayRead(pqr->xyz, &xyz);CHKERRQ(ierr);
  ierr = VecGetArrayRead(w, &weights);CHKERRQ(ierr);
  ierr = VecGetArrayRead(n, &normals);CHKERRQ(ierr);
  for (i = 0; i < Np; ++i) {
    for (j = 0; j < Nq; ++j) {
      PetscScalar G, dGdn;
      PetscReal   rvec[3];
      PetscReal   r = 0.0, dot = 0.0;

      for (d = 0; d < 3; ++d) {rvec[d] = coords[i*3+d] - xyz[j*3+d]; dot += rvec[d]*normals[i*3+d]; r += PetscSqr(rvec[d]);}
      r = PetscSqrtReal(r);

      if (r < 1e-10) {G = 0;     dGdn = 0;}
      else           {G = fac/r; dGdn = -dot*fac/PetscPowRealInt(r, 3);}

      if (potential)   {ierr = MatSetValue(*potential,   i, j, G,               INSERT_VALUES);CHKERRQ(ierr);}
      if (field)       {ierr = MatSetValue(*field,       i, j, dGdn,            INSERT_VALUES);CHKERRQ(ierr);}
      if (singleLayer) {ierr = MatSetValue(*singleLayer, j, i, G*weights[i],    INSERT_VALUES);CHKERRQ(ierr);}
      if (doubleLayer) {ierr = MatSetValue(*doubleLayer, j, i, dGdn*weights[i], INSERT_VALUES);CHKERRQ(ierr);}
    }
  }
  ierr = PetscLogFlops(16 * Np*Nq + 2);CHKERRQ(ierr);
  if (singleLayer) {ierr = PetscLogFlops(Np*Nq);CHKERRQ(ierr);}
  if (doubleLayer) {ierr = PetscLogFlops(Np*Nq);CHKERRQ(ierr);}
  ierr = VecRestoreArrayRead(coordinates, &coords);CHKERRQ(ierr);
  ierr = VecRestoreArrayRead(pqr->xyz, &xyz);CHKERRQ(ierr);
  ierr = VecRestoreArrayRead(w, &weights);CHKERRQ(ierr);
  ierr = VecRestoreArrayRead(n, &normals);CHKERRQ(ierr);
  if (potential)   {ierr = MatAssemblyBegin(*potential,   MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);ierr = MatAssemblyEnd(*potential,   MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);}
  if (field)       {ierr = MatAssemblyBegin(*field,       MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);ierr = MatAssemblyEnd(*field,       MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);}
  if (singleLayer) {ierr = MatAssemblyBegin(*singleLayer, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);ierr = MatAssemblyEnd(*singleLayer, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);}
  if (doubleLayer) {ierr = MatAssemblyBegin(*doubleLayer, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);ierr = MatAssemblyEnd(*doubleLayer, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);}
  ierr = PetscLogEventEnd(CalcStoQ_Event, 0, 0, 0, 0);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "makeSurfaceToSurfacePointOperators_Laplace"
/*@
  makeSurfaceToSurfacePointOperators_Laplace - Make V and K matrices mapping the surface to itself

  Input Parameters:
+ epsIn - the dielectric constant inside the protein
. epsOut - the dielectric constant outside the protein
. pqrData - the PQRData context
. coordinates - The vertex coordinates
. w - The vertex weights
- n - The vertex normals

  Output Parameters:
+ V - The single layer surface operator
- K - The double layer surface operator

  Level: beginner

.seealso: doAnalytical()
@*/
PetscErrorCode makeSurfaceToSurfacePointOperators_Laplace(Vec coordinates, Vec w, Vec n, Mat *V, Mat *K)
{
  const PetscReal    fac = 1.0/4.0/PETSC_PI;
  const PetscScalar *coords, *weights, *normals;
  PetscInt           Np, i, j, d;
  PetscErrorCode     ierr;

  PetscFunctionBeginUser;
  ierr = PetscLogEventBegin(CalcStoS_Event, 0, 0, 0, 0);CHKERRQ(ierr);
  ierr = VecGetLocalSize(w, &Np);CHKERRQ(ierr);
  if (V) {ierr = MatCreateSeqDense(PETSC_COMM_SELF, Np, Np, NULL, V);CHKERRQ(ierr);}
  if (K) {ierr = MatCreateSeqDense(PETSC_COMM_SELF, Np, Np, NULL, K);CHKERRQ(ierr);}
  ierr = VecGetArrayRead(coordinates, &coords);CHKERRQ(ierr);
  ierr = VecGetArrayRead(w, &weights);CHKERRQ(ierr);
  ierr = VecGetArrayRead(n, &normals);CHKERRQ(ierr);
  for (i = 0; i < Np; ++i) { /* Target points */
    for (j = 0; j < Np; ++j) { /* Source points */
      PetscReal   rvec[3];
      PetscReal   r = 0.0, dot = 0.0;

      for (d = 0; d < 3; ++d) {rvec[d] = coords[i*3+d] - coords[j*3+d]; dot += rvec[d]*normals[j*3+d]; r += PetscSqr(rvec[d]);}
      r = PetscSqrtReal(r);
      if (r > 1e-6) {
        if (V) {ierr = MatSetValue(*V, i, j, weights[j]* fac/r, INSERT_VALUES);CHKERRQ(ierr);}
        if (K) {ierr = MatSetValue(*K, i, j, weights[j]* dot*fac/PetscPowRealInt(r,3), INSERT_VALUES);CHKERRQ(ierr);}
      } else {
        const PetscReal R0 = PetscSqrtReal(weights[j]/PETSC_PI); /* radius of a circle with the area assoc with surfpt j */
        if (V) {ierr = MatSetValue(*V, i, j, (2 * PETSC_PI * R0) /4/PETSC_PI, INSERT_VALUES);CHKERRQ(ierr);}
      }
    }
  }
  ierr = PetscLogFlops(16 * Np*Np + 2);CHKERRQ(ierr);
  if (V) {ierr = PetscLogFlops(2 * Np*Np);CHKERRQ(ierr);}
  if (K) {ierr = PetscLogFlops(5 * Np*Np);CHKERRQ(ierr);}
  ierr = VecRestoreArrayRead(coordinates, &coords);CHKERRQ(ierr);
  ierr = VecRestoreArrayRead(w, &weights);CHKERRQ(ierr);
  ierr = VecRestoreArrayRead(n, &normals);CHKERRQ(ierr);
  if (V) {ierr = MatAssemblyBegin(*V, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);ierr = MatAssemblyEnd(*V, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);}
  if (K) {ierr = MatAssemblyBegin(*K, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);ierr = MatAssemblyEnd(*K, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);}
  ierr = PetscLogEventEnd(CalcStoS_Event, 0, 0, 0, 0);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "makeBEMPcmQualMatrices"
/*@
  makeBEMPcmQualMatrices - Make solvation matrix, L = C A^{-1} B in the Polarizable Continuum Model

  Input Parameters:
+ epsIn - the dielectric constant inside the protein
. epsOut - the dielectric constant outside the protein
. pqrData - the PQRData context
. coordinates - The vertex coordinates
. w - The vertex weights
- n - The vertex normals

  Output Parameters:
. L - The solvation matrix

  Level: beginner

.seealso: doAnalytical()
@*/
PetscErrorCode makeBEMPcmQualMatrices(DM dm, BEMType bem, PetscReal epsIn, PetscReal epsOut, PQRData *pqr, Vec coordinates, Vec w, Vec n, Mat *L)
{
  const PetscReal epsHat = (epsIn + epsOut)/(epsIn - epsOut);
  KSP             ksp;
  PC              pc;
  Mat             A, Bp, B, C, S, fact;
  Vec             d;
  PetscErrorCode  ierr;

  PetscFunctionBeginUser;
  switch (bem) {
  case BEM_POINT:
    ierr = makeSurfaceToSurfacePointOperators_Laplace(coordinates, w, n, NULL, &A);CHKERRQ(ierr);
    ierr = makeSurfaceToChargePointOperators(coordinates, w, n, pqr, NULL, &B, &C, NULL);CHKERRQ(ierr);
    ierr = PetscLogEventBegin(CalcL_Event, 0, 0, 0, 0);CHKERRQ(ierr);
    /* B = chargesurfop.dphidnCoul */
    ierr = MatDiagonalScale(B, w, NULL);CHKERRQ(ierr);
    ierr = MatScale(B, -1/epsIn);CHKERRQ(ierr);
    break;
  case BEM_PANEL:
    ierr = makeSurfaceToSurfacePanelOperators_Laplace(dm, w, NULL /*n*/, NULL, &A);CHKERRQ(ierr);
    ierr = makeSurfaceToChargePanelOperators(dm, w, NULL /*n*/, pqr, NULL, NULL, &C, &Bp);CHKERRQ(ierr);
    ierr = PetscLogEventBegin(CalcL_Event, 0, 0, 0, 0);CHKERRQ(ierr);
    /* Bp = chargesurfop.dlpToCharges */
    ierr = MatTranspose(Bp, MAT_INITIAL_MATRIX, &B);CHKERRQ(ierr);
    ierr = MatDestroy(&Bp);CHKERRQ(ierr);
    ierr = MatScale(B, -1/epsIn);CHKERRQ(ierr);
    break;
  default:
    SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONG, "Invalid BEM type: %d", bem);
  }
  /* C = chargesurfop.slpToCharges */
  ierr = MatScale(C, 4.0*PETSC_PI);CHKERRQ(ierr);
  /* A = surfsurfop.K */
  ierr = MatTranspose(A, MAT_REUSE_MATRIX, &A);CHKERRQ(ierr);
  ierr = MatDiagonalScale(A, NULL, w);CHKERRQ(ierr);
  ierr = VecDuplicate(w, &d);CHKERRQ(ierr);
  ierr = VecCopy(w, d);CHKERRQ(ierr);
  ierr = VecScale(d, epsHat/2.0);CHKERRQ(ierr);
  ierr = MatDiagonalSet(A, d, ADD_VALUES);CHKERRQ(ierr);
  ierr = VecDestroy(&d);CHKERRQ(ierr);

  ierr = KSPCreate(PetscObjectComm((PetscObject) A), &ksp);CHKERRQ(ierr);
  ierr = KSPSetType(ksp, KSPPREONLY);CHKERRQ(ierr);
  ierr = KSPGetPC(ksp, &pc);CHKERRQ(ierr);
  ierr = PCSetType(pc, PCLU);CHKERRQ(ierr);
  ierr = KSPSetOperators(ksp, A, A);CHKERRQ(ierr);
  ierr = KSPSetUp(ksp);CHKERRQ(ierr);
  ierr = PCFactorGetMatrix(pc, &fact);CHKERRQ(ierr);
  ierr = MatDuplicate(B, MAT_DO_NOT_COPY_VALUES, &S);CHKERRQ(ierr);
  ierr = MatMatSolve(fact, B, S);CHKERRQ(ierr);
  ierr = KSPDestroy(&ksp);CHKERRQ(ierr);
  ierr = MatDestroy(&A);CHKERRQ(ierr);
  ierr = MatDestroy(&B);CHKERRQ(ierr);
  ierr = MatMatMult(C, S, MAT_INITIAL_MATRIX, PETSC_DEFAULT, L);CHKERRQ(ierr);
  ierr = MatDestroy(&S);CHKERRQ(ierr);
  ierr = MatDestroy(&C);CHKERRQ(ierr);
  ierr = PetscLogEventEnd(CalcL_Event, 0, 0, 0, 0);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ComputeBEMResidual"
PetscErrorCode ComputeBEMResidual(SNES snes, Vec x, Vec r, void *ctx)
{
  Mat *A = (Mat *) ctx;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = MatMult(*A, x, r);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ComputeBEMJacobian"
PetscErrorCode ComputeBEMJacobian(SNES snes, Vec x, Mat J, Mat P, void *ctx)
{
  Mat *A = (Mat *) ctx;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = MatCopy(*A, P, SAME_NONZERO_PATTERN);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "nonlinearH"
/*@
  nonlinearH - the function h(En) defining the nonlinearity
  alpha*tanh(beta*En - gamma) + mu

  Input Parameters:
+ E - vector containing values of En
- ctx - context containing parameters alpha, beta, gamma, mu

  Output Parameters:
. h - vector containing values of h(En)

  Level: beginner

.seealso: makeBEMPcmQualReactionPotential
@*/
PetscErrorCode nonlinearH(Vec E, HContext *ctx, Vec *hEn)
{
  PetscInt length;
  PetscReal alpha, beta, gamma, mu;
  PetscErrorCode ierr;

  PetscFunctionBeginUser;
  alpha = ctx->alpha;
  beta  = ctx->beta;
  gamma = ctx->gamma;
  mu    = -alpha*PetscTanhReal(-gamma);
  ierr = VecGetSize(E, &length); CHKERRQ(ierr);
  //ierr = VecDuplicate(E, &h); CHKERRQ(ierr);
  
  PetscInt i;
  PetscReal val;
  for(i=0; i<length; ++i) {
    ierr = VecGetValues(E, 1, &i, &val); CHKERRQ(ierr);
    val = alpha*PetscTanhReal(beta*val - gamma) + mu;
    ierr = VecSetValue(*hEn, i, val, INSERT_VALUES); CHKERRQ(ierr);
  }
  ierr = VecAssemblyBegin(*hEn); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(*hEn); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ASCBq"
/*@ 
  ASCBq - forms right hand side Bq of ASC model
  
.seealso 
@*/
PetscErrorCode ASCBq(Vec sigma, Vec *Bq, NonlinearContext *ctx)
{
  PetscReal epsOut, epsIn, epsHat;
  Mat*      B;
  PetscErrorCode ierr;

  PetscFunctionBeginUser;
  epsOut = ctx->epsOut;
  epsIn  = ctx->epsIn;
  epsHat = (epsIn - epsOut)/epsIn;

  B = ctx->B;

  ierr = MatMult(*B, ctx->pqr->q, *Bq); CHKERRQ(ierr);
  ierr = VecScale(*Bq, epsHat); CHKERRQ(ierr);
  //ierr = VecView(ctx->pqr->q, PETSC_VIEWER_STDOUT_SELF);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "FormASCNonlinearMatrix"
/*@
  FormASCNonlinearMatrix - forms the left hand matrix for the ASC nonlinear BIE
  

@*/
PetscErrorCode FormASCNonlinearMatrix(Vec sigma, Mat *A, NonlinearContext *ctx)
{
  PetscReal epsOut, epsIn, epsHat, epsHat2;
  PetscInt  dim;
  Mat*      B;
  Mat*      K;
  Vec*      Bq;
  Vec*      w;
  Vec*      q;
  //Vec       iden;
  Vec       En;
  Vec       hEn;
  Vec       v1, v2, v3;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  epsOut  = ctx->epsOut;
  epsIn   = ctx->epsIn;
  B       = ctx->B;
  K       = ctx->K;
  Bq      = ctx->Bq;
  w       = ctx->w;
  q       = &(ctx->pqr->q);
  epsHat  = (epsIn - epsOut)/epsIn;
  epsHat2 = (epsOut + epsIn)/(2*epsOut);

  //get the dimension of sigma
  ierr = VecGetLocalSize(sigma, &dim); CHKERRQ(ierr);
  //initialize dense A matrix
  //if(A) {ierr = MatCreateSeqDense(PETSC_COMM_SELF, dim, dim, NULL, A);CHKERRQ(ierr);}

  //A = I + epsHat*(K - (1/2)*I)
  ierr = MatTranspose(*K, MAT_REUSE_MATRIX, A); CHKERRQ(ierr); //gives qualocated normal field operator
  //ierr = MatCopy(*K, *A, SAME_NONZERO_PATTERN); CHKERRQ(ierr);
  //ierr = MatScale(A, epsHat); CHKERRQ(ierr);
  //ierr = VecDuplicate(sigma, &iden); CHKERRQ(ierr);
  //ierr = VecSet(iden, epsHat2); CHKERRQ(ierr);
  //ierr = MatDiagonalSet(A, iden, ADD_VALUES); CHKERRQ(ierr);

  ierr = MatShift(*A, -0.5); CHKERRQ(ierr);
  ierr = MatScale(*A, epsHat); CHKERRQ(ierr);
  ierr = MatShift(*A, 1.0); CHKERRQ(ierr);

  /* Use sigma to solve for En. En = rhs - K*sigma */
  //Calculate En = -B*q - K*sigma
  ierr = VecDuplicate(sigma, &v1); CHKERRQ(ierr);
  ierr = VecDuplicate(*w   , &v2); CHKERRQ(ierr);
  ierr = MatCreateVecs(*B, NULL, &v3); CHKERRQ(ierr);

  //ierr = VecPointwiseMult(v3, sigma, *w); CHKERRQ(ierr);
  ierr = MatMult(*K, v3, v1); CHKERRQ(ierr);

  ierr = MatMult(*B, *q, v2); CHKERRQ(ierr);
  //ierr = VecPointwiseMult(v2, v3, *w); CHKERRQ(ierr);
  ierr = VecAXPY(v1, 1.0, v2); CHKERRQ(ierr);
  ierr = VecScale(v1, -1.0); CHKERRQ(ierr);
  ierr = VecDuplicate(v1, &En); CHKERRQ(ierr);
  ierr = VecCopy(v1, En); CHKERRQ(ierr);

  //compute h(En)
  ierr = VecDuplicate(En, &hEn); CHKERRQ(ierr);
  ierr = nonlinearH(En, ctx->hctx, &hEn); CHKERRQ(ierr);
  //printf("hEn:\n");
  //ierr = VecView(hEn, PETSC_VIEWER_STDOUT_SELF);
  
  //add h(En) to A
  //ierr = MatDiagonalSet(*A, hEn, ADD_VALUES); CHKERRQ(ierr);
  
  //printf("val should be %15.15f\n", 1+epsHat/2.);
  //ierr = MatView(*A, PETSC_VIEWER_STDOUT_SELF); CHKERRQ(ierr);

  //Scale columns by entries of w
  ierr = MatDiagonalScale(*A, NULL, *w); CHKERRQ(ierr);


  ierr = MatAssemblyBegin(*A, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);ierr = MatAssemblyEnd(*A, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "NonlinearPicard"
/*@
  NonlinearPicard - uses a picard iteration to solve for the charge layer sigma.
  
  Input Parameters:
  + lhs - function used to evaluate matrix A(x)
  . rhs - function evaluating the right hand side b(x)
  - guess - initial guess
  - ctx - additional user supplied data passed to lhs and rhs function evaluations (optional)

  Output Parameters:
  . sol - calculated solution

  Level: intermediate

.seealso: makeBEMPcmQualReactionPotential
@*/
PetscErrorCode NonlinearPicard(PetscErrorCode (*lhs)(Vec, Mat*, void*), PetscErrorCode (*rhs)(Vec, Vec*, void*), Vec guess, void *ctx, Vec *sol)
{
  PetscErrorCode ierr;
  PetscInt dim;
  Mat A;
  Vec b;
  Vec errvec;
  KSP ksp;
  PetscFunctionBegin;

  //get dimension of problem
  ierr = VecGetSize(guess, &dim); CHKERRQ(ierr);

  //initialize dense A matrix and b vector
  ierr = VecDuplicate(guess, &b); CHKERRQ(ierr);
  ierr = MatCreateSeqDense(PETSC_COMM_WORLD, dim, dim, NULL, &A); CHKERRQ(ierr);

  //make sol vector to be same length as guess and copy values
  //VecDuplicate(guess, sol);
  ierr = VecCopy(guess, *sol); CHKERRQ(ierr);
  
  //create linear solver context
  ierr = KSPCreate(PETSC_COMM_WORLD, &ksp); CHKERRQ(ierr);
  
  //initialize errvec for error calc
  ierr = VecDuplicate(guess, &errvec);
  ierr = VecSet(errvec, -1.2);

  //ierr = VecView(errvec, PETSC_VIEWER_STDOUT_SELF);
  PetscReal err = 1; 
  for(int iter=1; iter<=10; ++iter)
  {
    //printf("\n\nITERATION NUMBER %d\n", iter);
    ierr = (*lhs)(*sol, &A, ctx); CHKERRQ(ierr);
    ierr = (*rhs)(*sol, &b, ctx); CHKERRQ(ierr);
    //ierr = MatView(A, PETSC_VIEWER_STDOUT_SELF);
    ierr = KSPSetOperators(ksp, A, A); CHKERRQ(ierr);
    ierr = KSPSetTolerances(ksp, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT); CHKERRQ(ierr);
    ierr = KSPSetFromOptions(ksp); CHKERRQ(ierr);
    //copy current iteration value to prev
    ierr = VecCopy(*sol, errvec); CHKERRQ(ierr);

    //solve system
    ierr = KSPSolve(ksp, b, *sol); CHKERRQ(ierr);

    //calculate current error, will be inaccurate on first iteration
    ierr = VecAXPY(errvec, -1.0, *sol); CHKERRQ(ierr);
    ierr = VecNorm(errvec, NORM_2, &err); CHKERRQ(ierr);
    
    //printf("THE ERROR IS: %15.15f\n", err);
    //printf("sol:\n");
    //ierr = VecView(*sol, PETSC_VIEWER_STDOUT_SELF); CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "makeBEMPcmQualReactionPotential2"
/*@
  makeBEMPcmQualReactionPotential - Make the reaction potential, phi_react = Lq = C A^{-1} Bq in the Polarizable Continuum Model

  Input Parameters:
+ epsIn - the dielectric constant inside the protein
. epsOut - the dielectric constant outside the protein
. pqrData - the PQRData context
. coordinates - The vertex coordinates
. w - The vertex weights
- n - The vertex normals

  Output Parameters:
. react - The reaction potential

  Level: beginner

.seealso: doAnalytical()
@*/
PetscErrorCode makeBEMPcmQualReactionPotential2(DM dm, BEMType bem, PetscReal epsIn, PetscReal epsOut, PQRData *pqr, Vec coordinates, Vec w, Vec n, Vec react)
{
  //const PetscReal epsHat = (epsIn + epsOut)/(epsIn - epsOut);
  //SNES            snes;
  Mat             K, Bp, B, C;
  Vec             t0, t1;
  Vec             guess;
  PetscErrorCode  ierr;

  PetscFunctionBegin;
  switch (bem) {
  case BEM_POINT_MF:
    ierr = makeSurfaceToSurfacePointOperators_Laplace(coordinates, w, n, NULL, &K);CHKERRQ(ierr);
    ierr = makeSurfaceToChargePointOperators(coordinates, w, n, pqr, NULL, &B, &C, NULL);CHKERRQ(ierr);
    ierr = PetscLogEventBegin(CalcR_Event, 0, 0, 0, 0);CHKERRQ(ierr);
    /* B = chargesurfop.dphidnCoul */
    ierr = MatDiagonalScale(B, w, NULL);CHKERRQ(ierr);
    ierr = MatScale(B, -1/epsIn);CHKERRQ(ierr);
    break;
  case BEM_PANEL_MF:
    ierr = makeSurfaceToSurfacePanelOperators_Laplace(dm, w, NULL /*n*/, NULL, &K);CHKERRQ(ierr);
    ierr = makeSurfaceToChargePanelOperators(dm, w, NULL /*n*/, pqr, NULL, NULL, &C, &Bp);CHKERRQ(ierr);
    ierr = PetscLogEventBegin(CalcR_Event, 0, 0, 0, 0);CHKERRQ(ierr);
    /* Bp = chargesurfop.dlpToCharges */
    ierr = MatTranspose(Bp, MAT_INITIAL_MATRIX, &B);CHKERRQ(ierr);
    ierr = MatDestroy(&Bp);CHKERRQ(ierr);
    ierr = MatScale(B, -1/epsIn);CHKERRQ(ierr);
    break;
  default:
    SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONG, "Invalid BEM type: %d", bem);
  }
  /* C = chargesurfop.slpToCharges */
  ierr = MatScale(C, 4.0*PETSC_PI);CHKERRQ(ierr);
  /* A = surfsurfop.K */
  //ierr = MatTranspose(K, MAT_REUSE_MATRIX, &K);CHKERRQ(ierr);
  //ierr = MatDiagonalScale(A, NULL, w);CHKERRQ(ierr);
  //ierr = VecDuplicate(w, &d);CHKERRQ(ierr);
  //ierr = VecCopy(w, d);CHKERRQ(ierr);
  //ierr = VecScale(d, epsHat/2.0);CHKERRQ(ierr);
  //ierr = MatDiagonalSet(A, d, ADD_VALUES);CHKERRQ(ierr);
  //ierr = VecDestroy(&d);CHKERRQ(ierr);

  ierr = MatCreateVecs(B, NULL, &t0);CHKERRQ(ierr);
  ierr = VecDuplicate(t0, &t1);CHKERRQ(ierr);
  //ierr = MatMult(B, pqr->q, t0);CHKERRQ(ierr);


  
  /* Can do Picard by using the Jacobian that gets made, the rhs that is passed in, and NEWTONLS
       F(x) = A x - b,   J(x) = A,   J dx = F(0)  ==>  A dx = -b,   x = 0 - dx

       A(0) dx_1 = F(0) - b = A(0) 0 - b = -b
         x_1 = 0 - dx_1 = 0 - (-p_1) = p_1
       A(x_1) dx_2 = F(x_1) - b = A(x_1) x_1 - b  ==>  A(x_1) (dx_2 - x_1) = -b
         dx_2 - x_1 = -p_2
         x_2 = x_1 - dx_2 = x_1 - (-p_2 + x_1) = p_2
  */

  //ierr = MatDuplicate(A, MAT_DO_NOT_COPY_VALUES, &J);CHKERRQ(ierr);
  //ierr = VecDuplicate(t0, &t2);CHKERRQ(ierr);
  //PetscErrorCode NonlinearPicard(PetscErrorCode (*lhs)(Vec, Mat, void*), PetscErrorCode (*rhs)(Vec, Vec, void*), Vec guess, void *ctx, Vec sol)
  //PetscErrorCode FormASCNonlinearMatrix(Vec sigma, Mat *A, NonlinearContext *ctx);
  ierr = VecDuplicate(t0, &guess); CHKERRQ(ierr);
  ierr = VecZeroEntries(guess); CHKERRQ(ierr);
  NonlinearContext nctx;
  HContext         hctx = {.alpha = 0.0, .beta = -60, .gamma = -0.5};
  nctx.pqr    = pqr;
  nctx.epsIn  = epsIn;
  nctx.epsOut = epsOut;
  nctx.B      = &B;
  nctx.K      = &K;
  nctx.Bq     = &t0;
  nctx.w      = &w;
  nctx.hctx   = &hctx;
  ierr = NonlinearPicard((PetscErrorCode (*)(Vec, Mat*, void*))&FormASCNonlinearMatrix, (PetscErrorCode (*)(Vec, Vec*, void*))&ASCBq, guess, &nctx, &t1); CHKERRQ(ierr);
  /*
  ierr = SNESCreate(PetscObjectComm((PetscObject) dm), &snes);CHKERRQ(ierr);
  ierr = SNESSetFunction(snes, t2, ComputeBEMResidual, &A);CHKERRQ(ierr);
  ierr = SNESSetJacobian(snes, J, J, ComputeBEMJacobian, &A);CHKERRQ(ierr);
  ierr = SNESSetFromOptions(snes);CHKERRQ(ierr);
  ierr = SNESSolve(snes, t0, t1);CHKERRQ(ierr);
  ierr = SNESDestroy(&snes);CHKERRQ(ierr);
  */


  ierr = MatMult(C, t1, react);CHKERRQ(ierr);
  ierr = VecDestroy(&t0);CHKERRQ(ierr);
  ierr = VecDestroy(&t1);CHKERRQ(ierr);
  //ierr = VecDestroy(&t2);CHKERRQ(ierr);
  //ierr = MatDestroy(&J);CHKERRQ(ierr);
  ierr = MatDestroy(&K);CHKERRQ(ierr);
  ierr = MatDestroy(&B);CHKERRQ(ierr);
  ierr = MatDestroy(&C);CHKERRQ(ierr);
  ierr = PetscLogEventEnd(CalcR_Event, 0, 0, 0, 0);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "makeBEMPcmQualReactionPotential"
/*@
  makeBEMPcmQualReactionPotential - Make the reaction potential, phi_react = Lq = C A^{-1} Bq in the Polarizable Continuum Model

  Input Parameters:
+ epsIn - the dielectric constant inside the protein
. epsOut - the dielectric constant outside the protein
. pqrData - the PQRData context
. coordinates - The vertex coordinates
. w - The vertex weights
- n - The vertex normals

  Output Parameters:
. react - The reaction potential

  Level: beginner

.seealso: doAnalytical()
@*/
PetscErrorCode makeBEMPcmQualReactionPotential(DM dm, BEMType bem, PetscReal epsIn, PetscReal epsOut, PQRData *pqr, Vec coordinates, Vec w, Vec n, Vec react)
{
  const PetscReal epsHat = (epsIn + epsOut)/(epsIn - epsOut);
  SNES            snes;
  Mat             J, A, Bp, B, C;
  Vec             d, t0, t1, t2;
  PetscErrorCode  ierr;

  PetscFunctionBeginUser;
  switch (bem) {
  case BEM_POINT_MF:
    ierr = makeSurfaceToSurfacePointOperators_Laplace(coordinates, w, n, NULL, &A);CHKERRQ(ierr);
    ierr = makeSurfaceToChargePointOperators(coordinates, w, n, pqr, NULL, &B, &C, NULL);CHKERRQ(ierr);
    ierr = PetscLogEventBegin(CalcR_Event, 0, 0, 0, 0);CHKERRQ(ierr);
    /* B = chargesurfop.dphidnCoul */
    ierr = MatDiagonalScale(B, w, NULL);CHKERRQ(ierr);
    ierr = MatScale(B, -1/epsIn);CHKERRQ(ierr);
    break;
  case BEM_PANEL_MF:
    ierr = makeSurfaceToSurfacePanelOperators_Laplace(dm, w, NULL /*n*/, NULL, &A);CHKERRQ(ierr);
    ierr = makeSurfaceToChargePanelOperators(dm, w, NULL /*n*/, pqr, NULL, NULL, &C, &Bp);CHKERRQ(ierr);
    ierr = PetscLogEventBegin(CalcR_Event, 0, 0, 0, 0);CHKERRQ(ierr);
    /* Bp = chargesurfop.dlpToCharges */
    ierr = MatTranspose(Bp, MAT_INITIAL_MATRIX, &B);CHKERRQ(ierr);
    ierr = MatDestroy(&Bp);CHKERRQ(ierr);
    ierr = MatScale(B, -1/epsIn);CHKERRQ(ierr);
    break;
  default:
    SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONG, "Invalid BEM type: %d", bem);
  }
  /* C = chargesurfop.slpToCharges */
  ierr = MatScale(C, 4.0*PETSC_PI);CHKERRQ(ierr);
  /* A = surfsurfop.K */
  ierr = MatTranspose(A, MAT_REUSE_MATRIX, &A);CHKERRQ(ierr);
  ierr = MatDiagonalScale(A, NULL, w);CHKERRQ(ierr);
  ierr = VecDuplicate(w, &d);CHKERRQ(ierr);
  ierr = VecCopy(w, d);CHKERRQ(ierr);
  ierr = VecScale(d, epsHat/2.0);CHKERRQ(ierr);
  ierr = MatDiagonalSet(A, d, ADD_VALUES);CHKERRQ(ierr);
  ierr = VecDestroy(&d);CHKERRQ(ierr);

  ierr = MatCreateVecs(B, NULL, &t0);CHKERRQ(ierr);
  ierr = VecDuplicate(t0, &t1);CHKERRQ(ierr);
  ierr = MatMult(B, pqr->q, t0);CHKERRQ(ierr);
  
  /* Can do Picard by using the Jacobian that gets made, the rhs that is passed in, and NEWTONLS
       F(x) = A x - b,   J(x) = A,   J dx = F(0)  ==>  A dx = -b,   x = 0 - dx

       A(0) dx_1 = F(0) - b = A(0) 0 - b = -b
         x_1 = 0 - dx_1 = 0 - (-p_1) = p_1
       A(x_1) dx_2 = F(x_1) - b = A(x_1) x_1 - b  ==>  A(x_1) (dx_2 - x_1) = -b
         dx_2 - x_1 = -p_2
         x_2 = x_1 - dx_2 = x_1 - (-p_2 + x_1) = p_2
  */

  ierr = MatDuplicate(A, MAT_DO_NOT_COPY_VALUES, &J);CHKERRQ(ierr);
  ierr = VecDuplicate(t0, &t2);CHKERRQ(ierr);
  ierr = SNESCreate(PetscObjectComm((PetscObject) dm), &snes);CHKERRQ(ierr);
  ierr = SNESSetFunction(snes, t2, ComputeBEMResidual, &A);CHKERRQ(ierr);
  ierr = SNESSetJacobian(snes, J, J, ComputeBEMJacobian, &A);CHKERRQ(ierr);
  ierr = SNESSetFromOptions(snes);CHKERRQ(ierr);
  ierr = SNESSolve(snes, t0, t1);CHKERRQ(ierr);
  ierr = SNESDestroy(&snes);CHKERRQ(ierr);

  ierr = MatMult(C, t1, react);CHKERRQ(ierr);
  ierr = VecDestroy(&t0);CHKERRQ(ierr);
  ierr = VecDestroy(&t1);CHKERRQ(ierr);
  ierr = VecDestroy(&t2);CHKERRQ(ierr);
  ierr = MatDestroy(&J);CHKERRQ(ierr);
  ierr = MatDestroy(&A);CHKERRQ(ierr);
  ierr = MatDestroy(&B);CHKERRQ(ierr);
  ierr = MatDestroy(&C);CHKERRQ(ierr);
  ierr = PetscLogEventEnd(CalcR_Event, 0, 0, 0, 0);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "CalculateAnalyticSolvationEnergy"
/*@
  CalculateAnalyticSolvationEnergy - Calculate the solvation energy 1/2 q^T L q

  Input Parameters:
+ epsIn - the dielectric constant inside the protein
. epsOut - the dielectric constant outside the protein
. pqrData - the PQRData context
. R - The sphere radius
- Nmax - The multipole order

  Output Parameters:
+ react - The reaction potential, L q
- E - The solvation energy

  Level: beginner

.seealso: doAnalytical()
@*/
PetscErrorCode CalculateAnalyticSolvationEnergy(PetscReal epsIn, PetscReal epsOut, PQRData *pqr, PetscReal R, PetscInt Nmax, Vec react, PetscReal *E)
{
  const PetscReal q     = ELECTRON_CHARGE;
  const PetscReal Na    = AVOGADRO_NUMBER;
  const PetscReal JperC = 4.184; /* Jouled/Calorie */
  const PetscReal cf    = Na * (q*q/EPSILON_0)/JperC * (1e10/1000) * 1/4/PETSC_PI; /* kcal ang/mol */
  Mat             L;
  PetscErrorCode  ierr;

  PetscFunctionBeginUser;
  PetscValidPointer(pqr, 3);
  PetscValidPointer(E, 6);
  ierr = doAnalytical(R, epsIn, epsOut, pqr, Nmax, &L);CHKERRQ(ierr);
  ierr = MatMult(L, pqr->q, react);CHKERRQ(ierr);
  ierr = MatDestroy(&L);CHKERRQ(ierr);
  ierr = VecDot(pqr->q, react, E);CHKERRQ(ierr);
  *E  *= cf * 0.5;
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "CalculateBEMSolvationEnergy"
/*@
  CalculateBEMSolvationEnergy - Calculate the solvation energy 1/2 q^T L q

  Input Parameters:
+ dm - The DM
. prefix - A prefix to use for the objects created
. bem - The type BEM method
. epsIn - the dielectric constant inside the protein
. epsOut - the dielectric constant outside the protein
. pqrData - the PQRData context
. w - The weights
- n - The normals

  Output Parameters:
+ react - The reaction potential, L q
- E - The solvation energy

  Level: beginner

.seealso: doAnalytical()
@*/
PetscErrorCode CalculateBEMSolvationEnergy(DM dm, const char prefix[], BEMType bem, PetscReal epsIn, PetscReal epsOut, PQRData *pqr, Vec w, Vec n, Vec react, PetscReal *E)
{
  const PetscReal q     = ELECTRON_CHARGE;
  const PetscReal Na    = AVOGADRO_NUMBER;
  const PetscReal JperC = 4.184; /* Jouled/Calorie */
  const PetscReal cf    = Na * (q*q/EPSILON_0)/JperC * (1e10/1000) * 1/4/PETSC_PI; /* kcal ang/mol */
  Mat             L;
  Vec             coords;
  PetscErrorCode  ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(dm, DM_CLASSID, 1);
  PetscValidPointer(pqr, 6);
  PetscValidPointer(E, 10);
  switch (bem) {
  case BEM_POINT:
  case BEM_PANEL:
    ierr = DMGetCoordinatesLocal(dm, &coords);CHKERRQ(ierr);
    ierr = makeBEMPcmQualMatrices(dm, bem, epsIn, epsOut, pqr, coords, w, n, &L);CHKERRQ(ierr);
    ierr = PetscObjectSetName((PetscObject) L, "L");CHKERRQ(ierr);
    ierr = PetscObjectSetOptionsPrefix((PetscObject) L, prefix);CHKERRQ(ierr);
    ierr = MatViewFromOptions(L, NULL, "-mat_view");CHKERRQ(ierr);

    ierr = PetscLogEventBegin(CalcE_Event, L, react, pqr->q, 0);CHKERRQ(ierr);
    ierr = MatMult(L, pqr->q, react);CHKERRQ(ierr);
    ierr = MatDestroy(&L);CHKERRQ(ierr);
    break;
  case BEM_POINT_MF:
  case BEM_PANEL_MF:
    ierr = DMGetCoordinatesLocal(dm, &coords);CHKERRQ(ierr);
    ierr = makeBEMPcmQualReactionPotential2(dm, bem, epsIn, epsOut, pqr, coords, w, n, react);CHKERRQ(ierr);
    L = NULL;
    ierr = PetscLogEventBegin(CalcE_Event, L, react, pqr->q, 0);CHKERRQ(ierr);
    break;
  }
  ierr = VecDot(pqr->q, react, E);CHKERRQ(ierr);
  *E  *= cf * 0.5;
  ierr = PetscLogEventEnd(CalcE_Event, L, react, pqr->q, 0);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc, char **argv)
{
  /* Constants */
  //const PetscReal  q     = ELECTRON_CHARGE;
  //const PetscReal  Na    = AVOGADRO_NUMBER;
  //const PetscReal  JperC = 4.184; /* Jouled/Calorie */
  //const PetscReal  kB    = Na * BOLTZMANN_K/4.184/1000.0; /* Now in kcal/K/mol */
  //const PetscReal  cf    = Na * (q*q/EPSILON_0)/JperC * (1e10/1000) * 1/4/PETSC_PI; /* kcal ang/mol */
  /* Problem data */
  DM               dm;//, dmSimple;
  PQRData          pqr;
  PetscSurface     msp;
  Vec              panelAreas, vertWeights, vertNormals, react;
  PetscReal        totalArea;
  //PetscInt         Np;
  SolvationContext ctx;
  /* Solvation Energies */
  PetscScalar      Eref = 0.0, ESimple = 0.0, ESurf = 0.0, ESurfMF = 0.0, EPanel = 0.0, EMSP = 0.0;
  PetscLogStage    stageSimple, stageSurf, stageSurfMF, stagePanel, stageMSP;
  PetscErrorCode   ierr;

  ierr = PetscInitialize(&argc, &argv, NULL, NULL);CHKERRQ(ierr);
  ierr = PetscLogDefaultBegin();CHKERRQ(ierr);
  ierr = ProcessOptions(PETSC_COMM_WORLD, &ctx);CHKERRQ(ierr);
  /* Make PQR */
  if (ctx.isSphere) {
    ierr = makeSphereChargeDistribution(ctx.R, ctx.numCharges, ctx.h, PETSC_DETERMINE, &pqr);CHKERRQ(ierr);
    ierr = PQRViewFromOptions(&pqr);CHKERRQ(ierr);
  } else {
    ierr = PQRCreateFromPDB(PETSC_COMM_WORLD, ctx.pdbFile, ctx.crgFile, &pqr);CHKERRQ(ierr);
  }
  ierr = VecDuplicate(pqr.q, &react);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) react, "Reaction Potential");CHKERRQ(ierr);
  /* Make surface */
  ierr = loadSrfIntoSurfacePoints(PETSC_COMM_WORLD, ctx.srfFile, &vertNormals, &vertWeights, &panelAreas, &totalArea, &dm);CHKERRQ(ierr);
  {
    PetscInt cStart, cEnd, vStart, vEnd;

    ierr = DMPlexGetHeightStratum(dm, 0, &cStart, &cEnd);CHKERRQ(ierr);
    ierr = DMPlexGetDepthStratum(dm, 0, &vStart, &vEnd);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD, "SRF %D vertices %D cells\n", vEnd-vStart, cEnd-cStart);CHKERRQ(ierr);
  }
  ierr = PetscSurfaceCreateMSP(PETSC_COMM_WORLD, ctx.pntFile, &msp);CHKERRQ(ierr);
  if (msp.weights) {
    PetscInt Nv;

    ierr = VecGetSize(msp.weights, &Nv);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD, "MSP %D vertices\n", Nv);CHKERRQ(ierr);
  }
  /* Calculate solvation energy */
  ierr = PetscLogStageRegister("Point Surface", &stageSurf);CHKERRQ(ierr);
  ierr = PetscLogStageRegister("Point Surface MF", &stageSurfMF);CHKERRQ(ierr);
  ierr = PetscLogStageRegister("Panel Surface", &stagePanel);CHKERRQ(ierr);
  ierr = PetscLogStagePush(stageSurf);CHKERRQ(ierr);
  ierr = CalculateBEMSolvationEnergy(dm, "lsrf_", BEM_POINT, ctx.epsIn, ctx.epsOut, &pqr, vertWeights, vertNormals, react, &ESurf);CHKERRQ(ierr);
  ierr = PetscLogStagePop();CHKERRQ(ierr);
  ierr = PetscLogStagePush(stageSurfMF);CHKERRQ(ierr);
  ierr = CalculateBEMSolvationEnergy(dm, "lsrf_mf_", BEM_POINT_MF, ctx.epsIn, ctx.epsOut, &pqr, vertWeights, vertNormals, react, &ESurfMF);CHKERRQ(ierr);
  ierr = PetscLogStagePop();CHKERRQ(ierr);
  ierr = PetscLogStagePush(stagePanel);CHKERRQ(ierr);
  ierr = CalculateBEMSolvationEnergy(dm, "lpanel_", BEM_PANEL, ctx.epsIn, ctx.epsOut, &pqr, panelAreas, vertNormals, react, &EPanel);CHKERRQ(ierr);
  ierr = PetscLogStagePop();CHKERRQ(ierr);
  if (msp.weights) {
    ierr = PetscLogStageRegister("MSP Surface", &stageMSP);CHKERRQ(ierr);
    ierr = PetscLogStagePush(stageMSP);CHKERRQ(ierr);
    ierr = CalculateBEMSolvationEnergy(msp.dm, "lmsp_", BEM_POINT, ctx.epsIn, ctx.epsOut, &pqr, msp.weights, msp.normals, react, &EMSP);CHKERRQ(ierr);
    ierr = PetscLogStagePop();CHKERRQ(ierr);
  }
  /* Verification */
  if (ctx.isSphere) {
    const PetscInt Np = PetscCeilReal(4.0 * PETSC_PI * PetscSqr(ctx.R))*ctx.density;
    DM             dmSimple;
    Vec            vertWeightsSimple, vertNormalsSimple;

    ierr = PetscPrintf(PETSC_COMM_WORLD, "Total area: %g Sphere area: %g\n", totalArea, 4*PETSC_PI*PetscPowRealInt(ctx.R, 2));CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD, "Simple %D vertices\n", Np);CHKERRQ(ierr);
    ierr = makeSphereSurface(PETSC_COMM_WORLD, ctx.origin, ctx.R, Np, &vertWeightsSimple, &vertNormalsSimple, NULL, &dmSimple);CHKERRQ(ierr);

    ierr = PetscLogStageRegister("Simple Surface", &stageSimple);CHKERRQ(ierr);
    ierr = PetscLogStagePush(stageSimple);CHKERRQ(ierr);
    ierr = CalculateBEMSolvationEnergy(dmSimple, "lsimple_", BEM_POINT, ctx.epsIn, ctx.epsOut, &pqr, vertWeightsSimple, vertNormalsSimple, react, &ESimple);CHKERRQ(ierr);
    ierr = PetscLogStagePop();CHKERRQ(ierr);
    ierr = CalculateAnalyticSolvationEnergy(ctx.epsIn, ctx.epsOut, &pqr, ctx.R, ctx.Nmax, react, &Eref);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD, "Eref = %.6f ESurf   = %.6f Error = %.6f Rel. error = %.4f\n", Eref, ESurf,   Eref-ESurf,   (Eref-ESurf)/Eref);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD, "Eref = %.6f ESurfMF = %.6f Error = %.6f Rel. error = %.4f\n", Eref, ESurfMF, Eref-ESurfMF, (Eref-ESurfMF)/Eref);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD, "Eref = %.6f ESimple = %.6f Error = %.6f Rel. error = %.4f\n", Eref, ESimple, Eref-ESimple, (Eref-ESimple)/Eref);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD, "Eref = %.6f EPanel  = %.6f Error = %.6f Rel. error = %.4f\n", Eref, EPanel,  Eref-EPanel,  (Eref-EPanel)/Eref);CHKERRQ(ierr);

    ierr = VecDestroy(&vertWeightsSimple);CHKERRQ(ierr);
    ierr = VecDestroy(&vertNormalsSimple);CHKERRQ(ierr);
    ierr = DMDestroy(&dmSimple);CHKERRQ(ierr);
  } else {
    Eref = 1.0; /* TODO This should be higher resolution BEM */
    ierr = PetscPrintf(PETSC_COMM_WORLD, "Eref = %.6f ESurf  = %.6f Error = %.6f Rel. error = %.4f\n", Eref, ESurf,  Eref-ESurf,  (Eref-ESurf)/Eref);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD, "Eref = %.6f ESurfMF = %.6f Error = %.6f Rel. error = %.4f\n", Eref, ESurfMF, Eref-ESurfMF, (Eref-ESurfMF)/Eref);CHKERRQ(ierr);
    //ierr = PetscPrintf(PETSC_COMM_WORLD, "Eref = %.6f ESimple = %.6f Error = %.6f Rel. error = %.4f\n", Eref, ESimple, Eref-ESimple, (Eref-ESimple)/Eref);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD, "Eref = %.6f EPanel = %.6f Error = %.6f Rel. error = %.4f\n", Eref, EPanel, Eref-EPanel, (Eref-EPanel)/Eref);CHKERRQ(ierr);
    if (msp.weights) {ierr = PetscPrintf(PETSC_COMM_WORLD, "Eref = %.6f EMSP   = %.6f Error = %.6f Rel. error = %.4f\n", Eref, EMSP,   Eref-EMSP,   (Eref-EMSP)/Eref);CHKERRQ(ierr);}
  }
  /* Output flops */
  {
    PetscStageLog stageLog;

    ierr = PetscLogGetStageLog(&stageLog);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD, "Flops_Surf  = %.4e Flops_S2S_Surf  = %.4e\n",
                       stageLog->stageInfo[stageSurf].perfInfo.flops, stageLog->stageInfo[stageSurf].eventLog->eventInfo[CalcStoS_Event].flops);CHKERRQ(ierr);
    if (ctx.isSphere) {
      ierr = PetscPrintf(PETSC_COMM_WORLD, "Flops_SurfMF= %.4e Flops_S2S_Surf  = %.4e\n",
                         stageLog->stageInfo[stageSurfMF].perfInfo.flops, stageLog->stageInfo[stageSurfMF].eventLog->eventInfo[CalcStoS_Event].flops);CHKERRQ(ierr);
    }
    ierr = PetscPrintf(PETSC_COMM_WORLD, "Flops_Panel = %.4e Flops_S2S_Panel = %.4e\n",
                       stageLog->stageInfo[stagePanel].perfInfo.flops, stageLog->stageInfo[stagePanel].eventLog->eventInfo[CalcStoS_Event].flops);CHKERRQ(ierr);
    if (msp.weights) {ierr = PetscPrintf(PETSC_COMM_WORLD, "Flops_MSP   = %.4e Flops_S2S_MSP   = %.4e\n",
                                         stageLog->stageInfo[stageMSP].perfInfo.flops, stageLog->stageInfo[stageMSP].eventLog->eventInfo[CalcStoS_Event].flops);CHKERRQ(ierr);
    }
    if (ctx.isSphere) {ierr = PetscPrintf(PETSC_COMM_WORLD, "Flops_Simple = %.4e Flops_S2S_Simple = %.4e\n",
                                          stageLog->stageInfo[stageSimple].perfInfo.flops, stageLog->stageInfo[stageSimple].eventLog->eventInfo[CalcStoS_Event].flops);CHKERRQ(ierr);
    }
  }
  /* Cleanup */
  ierr = VecDestroy(&vertWeights);CHKERRQ(ierr);
  ierr = VecDestroy(&vertNormals);CHKERRQ(ierr);
  ierr = VecDestroy(&panelAreas);CHKERRQ(ierr);
  ierr = VecDestroy(&react);CHKERRQ(ierr);
  ierr = DMDestroy(&dm);CHKERRQ(ierr);
  ierr = PetscSurfaceDestroy(&msp);CHKERRQ(ierr);
  ierr = PQRDestroy(&pqr);CHKERRQ(ierr);
  ierr = PetscFinalize();
  return 0;
}
