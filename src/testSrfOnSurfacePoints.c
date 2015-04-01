#include <petsc.h>
#include <constants.h>
#include <surface.h>

typedef struct {
  Vec q;   /* Charge values */
  Vec xyz; /* Charge coordinates, always 3D */
  Vec R;   /* Charge radii */
} PQRData;

typedef struct {
  PetscInt numCharges; /* Number of atomic charges in the protein */
  PetscInt Nmax;       /* Order of the multipole expansion */
} SolvationContext;

#undef __FUNCT__
#define __FUNCT__ "ProcessOptions"
PetscErrorCode ProcessOptions(MPI_Comm comm, SolvationContext *ctx)
{
  PetscErrorCode ierr;

  PetscFunctionBeginUser;
  ctx->numCharges = 100;
  ctx->Nmax       = 100;

  ierr = PetscOptionsBegin(comm, "", "Solvation Problem Options", "BIBEE");CHKERRQ(ierr);
    ierr = PetscOptionsInt("-num_charges", "The number of atomic charges in the protein", "testSrfOnSurfacePoints", ctx->numCharges, &ctx->numCharges, PETSC_NULL);CHKERRQ(ierr);
    ierr = PetscOptionsInt("-nmax",        "The order of the multipole expansion", "testSrfOnSurfacePoints", ctx->Nmax, &ctx->Nmax, PETSC_NULL);CHKERRQ(ierr);
  ierr = PetscOptionsEnd();
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "PQRView"
PetscErrorCode PQRView(PQRData *pqr, PetscViewer viewer)
{
  PetscErrorCode ierr;

  PetscFunctionBeginUser;
  ierr = VecView(pqr->xyz, viewer);CHKERRQ(ierr);
  ierr = VecView(pqr->q,   viewer);CHKERRQ(ierr);
  ierr = VecView(pqr->R,   viewer);CHKERRQ(ierr);
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
  PetscReal      sqz2   = PetscSqrtReal(1.0 - PetscSqr(PetscRealPart(z)));
  PetscReal      hsqz2  = 0.5*sqz2;
  PetscReal      ihsqz2 = PetscRealPart(z)/hsqz2;
  PetscReal      fac    = 1.0;
  PetscInt       pre    = l % 2 ? -1 : 1;
  PetscInt       m;
  PetscErrorCode ierr;

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

    ierr = PetscOptionsGetRealArray(NULL, "-test", vals, &nmax, &flg);CHKERRQ(ierr);
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
  PetscReal      epsHat = 2.0*(epsIn - epsOut)/(epsIn + epsOut);
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

    ierr = PetscPrintf(PETSC_COMM_SELF, "Computing charge %d\n", q);CHKERRQ(ierr);
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
#define __FUNCT__ "makeSurfaceToChargeOperators"
/*@
  makeSurfaceToChargeOperators - Make matrices???

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
PetscErrorCode makeSurfaceToChargeOperators(Vec coordinates, Vec w, Vec n, PQRData *pqr, Mat *potential, Mat *field, Mat *singleLayer, Mat *doubleLayer)
{
  const PetscScalar *coords, *xyz, *weights, *normals;
  PetscInt           Nq, Np;
  PetscInt           i, j, d;
  PetscErrorCode     ierr;

  PetscFunctionBeginUser;
  ierr = VecGetLocalSize(pqr->q, &Nq);CHKERRQ(ierr);
  ierr = VecGetLocalSize(w, &Np);CHKERRQ(ierr);
  if (potential)   {ierr = MatCreateSeqDense(PETSC_COMM_SELF, Np, Nq, NULL, potential);CHKERRQ(ierr);}
  if (field)       {ierr = MatCreateSeqDense(PETSC_COMM_SELF, Np, Nq, NULL, field);CHKERRQ(ierr);}
  if (singleLayer) {ierr = MatCreateSeqDense(PETSC_COMM_SELF, Nq, Np, NULL, singleLayer);CHKERRQ(ierr);}
  if (doubleLayer) {ierr = MatCreateSeqDense(PETSC_COMM_SELF, Nq, Np, NULL, doubleLayer);CHKERRQ(ierr);}
  PetscPrintf(PETSC_COMM_SELF, "Starting\n");
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

      if (r < 1e-10) {G = 0;                dGdn = 0;}
      else           {G = 1.0/4/PETSC_PI/r; dGdn = -dot/4/PETSC_PI/PetscPowRealInt(r, 3);}

      if (potential)   {ierr = MatSetValue(*potential,   i, j, G,               INSERT_VALUES);CHKERRQ(ierr);}
      if (field)       {ierr = MatSetValue(*field,       i, j, dGdn,            INSERT_VALUES);CHKERRQ(ierr);}
      if (singleLayer) {ierr = MatSetValue(*singleLayer, j, i, G*weights[i],    INSERT_VALUES);CHKERRQ(ierr);}
      if (doubleLayer) {ierr = MatSetValue(*doubleLayer, j, i, dGdn*weights[i], INSERT_VALUES);CHKERRQ(ierr);}
    }
  }
  ierr = VecRestoreArrayRead(coordinates, &coords);CHKERRQ(ierr);
  ierr = VecRestoreArrayRead(pqr->xyz, &xyz);CHKERRQ(ierr);
  ierr = VecRestoreArrayRead(w, &weights);CHKERRQ(ierr);
  ierr = VecRestoreArrayRead(n, &normals);CHKERRQ(ierr);
  if (potential)   {ierr = MatAssemblyBegin(*potential,   MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);ierr = MatAssemblyEnd(*potential,   MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);}
  if (field)       {ierr = MatAssemblyBegin(*field,       MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);ierr = MatAssemblyEnd(*field,       MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);}
  if (singleLayer) {ierr = MatAssemblyBegin(*singleLayer, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);ierr = MatAssemblyEnd(*singleLayer, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);}
  if (doubleLayer) {ierr = MatAssemblyBegin(*doubleLayer, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);ierr = MatAssemblyEnd(*doubleLayer, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);}
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "makeSurfaceToSurfaceLaplaceOperators"
/*@
  makeSurfaceToSurfaceLaplaceOperators - Make matrices???

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
PetscErrorCode makeSurfaceToSurfaceLaplaceOperators(Vec coordinates, Vec w, Vec n, Mat *V, Mat *K)
{
  const PetscScalar *coords, *weights, *normals;
  PetscInt           Np, i, j, d;
  PetscErrorCode     ierr;

  PetscFunctionBeginUser;
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
        if (V) {ierr = MatSetValue(*V, i, j, weights[j]* 1/4/PETSC_PI/r, INSERT_VALUES);CHKERRQ(ierr);}
        if (K) {ierr = MatSetValue(*K, i, j, weights[j]* dot/4/PETSC_PI/PetscPowRealInt(r,3), INSERT_VALUES);CHKERRQ(ierr);}
      } else {
        const PetscReal R0 = PetscSqrtReal(weights[j]/PETSC_PI); /* radius of a circle with the area assoc with surfpt j */
        if (V) {ierr = MatSetValue(*V, i, j, (2 * PETSC_PI * R0) /4/PETSC_PI, INSERT_VALUES);CHKERRQ(ierr);}
      }
    }
  }
  ierr = VecRestoreArrayRead(coordinates, &coords);CHKERRQ(ierr);
  ierr = VecRestoreArrayRead(w, &weights);CHKERRQ(ierr);
  ierr = VecRestoreArrayRead(n, &normals);CHKERRQ(ierr);
  if (V) {ierr = MatAssemblyBegin(*V, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);ierr = MatAssemblyEnd(*V, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);}
  if (K) {ierr = MatAssemblyBegin(*K, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);ierr = MatAssemblyEnd(*K, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);}
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "makeBEMEcfQualMatrices"
/*@
  makeBEMEcfQualMatrices - Make matrices???

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
PetscErrorCode makeBEMEcfQualMatrices(PetscReal epsIn, PetscReal epsOut, PQRData *pqr, Vec coordinates, Vec w, Vec n, Mat *L)
{
  const PetscReal epsHat = (epsIn + epsOut)/(epsIn - epsOut);
  KSP             ksp;
  PC              pc;
  Mat             A, B, C, S, fact;
  Vec             d;
  PetscErrorCode  ierr;

  PetscFunctionBeginUser;
  ierr = makeSurfaceToSurfaceLaplaceOperators(coordinates, w, n, NULL, &A);CHKERRQ(ierr);
  ierr = makeSurfaceToChargeOperators(coordinates, w, n, pqr, NULL, &B, &C, NULL);CHKERRQ(ierr);

  /* B = chargesurfop.dphidnCoul */
  ierr = MatDiagonalScale(B, w, NULL);CHKERRQ(ierr);
  ierr = MatScale(B, -1/epsIn);CHKERRQ(ierr);

  /* C = chargesurfop.slpToCharges */
  ierr = MatScale(C, 4.0*PETSC_PI);CHKERRQ(ierr);
  /* A = surfsurfop.K */
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
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc, char **argv)
{
  DM               dm;
  PQRData          pqr;
  Vec              vertCoords, vertWeights, vertNormals, vertCoordsSimple, vertWeightsSimple, vertNormalsSimple, vertAnglesSimple;
  SolvationContext ctx;
  PetscReal      q  = ELECTRON_CHARGE;
  PetscReal      Na = 6.0221415e23;
  PetscReal      kB = Na * BOLTZMANN_K/4.184/1000.0; /* Now in kcal/K/mol */
  PetscReal      R           = 6.0;
  PetscReal      epsIn       =  4;
  PetscReal      epsOut      = 80;
  PetscReal      conv_factor = 332.112;
  PetscReal      origin[3]   = {0.0, 0.0, 0.0};
  PetscReal      density     = 1.0;
  PetscReal      h           = 1.0;
  PetscInt       numPoints   = ceil(4.0 * PETSC_PI * R*R);
  char           srfFile[PETSC_MAX_PATH_LEN];
  PetscScalar    Eref, ESimple = 0.0, ESRF = 0.0;
  Vec            react;
  Mat            Lref, LSRF, LSimple;
  PetscErrorCode ierr;

  ierr = PetscInitialize(&argc, &argv, NULL, NULL);CHKERRQ(ierr);
  ierr = ProcessOptions(PETSC_COMM_WORLD, &ctx);CHKERRQ(ierr);
  ierr = PetscStrcpy(srfFile, "../../jay-pointbem/geometry/sphere_R6_vdens1.srf");CHKERRQ(ierr);
  ierr = makeSphereChargeDistribution(R, ctx.numCharges, h, PETSC_DETERMINE, &pqr);
  ierr = PQRView(&pqr, PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);

  ierr = loadSrfIntoSurfacePoints(PETSC_COMM_WORLD, srfFile, &vertNormals, &vertWeights, &dm);CHKERRQ(ierr);
  ierr = makeSphereSurface(PETSC_COMM_WORLD, origin, R, numPoints, &vertCoordsSimple, &vertWeightsSimple, &vertNormalsSimple, &vertAnglesSimple);CHKERRQ(ierr);

  ierr = DMGetCoordinatesLocal(dm, &vertCoords);CHKERRQ(ierr);
  ierr = makeBEMEcfQualMatrices(epsIn, epsOut, &pqr, vertCoords,       vertWeights,       vertNormals,       &LSRF);CHKERRQ(ierr);
  ierr = makeBEMEcfQualMatrices(epsIn, epsOut, &pqr, vertCoordsSimple, vertWeightsSimple, vertNormalsSimple, &LSimple);CHKERRQ(ierr);
  ierr = doAnalytical(R, epsIn, epsOut, &pqr, ctx.Nmax, &Lref);CHKERRQ(ierr);

  ierr = PetscObjectSetName((PetscObject) LSRF, "LSRF");CHKERRQ(ierr);
  ierr = PetscObjectSetOptionsPrefix((PetscObject) LSRF, "lsrf_");CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) LSimple, "LSimple");CHKERRQ(ierr);
  ierr = PetscObjectSetOptionsPrefix((PetscObject) LSimple, "lsimple_");CHKERRQ(ierr);
  ierr = MatViewFromOptions(LSRF, NULL, "-mat_view");CHKERRQ(ierr);

  ierr = VecDuplicate(pqr.q, &react);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) react, "Reaction Potential");CHKERRQ(ierr);
  ierr = MatMult(Lref, pqr.q, react);CHKERRQ(ierr);
  ierr = VecDot(pqr.q, react, &Eref);CHKERRQ(ierr);
  Eref    *= conv_factor * 0.5;
  ierr = MatMult(LSRF, pqr.q, react);CHKERRQ(ierr);
  ierr = VecDot(pqr.q, react, &ESRF);CHKERRQ(ierr);
  ESRF    *= conv_factor * 0.5;
  ierr = MatMult(LSimple, pqr.q, react);CHKERRQ(ierr);
  ierr = VecDot(pqr.q, react, &ESimple);CHKERRQ(ierr);
  ESimple *= conv_factor * 0.5;

  ierr = PetscPrintf(PETSC_COMM_WORLD, "Eref = %.6f ESRF    = %.6f Error = %.6f Rel. error = %.4f\n", Eref, ESRF, Eref-ESRF, (Eref-ESRF)/Eref);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "Eref = %.6f ESimple = %.6f Error = %.6f Rel. error = %.4f\n", Eref, ESimple, Eref-ESimple, (Eref-ESimple)/Eref);

  ierr = VecDestroy(&vertCoordsSimple);CHKERRQ(ierr);
  ierr = VecDestroy(&vertWeightsSimple);CHKERRQ(ierr);
  ierr = VecDestroy(&vertNormalsSimple);CHKERRQ(ierr);
  ierr = VecDestroy(&vertAnglesSimple);CHKERRQ(ierr);
  ierr = VecDestroy(&vertWeights);CHKERRQ(ierr);
  ierr = VecDestroy(&vertNormals);CHKERRQ(ierr);
  ierr = DMDestroy(&dm);CHKERRQ(ierr);
  ierr = MatDestroy(&Lref);CHKERRQ(ierr);
  ierr = MatDestroy(&LSRF);CHKERRQ(ierr);
  ierr = MatDestroy(&LSimple);CHKERRQ(ierr);
  ierr = VecDestroy(&react);CHKERRQ(ierr);
  ierr = PQRDestroy(&pqr);CHKERRQ(ierr);
  ierr = PetscFinalize();
  return 0;
}
