#include <petsc.h>
#include "constants.h"
#include "molecule.h"
#include "sphere.h"

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

//#undef __FUNCT__
//#define __FUNCT__ "legendre"
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
//PetscErrorCode legendre(PetscInt l, PetscInt m, PetscScalar x, PetscScalar *leg, PetscScalar *err)
//{
//  /* The error estimate is based on worst case scenario and the significant digits, and thus
//     based on the largest polynomial coefficient and machine error, "eps" */
//  PetscInt       maxcf;    /* largest polynomial coefficient */
//  PetscReal      cfnm = 1; /* proportionality constant for m < 0 polynomials compared to m > 0 */
//  PetscReal      cl, x2, p, f1, f2, f3;
//  PetscInt       px, j;
//  PetscErrorCode ierr;
//
//  PetscFunctionBeginUser;
//  /* The polynomials are not defined for |x| > 1 */
//  if (PetscAbsScalar(x) > 1) SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE, "Invalid input |x %g| > 1", PetscRealPart(x));
//  /* Could also define this to be 0 */
//  if (PetscAbsInt(m) > PetscAbsInt(l)) SETERRQ2(PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE, "Invalid input m %d > l %d", m, l);
//  if (l < 0) l = -(l+1);
//  if (m < 0) {
//    PetscReal num, den;

//    m    = -m;
//    ierr = factorial(l-m, &num);CHKERRQ(ierr);
//    ierr = factorial(l+m, &den);CHKERRQ(ierr);
//    cfnm = PetscPowInt(-1, m)*num/den;
//  }
//  /* Calculate coef of maximum degree in x from the explicit analytical formula */
//  ierr  = factorial(2*l, &f1);CHKERRQ(ierr);
//  ierr  = factorial(l,   &f2);CHKERRQ(ierr);
//  ierr  = factorial(l-m, &f3);CHKERRQ(ierr);
//  cl    = PetscPowInt(-1, m) * cfnm * f1/((1 << l)*f2*f3);
//  maxcf = PetscAbsInt(cl);
//  px    = l-m;
//  /* Power of x changes from one term to the next by 2. Also needed for sqrt(1-x^2). */
//  x2    = x*x; /* TODO make pointwise square */
//  /* Calculate efficiently P_l^m (x)/sqrt(1-x^2)^(m/2) - that is, only the polynomial part.
//     At least one coefficient is guaranteed to exist - there is no null Legendre polynomial. */
//  p     = cl; /* TODO make an array of cl */
// for (j = l-1; j >= 0; --j) {
//    /* Check the exponent of x for current coefficient, px. If it is 0 or 1, just exit the loop */
//    if (px < 2) break;
//    /* If current exponent is >=2, there is a "next" coefficient; multiply p by x2 and add it. Calculate the current coefficient */
//    cl = -(j+j+2-l-m)*(j+j+1-l-m)/(2*(j+j+1)*(l-j))*cl;
//    
//    if (maxcf < PetscAbsReal(cl)) maxcf = PetscAbsReal(cl);
//    /* ...and add to the polynomial */
//    p = p*x2 + cl; /* TODO make this pointwise multiply */
//    /* Decrease the exponent of x - this is the exponent of x corresponding to the newly added coefficient */
//    px -= 2;
//  }
//  /* Estimate the error */
//  if (err) *err = maxcf*PETSC_MACHINE_EPSILON;

  /* Now we're done adding coefficients. However, if the exponent of x
     corresponding to the last added coefficient is 1 (polynomial is odd),
     multiply the polynomial by x */
//  if (px == 1) p = p*x;

  /* All that's left is to multiply the whole thing with sqrt(1-x^2)^(m/2). No further calculations are needed if m = 0. */
//  if (m == 0) {*leg = p; PetscFunctionReturn(0);}

//  x2 = 1-x2;
//  /* First, multiply by the integer part of m/2 */
//  for (j = 1; j < PetscFloorReal(m/2.0); ++j) p = p*x2; /* TODO make this pointwise multiply */
//  /* If m is odd, there is an additional factor sqrt(1-x^2) */
//  if (m % 2) p = p*PetscSqrtReal(x2); /* TODO make this pointwise multiply */
//  *leg = p;
//  PetscFunctionReturn(0);
//}


#undef __FUNCT__
#define __FUNCT__ "convertToSpherical"
PetscErrorCode convertToSpherical(const PetscScalar xyz[], PetscScalar r[])
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
#define __FUNCT__ "computePotentialSpherical"
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
PetscErrorCode computePotentialSpherical(PQRData *pqr, PetscInt Nmax, Vec Bnm, Vec phi)
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

