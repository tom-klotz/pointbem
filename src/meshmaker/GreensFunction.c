#include "meshmaker.h"

// GRYCUK_KERNEL evaluated using LJ6!!
real GreensFunction(Vector3D dest, Vector3D src, BEMKernelType kerneltype, void* parameters) {
   switch (kerneltype) {
   case CONSTANT_KERNEL: return 1.0;
   case X_KERNEL: return src->x + 1.0;
   case Y_KERNEL: return src->y + 1.0;
   case Z_KERNEL: return src->z + 1.0;
   case POISSON_KERNEL: return GreensFunction_oneoverr(dest, src, parameters);
   case HELMHOLTZ_KERNEL: return GreensFunction_ekroverr(dest, src, parameters);
   case DESINGULARIZED_HELMHOLTZ_KERNEL: return GreensFunction_ekroverr_desingularized(dest, src, parameters);
   case LJ_KERNEL: return GreensFunction_LJ(dest, src, parameters);
   case LJ12_KERNEL: return GreensFunction_LJ12(dest, src, parameters);
   case LJ6_KERNEL: return GreensFunction_LJ6(dest, src, parameters);
   case MONOMIAL_KERNEL: return GreensFunction_monomial(dest, src, parameters);
   case INVERSEPOWER_KERNEL: return GreensFunction_inversepower(dest, src, parameters);
   case GHOSH_KERNEL: return GreensFunction_Ghosh(dest, src, parameters);
   case GRYCUK_KERNEL: return GreensFunction_LJ6(dest, src, parameters);
   case LESYNG_KERNEL: return GreensFunction_Lesyng(dest, src, parameters);
	case XGB_KERNEL: return GreensFunction_XGB(dest, src, parameters);
   default: printf("No Green's Function!\n"); return 1.0;
   }
}

// GRYCUK_KERNEL evaluated using LJ6 functions!!
real GreensFunction_deriv(Vector3D dest, Vector3D src, BEMKernelType kerneltype, void* parameters, Vector3D direction) {
   switch (kerneltype) {
   case POISSON_KERNEL: return GreensFunction_oneoverr_deriv(dest, src, parameters, direction);
   case HELMHOLTZ_KERNEL: return GreensFunction_ekroverr_deriv(dest, src, parameters, direction);
   case DESINGULARIZED_HELMHOLTZ_KERNEL: return GreensFunction_ekroverr_desingularized_deriv(dest, src, parameters, direction);
   case LJ12_KERNEL: return GreensFunction_LJ12_deriv(dest, src, parameters, direction);
   case LJ6_KERNEL: return GreensFunction_LJ6_deriv(dest, src, parameters, direction);
   case MONOMIAL_KERNEL: return GreensFunction_monomial_deriv(dest, src, parameters, direction);
   case INVERSEPOWER_KERNEL: return GreensFunction_inversepower_deriv(dest, src, parameters, direction);
   case GHOSH_KERNEL: return GreensFunction_Ghosh_deriv(dest, src, parameters, direction);
   case GRYCUK_KERNEL: return GreensFunction_LJ6_deriv(dest, src, parameters, direction);
   case LESYNG_KERNEL: return GreensFunction_Lesyng_deriv(dest, src, parameters, direction);
   case VOLUME_KERNEL: return Vector3D_dot(src, direction) / 3.0;
	case XGB_KERNEL: return GreensFunction_XGB_deriv(dest, src, parameters, direction);
   default: printf("No Derivative Greens's Function!\n"); return 1.0;
   }
}
real GreensFunction_doublederiv(Vector3D dest, Vector3D src, BEMKernelType kerneltype,
										  void* parameters, Vector3D direction, Vector3D direction2)
{
  real G = 0.0;
  Vector3D rvector = Vector3D_allocate();
  real r;
  if (kerneltype != JUFFER_KERNEL) {
	 printf("GreensFunction_doublederiv is only defined for JUFFER_KERNEL, because it has no hypersingularity.\n");
	 exit(-1);
  }
  Vector3D_sub(rvector, dest, src);
  r = Vector3D_length(rvector);
  Vector3D_normalize(rvector);
  
  real kappa = ((real *)parameters)[0];
  real cosTheta = Vector3D_dot(direction, rvector) / r;
  real cosTheta0 = Vector3D_dot(direction2, rvector) / r;
  if (r < 1e-6) {
	 G += kappa*kappa*(cosTheta * cosTheta0 - Vector3D_dot(direction, direction2))/(8.0 * M_PI * r);
	 G += PetscPowRealInt(kappa, 3) * Vector3D_dot(direction, direction2) / (12.0 * M_PI);
  } else {
	 G += kappa * r * exp(-kappa * r) *
		((Vector3D_dot(direction, direction2) - 3 * cosTheta * cosTheta0)/ (4.0 * M_PI * PetscPowRealInt(r, 3)));
	 G += -kappa * kappa * exp(-kappa * r) * cosTheta * cosTheta0 / (4.0 * M_PI * r);
  }
  
  Vector3D_free(rvector);
  return G;
}

real GreensFunction_oneoverr(Vector3D dest, Vector3D src, void* parameters) {
   real r = Vector3D_distance(dest, src);
   if (r == 0.0)
      return 0.0;
   else
      return 1.0 / r;
}

/*
real GreensFunction_oneoverr_deriv(Vector3D dest, Vector3D src, void* parameters, Vector3D direction) {
   Vector3D rvector = Vector3D_allocate();
   Vector3D_sub(rvector, dest, src);
   real r = Vector3D_length(rvector);
   Vector3D_normalize(rvector);

   real g = -1.0 / (r * r) * Vector3D_dot(rvector, direction);
   Vector3D_free(rvector);

   return g;
}
*/

real GreensFunction_oneoverr_deriv(Vector3D dest, Vector3D src, void* parameters, Vector3D direction) {
   struct _Vector3D rvector;

   rvector.x = dest->x - src->x;
   rvector.y = dest->y - src->y;
   rvector.z = dest->z - src->z;

   real r2 = rvector.x*rvector.x + rvector.y*rvector.y + rvector.z*rvector.z;
   real r = sqrt(r2);

   rvector.x /= r;
   rvector.y /= r;
   rvector.z /= r;

   real dot = rvector.x*direction->x + rvector.y*direction->y + rvector.z*direction->z;

   return -dot / r2;
}

real GreensFunction_ekroverr(Vector3D dest, Vector3D src, void* parameters) {
   real r = Vector3D_distance(dest, src);
   return exp(-((real*)parameters)[0] * r) / r;
}

real GreensFunction_ekroverr_deriv(Vector3D dest, Vector3D src, void* parameters, Vector3D direction) {
   struct _Vector3D rvector;

   rvector.x = dest->x - src->x;
   rvector.y = dest->y - src->y;
   rvector.z = dest->z - src->z;

   real r2 = rvector.x*rvector.x + rvector.y*rvector.y + rvector.z*rvector.z;
   real r = sqrt(r2);

   rvector.x /= r;
   rvector.y /= r;
   rvector.z /= r;

   real dot = rvector.x*direction->x + rvector.y*direction->y + rvector.z*direction->z;

   real kappa = ((real*)parameters)[0];

   return -dot * exp(-kappa * r) * (kappa * r + 1.0) / r2;
}

real GreensFunction_ekroverr_desingularized(Vector3D dest, Vector3D src, void* parameters) {
   real r = Vector3D_distance(dest, src);
   return (exp(-((real*)parameters)[0] * r) - 1.0)/ r;
}

real GreensFunction_ekroverr_desingularized_deriv(Vector3D dest, Vector3D src, void* parameters, Vector3D direction) {
   struct _Vector3D rvector;

   rvector.x = dest->x - src->x;
   rvector.y = dest->y - src->y;
   rvector.z = dest->z - src->z;

   real r2 = rvector.x*rvector.x + rvector.y*rvector.y + rvector.z*rvector.z;
   real r = sqrt(r2);

   rvector.x /= r;
   rvector.y /= r;
   rvector.z /= r;

   real dot = rvector.x*direction->x + rvector.y*direction->y + rvector.z*direction->z;

   real kappa = ((real*)parameters)[0];
   real ekr = exp(-kappa * r);
   return dot * (-kappa * ekr / r - (ekr - 1.0) / r2);
}

real GreensFunction_LJ(Vector3D dest, Vector3D src, void* parameters) {
   real r = Vector3D_distance(dest, src);
   return -((real*)parameters)[1] / (90.0 * pow(r, 10.0)) +
      ((real*)parameters)[0] / (12.0 * pow(r, 4.0));
}

real GreensFunction_LJ12(Vector3D dest, Vector3D src, void* parameters) {
   real r = Vector3D_distance(dest, src);
   return -((real*)parameters)[0] / (90.0 * pow(r, 10.0));
}

real GreensFunction_LJ6(Vector3D dest, Vector3D src, void* parameters) {
   real r = Vector3D_distance(dest, src);
   return ((real*)parameters)[0] / (12.0 * pow(r, 4.0));
}

real GreensFunction_LJ12_deriv(Vector3D dest, Vector3D src, void* parameters, Vector3D direction) {
   Vector3D rvector = Vector3D_allocate();
   Vector3D_sub(rvector, dest, src);
   real r = Vector3D_length(rvector);
   Vector3D_normalize(rvector);

   real g = ((real*)parameters)[0] / (9.0 * pow(r, 11.0)) * Vector3D_dot(rvector, direction);
   Vector3D_free(rvector);
   return g;
}

real GreensFunction_LJ6_deriv(Vector3D dest, Vector3D src, void* parameters, Vector3D direction) {
   Vector3D rvector = Vector3D_allocate();
   Vector3D_sub(rvector, dest, src);
   real r = Vector3D_length(rvector);
   Vector3D_normalize(rvector);

   real g = -((real*)parameters)[0] / (3.0 * pow(r, 5.0)) * Vector3D_dot(rvector, direction);
   //   printf("%lf  %lf  %lf  \n",r, g, (3.0 * pow(r, 5.0)));
   Vector3D_free(rvector);
   return g;
}

real GreensFunction_Lesyng_deriv(Vector3D dest, Vector3D src, void* parameters, Vector3D direction) {
   real n = ((real*)parameters)[0];
   Vector3D rvector = Vector3D_allocate();
   Vector3D_sub(rvector, dest, src);
   real r = Vector3D_length(rvector);

   Vector3D_normalize(rvector);
   real g = -1. * Vector3D_dot(rvector, direction)/((n-3) * pow(r, (n-1)));
   //   printf("%lf  %lf  %lf\n",r,g,((n-3) * pow(r, (n-1))));
   Vector3D_free(rvector);
   return g;

}

real GreensFunction_Ghosh(Vector3D dest, Vector3D src, void* parameters) {
   printf("GreensFunction Ghosh is not implemented yet, sorry!\n");
   exit(-1);
}

real GreensFunction_Ghosh_deriv(Vector3D dest, Vector3D src, void* parameters, Vector3D direction) {
   Vector3D rvector = Vector3D_allocate();
   Vector3D_sub(rvector, dest, src);
   real r = Vector3D_length(rvector);
   // treacherous: we are using explicitly the formula from Ghosh et
   // al (1998).  NOT what would be derived from r^-n volume calculation!!
   real g = ((real*)parameters)[0] * Vector3D_dot(rvector, direction)  / pow(r, 4.0);
   Vector3D_free(rvector);
   return g;
}

real GreensFunction_XGB(Vector3D dest, Vector3D src, void* parameters) {
  printf("GreensFunction_XGB is not implemented... it shouldn't ever be needed!\n");
  exit(-1);
}

real GreensFunction_XGB_deriv(Vector3D src, Vector3D dest, void* parameters, Vector3D direction) {
  Vector3D ultDest = Vector3D_allocate();
  Vector3D rvector_srcToQuad = Vector3D_allocate();
  Vector3D rvector_quadToDdest = Vector3D_allocate();

  ultDest->x = ((real*)parameters)[0];
  ultDest->y = ((real*)parameters)[1];
  ultDest->z = ((real*)parameters)[2];

  Vector3D_sub(rvector_srcToQuad, dest, src);
  real r_srcQuad = Vector3D_length(rvector_srcToQuad);
  Vector3D_sub(rvector_quadToDdest, ultDest, dest);
  real r_quadUlt = Vector3D_length(rvector_quadToDdest);

/*   printf("ultDest=(%lf, %lf, %lf)\n", ultDest->x, ultDest->y, ultDest->z); */
/*   printf("Dest=(%lf, %lf, %lf)\n", dest->x, dest->y, dest->z); */
/*   printf("r1 = %lf\nr2 = %lf\n",r_srcQuad, r_quadUlt); */
  
  real g = Vector3D_dot(rvector_srcToQuad, direction) / (pow(r_srcQuad, 3.0) * r_quadUlt);
  
  Vector3D_free(rvector_quadToDdest);
  Vector3D_free(rvector_srcToQuad);
  Vector3D_free(ultDest);
  return g;
}


real GreensFunction_Lesyng(Vector3D dest, Vector3D src, void* parameters) {
   real n = ((real*)parameters)[0];
   Vector3D rvector = Vector3D_allocate();
   Vector3D_sub(rvector, dest, src);
   real r = Vector3D_length(rvector);
   real g = 1./(n-3) * 1./(n-2) * pow(r, -(n-2));
   Vector3D_free(rvector);
   return g;
}


real GreensFunction_monomial(Vector3D dest, Vector3D src, void* parameters) {
   real* params = (real*)parameters;

   unsigned int x = (unsigned int)params[0];
   unsigned int y = (unsigned int)params[1];
   unsigned int z = (unsigned int)params[2];
   real centerx = params[3];
   real centery = params[4];
   real centerz = params[5];
   real gridspacing = params[6];

   real px = (src->x - centerx) / gridspacing;
   real py = (src->y - centery) / gridspacing;
   real pz = (src->z - centerz) / gridspacing;

   return PetscPowRealInt(px, x) * PetscPowRealInt(py, y) * PetscPowRealInt(pz, z);
}

real GreensFunction_monomial_deriv(Vector3D dest, Vector3D src, void* parameters, Vector3D direction) {
   real* params = (real*)parameters;

   unsigned int x = (unsigned int)params[0];
   unsigned int y = (unsigned int)params[1];
   unsigned int z = (unsigned int)params[2];
   real centerx = params[3];
   real centery = params[4];
   real centerz = params[5];
   real gridspacing = params[6];

   real px = (src->x - centerx) / gridspacing;
   real py = (src->y - centery) / gridspacing;
   real pz = (src->z - centerz) / gridspacing;

   real dx = 0.0, dy = 0.0, dz = 0.0;

   if (x > 0) dx = x * PetscPowRealInt(px, x-1);
   if (y > 0) dy = y * PetscPowRealInt(py, y-1);
   if (z > 0) dz = z * PetscPowRealInt(pz, z-1);

   /* Direction is negative only for curved panels!  Arg! */

   return (dx * PetscPowRealInt(py, y) * PetscPowRealInt(pz, z) * -direction->x +
           PetscPowRealInt(px, x) * dy * PetscPowRealInt(pz, z) * -direction->y +
           PetscPowRealInt(px, x) * PetscPowRealInt(py, y) * dz * -direction->z) / gridspacing;
}

real GreensFunction_inversepower(Vector3D dest, Vector3D src, void* parameters) {
   real r = Vector3D_distance(dest, src);
   real power = ((real*)parameters)[0];

   if (r == 0.0)
      return 0.0;
   else
      return 1.0 / (pow(r, power));
}

real GreensFunction_inversepower_deriv(Vector3D dest, Vector3D src, void* parameters, Vector3D direction) {
   Vector3D rvector = Vector3D_allocate();
   Vector3D_sub(rvector, dest, src);
   real r = Vector3D_length(rvector);
   Vector3D_normalize(rvector);
   real power = ((real*)parameters)[0];
   real g;

   if (r == 0.0)
      g = 0.0;
   else
      g = -power / (pow(r, power+1.0)) * Vector3D_dot(rvector, direction);

   Vector3D_free(rvector);
   return g;
}
