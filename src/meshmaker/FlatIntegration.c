#include "meshmaker.h"

#define FIVE3 1.666666666667
#define SEVEN3 2.3333333333333
#define ONE6 0.16666666666667
#define ONE3 0.3333333333333
#define FT3 4.666666666667
#define LIMITFOURTH 9.0
#define LIMITSECOND 36.0
#define EQUIV_TOL 1.0e-9

// return gradient in parameters

void FlatIntegration_oneoverr_grad(Vector3D point, FlatPanel panel, void* parameters, real* slp, real* dlp) {
   Vector3D grad_dir;
   unsigned int i;
   grad_dir = Vector3D_allocate();
   for (i = 0; i < 3; i++) {
      grad_dir->x = grad_dir->y = grad_dir->z = 0.0;
      if (i == 0) 
         grad_dir->x = 1.0;
      else if (i == 1)
         grad_dir->y = 1.0;
      else if (i == 2) 
         grad_dir->z = 1.0;
      FlatIntegration_oneoverr_deriv(point, panel, (void *)&grad_dir, slp, dlp);
      ((real *)parameters)[i] = *slp;
   }

   Vector3D_free(grad_dir);
}

void FlatIntegration_oneoverr_deriv_qual(FlatPanel point, FlatPanel panel, void* parameters, real* slp, real* dlp) {
   real myslp, mydlp;

   if (point == panel) {
      *dlp = 0.0;
      return;
   }

   FlatIntegration_oneoverr(panel->centroid, point, parameters, &myslp, &mydlp);

   *dlp = mydlp * panel->area;
}

void FlatIntegration_oneoverr_deriv_qual_point(Vector3D point, FlatPanel panel, void* parameters, real* slp, real* dlp) {
   real myslp, mydlp;

   if (Vector3D_equal(panel->centroid, point)) {
      *dlp = 0.0;
      return;
   }

   FlatIntegration_oneoverr(point, panel, parameters, &myslp, &mydlp);

   *dlp = mydlp;
}

void FlatIntegration_oneoverr_deriv(Vector3D point, FlatPanel panel, void* parameters, real* slp, real* dlp) {
   Vector3D pmc, nrm;
   real xn, yn, zn, xsq, ysq, zsq, rsq, diagsq, dtol;
   unsigned int OK = 1, i;
   real znabs, xmxv[3], ymyv[3], fe[3], r[3], xri[3], yri[3];
   real fs = 0.0, fd = 0.0;
   real fsx = 0.0, fsy = 0.0;
   real fdx = 0.0, fdy = 0.0, fdz = 0.0;

   pmc = Vector3D_allocate();

   Vector3D_sub(pmc, point, panel->centroid);

   xn = Vector3D_dot(panel->panelaxis[0], pmc);
   yn = Vector3D_dot(panel->panelaxis[1], pmc);
   zn = Vector3D_dot(panel->panelaxis[2], pmc);
   
   Vector3D_free(pmc);

   if (parameters) {
      nrm = Vector3D_allocate();

      nrm->x = Vector3D_dot(panel->panelaxis[0], (Vector3D)parameters);
      nrm->y = Vector3D_dot(panel->panelaxis[1], (Vector3D)parameters);
      nrm->z = Vector3D_dot(panel->panelaxis[2], (Vector3D)parameters);
   }

   xsq = xn * xn;
   ysq = yn * yn;
   zsq = zn * zn;
   rsq = xsq + ysq + zsq;
   dtol = EQUIV_TOL * panel->min_diag;
   diagsq = panel->max_diag * panel->max_diag;
   znabs = fabs(zn);

   /* If the evaluation point is far enough away from the panel,
      compute the influence approximately using moments */

   if (rsq > (LIMITFOURTH * diagsq)) {
      real* s = panel->moments;
      /* First, second moments. */
      real r2Inv = 1.0 / rsq;
      real rInv = sqrt(r2Inv);
      real r3Inv = r2Inv * rInv;
      real r5Inv = r3Inv * r2Inv;
      real zr2Inv = zn * r2Inv;
      real ss1 = s[1] * rInv;
      real ss3 = -(s[3] + s[10]) * r3Inv;
      real ss5 = (xsq * s[10] + (xn * yn * s[7]) + ysq * s[3]) * r5Inv;
      fs = ss1 + ONE3 * ss3 + ss5;
      real fdsum = ss1 + ss3 + 5.0 * ss5;
      fd = zr2Inv * fdsum;
      real rss3 = r2Inv*ss1;
      real ssx3 = -xn*rss3;
      real ssy3 = -yn*rss3;
      real ssx5 = (xn*(s[3]+3.0*s[10])+yn*s[7])*r5Inv;
      real ssy5 = (yn*(s[10]+3.0*s[3])+xn*s[7])*r5Inv;
      real rss7 = -5.0*r2Inv*ss5;
      real ssx7 = xn*rss7;
      real ssy7 = yn*rss7;
      fsx = ssx3 + ssx5 + ssx7;
      fsy = ssy3 + ssy5 + ssy7;
      fdx = zr2Inv*(3.0*ssx3+5.0*ssx5+7.0*ssx7);
      fdy = zr2Inv*(3.0*ssy3+5.0*ssy5+7.0*ssy7);
      fdz = r2Inv*fdsum - zr2Inv*zr2Inv*(3.0*ss1 + 5.0*ss3 + 35.0*ss5);

      if (rsq < (LIMITSECOND * diagsq)) {
         /* Third and fourth moments added for diagsq/r2 between 40 and 150. */
         real s914 = s[9] + s[14];
         real s813 = s[8] + s[13];
         real s411 = s[4] + s[11];
         real s512 = s[5] + s[12];
         real s1215 = s[12] + s[15];
         real r7Inv = r5Inv * r2Inv;
         real r9Inv = r7Inv * r2Inv;
         real ss5 = (-xn * s813 - yn * s411 + 0.1 * (s512 + s1215)) * r5Inv;
         real ss7 = (FIVE3 *((xn * xsq * s[13] + yn * ysq * s[4])
                             + 3.0 * xn * yn * (xn * s[11]  +  yn * s[8]))
                     - xsq * s1215 - ysq * s512 - xn * yn * s914) * r7Inv;
         real ss9 = (7.0 * (ONE6 * (xsq * xsq * s[15] + ysq * ysq * s[5])
                            + xsq * ysq * s[12])
                     + SEVEN3 * xn * yn * (xsq * s[14] + ysq * s[9])) * r9Inv;
         fs += ss5 + ss7 + ss9;
         fdsum = 5.0 * ss5 + 7.0 * ss7 + 9.0 * ss9;
         fd += zr2Inv * fdsum;
         real txy = 2*xn*yn;
         ssx5 = -s813*r5Inv;
         ssy5 = -s411*r5Inv;
         rss7 =  5.0*r2Inv*ss5;
         ssx7 = (5.0*(xn*xn*s[13] + txy*s[11] + yn*yn*s[8]) - s1215*(xn+xn) -
                 yn*s914)*r7Inv - xn*rss7;
         ssy7 = (5.0*(yn*yn*s[4] + xn*xn*s[11] + txy*s[8]) - s512*(yn+yn) -
                 xn*s914)*r7Inv - yn*rss7;
         real rss9 = 7.0*ss7*r2Inv;
         real ssx9 = (FT3*xn*xsq*s[15] + 14.0*xn*ysq*s[12] + 49.0*yn*(xsq*s[14] +
                                                                      ONE3*ysq*s[9]))*r9Inv - xn*rss9;
         real ssy9 = (FT3*yn*ysq*s[5] + 14.0*yn*xsq*s[12] + 49.0*xn*(ysq*s[9] +
                                                                     ONE3*xsq*s[14]))*r9Inv - yn*rss9;
         real rss11 = 9.0*ss9*r2Inv;
         real ssx11 = -xn*rss11;
         real ssy11 = -yn*rss11;

         fsx += ssx5+ssx7+ssx9+ssx11;
         fsy += ssy5+ssy7+ssy9+ssy11;
         fdx += zr2Inv*(5.0*ssx5 + 7.0*ssx7 + 9.0*ssx9 + 11.0*ssx11);
         fdy += zr2Inv*(5.0*ssy5 + 7.0*ssy7 + 9.0*ssy9 + 11.0*ssy11);
         fdz += r2Inv*fdsum - zr2Inv*zr2Inv*(35.0*ss5 + 63.0*ss7 + 99.0*ss9);
      }
   }
   else {
      /* Otherwise, compute the influence analytically */

      if (znabs < dtol) {
         zn = 0.5 * dtol;
         znabs = 0.5 * dtol;
      }
   
      for (i = 0; i < 3; i++) {
         real xc = xn - panel->panelvertex[i]->x;
         real yc = yn - panel->panelvertex[i]->y;
         real zc = zn - panel->panelvertex[i]->z;
         xmxv[i] = xc;
         ymyv[i] = yc;
         fe[i] = xc*xc + zc*zc;
         r[i] = sqrt(yc*yc + fe[i]);
         if (r[i] < (1.005 * znabs))
            OK = 0;
         xri[i] = xmxv[i]/r[i];
         yri[i] = ymyv[i]/r[i];
      }
   
      for (i = 0; i < 3; i++) {
         unsigned int next;
         real v, fln, arg, s1, c1, s2, c2, s12, c12, val;

         if (i == 2)
            next = 0;
         else
            next = i + 1;
   
         v = xmxv[i]*panel->contributionS[i] - ymyv[i]*panel->contributionC[i];
   
         arg = (r[i] + r[next] - panel->edgelength[i]) / (r[i] + r[next] + panel->edgelength[i]);
   
         fln = -log(arg);

         if (arg > 0.0)
            fs += v * fln;

         if (arg > 0.0) {
            real fac = (r[i] + r[next] - panel->edgelength[i])*(r[i] + r[next] + panel->edgelength[i]);
            fac = v*(panel->edgelength[i]+panel->edgelength[i])/fac;
            fsx += fln*panel->contributionS[i] - fac*(xri[i] + xri[next]);
            fsy -= fln*panel->contributionC[i] + fac*(yri[i] + yri[next]);
            fdz -= fac*( 1.0/r[i] + 1.0/r[next] );
         }
   
         if (OK) {
            s1 = v * r[i];
            c1 = znabs*(xmxv[i]*panel->contributionC[i] + ymyv[i]*panel->contributionS[i]);
            s2 = v * r[next];
            c2 = znabs*(xmxv[next]*panel->contributionC[i] + ymyv[next]*panel->contributionS[i]);
         } else {
            s1 = (fe[i]*panel->contributionS[i]) - (xmxv[i]*ymyv[i]*panel->contributionC[i]);
            c1 = znabs*r[i]*panel->contributionC[i];
            s2 = (fe[next]*panel->contributionS[i]) - (xmxv[next]*ymyv[next]*panel->contributionC[i]);
            c2 = znabs*r[next]*panel->contributionC[i];
         }
   
         s12 = (s1*c2) - (s2*c1);
         c12 = (c1*c2) + (s1*s2);
   
         val = atan2(s12, c12);
   
         fd += val;

         real u1 = xmxv[i]*panel->contributionC[i] + ymyv[i]*panel->contributionS[i];
         real u2 = xmxv[next]*panel->contributionC[i] + ymyv[next]*panel->contributionS[i];
         if (!OK) {
            real rr = r[i]*r[i];
            real fh1  = xmxv[i]*ymyv[i];
            real fh2  = xmxv[next]*ymyv[next];
            real fac  = c1/((c1*c1+s1*s1)*rr );
            fdx += (rr*v+fh1*u1)*fac;
            fdy -= fe[i]*u1*fac;
            rr   = r[next]*r[next];
            fac  = c2/((c2*c2+s2*s2)*rr);
            fdx -= (rr*v+fh2*u2)*fac;
            fdy += fe[next]*u2*fac;
         }
         else {
            real fac  = zn/(c1*c1+s1*s1);
            fdx += (u1*v*xri[i]+r[i]*ymyv[i])*fac;
            fdy += (u1*v*yri[i]-r[i]*xmxv[i])*fac;
            fac  = zn/(c2*c2+s2*s2);
            fdx -= (u2*v*xri[next]+r[next]*ymyv[next])*fac;
            fdy -= (u2*v*yri[next]-r[next]*xmxv[next])*fac;
         }
      }      

      if (fd < 0.0)
         fd += 2.0 * M_PI;
   
      if (zn < 0.0)
         fd = -fd;
   
      fs -= zn * fd;

      if (Vector3D_equal(panel->centroid, point))
         fd = -2.0 * M_PI;  // this is internally inconsistent

      fsx -= zn*fdx;
      fsy -= zn*fdy;
   }

   if (parameters) {
      real nDrvSrc, nDrvDip;

      if (rsq < (dtol * dtol))
         nDrvSrc = 0.0;
      else
         nDrvSrc = nrm->x*fsx + nrm->y*fsy - nrm->z*fd;

      nDrvDip = nrm->x*fdx + nrm->y*fdy + nrm->z*fdz;

      *slp = nDrvSrc;
      *dlp = nDrvDip;

      Vector3D_free(nrm);
   }
   else {
      *slp = fs;
      *dlp = fd;
   }
}

void FlatIntegration_oneoverr(Vector3D point, FlatPanel panel, void* parameters, real* slp, real* dlp) {
   Vector3D pmc;
   real xn, yn, zn, xsq, ysq, zsq, rsq, diagsq;
   unsigned int OK = 1, i;
   real znabs, xmxv[3], ymyv[3], fe[3], r[3];
   real fs = 0.0, fd = 0.0;

   pmc = Vector3D_allocate();

   Vector3D_sub(pmc, point, panel->centroid);

   xn = Vector3D_dot(panel->panelaxis[0], pmc);
   yn = Vector3D_dot(panel->panelaxis[1], pmc);
   zn = Vector3D_dot(panel->panelaxis[2], pmc);
   
   Vector3D_free(pmc);
   
   xsq = xn * xn;
   ysq = yn * yn;
   zsq = zn * zn;
   rsq = xsq + ysq + zsq;
   diagsq = panel->max_diag * panel->max_diag;

   /* If the evaluation point is far enough away from the panel,
      compute the influence approximately using moments */

   if (rsq > (LIMITFOURTH * diagsq)) {
      real fs = (real)0, fd = (real)0;
      real* s = panel->moments;
      /* First, second moments. */
      real r2Inv = 1.0 / rsq;
      real rInv = sqrt(r2Inv);
      real r3Inv = r2Inv * rInv;
      real r5Inv = r3Inv * r2Inv;
      real zr2Inv = zn * r2Inv;
      real ss1 = s[1] * rInv;
      real ss3 = -(s[3] + s[10]) * r3Inv;
      real ss5 = (xsq * s[10] + (xn * yn * s[7]) + ysq * s[3]) * r5Inv;
      real fdsum;
      fs = ss1 + ONE3 * ss3 + ss5;
      fdsum = ss1 + ss3 + 5.0 * ss5;
      fd = zr2Inv * fdsum;
      if (rsq < (LIMITSECOND * diagsq)) {
         /* Third and fourth moments added for diagsq/r2 between 40 and 150. */
         real s914 = s[9] + s[14];
         real s813 = s[8] + s[13];
         real s411 = s[4] + s[11];
         real s512 = s[5] + s[12];
         real s1215 = s[12] + s[15];
         real r7Inv = r5Inv * r2Inv;
         real r9Inv = r7Inv * r2Inv;
         real ss5 = (-xn * s813 - yn * s411 + 0.1 * (s512 + s1215)) * r5Inv;
         real ss7 = (FIVE3 *((xn * xsq * s[13] + yn * ysq * s[4])
                             + 3.0 * xn * yn * (xn * s[11]  +  yn * s[8]))
                     - xsq * s1215 - ysq * s512 - xn * yn * s914) * r7Inv;
         real ss9 = (7.0 * (ONE6 * (xsq * xsq * s[15] + ysq * ysq * s[5])
                            + xsq * ysq * s[12])
                     + SEVEN3 * xn * yn * (xsq * s[14] + ysq * s[9])) * r9Inv;
         fs += ss5 + ss7 + ss9;
         fdsum = 5.0 * ss5 + 7.0 * ss7 + 9.0 * ss9;
         fd += zr2Inv * fdsum;
      }

      *slp = fs;
      *dlp = fd;
      
      return;
   }

   /* Otherwise, compute the influence analytically */

   znabs = fabs(zn);
   
   for (i = 0; i < 3; i++) {
      real xc = xn - panel->panelvertex[i]->x;
      real yc = yn - panel->panelvertex[i]->y;
      real zc = zn - panel->panelvertex[i]->z;
      xmxv[i] = xc;
      ymyv[i] = yc;
      fe[i] = xc*xc + zc*zc;
      r[i] = sqrt(yc*yc + fe[i]);
      if (r[i] < (1.005 * znabs))
         OK = 0;
   }
   
   for (i = 0; i < 3; i++) {
      unsigned int next;
      real v, arg, s1, c1, s2, c2, s12, c12, val;

      if (i == 2)
         next = 0;
      else
         next = i + 1;
   
      v = xmxv[i]*panel->contributionS[i] - ymyv[i]*panel->contributionC[i];
   
      arg = (r[i] + r[next] - panel->edgelength[i]) / (r[i] + r[next] + panel->edgelength[i]);
   
      if (arg > 0.0)
         fs -= v * log(arg);
   
      if (OK) {
         s1 = v * r[i];
         c1 = znabs*(xmxv[i]*panel->contributionC[i] + ymyv[i]*panel->contributionS[i]);
         s2 = v * r[next];
         c2 = znabs*(xmxv[next]*panel->contributionC[i] + ymyv[next]*panel->contributionS[i]);
      } else {
         s1 = (fe[i]*panel->contributionS[i]) - (xmxv[i]*ymyv[i]*panel->contributionC[i]);
         c1 = znabs*r[i]*panel->contributionC[i];
         s2 = (fe[next]*panel->contributionS[i]) - (xmxv[next]*ymyv[next]*panel->contributionC[i]);
         c2 = znabs*r[next]*panel->contributionC[i];
      }
   
      s12 = (s1*c2) - (s2*c1);
      c12 = (c1*c2) + (s1*s2);
   
      val = atan2(s12, c12);
   
      fd += val;
   }
   
   if (fd < 0.0)
      fd += 2.0 * M_PI;
   
   if (zn < 0.0)
      fd = -fd;
   
   fs -= zn * fd;

   if (Vector3D_equal(panel->centroid, point))
      fd = 2.0 * M_PI;

   *slp = fs;
   *dlp = fd;
}

void calcp_rotatePoints(Vector3D v1, Vector3D v2, Vector3D v3, Vector3D normal) {
   Vector3D v1o, v2o, v3o;
   real theta, costheta, sintheta;
   Vector3D Np;
   real Zm[9],Ym[9];

   v1o = Vector3D_allocate(); v2o = Vector3D_allocate(); v3o = Vector3D_allocate();
   Np = Vector3D_allocate();
   Vector3D_copy(v1o, v1); Vector3D_copy(v2o, v2); Vector3D_copy(v3o, v3);
  
   theta = -atan2(normal->y, normal->x);

   costheta = cos(theta); sintheta = sin(theta);

   Zm[0] = costheta; Zm[1] = -sintheta; Zm[2] = 0.0;
   Zm[3] = sintheta; Zm[4] =  costheta; Zm[5] = 0.0;
   Zm[6] = 0.0;  Zm[7] = 0.0; Zm[8] = 1.0;
  
   Vector3D_transform(v1, Zm, v1o);
   Vector3D_transform(v2, Zm, v2o);
   Vector3D_transform(v3, Zm, v3o);

   Vector3D_cross(Np, v1, v2);
   theta = M_PI/2 - atan2(Np->z, Np->x);

   costheta = cos(theta); sintheta = sin(theta);
   Ym[0] = costheta; Ym[1] = 0.0; Ym[2] = -sintheta;
   Ym[3] = 0.0; Ym[4] = 1.0; Ym[5] = 0.0;
   Ym[6] = sintheta; Ym[7] = 0.0; Ym[8] = costheta;

   Vector3D_transform(v1o, Ym, v1);
   Vector3D_transform(v2o, Ym, v2);
   Vector3D_transform(v3o, Ym, v3);

   Vector3D_copy(v1, v1o); Vector3D_copy(v2, v2o); Vector3D_copy(v3,v3o);

   Vector3D_free(v1o);
   Vector3D_free(v2o);
   Vector3D_free(v3o);
   Vector3D_free(Np);
}

void calcp_LogIntDoLine(real z, real x, real y1, real y2, BEMKernelType kernel, BEMLayerType layer, QuadratureRule qr, void* parameters, real* integral) {
   real R4term1, R4term2, R4term, R10term1, R10term2, R10term, BCAterm1, BCAterm2, fcnval;
   real startY, endY;
   real lineInt;
   real phi1, phi2, dPhi;
   real theta, RTheta, signZ, absZ;
   unsigned int i;

   if (z > 0.0)
      absZ = z;
   else 
      absZ = -z;

   if (x < 0)
      x = -x;

   startY = y1;
  
   if (x < y1) {
      if (y2 > 10.0 * y1)
         endY = 10.0 * y1;
      else
         endY = y2;
   } else {
      if (x > y2)
         endY = y2;
      else
         endY = x;
   }

   if (z > 0.0)
      signZ = 1.0;
   else if (z < 0.0)
      signZ = -1.0;
   else
      signZ = 0.0;
  
   if ((kernel == HELMHOLTZ_KERNEL) && (((real*)parameters)[0] < 1e-6))
      kernel = POISSON_KERNEL;
  
   lineInt = 0.0;  /* initialize the total line integral */
  
   /* now loop over all Y values to do complete line FlatIntegration */

   while (startY < y2) {
      phi1 = atan2(startY, x);  /* find the two angles between which we do 
                                   our integrating */
      phi2 = atan2(endY, x);
      //    printf("phi1 = %f\nphi2 = %f\n",phi1, phi2);
      dPhi = phi2 - phi1;
    
      /* now integrate over this angle range, using our given quadrature rule */
      for (i = 0; i < qr->order; i++) {
         theta = phi1 + (phi2-phi1) * qr->x[i];
         /* these are our variables in the loop */
         RTheta = sqrt(PetscPowRealInt(x * (1.0/cos(theta)),2) + PetscPowRealInt(z,2));
         /* our current angle and dist from point to line at that angle */
         if (kernel == POISSON_KERNEL) {
            if (layer == SINGLE_LAYER_INT)
               lineInt = lineInt + qr->w[i] * dPhi * (RTheta - absZ);
            else if (layer == DOUBLE_LAYER_INT) {
               /*                          printf("try this: %f\n", z/RTheta - signZ); */
               lineInt = lineInt + qr->w[i] * dPhi * (z/RTheta - signZ);
            } else 
               printf("LogIntDoLine: argument 'layer' must be SINGLE_LAYER_INT or DOUBLE_LAYER_INT!\n");
         }
         else if (kernel == HELMHOLTZ_KERNEL) {
            real KAPPA = ((real*)parameters)[0];
            if (layer == SINGLE_LAYER_INT)
               lineInt = lineInt - qr->w[i] * dPhi * (exp(-KAPPA * RTheta) - exp(-KAPPA * absZ))/KAPPA;
            else if (layer == DOUBLE_LAYER_INT)
               lineInt = lineInt - qr->w[i] * dPhi * (signZ * exp(-KAPPA * absZ) - (z/RTheta) * exp(-KAPPA * RTheta));
            else
               printf( "LogIntDoLine: argument 'layer' must be SINGLE_LAYER_INT or DOUBLE_LAYER_INT!\n");
         }
         else if (kernel == LJ_KERNEL) {
            if (layer == SINGLE_LAYER_INT) {
               printf("LJ integral only valid for double layer integration!\n");
               exit(-5);
            }
            else if (layer == DOUBLE_LAYER_INT) {
               real R4coeff = ((real*)parameters)[0];
               real R10coeff = ((real*)parameters)[1];
             
               R4term1 = -1. / (2 * 2.0 * PetscPowRealInt(RTheta, 4));
               R4term2 = -1. / (2 * 2.0 * PetscPowRealInt(fabs(z), 4));
               R4term = R4coeff/3. * (R4term1 - R4term2);
               R10term1 = -1. / (2 * 5.0 * PetscPowRealInt(RTheta, 10));
               R10term2 = -1. / (2 * 5.0 * PetscPowRealInt(fabs(z), 10));
               R10term = R10coeff / 9. * (R10term1 - R10term2);
               fcnval = z * (R10term - R4term);
               lineInt = lineInt + qr->w[i] * dPhi * fcnval;
            }
            else
               printf("LogIntDoLine: argument 'layer must be SINGLE_LAYER_INT or DOUBLE_LAYER_INT!\n");
         }
         else if (kernel == LJ12_KERNEL) {
            if (layer == SINGLE_LAYER_INT) {
               printf("LJ12 integral only valid for double layer integration!\n");
               exit(-5);
            }
            else if (layer == DOUBLE_LAYER_INT) {
               real R10coeff = ((real*)parameters)[0];
             
               R10term1 = -1. / (2 * 5.0 * PetscPowRealInt(RTheta, 10));
               R10term2 = -1. / (2 * 5.0 * PetscPowRealInt(fabs(z),10));
               R10term = R10coeff / 9. * (R10term1 - R10term2);
               fcnval = z * R10term;
               lineInt = lineInt + qr->w[i] * dPhi * fcnval;
            }
            else
               printf("LogIntDoLine: argument 'layer must be SINGLE_LAYER_INT or DOUBLE_LAYER_INT!\n");
         }
         else if (kernel == LJ6_KERNEL) {
			  if (layer == SINGLE_LAYER_INT) {
				 printf("LJ6 integral only valid for double layer integration!\n");
				 exit(-5);
			  }
			  else if (layer == DOUBLE_LAYER_INT) {
				 real R4coeff = ((real*)parameters)[0];
				 
				 R4term1 = -1. / (2 * 2.0 * PetscPowRealInt(RTheta, 4));
				 R4term2 = -1. / (2 * 2.0 * PetscPowRealInt(fabs(z), 4));
				 R4term = R4coeff/3. * (R4term1 - R4term2);
				 fcnval = z * (-R4term);
				 lineInt = lineInt + qr->w[i] * dPhi * fcnval;
			  }
			  else
				 printf("LogIntDoLine: argument 'layer must be SINGLE_LAYER_INT or DOUBLE_LAYER_INT!\n");
			}
         else if (kernel == LESYNG_KERNEL) {
			  if (layer == SINGLE_LAYER_INT) {
				 printf("LESYNG integral only valid for double layer integration!\n");
				 exit(-5);
			  }
			  else if (layer == DOUBLE_LAYER_INT) {
				 real n = ((real*)parameters)[0];
				 BCAterm1 = -1.0 / ((n - 2) * pow(RTheta, n - 2));
				 BCAterm2 = -1.0 / ((n - 2) * pow(fabs(z), n - 2));
				 fcnval  =  z * 1.0 / (n - 3) * (BCAterm1 - BCAterm2);
				 lineInt = lineInt + qr->w[i] * dPhi * fcnval;
			  }
			  else
				 printf("LogIntDoLine: argument 'layer must be SINGLE_LAYER_INT or DOUBLE_LAYER_INT!\n");
			}
         else
            printf( "LogIntDoLine: argument 'kernel' must be POISSON_KERNEL or HELMHOLTZ_KERNEL!\n");
      }
      
      /* reset startY, endY for the next set of angles, or set startY to exit loop */

      startY = endY;
 
      if (endY == y2)
         endY = 10.0 * y2;
      else if (endY * 10.0 > y2)
         endY = y2;
      else
         endY = 10.0 * endY;   
   }
    
   *integral = lineInt;
}

void calcp_LogInt(Vector3D v1, Vector3D v2, Vector3D p, BEMKernelType kernel, BEMLayerType layer, QuadratureRule qr, void* parameters, real* integral) {
   int crap = 0;
   Vector3D v1o, v2o, po, n;
   Vector3D v1n, v2n, pn; // on the plane stuff
   Vector3D dv1v2 = Vector3D_allocate();
   Vector3D dv1p = Vector3D_allocate();
   Vector3D dv2p = Vector3D_allocate();
   Vector3D intPoint = Vector3D_allocate();
   Vector3D newXaxis = Vector3D_allocate();
   Vector3D newYaxis = Vector3D_allocate();
   real d, dz;
   real startY, endY;
   real x1, x2, y1, y2, testNorm;
   real edgeLength, c, s, z;
   real nxy, alpha;
   real dLineInt = 0.0;
   Vector3D translation; translation = Vector3D_allocate();
   v1o = Vector3D_allocate(); v2o = Vector3D_allocate();
   po = Vector3D_allocate(); n = Vector3D_allocate();


   Vector3D_copy(v1o, v1);
   Vector3D_copy(v2o, v2);
   Vector3D_copy(po, p);

   /* set the line up so that the normal to the triangle (v1, v2, p) is
      parallel to the z axis, ie rotate the coord system so that the triangle
      is parallel to the xy plane */

   Vector3D_cross(n, v1, v2);
   nxy = sqrt(PetscPowRealInt(n->x,2) + PetscPowRealInt(n->y,2));

#ifdef REAL_IS_DOUBLE
   if (nxy > 10 * DBL_EPSILON)  // machine precision tolerance for normal to not be straight z;
      calcp_rotatePoints(v1, v2, p, n);
#else
#ifdef REAL_IS_FLOAT
   if (nxy > 10 * FLT_EPSILON) { // machine precision tolerance for normal to not be straight z;
      calcp_rotatePoints(v1, v2, p, n);
   }
#endif
#endif

   /* translate so that the line is at z = 0, and the observation point is at (0,0,pz) */

   translation->x = -p->x;
   translation->y = -p->y;
   translation->z = -(v1->z);
   Vector3D_add(v1, v1, translation);
   Vector3D_add(v2, v2, translation);
   Vector3D_add( p,  p, translation);
   Vector3D_free(translation);
   /* right up to here (11:00 AM 7/30) */

   dz = v1->z - v2->z;

   if (dz < 0)
      dz = -dz;

#ifdef REAL_IS_DOUBLE
   if (dz > 1e-6) {
      printf("LogInt: transforming line onto x-y plane failed! (dz = %f)\n", dz);
      exit(-4);
   }
#else
#ifdef REAL_IS_FLOAT
   if (dz > 1e-6) {
      printf("LogInt: transforming line onto x-y plane failed! (dz = %f)\n", dz);
      exit(-4);
   }
#endif
#endif

   v1n = Vector3D_allocate();
   v2n = Vector3D_allocate();
   pn = Vector3D_allocate();
   v1n->x = v1->x; v1n->y = v1->y;
   v2n->x = v2->x; v2n->y = v2->y;
   pn->x = p->x; pn->y = p->y;
   Vector3D_scale(pn, -1.);
   Vector3D_add(dv1p, v1n, pn);
   Vector3D_add(dv2p, v2n, pn);
   Vector3D_scale(v2n, -1.);

   Vector3D_add(dv1v2, v1n, v2n);
   Vector3D_normalize(dv1v2);
   Vector3D_scale(v2n, -1.);
   c = dv1v2->x;  // see notes 11/19/04--this is a more stable rotation
   s = dv1v2->y;
   Vector3D_copy(intPoint, v1n);  // little hardcoded rotations
   v1n->x = s * intPoint->x - c * intPoint->y;
   v1n->y = c * intPoint->x + s * intPoint->y;
   Vector3D_copy(intPoint, v2n);
   v2n->x = s * intPoint->x - c * intPoint->y;
   v2n->y = c * intPoint->x + s * intPoint->y;
   d = v1n->x;

   // hardcoded (relative) constants
#ifdef REAL_IS_DOUBLE
   if (fabs(d) < 10 * DBL_EPSILON) {
#else
#ifdef REAL_IS_FLOAT
   if (fabsf(d) < 10 * FLT_EPSILON) {
#endif
#endif
      *integral = 0.0;
   } else {
      z  = p->z;
      x1 = v1n->x;
      y1 = v2n->y;
      y2 = v1n->y;
      if ( ( y2 < 0  ) && (y1 < 0)) {
         y1 = - y1;
         y2 = - y2;
      }
      
      if (y2 < y1) {
         testNorm = y1;
         y1 = y2;
         y2 = testNorm;
      }
      
      if (y1 >= 0)
         startY = y1;
      else
         startY = 0;
      
      if (x1 <= y2)
         endY = x1;
      else
         endY = y2;
      
      calcp_LogIntDoLine(z, x1, startY, y2, kernel, layer, qr, parameters, integral);
      
      if (y1 < 0.0) {
         y1 = -y1;
         calcp_LogIntDoLine(z, x1, 0.0, y1, kernel, layer, qr, parameters, &dLineInt);
         *integral = *integral + dLineInt;
      }
   }
   Vector3D_copy(p, po);
   Vector3D_copy(v1, v1o);
   Vector3D_copy(v2, v2o);
   
   /* make sure we clean up acceptably here */
   
   Vector3D_free(newXaxis);
   Vector3D_free(newYaxis);
   Vector3D_free(intPoint);
   Vector3D_free(dv1v2);
   Vector3D_free(dv1p);
   Vector3D_free(dv2p);
   Vector3D_free(v1n);
   Vector3D_free(v2n);
   Vector3D_free(pn);
   Vector3D_free(n);
   Vector3D_free(po);
   Vector3D_free(v1o);
   Vector3D_free(v2o);
#ifdef I_HATE_IFDEFS_BECAUSE_THEY_MESS_UP_EMACS_TABS
   }
#endif
}  

void calcp_directQuadrature(Vector3D point, FlatPanel panel, BEMKernelType kernel, BEMLayerType layer, QuadratureRule qr, void* parameters, real* integral) {
   unsigned int i, j;
   real minY, maxY, dY, Xfinal;
   real Xcurrent, Ycurrent, currentHeight, r;
   real curWeight, curFcnValue = 0.0;
   Vector3D v1;
   Vector3D v2poop;
   Vector3D v3;
   Vector3D pnt;
   Vector3D Dv2;
   Vector3D curPoint;
   real A[9];
   v1 = Vector3D_allocate();   
   v2poop = Vector3D_allocate();
   v3 = Vector3D_allocate();   
   pnt = Vector3D_allocate();
   Dv2 = Vector3D_allocate();
   curPoint = Vector3D_allocate();
   Vector3D_copy(Dv2, panel->panelvertex[1]);
   Vector3D_scale(Dv2,-1.0);
   Vector3D_add(panel->panelvertex[0], panel->panelvertex[0], Dv2);
   Vector3D_add(panel->panelvertex[1], panel->panelvertex[1], Dv2);
   Vector3D_add(panel->panelvertex[2], panel->panelvertex[2], Dv2);
   Vector3D_add(point, point, Dv2);
   Vector3D_scale(Dv2,-1.0);   
   for (i = 0; i < 9; i++)   
      A[i] = 0;
   A[1] = -1; A[3] = 1; A[8] = 1; 
   
   Vector3D_transform(v1, A, panel->panelvertex[0]);
   Vector3D_transform(v2poop, A, panel->panelvertex[1]);
   Vector3D_transform(v3, A, panel->panelvertex[2]);
   Vector3D_transform(pnt, A, point);
   Vector3D_add(panel->panelvertex[0], panel->panelvertex[0], Dv2);
   Vector3D_add(panel->panelvertex[1], panel->panelvertex[1], Dv2);
   Vector3D_add(panel->panelvertex[2], panel->panelvertex[2], Dv2);
   Vector3D_add(point, point, Dv2);
   
   *integral = 0.0;
   minY = v1->y; if (v2poop->y < minY) minY = v2poop->y; if (v3->y < minY) minY = v3->y;
   maxY = v1->y; if (v2poop->y > maxY) maxY = v2poop->y; if (v3->y > maxY) maxY = v3->y;
   dY = maxY - minY;
   Xfinal = v1->x;
   for (i = 0; i < qr->order; i++) {
      Xcurrent = Xfinal * qr->x[i];
      for (j = 0; j < qr->order; j++) {
         currentHeight = dY * qr->x[i];
         Ycurrent = minY * qr->x[i] + currentHeight * qr->x[j];
         curPoint->x = Xcurrent;
         curPoint->z = 0;
         curPoint->y = Ycurrent;
         /*Vector3D_sub(curPoint, pnt, curPoint);*/
         /*Vector3D_addscaled(-1.0, curPoint, 1.0, pnt);*/
         r = Vector3D_distance(curPoint, pnt);
         if (kernel == LJ_KERNEL)
            curFcnValue = (pnt->z / r) * (-((real*)parameters)[0]/(3 * PetscPowRealInt(r,5)) + ((real*)parameters)[1]/(9 * PetscPowRealInt(r,11)));
         else if (kernel == LJ12_KERNEL)
            curFcnValue = (pnt->z / r) * ((real*)parameters)[0]/(9 * PetscPowRealInt(r,11));
         else if (kernel == LJ6_KERNEL)
            curFcnValue = (pnt->z / r) * (-((real*)parameters)[0]/(3 * PetscPowRealInt(r,5)));
         else if (kernel == POISSON_KERNEL) {
            if (layer == SINGLE_LAYER_INT) 
               curFcnValue = 1.0 / r;
            else if (layer == DOUBLE_LAYER_INT)
               curFcnValue = (pnt->z / r) * (1.0 / r);
         } 
         else if (kernel == HELMHOLTZ_KERNEL) {
            if (layer == SINGLE_LAYER_INT) 
               curFcnValue = exp(-((real*)parameters)[0] * r) / r;
            else if (layer == DOUBLE_LAYER_INT)
               curFcnValue = (pnt->z / r) * (exp(-((real*)parameters)[0] * r) / r);
         } 
         curWeight = Xfinal * qr->w[i] * currentHeight * qr->w[j];
         *integral += curWeight * curFcnValue;
      }
   }
   
   Vector3D_free(v1);  Vector3D_free(v2poop);  Vector3D_free(v3); Vector3D_free(pnt);
   Vector3D_free(Dv2); Vector3D_free(curPoint);
}

void calcp(Vector3D point, FlatPanel panel, BEMKernelType kernel, BEMLayerType layer, QuadratureRule qr, void* parameters, real *integral) {
   Vector3D pnt, pnttmp;
   Vector3D thispoint, nextpoint;
	unsigned int j, next, usedirect;
	unsigned int countLeft = 0;
   real normpoint, signD, dFlatPanelInt, d;
   real ztol = 1e-5;

   normpoint = Vector3D_length(point);
   thispoint = Vector3D_allocate();
   nextpoint = Vector3D_allocate();
	
   pnttmp = Vector3D_allocate();
   pnt = Vector3D_allocate();
   Vector3D_sub(pnttmp, point, panel->centroid);
   pnt->x = Vector3D_dot(panel->panelaxisnum[0], pnttmp);
   pnt->y = Vector3D_dot(panel->panelaxisnum[1], pnttmp);
   pnt->z = Vector3D_dot(panel->panelaxisnum[2], pnttmp);
   Vector3D_free(pnttmp);

   *integral = 0.0;

  //if (pnt->z < 0.2) { // this is probably crap
  if (0) { // this is probably crap
      usedirect = 1;
      calcp_directQuadrature(pnt, panel, kernel, layer, qr, parameters, integral);
      }
      else {
         for (j = 0; j < 3; j++) {
            /* figure out what line to integrate: goes between v_j and v_next */
            next = j + 1;

            if (j == 2)
               next = 0;

            d = Vector3D_dot(pnt, panel->edges[j]) - panel->edgeRHS[j];
            if ((d > 0.0 &&  panel->edgeOV[j] > 0.0) ||
                (d < 0.0 &&  panel->edgeOV[j] < 0.0)) {
               countLeft = countLeft + 1;
               signD = 1.0;
            }
            else {
               signD = -1.0;
            }
         
         
            Vector3D_copy(thispoint, panel->panelvertexnum[j]);
            Vector3D_copy(nextpoint, panel->panelvertexnum[next]);

            calcp_LogInt(thispoint, nextpoint, pnt, kernel, layer, qr, parameters, &dFlatPanelInt);
            *integral += signD * dFlatPanelInt;
         }
      }

      /* 2 little checks: 1 is double layer off the panel == 0*/
 
      if ((layer == DOUBLE_LAYER_INT) && (fabs(pnt->z) < 10 * FLT_EPSILON)) {
         if (countLeft == 3) {
            *integral = 2.0 * M_PI;
         } else {
            *integral = 0.0;
         }
      }
      
      Vector3D_free(pnt);
      Vector3D_free(thispoint);
      Vector3D_free(nextpoint);  
   }

   void FlatIntegration_LJ(Vector3D point, FlatPanel panel, void* parameters, real *slp, real *dlp)
      {
         static QuadratureRule qr = NULL;

#ifdef OMP
#pragma omp critical
#endif
         if (qr == NULL)
            qr = QuadratureRule_allocate(4);

         calcp(point, panel, LJ_KERNEL, DOUBLE_LAYER_INT, qr, parameters, dlp);
         *slp = *dlp;
      }

   void FlatIntegration_LJ12(Vector3D point, FlatPanel panel, void* parameters, real *slp, real *dlp) {
      static QuadratureRule qr = NULL;

#ifdef OMP
#pragma omp critical
#endif
      if (qr == NULL)
         qr = QuadratureRule_allocate(4);

      calcp(point, panel, LJ12_KERNEL, DOUBLE_LAYER_INT, qr, parameters, dlp);
      *slp = *dlp;
   }

   void FlatIntegration_LJ6(Vector3D point, FlatPanel panel, void* parameters, real *slp, real *dlp) {
      static QuadratureRule qr = NULL;

#ifdef OMP
#pragma omp critical
#endif
      if (qr == NULL)
         qr = QuadratureRule_allocate(4);

      calcp(point, panel, LJ6_KERNEL, DOUBLE_LAYER_INT, qr, parameters, dlp);
      *slp = *dlp;
   }

   void FlatIntegration_oneoverr_numerical(Vector3D point, FlatPanel panel, void* parameters, real *slp, real *dlp) {
      real slpa, dlpa;
      static QuadratureRule qr = NULL;

#ifdef OMP
#pragma omp critical
#endif
      if (qr == NULL)
         qr = QuadratureRule_allocate(4);

      calcp(point, panel, POISSON_KERNEL, SINGLE_LAYER_INT, qr, parameters, slp);
      calcp(point, panel, POISSON_KERNEL, DOUBLE_LAYER_INT, qr, parameters, dlp);

      /*   FlatIntegration_oneoverr(point, panel, parameters, &slpa, &dlpa); */
      /*   if ((fabs(*slp-slpa) > 1e-3 * fabs(slpa)) || */
      /*       (fabs(*dlp-dlpa) > 1e-3 * fabs(dlpa)) ) { */
      /*     printf("ERROR:\n"); */
      /*     printf("SLP: %f ASLP: %f DLP: %F ADLP: %f\n", *slp, slpa, *dlp, dlpa); */
      /*     printf("v1->x = %f;  v1->y = %f;  v1->z = %f;\n",panel->vertex[0]->x, */
      /*            panel->vertex[0]->y,panel->vertex[0]->z); */
      /*     printf("v2->x = %f;  v2->y = %f;  v2->z = %f;\n",panel->vertex[1]->x, */
      /*            panel->vertex[1]->y,panel->vertex[1]->z); */
      /*     printf("v3->x = %f;  v3->y = %f;  v3->z = %f;\n",panel->vertex[2]->x, */
      /*            panel->vertex[2]->y,panel->vertex[2]->z); */
      /*     printf("point->x = %f;  point->y = %f;  point->z = %f;\n", point->x, */
      /*            point->y, point->z); */
      /*     //exit(0); */
      /*  } */
   }

void FlatIntegration_ekroverr_deriv_numerical(Vector3D point, FlatPanel panel, Vector3D direction, void* parameters, real *slp, real *dlp) {
  FlatIntegration_oneoverr_deriv(point, panel, (void *)direction, slp, dlp);
  if (((real *)parameters)[0] < 1e-6) {
	 return;
  }
  real extraparam[4];
  extraparam[0] = ((real *)parameters)[0];
  extraparam[1] = direction->x;
  extraparam[2] = direction->y;
  extraparam[3] = direction->z;
  real junk= FlatIntegration_maltquad(point, panel, DESINGULARIZED_HELMHOLTZ_KERNEL, (void *)extraparam, NORMDERIV_SINGLE_LAYER_INT);
  *slp-=junk;
  //  *dlp += FlatIntegration_maltquad(point, panel, DESINGULARIZED_HELMHOLTZ_KERNEL, (void *)extraparam, NORMDERIV_DOUBLE_LAYER_INT);
}

void FlatIntegration_ekroverr_numerical(Vector3D point, FlatPanel panel, void* parameters, real *slp, real *dlp) {
      static QuadratureRule qr = NULL;
		if (((real *)parameters)[0] < 1e-6) {
		  FlatIntegration_oneoverr(point, panel, parameters, slp, dlp);
		  return;
		}
#ifdef OMP
#pragma omp critical
#endif
      if (qr == NULL)
         qr = QuadratureRule_allocate(4);

      calcp(point, panel, HELMHOLTZ_KERNEL, SINGLE_LAYER_INT, qr, parameters, slp);
      calcp(point, panel, HELMHOLTZ_KERNEL, DOUBLE_LAYER_INT, qr, parameters, dlp);
}

   void FlatIntegration_ekroverr_desingularized(Vector3D point, FlatPanel panel, void* parameters, real *slp, real *dlp) {
      FlatIntegration_oneoverr(point, panel, parameters, slp, dlp);
      *slp += FlatIntegration_maltquad(point, panel, DESINGULARIZED_HELMHOLTZ_KERNEL, parameters, SINGLE_LAYER_INT);
      *dlp += FlatIntegration_maltquad(point, panel, DESINGULARIZED_HELMHOLTZ_KERNEL, parameters, DOUBLE_LAYER_INT);
   }

   void FlatIntegration_general(Vector3D point, FlatPanel panel, BEMKernelType kernel, BEMLayerType layer, void* parameters, real *integral) {
      if (kernel == POISSON_KERNEL) {
         real slp, dlp;
         FlatIntegration_oneoverr(point, panel, parameters, &slp, &dlp);
         if (layer == SINGLE_LAYER_INT)
            *integral = slp;
         else if (layer == DOUBLE_LAYER_INT)
            *integral = dlp;
      }
      else if (kernel == HELMHOLTZ_KERNEL) {
         real slp, dlp;
         FlatIntegration_ekroverr_numerical(point, panel, parameters, &slp, &dlp);
         if (layer == SINGLE_LAYER_INT)
            *integral = slp;
         else if (layer == DOUBLE_LAYER_INT)
            *integral = dlp;
      }
   }

   void gen_Stroud_Rule(real *xtab, real *ytab, real *wtab) {
      real xtab1[4] = {
         -0.861136311594052575223946488893,
         -0.339981043584856264802665759103,
         0.339981043584856264802665759103,
         0.861136311594052575223946488893
      };
      real wtab1[4] = {
         0.347854845137453857373063949222E+00,
         0.652145154862546142626936050778E+00,
         0.652145154862546142626936050778E+00,
         0.347854845137453857373063949222E+00,
      };
      real xtab2[4] = {
         0.0571041961E+00,
         0.2768430136E+00,
         0.5835904324E+00,
         0.8602401357E+00
      };
      real wtab2[4] = {
         0.1355069134E+00,
         0.2034645680E+00,
         0.1298475476E+00,
         0.0311809709E+00
      };

      unsigned int k = 0,i,j;
      unsigned int norder2 = 4;
      for (i = 0; i < norder2; i++)
         xtab1[i] = 0.5 * (xtab1[i] + 1.0);
  
      for (i = 0; i < norder2; i++) {
         for (j = 0; j < norder2; j++) {
            xtab[k] = xtab2[j];
            ytab[k] = xtab1[i] * (1.00 - xtab2[j]);
            wtab[k] = 0.5 * wtab1[i] * wtab2[j];
            k = k + 1;
         }
      }
      /*   for (i=0; i < k; i++) */
      /*     printf("xtab[%d] = %f\n", i, xtab[i]); */
  
      /*   for (i=0; i < k; i++) */
      /*     printf("ytab[%d] = %f\n", i, ytab[i]); */

      /*   for (i=0; i < k; i++) */
      /*     printf("wtab[%d] = %f\n", i, wtab[i]); */
   }



 
void FlatIntegration_Ghosh(Vector3D point, FlatPanel panel, void* parameters, real* slp, real* dlp) {
  static QuadratureRule qr = NULL;

#ifdef OMP
#pragma omp critical
#endif
  if (qr == NULL)
     qr = QuadratureRule_allocate(4);

  calcp(point, panel, GHOSH_KERNEL, DOUBLE_LAYER_INT, qr, parameters, dlp);
  *slp = *dlp;
}

// the Grycuk 1/r^6 integral is accomplished using the LJ6!!
void FlatIntegration_Grycuk(Vector3D point, FlatPanel panel, void* parameters, real* slp, real* dlp) {
  static QuadratureRule qr = NULL;

#ifdef OMP
#pragma omp critical
#endif
  if (qr == NULL)
     qr = QuadratureRule_allocate(4);
  
  calcp(point, panel, LJ6_KERNEL, DOUBLE_LAYER_INT, qr, parameters, dlp);
  *slp = *dlp;
}

void FlatIntegration_Lesyng(Vector3D point, FlatPanel panel, void* parameters, real* slp, real* dlp) {
  static QuadratureRule qr = NULL;

#ifdef OMP
#pragma omp critical
#endif
  if (qr == NULL)
     qr = QuadratureRule_allocate(4);

  calcp(point, panel, LESYNG_KERNEL, DOUBLE_LAYER_INT, qr, parameters, dlp);
  *slp = *dlp;
}

real FlatIntegration_maltquad(Vector3D point, FlatPanel panel, BEMKernelType kerneltype, void* parameters, BEMLayerType layertype) {
   static QuadratureRule qr = NULL;

#ifdef OMP
#pragma omp critical
#endif
   if (qr == NULL)
      qr = QuadratureRule_allocate(4);

   unsigned int center, other1, other2;
   unsigned int c, o1, o2;

   center = 0; other1 = 0; other2 = 0;
	Vector3D axis1 = Vector3D_allocate();
   real axis1length;
   Vector3D axis1norm = Vector3D_allocate();
   Vector3D otherside = Vector3D_allocate();
   real axis1dot = 0.0; 

   for (c = 0; c < 3; c++) {
      for (o1 = 0; o1 < 3; o1++) {
         if (o1 == c) continue;

         o2 = 3 - o1 - c;

         Vector3D_sub(axis1, panel->vertex[o1], panel->vertex[c]);
         axis1length = Vector3D_length(axis1);
         Vector3D_copy(axis1norm, axis1);
         Vector3D_scale(axis1norm, 1.0 / axis1length);
         Vector3D_sub(otherside, panel->vertex[o2], panel->vertex[c]);
         axis1dot = Vector3D_dot(otherside, axis1) / Vector3D_dot(axis1, axis1);
         if ((axis1dot > 0.0) && (axis1dot < 1.0)) {
            center = c;
            other1 = o1;
            other2 = o2;
            break;
         }
      }

      if ((axis1dot > 0.0) && (axis1dot < 1.0)) break;
   }

   Vector3D axis1proj = Vector3D_allocate();
   Vector3D_addscaled(axis1proj, panel->vertex[center], axis1dot, axis1);
   Vector3D axis2 = Vector3D_allocate();
   Vector3D_sub(axis2, panel->vertex[other2], axis1proj);
   real height = Vector3D_length(axis2);
   real length1 = Vector3D_distance(axis1proj, panel->vertex[center]);
   real length2 = Vector3D_distance(panel->vertex[other1], axis1proj);
   real slope1 = -length1 / height;
   real slope2 = length2 / height;

   unsigned int i, j;
   Vector3D currentx = Vector3D_allocate();
   Vector3D currentpoint = Vector3D_allocate();

   real integral = 0.0;
   real weightsum = 0.0;

   for (i = 0; i < qr->order; i++) {
      Vector3D_addscaled(currentx, axis1proj, qr->x[i], axis2);
      real starty = slope2 * qr->x[i] * height - length2;
      real endy = slope1 * qr->x[i] * height + length1;
      real dy = endy - starty;
      for (j = 0; j < qr->order; j++) {
         Vector3D_addscaled(currentpoint, currentx, qr->x[j] * dy - endy, axis1norm);
         real currentweight = height * qr->w[i] * dy * qr->w[j];
         real value = 0.0;

         if (layertype == SINGLE_LAYER_INT)
            value = GreensFunction(point, currentpoint, kerneltype, parameters);
         else if (layertype == DOUBLE_LAYER_INT)
            value = GreensFunction_deriv(point, currentpoint, kerneltype, parameters, panel->panelaxis[2]);
			else if (layertype == NORMDERIV_SINGLE_LAYER_INT) {
			  Vector3D direction = Vector3D_allocate();
			  direction->x = ((real *)parameters)[1];
			  direction->y = ((real *)parameters)[2];
			  direction->z = ((real *)parameters)[3];
			  value = -GreensFunction_deriv(point, currentpoint, kerneltype, parameters, direction);
			  Vector3D_free(direction);
			} else if (layertype == NORMDERIV_DOUBLE_LAYER_INT) {
			  Vector3D direction = Vector3D_allocate();
			  Vector3D direction2 = Vector3D_allocate();
			  direction->x = ((real *)parameters)[1];
			  direction->y = ((real *)parameters)[2];
			  direction->z = ((real *)parameters)[3];
			  direction2->x = ((real *)parameters)[4];
			  direction2->y = ((real *)parameters)[5];
			  direction2->z = ((real *)parameters)[6];
			  //			  printf("normderiv double malt quad is not implemented yet!\n");
			  value = 4.0 * M_PI * GreensFunction_doublederiv(point, currentpoint, kerneltype, parameters, direction, direction2);
			  Vector3D_free(direction);
			  Vector3D_free(direction2);
			}
         integral += value * currentweight;
         weightsum += currentweight;
      }
   }

   if (fabs(weightsum - panel->area) > 1e-6)
      printf("MALTQUAD FUTZ\n");

   Vector3D_free(currentpoint);
   Vector3D_free(currentx);
   Vector3D_free(axis2);
   Vector3D_free(axis1proj);
   Vector3D_free(otherside);
   Vector3D_free(axis1norm);
   Vector3D_free(axis1);

   return integral;
}

