#include "meshmaker.h"

void Integration_oneoverr_TOR(Vector3D point, TOR tor, void* parameters, real* slp, real* dlp) {
   *slp = 0.0;
   *dlp = 0.0;
   Integration_oneoverr_TOR_single(point, tor, parameters, slp);
   Integration_oneoverr_TOR_double(point, tor, parameters, dlp);
}

void Integration_oneoverr_TOR_single(Vector3D point, TOR tor, void* parameters, real* slp) {
   *slp = 0.0;
   Integration_general_TOR_single(point, tor, POISSON_KERNEL, parameters, slp);
}

void Integration_general_TOR_single(Vector3D point, TOR tor,
                                    BEMKernelType kernel, void* parameters, real* slp) {
   real self_term_distance_tolerance = 1e-5;
   real near_term_distance_tolerance = 3.0;  // equivalent to threshRadius in matlab
   Vector3D localpnt = Vector3D_allocate();
   *slp = 0.0;

   // NOTE: this is not the same order of transformation as done in matlab version!!  there,
   // everything is done in A' * (x - y) to get into local coords.  here we do A * x - y
   Vector3D_transformVecs(localpnt, tor->Tmat[0], tor->Tmat[1], tor->Tmat[2], point);
   Vector3D_add(localpnt, localpnt, tor->Tvec);

   if (Vector3D_distance(point, tor->centroid) < self_term_distance_tolerance)
      Integration_general_TOR_single_self(localpnt, tor, kernel, parameters, slp);
   //else if (Vector3D_distance(point, tor->centroid) < near_term_distance_tolerance * tor->maxEdgelength)
   else if (Vector3D_distance(point, tor->centroid) < near_term_distance_tolerance * sqrt(tor->area))
      Integration_general_TOR_single_near(localpnt, tor, kernel, parameters, slp);
   else
      Integration_general_TOR_single_local(localpnt, tor, kernel, parameters, slp);

   Vector3D_free(localpnt);
}

void Integration_general_TOR_single_local(Vector3D point, TOR tor,
                                          BEMKernelType kernel, void* parameters, real* slp) {
  
   real startTheta, endTheta, startPsi, endPsi, c, a, dTheta, dPsi;
   real curTheta, curPsi;
   real GnFcnVal, detJ, curWeight;
   real near_tolerance = 3.0;
   Vector3D srcpoint = Vector3D_allocate();
   unsigned int i, j;

   startTheta = tor->startTheta;
   endTheta = tor->endTheta;
   startPsi =tor->startPsi;
   endPsi = tor->endPsi;
   c = tor->c;
   a = tor->a;
   dTheta = endTheta - startTheta;
   dPsi = endPsi - startPsi;
   static QuadratureRule qr = NULL;
   //if (Vector3D_distance(point, tor->centroid) < near_tolerance * tor->maxEdgelength) {
   if (Vector3D_distance(point, tor->centroid) < near_tolerance * sqrt(tor->area)) {
      Integration_general_TOR_single_near(point, tor, kernel, parameters, slp);
      Vector3D_free(srcpoint);
      return;
   }
#ifdef OMP
#pragma omp critical
#endif
   if (qr == NULL) {
      qr = QuadratureRule_allocate(4);
   }

   for (i = 0; i < qr->order; i++) {
      curTheta = startTheta + qr->x[i] * dTheta;
      for (j = 0; j < qr->order; j++) {
         curPsi = startPsi + qr->x[j] * dPsi;
         torToCart(srcpoint, NULL, c, a, curTheta, curPsi);
         GnFcnVal = GreensFunction(point, srcpoint, kernel, parameters);
         detJ = a * (c + a * cos(curPsi));
         curWeight = fabs(dTheta) * qr->w[i] * fabs(dPsi) * qr->w[j];
         *slp += GnFcnVal * detJ * curWeight;
      }
   }

   Vector3D_free(srcpoint);
}

void Integration_general_TOR_single_self(Vector3D point, TOR tor,
                                         BEMKernelType kernel, void* parameters, real* slp) {

   real thresholdForSmallPanelApprox = 1e-5;
   real maxRad = tor->maxEdgelength;

   real startTheta, startPsi, endTheta, endPsi, c, a;
   real midTheta, midPsi, one4thTheta, one4thPsi, three4thTheta, three4thPsi;
   real otherterms, curInt;
   Vector3D origin, Xaxis, Yaxis, Zaxis;
   TOR tmpTOR;

   if (maxRad < thresholdForSmallPanelApprox) {
      return;
   }

   origin = Vector3D_allocate();
   Xaxis = Vector3D_allocate();
   Yaxis = Vector3D_allocate();
   Zaxis = Vector3D_allocate();
   Xaxis->x = 1.0;
   Yaxis->y = 1.0;
   Zaxis->z = 1.0;

   startTheta = tor->startTheta;
   endTheta = tor->endTheta;
   startPsi = tor->startPsi;
   endPsi = tor->endPsi;
   c = tor->c;
   a = tor->a;
   midTheta = (startTheta+endTheta)/2.;
   midPsi = (startPsi+endPsi)/2.;
   one4thTheta = startTheta + (endTheta-startTheta)/4.0;
   one4thPsi = startPsi + (endPsi-startPsi)/4.0;
   three4thTheta = startTheta + 3.0 * (endTheta-startTheta)/4.0;
   three4thPsi = startPsi + 3.0 * (endPsi-startPsi)/4.0;

   otherterms = 0.0;
   // p1
   tmpTOR = TOR_allocate(c, a, one4thTheta, midTheta, three4thPsi, endPsi,
                         origin, Zaxis, Xaxis, Yaxis, 0, 0, tor->cavity, 0, 0, 0);
   Integration_general_TOR_single_local(point, tmpTOR, kernel, parameters, &otherterms);
   TOR_free(tmpTOR);

   // p2
   tmpTOR = TOR_allocate(c, a, startTheta, one4thTheta, three4thPsi, endPsi,
                         origin, Zaxis, Xaxis, Yaxis, 0, 0, tor->cavity, 0, 0, 0);
   Integration_general_TOR_single_local(point, tmpTOR, kernel, parameters, &otherterms);
   TOR_free(tmpTOR);
  
   // p3
   tmpTOR = TOR_allocate(c, a, startTheta, one4thTheta, midPsi, three4thPsi,
                         origin, Zaxis, Xaxis, Yaxis, 0 ,0, tor->cavity, 0, 0 ,0);
   Integration_general_TOR_single_local(point, tmpTOR, kernel, parameters, &otherterms);
   TOR_free(tmpTOR);
  
   // p4
   tmpTOR = TOR_allocate(c, a, startTheta, one4thTheta, one4thPsi, midPsi,
                         origin, Zaxis, Xaxis, Yaxis, 0, 0, tor->cavity, 0, 0, 0);
   Integration_general_TOR_single_local(point, tmpTOR, kernel, parameters, &otherterms);
   TOR_free(tmpTOR);

   // p5
   tmpTOR = TOR_allocate(c, a, startTheta, one4thTheta, startPsi, one4thPsi,
                         origin, Zaxis, Xaxis, Yaxis, 0 ,0, tor->cavity, 0 ,0 ,0);
   Integration_general_TOR_single_local(point, tmpTOR, kernel, parameters, &otherterms);
   TOR_free(tmpTOR);

   // p6
   tmpTOR = TOR_allocate(c, a, one4thTheta, midTheta, startPsi, one4thPsi,
                         origin, Zaxis, Xaxis, Yaxis, 0, 0, tor->cavity, 0, 0 ,0);
   Integration_general_TOR_single_local(point, tmpTOR, kernel, parameters, &otherterms);
   TOR_free(tmpTOR);

   *slp += 2 * otherterms;

   // remaining self-ness
   tmpTOR = TOR_allocate(c, a, one4thTheta, three4thTheta, one4thPsi, three4thPsi,
                         origin, Zaxis, Xaxis, Yaxis, 0,0 ,tor->cavity, 0,0,0);
   Integration_general_TOR_single_self(point, tmpTOR, kernel, parameters, slp);
   TOR_free(tmpTOR);

   // free everything else
   Vector3D_free(origin);
   Vector3D_free(Zaxis);
   Vector3D_free(Xaxis);
   Vector3D_free(Yaxis);
    
}

void Integration_general_TOR_single_near(Vector3D point, TOR tor,
                                         BEMKernelType kernel, void* parameters, real* slp) {
   real startTheta, startPsi, endTheta, endPsi, c, a;
   real midTheta, midPsi, one4thTheta, one4thPsi, three4thTheta, three4thPsi;
   real otherterms, curInt;
   Vector3D origin, Xaxis, Yaxis, Zaxis;
   TOR tmpTOR;
   // added to prevent single-precision futz
   real thresholdForSmallPanelApprox = 1e-5;
   real maxRad = tor->maxEdgelength;
   if (maxRad < thresholdForSmallPanelApprox) {
      return;
   }
   //end
   origin = Vector3D_allocate();
   Xaxis = Vector3D_allocate();
   Yaxis = Vector3D_allocate();
   Zaxis = Vector3D_allocate();
   Xaxis->x = 1.0;
   Yaxis->y = 1.0;
   Zaxis->z = 1.0;
  
   startTheta = tor->startTheta;
   endTheta = tor->endTheta;
   startPsi = tor->startPsi;
   endPsi = tor->endPsi;
   c = tor->c;
   a = tor->a;
   midTheta = (startTheta+endTheta)/2.;
   midPsi = (startPsi+endPsi)/2.;

   // p1
   tmpTOR = TOR_allocate(c, a, startTheta, midTheta, startPsi, midPsi,
                         origin, Zaxis, Xaxis, Yaxis, 0, 0, tor->cavity, 0, 0, 0);
   Integration_general_TOR_single_local(point, tmpTOR, kernel, parameters, slp);
   TOR_free(tmpTOR);

   // p2
   tmpTOR = TOR_allocate(c, a, midTheta, endTheta, startPsi, midPsi,
                         origin, Zaxis, Xaxis, Yaxis, 0, 0, tor->cavity, 0, 0, 0);
   Integration_general_TOR_single_local(point, tmpTOR, kernel, parameters, slp);
   TOR_free(tmpTOR);

   // p3
   tmpTOR = TOR_allocate(c, a, startTheta, midTheta,  midPsi, endPsi,
                         origin, Zaxis, Xaxis, Yaxis, 0, 0, tor->cavity, 0, 0, 0);
   Integration_general_TOR_single_local(point, tmpTOR, kernel, parameters, slp);
   TOR_free(tmpTOR);

   // p4
   tmpTOR = TOR_allocate(c, a, midTheta, endTheta, midPsi, endPsi,
                         origin, Zaxis, Xaxis, Yaxis, 0, 0, tor->cavity, 0, 0, 0);
   Integration_general_TOR_single_local(point, tmpTOR, kernel, parameters, slp);
   TOR_free(tmpTOR);

   Vector3D_free(origin);
   Vector3D_free(Xaxis);
   Vector3D_free(Yaxis);
   Vector3D_free(Zaxis);
}

void Integration_oneoverr_TOR_double(Vector3D point, TOR tor, void* parameters, real* dlp) {
   real self_term_distance_tolerance = 1e-5;
   real near_term_distance_tolerance = 3.0;  // equivalent to threshRadius in matlab

   //if (Vector3D_distance(point, tor->centroid) < near_term_distance_tolerance * tor->maxEdgelength)
   //if (Vector3D_distance(point, tor->centroid) < near_term_distance_tolerance * sqrt(tor->area))
      Integration_oneoverr_TOR_double_DJW(point, tor, parameters, dlp);

   // should no longer be needed, handled by Panel_quadrature
/*
   else {
      Vector3D localpnt = Vector3D_allocate();

      // NOTE: this is not the same order of transformation as done in matlab version!!  there,
      // everything is done in A' * (x - y) to get into local coords.  here we do A * x - y
      Vector3D_transformVecs(localpnt, tor->Tmat[0], tor->Tmat[1], tor->Tmat[2], point);
      Vector3D_add(localpnt, localpnt, tor->Tvec);

      Integration_oneoverr_TOR_double_local(localpnt, tor, parameters, dlp);

      Vector3D_free(localpnt);
   }
*/
}

// implements calcpcd_TOR.m
void Integration_oneoverr_TOR_double_DJW(Vector3D point, TOR tor, void* parameters, real* dlp) {
   Vector3D localpnt; // is our observation point: will be the origin

   Vector3D centroid;
   Vector3D normalAtCentroid;

   static QuadratureRule DJW_dipole = NULL;
   static QuadratureRule DJW_dipole_self = NULL;

   Vector3D center, normal,localX, localY;
   Vector3D X, Y, Z;
   real sphRadius;
   real v[4][2];
   unsigned int i, next, j;
   real self_term_tolerance, gen_axis_tolerance, off_axis_tolerance;
   real theta;
   Vector3D circleCenter;
   Vector3D X_l, Y_l;
   real circleRadius;
   real t1, t2;
   real psi;

   real t;
   Vector3D p, dpdt, pProj;
   real ProjPhi, ddu_atan_u, dudt, detJ;

#ifdef OMP
#pragma omp critical
#endif
   if (DJW_dipole == NULL) {
      DJW_dipole = QuadratureRule_allocate(64);
   }
#ifdef OMP
#pragma omp critical
#endif
   if (DJW_dipole_self == NULL) {
      DJW_dipole_self = QuadratureRule_allocate(512);
   }
   QuadratureRule qr = DJW_dipole;

   *dlp = 0.0;

   localpnt = Vector3D_allocate();
   centroid = Vector3D_allocate();
   normalAtCentroid = Vector3D_allocate();
   X = Vector3D_allocate();
   Y = Vector3D_allocate();
   Z = Vector3D_allocate();
   center = Vector3D_allocate();
   normal = Vector3D_allocate();
   localX = Vector3D_allocate();
   localY = Vector3D_allocate();
  
   circleCenter = Vector3D_allocate();
   X_l = Vector3D_allocate();
   Y_l = Vector3D_allocate();
   p = Vector3D_allocate();
   dpdt = Vector3D_allocate();
   pProj = Vector3D_allocate();

   Vector3D_copy(localpnt, point);
   Vector3D_sub(localpnt, localpnt, point);
   Vector3D_sub(centroid, tor->centroid, point);
   Vector3D_sub(center, tor->center, point);
   Vector3D_copy(normalAtCentroid, tor->normalAtCentroid);
   Vector3D_copy(normal, tor->normal);
   Vector3D_copy(localX, tor->localX);
   Vector3D_copy(localY, tor->localY);

   Vector3D_transformVecs_inverse(centroid, tor->localX, tor->localY, tor->normal, centroid);
   Vector3D_transformVecs_inverse(center, tor->localX, tor->localY, tor->normal, center);
   Vector3D_transformVecs_inverse(normalAtCentroid, tor->localX, tor->localY, tor->normal, normalAtCentroid);
   Vector3D_transformVecs_inverse(normal, tor->localX, tor->localY, tor->normal, normal);
   Vector3D_transformVecs_inverse(localX, tor->localX, tor->localY, tor->normal, localX);
   Vector3D_transformVecs_inverse(localY, tor->localX, tor->localY, tor->normal, localY);
  
   sphRadius = 8;  // it don't matter to jesus
   self_term_tolerance = 1e-5;
   off_axis_tolerance = 1e-5;
   gen_axis_tolerance = 1e-5;
   if (Vector3D_distance(centroid, localpnt) < self_term_tolerance) {
      // self term
      Vector3D_copy(Z, normalAtCentroid);
      X->z = 1.0;
      Vector3D_addscaled(X, X, -Vector3D_dot(X,Z), Z);
      if (Vector3D_length(X) < gen_axis_tolerance) {
         X->x = 1.0; X->y = 0.0; X->z = 0.0;
         Vector3D_addscaled(X, X, -Vector3D_dot(X,Z), Z);
      }
      Vector3D_normalize(X);
      Vector3D_cross(Y, Z, X);
      qr = DJW_dipole_self;
   } else if (sqrt(centroid->x*centroid->x + centroid->y * centroid->y) < off_axis_tolerance) {
      Vector3D_copy(Z, centroid);
      Vector3D_normalize(Z);
      X->z = 1.0;
      Vector3D_addscaled(X, X, -Vector3D_dot(X,Z), Z);
      if (Vector3D_length(X) < gen_axis_tolerance) {
         X->x = 1.0; X->y = 0.0; X->z = 0.0;
         Vector3D_addscaled(X, X, -Vector3D_dot(X,Z), Z);
      }
      Vector3D_normalize(X);
      Vector3D_cross(Y, Z, X);

   } else {
      X->x = 1.0;
      Z->z = (centroid->z > 0.0)?1.0:-1.0;
      Vector3D_cross(Y, Z, X);
   }

   Vector3D_transformVecs_inverse(localpnt, X, Y, Z, localpnt);
   Vector3D_transformVecs_inverse(centroid, X, Y, Z, centroid);
   Vector3D_transformVecs_inverse(normalAtCentroid, X, Y, Z, normalAtCentroid);
   Vector3D_transformVecs_inverse(normal, X, Y, Z, normal);
   Vector3D_transformVecs_inverse(center, X, Y, Z, center);
   Vector3D_transformVecs_inverse(localX, X, Y, Z, localX);
   Vector3D_transformVecs_inverse(localY, X, Y, Z, localY);

   v[0][0] = tor->startTheta;  v[0][1] = tor->startPsi;
   v[1][0] = tor->startTheta;  v[1][1] = tor->endPsi;
   v[2][0] = tor->endTheta;    v[2][1] = tor->endPsi;
   v[3][0] = tor->endTheta;    v[3][1] = tor->startPsi;

   for (i = 0; i < 4; i++) {
      next = i + 1;
      if (next == 4)
         next = 0;

      if (i & 1) { // theta edge (looks backwards from matlab code b/c we're 0 based)
         psi = v[i][1];
         Vector3D_addscaled(circleCenter, center, tor->a * sin(psi), normal);
         circleRadius = tor->c + tor->a * cos(psi);
         Vector3D_copy(X_l, localX);
         Vector3D_copy(Y_l, localY);
         t1 = v[i][0];
         t2 = v[next][0];
      } else {
         theta = v[i][0];
         Vector3D_copy(circleCenter, center);
         Vector3D_addscaled(circleCenter, circleCenter, tor->c * cos(theta), localX);
         Vector3D_addscaled(circleCenter, circleCenter, tor->c * sin(theta), localY);
         Vector3D_copy(X_l, localX);
         Vector3D_scale(X_l, cos(theta));
         Vector3D_addscaled(X_l, X_l, sin(theta), localY);
         Vector3D_normalize(X_l);
         Vector3D_copy(Y_l, normal);
         circleRadius = tor->a;
         t1 = v[i][1];
         t2 = v[next][1];
      }

      // original loop

/*
      for (j = 0; j < qr->order; j++) {
         t = t1 + (t2 - t1) * qr->x[j];
         Vector3D_copy(p, circleCenter);
         Vector3D_addscaled(p, p, circleRadius * cos(t), X_l);
         Vector3D_addscaled(p, p, circleRadius * sin(t), Y_l);
         Vector3D_copy(dpdt, X_l);
         Vector3D_scale(dpdt, -sin(t) * circleRadius);
         Vector3D_addscaled(dpdt, dpdt, circleRadius * cos(t), Y_l);
         findSpherePoint(pProj, localpnt, sphRadius, p);
         //ProjPhi = acos(pProj->z / sphRadius);
         ddu_atan_u = (p->x*p->x) / (p->x*p->x + p->y * p->y);
         dudt = (p->y*dpdt->x - p->x *dpdt->y) / (p->x*p->x);
         detJ = ddu_atan_u * dudt;

         *dlp += (t2-t1) * qr->w[j] * sphRadius*sphRadius
            * detJ * (1 - pProj->z/sphRadius);
      }
*/

      // new super optimized loop

      real rsint, rcost, sphRadius2 = sphRadius * sphRadius, length;

      for (j = 0; j < qr->order; j++) {
         t = t1 + (t2 - t1) * qr->x[j];

         // collect these often used expensive terms
         rsint = circleRadius * sin(t);
         rcost = circleRadius * cos(t);

         // remove function calls from setting up p
         p->x = circleCenter->x + rcost * X_l->x + rsint * Y_l->x;
         p->y = circleCenter->y + rcost * X_l->y + rsint * Y_l->y;
         p->z = circleCenter->z + rcost * X_l->z + rsint * Y_l->z;

         // findSpherePoint replacement
         pProj->x = p->x - localpnt->x;
         pProj->y = p->y - localpnt->y;
         pProj->z = p->z - localpnt->z;
         length = sqrt(pProj->x*pProj->x+pProj->y*pProj->y+pProj->z*pProj->z);
         pProj->z = localpnt->z + sphRadius / length * pProj->z;

         // remove function calls from setting up dpdt
         dpdt->x = -rsint * X_l->x + rcost * Y_l->x;
         dpdt->y = -rsint * X_l->y + rcost * Y_l->y;

         detJ = (p->y*dpdt->x - p->x*dpdt->y) / (p->x*p->x + p->y*p->y);

         *dlp += (t2-t1) * qr->w[j] * sphRadius2
            * detJ * (1.0 - pProj->z/sphRadius);
      }
    
   }

   *dlp /= sphRadius * sphRadius;

   // JAY BOO (removed to allow flipping of torus panels)

   //if (tor->endTheta < tor->startTheta)
   //   *dlp = -(*dlp);

   // MALTY BOO
                                                                
   if (Vector3D_equal(point, tor->centroid))
      *dlp += 4.0 * M_PI;
  
   Vector3D_free(localpnt);
   Vector3D_free(centroid);
   Vector3D_free(normalAtCentroid);
   Vector3D_free(X);
   Vector3D_free(Y);
   Vector3D_free(Z);
   Vector3D_free(circleCenter);
   Vector3D_free(X_l);
   Vector3D_free(Y_l);
   Vector3D_free(p);
   Vector3D_free(dpdt);
   Vector3D_free(pProj);
   Vector3D_free(normal);
   Vector3D_free(center);
   Vector3D_free(localX);
   Vector3D_free(localY);
}

void Integration_oneoverr_TOR_double_local(Vector3D point, TOR tor,
                                           void* parameters, real* dlp) {
  
   real startTheta, endTheta, startPsi, endPsi, c, a, dTheta, dPsi;
   real curTheta, curPsi;
   real GnFcnVal, detJ, curWeight;
   real near_tolerance = 2.0;
   Vector3D srcpoint = Vector3D_allocate();
   Vector3D normal = Vector3D_allocate();
   unsigned int i, j;

   startTheta = tor->startTheta;
   endTheta = tor->endTheta;
   startPsi =tor->startPsi;
   endPsi = tor->endPsi;
   c = tor->c;
   a = tor->a;
   dTheta = endTheta - startTheta;
   dPsi = endPsi - startPsi;
   static QuadratureRule qr = NULL;

#ifdef OMP
#pragma omp critical
#endif
   if (qr == NULL) {
      qr = QuadratureRule_allocate(4);
   }

   *dlp = 0.0;

   for (i = 0; i < qr->order; i++) {
      curTheta = startTheta + qr->x[i] * dTheta;
      for (j = 0; j < qr->order; j++) {
         curPsi = startPsi + qr->x[j] * dPsi;
         torToCart(srcpoint, normal, c, a, curTheta, curPsi);
         if (tor->cavity)
            Vector3D_scale(normal, -1.0);
         GnFcnVal = GreensFunction_deriv(point, srcpoint, POISSON_KERNEL, parameters, normal);
         detJ = a * (c + a * cos(curPsi));
         curWeight = fabs(dTheta) * qr->w[i] * fabs(dPsi) * qr->w[j];
         *dlp += GnFcnVal * detJ * curWeight;
      }
   }

   Vector3D_free(srcpoint);
   Vector3D_free(normal);

   // malty BOO, answer seems to be negated

   *dlp = -(*dlp);
}
