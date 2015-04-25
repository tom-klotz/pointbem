#include "meshmaker.h"

real realabs(real a) {
   if (a > 0.0)
      return a;
   else
      return -a;
}

void Integration_oneoverr_GST(Vector3D point, GST gst, void* parameters, real* slp, real *dlp) {
   Integration_oneoverr_GST_single(point, gst, parameters, slp);
   Integration_oneoverr_GST_double(point, gst, parameters, dlp);
}

void Integration_oneoverr_GST_single(Vector3D point, GST gst, void* parameters, real* slp) {
   Integration_general_GST_single(point, gst, POISSON_KERNEL, parameters, slp);
}

void Integration_general_GST_single(Vector3D point, GST gst,
                                    BEMKernelType kernel, void* parameters, real *slp) {
   unsigned int order = 4;
   unsigned int numcoeffs, i;
   real *IntegralsOfPolynomialTerms;
   unsigned int *kpvec;
   unsigned int *epvec;
   unsigned int numquadpoints;
   Vector coeffs;
   Vector Pvec;
   numquadpoints = gst->flatpanel->numquadpoints;
   GST_computePolynomialIntegrals(&numcoeffs, &IntegralsOfPolynomialTerms, &kpvec, &epvec,
                                  order, gst, point, kernel, SINGLE_LAYER_INT, parameters);

   Pvec = Vector_allocate(numquadpoints);
   GST_computePolynomialValues(Pvec, gst, point, kernel, SINGLE_LAYER_INT, parameters);
  
   coeffs = Vector_allocate(numcoeffs);

   /*
     Matrix vdmInverse = Matrix_allocate(numcoeffs, numquadpoints);
     Matrix_pseudoinverse(vdmInverse, vandermondeMat, numquadpoints, numcoeffs);
     Matrix_multiplyvector(coeffs, vdmInverse, Pvec, numcoeffs, numquadpoints);
     Matrix_free(vdmInverse);
   */

   //Matrix_lsq_solve(coeffs, vandermondeMat, Pvec, numquadpoints, numcoeffs);

#ifdef OMP
#pragma omp critical
#endif
   if (gst->vandermondeInverse == NULL) {
      Matrix vandermondeMat = Matrix_allocate(numquadpoints, numcoeffs);
      GST_computeVandermondeMatrix(vandermondeMat, numcoeffs, gst, kpvec, epvec);
      gst->vandermondeInverse = Matrix_allocate(numcoeffs, numquadpoints);
      Matrix_pseudoinverse_droptol(gst->vandermondeInverse, vandermondeMat, numquadpoints, numcoeffs, 1e-6);
      Matrix_free(vandermondeMat);
   }

   Matrix_multiplyvector(coeffs, gst->vandermondeInverse, Pvec, numcoeffs, numquadpoints);

   *slp = 0.0;
   for (i = 0; i < numcoeffs; i++) {
      *slp += coeffs[i] * IntegralsOfPolynomialTerms[i];
   }

   Vector_free(coeffs);
   Vector_free(Pvec);
   free(IntegralsOfPolynomialTerms);
   free(kpvec);
   free(epvec);
}

void Integration_general_GST_double(Vector3D point, GST gst,
                                    BEMKernelType kernel, void* parameters, real *dlp) {
   unsigned int order = 4;
   unsigned int numcoeffs, i;
   real *IntegralsOfPolynomialTerms;
   unsigned int *kpvec;
   unsigned int *epvec;
   unsigned int numquadpoints;
   Vector coeffs;
   Vector Pvec;
   numquadpoints = gst->flatpanel->numquadpoints;
   GST_computePolynomialIntegrals(&numcoeffs, &IntegralsOfPolynomialTerms, &kpvec, &epvec,
                                  order, gst, point, kernel, DOUBLE_LAYER_INT, parameters);

   Pvec = Vector_allocate(numquadpoints);
   GST_computePolynomialValues(Pvec, gst, point, kernel, DOUBLE_LAYER_INT, parameters);
  
   coeffs = Vector_allocate(numcoeffs);

   /*
     Matrix vdmInverse = Matrix_allocate(numcoeffs, numquadpoints);
     Matrix_pseudoinverse(vdmInverse, vandermondeMat, numquadpoints, numcoeffs);
     Matrix_multiplyvector(coeffs, vdmInverse, Pvec, numcoeffs, numquadpoints);
     Matrix_free(vdmInverse);
   */

   //Matrix_lsq_solve(coeffs, vandermondeMat, Pvec, numquadpoints, numcoeffs);

#ifdef OMP
#pragma omp critical
#endif
   if (gst->vandermondeInverse == NULL) {
      Matrix vandermondeMat = Matrix_allocate(numquadpoints, numcoeffs);
      GST_computeVandermondeMatrix(vandermondeMat, numcoeffs, gst, kpvec, epvec);
      gst->vandermondeInverse = Matrix_allocate(numcoeffs, numquadpoints);
      Matrix_pseudoinverse_droptol(gst->vandermondeInverse, vandermondeMat, numquadpoints, numcoeffs, 1e-6);
      Matrix_free(vandermondeMat);
   }

   Matrix_multiplyvector(coeffs, gst->vandermondeInverse, Pvec, numcoeffs, numquadpoints);

   *dlp = 0.0;
   for (i = 0; i < numcoeffs; i++) {
      *dlp += coeffs[i] * IntegralsOfPolynomialTerms[i];
   }

   Vector_free(coeffs);
   Vector_free(Pvec);
   free(IntegralsOfPolynomialTerms);
   free(kpvec);
   free(epvec);
}

void GST_computePolynomialIntegrals(unsigned int* numcoeffs, real** Integrals,
                                    unsigned int** kpvec, unsigned int** epvec,
                                    unsigned int order, GST gst, Vector3D point,                                  BEMKernelType kernel, BEMLayerType layertype, void* parameters) {
   Vector3D localpnt, quadpoint, normal;
   real checkRadius = 4.0;
   *numcoeffs = (unsigned int) ((order * (order + 1))/2);
   *Integrals = (real *)calloc(*numcoeffs, sizeof(real));
   *kpvec = (unsigned int *)calloc(*numcoeffs, sizeof(unsigned int));
   *epvec = (unsigned int *)calloc(*numcoeffs, sizeof(unsigned int));
   unsigned int ksiPower, etaPower, curorder;
   unsigned int i, j;
   real polyVal, GreensFcnVal, fcnVal;
   real *ellipse_integrals;
   ellipse_integrals = (real *)calloc(*numcoeffs, sizeof(real));
  
   localpnt = Vector3D_allocate();
   quadpoint = Vector3D_allocate();
   normal = Vector3D_allocate();
   normal->z = 1.0;  // flat panel points straight up in z ; can probably just use flatpanel->normal
   // convert point to flat panel local coords
   Vector3D_transformVecs(localpnt, gst->Tmat[0], gst->Tmat[1], gst->Tmat[2], point);
   Vector3D_add(localpnt, localpnt, gst->Tvec);
  
   // do panel moment calculations: if nearfield, do via xin wang (thanks xin!) recurision code
   curorder = 0; ksiPower = 0; etaPower = 0;
   if ( (Vector3D_distance(localpnt, gst->flatpanel->centroid) < gst->flatpanel->maxEdgelength * checkRadius)
        &&
        (kernel == POISSON_KERNEL) ) {
      greenInt_XinWang(*Integrals, gst->flatpanel, localpnt, order);
   } else { // calculate far field moments via numquad
      for (i = 0; i < *numcoeffs; i++) {
         for (j = 0; j < gst->flatpanel->numquadpoints; j++) {
            polyVal = PetscPowRealInt(gst->flatpanel->ksi[j], ksiPower) * PetscPowRealInt(gst->flatpanel->eta[j],  etaPower);
            quadpoint->x = gst->flatpanel->ksi[j];
            quadpoint->y = gst->flatpanel->eta[j];
            if (layertype == SINGLE_LAYER_INT)
               GreensFcnVal = GreensFunction(localpnt, quadpoint, kernel, parameters);
            else
               GreensFcnVal = GreensFunction_deriv(localpnt, quadpoint, kernel, parameters, gst->flatpanel->normal);
            fcnVal = polyVal * GreensFcnVal;
            (*Integrals)[i] += gst->flatpanel->weights[j] * fcnVal;
         }
         // update ksipower, etapower
         if (ksiPower == 0) {
            curorder = curorder + 1;
            ksiPower = curorder;
            etaPower = 0;
         } else {
            etaPower = etaPower + 1;
            ksiPower = ksiPower - 1;
         }
      }
   }

   // assemble kpvec, epvec
   curorder = 0;
   ksiPower = 0;
   etaPower = 0;
   for (i = 0; i < *numcoeffs; i++) {
      (*kpvec)[i] = ksiPower;
      (*epvec)[i] = etaPower;
      if (ksiPower == 0) {
         curorder = curorder + 1;
         ksiPower = curorder;
         etaPower = 0;
      } else {
         etaPower = etaPower + 1;
         ksiPower = ksiPower - 1;
      }
   }

   // find elliptical corrections of polynomial terms
   ellipse_int2(ellipse_integrals, gst, point, *numcoeffs, *kpvec, *epvec, kernel, layertype, parameters);

   for (i = 0; i < *numcoeffs; i++)
      (*Integrals)[i] += ellipse_integrals[i];

   free(ellipse_integrals);
   Vector3D_free(localpnt);
   Vector3D_free(quadpoint);
   Vector3D_free(normal);
}

void GST_computeVandermondeMatrix(Matrix vandermondeMat, unsigned int numcoeffs, GST gst, unsigned int* kpvec, unsigned int* epvec) {
   unsigned int i, j;
   for (i = 0; i < gst->flatpanel->numquadpoints; i++) {
      for (j = 0; j < numcoeffs; j++) {
         vandermondeMat[i][j] = PetscPowRealInt(gst->flatpanel->ksi[i], kpvec[j]) * PetscPowRealInt(gst->flatpanel->eta[i], epvec[j]);
      }
   }
}

void GST_computePolynomialValues(Vector Pvec, GST gst, Vector3D point, BEMKernelType kernel, BEMLayerType layertype, void* parameters) {
   unsigned int i;
   Vector3D Xflat, Xcurved, normal, normalSph, localpnt, transformedGSTcenter;
   real G_flat, G_curve, Jdet;

   Xflat = Vector3D_allocate();
   Xcurved = Vector3D_allocate();
   normal = Vector3D_allocate();
   normalSph = Vector3D_allocate();
   localpnt = Vector3D_allocate();
   transformedGSTcenter = Vector3D_allocate();
  
   // transform obspnt, gst center
   Vector3D_transformVecs(localpnt, gst->Tmat[0], gst->Tmat[1], gst->Tmat[2], point);
   Vector3D_add(localpnt, localpnt, gst->Tvec);
   Vector3D_transformVecs(transformedGSTcenter, gst->Tmat[0], gst->Tmat[1], gst->Tmat[2], gst->center);
   Vector3D_add(transformedGSTcenter, transformedGSTcenter, gst->Tvec);
   normal->z = 1.0;

   for (i = 0; i < gst->flatpanel->numquadpoints; i++) {
      Xflat->x = gst->flatpanel->ksi[i];
      Xflat->y = gst->flatpanel->eta[i];
      if (layertype == SINGLE_LAYER_INT)
         G_flat = GreensFunction(localpnt, Xflat, kernel, parameters);
      else
         G_flat = GreensFunction_deriv(localpnt, Xflat, kernel, parameters, normal);
    
      findSpherePoint(Xcurved, transformedGSTcenter, gst->radius, Xflat);
      Vector3D_sub(normalSph, Xcurved, transformedGSTcenter);
      Vector3D_normalize(normalSph);

      if (layertype == SINGLE_LAYER_INT)
         G_curve = GreensFunction(localpnt, Xcurved, kernel, parameters);
      else
         G_curve = GreensFunction_deriv(localpnt, Xcurved, kernel, parameters, normalSph);

      Jdet = findJacobianDet(transformedGSTcenter, gst->radius, Xflat);

      Pvec[i] = (G_curve / G_flat) * Jdet;
   }

   Vector3D_free(localpnt);
   Vector3D_free(normalSph);
   Vector3D_free(transformedGSTcenter);
   Vector3D_free(Xflat);
   Vector3D_free(Xcurved);
   Vector3D_free(normal);
}

real findJacobianDet(Vector3D center, real radius, Vector3D pnt) {
   real u, v, r, r3;
   real dxdksiA, dxdetaA, dydksiA, dydetaA, dzdksiA, dzdetaA, Jdet;
   Vector3D planepnt = Vector3D_allocate();
   planepnt->x = pnt->x;
   planepnt->y = pnt->y;
   planepnt->z = 0.0;

   u = planepnt->x - center->x;
   v = planepnt->y - center->y;
   r = Vector3D_distance(center, planepnt);
   r3 = r * r * r;

   dxdksiA = radius * (v*v + center->z * center->z) / r3;
   dxdetaA = radius * -u * v / r3;
   dydksiA = radius * -u * v / r3;
   dydetaA = radius * (u*u + center->z * center->z) / r3;
   dzdksiA = radius * -center->z * u / r3;
   dzdetaA = radius * -center->z * v / r3;

   Jdet = sqrt( PetscPowRealInt(dxdksiA*dydetaA - dydksiA*dxdetaA,2) +
                PetscPowRealInt(dydksiA*dzdetaA - dzdksiA*dydetaA,2) +
                PetscPowRealInt(dzdksiA*dxdetaA - dxdksiA*dzdetaA,2));
  
   Vector3D_free(planepnt);
   return Jdet;
}

void ellipse_int2(real *ellipse_integrals, GST gst, Vector3D point, unsigned int numcoeffs,
                  unsigned int* kpvec, unsigned int* epvec, BEMKernelType kernel, BEMLayerType layertype, void* parameters) {
   real nma_tolerance = 1e-6;  // SINGLE PRECISION WARNING BLAH DE BLAH
   unsigned int i;
   for (i = 0; i < numcoeffs; i++)
      ellipse_integrals[i] = 0.0;
  
   for (i = 0; i < 3; i++) {
      if (gst->factor[i] != 0) {
         conic_section_integrate(ellipse_integrals, i, gst, point, numcoeffs, kpvec, epvec, kernel, layertype, parameters);
      }
   }
  
}

// implements conic_section_integrate3.m
void conic_section_integrate(real *ellipse_integrals, unsigned int arcnumber, GST gst, Vector3D point,
                             unsigned int numcoeffs, unsigned int* kpvec, unsigned int* epvec,
                             BEMKernelType kernel, BEMLayerType layertype, void *parameters) {
   Vector3D pntInFlatpanel,vert1,vert2, edgePointAtTheta, curQuadPoint, curQuadPointInFlatpanel;
   Vector3D normal;
   real theta1, theta2, asquared, bsquared;
   real radius1,radius2, deltaTheta, curtheta;
   real rconic, rline, startR, endR, deltaR, curR;
   real curGreensFunc, curWeight, ksi, eta, curFcnValue;
  
   unsigned int i, j, k;
   static QuadratureRule qr = NULL;
#ifdef OMP
#pragma omp critical
#endif
   if (qr == NULL) {
      qr = QuadratureRule_allocate(4);
   }
   pntInFlatpanel = Vector3D_allocate();
   vert1 = Vector3D_allocate();
   vert2 = Vector3D_allocate();
   edgePointAtTheta = Vector3D_allocate();
   curQuadPoint = Vector3D_allocate();
   curQuadPointInFlatpanel = Vector3D_allocate();
   normal = Vector3D_allocate();
  
   Vector3D_transformVecs(pntInFlatpanel, gst->Tmat[0], gst->Tmat[1], gst->Tmat[2], point);
   Vector3D_add(pntInFlatpanel, pntInFlatpanel, gst->Tvec);

   theta1 = gst->conic[arcnumber]->theta1Int;
   theta2 = gst->conic[arcnumber]->theta2Int;
   asquared = gst->conic[arcnumber]->asquared;
   bsquared = gst->conic[arcnumber]->bsquared;

   radius1 = (asquared * bsquared) / (bsquared * cos(theta1) * cos(theta1) - asquared * sin(theta1) * sin(theta1));
   vert1->x = sqrt(radius1) * cos(theta1);
   vert1->y = sqrt(radius1) * sin(theta1);
   radius2 = (asquared * bsquared) / (bsquared * cos(theta2) * cos(theta2) - asquared * sin(theta2) * sin(theta2));
   vert2->x = sqrt(radius2) * cos(theta2);
   vert2->y = sqrt(radius2) * sin(theta2);

   normal->z = 1.0;

   deltaTheta = theta2 - theta1;
   for (i = 0; i < qr->order; i++) {
      curtheta = theta1 + qr->x[i] * deltaTheta;
      rconic = sqrt((asquared * bsquared) / (bsquared * cos(curtheta) * cos(curtheta) - asquared * sin(curtheta) * sin(curtheta)));
      findIntersection(edgePointAtTheta, vert1, vert2, curtheta); // malty i trust this gets inlined??
      rline = Vector3D_length(edgePointAtTheta);

      if (rline < rconic) {
         startR = rline;
         endR = rconic;
      } else {
         startR = rconic;
         endR = rline;
      }
      deltaR = endR - startR;
    
      for (j = 0; j < qr->order; j++) {
         curR = startR + qr->x[j] * deltaR;
         curQuadPoint->x = curR * cos(curtheta);
         curQuadPoint->y = curR * sin(curtheta);
         Vector3D_sub(curQuadPointInFlatpanel, curQuadPoint, gst->conic[arcnumber]->Tvec);
         Vector3D_transformVecs_inverse(curQuadPointInFlatpanel, gst->conic[arcnumber]->Tmat[0],
                                        gst->conic[arcnumber]->Tmat[1],
                                        gst->conic[arcnumber]->Tmat[2], curQuadPointInFlatpanel);
         if (layertype == SINGLE_LAYER_INT)
            curGreensFunc = GreensFunction(pntInFlatpanel, curQuadPointInFlatpanel, kernel, parameters);
         else
            curGreensFunc = GreensFunction_deriv(pntInFlatpanel, curQuadPointInFlatpanel, kernel, parameters, normal);
         curWeight = realabs(deltaTheta) * qr->w[i] * realabs(deltaR) * qr->w[j] * curR;
         ksi = curQuadPointInFlatpanel->x;
         eta = curQuadPointInFlatpanel->y;
         for (k = 0; k < numcoeffs; k++) {
            curFcnValue = PetscPowRealInt(ksi, kpvec[k]) * PetscPowRealInt(eta, epvec[k]) * curGreensFunc;
            ellipse_integrals[k] += gst->factor[arcnumber] * curFcnValue * curWeight;
         }
      }
   }
  
   Vector3D_free(vert1);
   Vector3D_free(vert2);
   Vector3D_free(pntInFlatpanel);
   Vector3D_free(edgePointAtTheta);
   Vector3D_free(curQuadPoint);
   Vector3D_free(curQuadPointInFlatpanel);
   Vector3D_free(normal);
}

void findIntersection(Vector3D point, Vector3D v1, Vector3D v2, real theta) {
   real angleEdge, rhsEdge, invdet;
   Vector3D edge;

   edge = Vector3D_allocate();
   Vector3D_sub(edge, v1, v2);
   angleEdge = atan2(edge->y, edge->x);
   rhsEdge = -sin(angleEdge) * v2->x + cos(angleEdge) * v2->y;
   invdet = 1.0 / (-sin(angleEdge) * cos(theta) - cos(angleEdge) * -sin(theta));
   point->x = invdet * cos(theta) * rhsEdge;
   point->y = invdet * sin(theta) * rhsEdge;
   point->z = 0;
   Vector3D_free(edge);
}


// implements calcpcd_GST.m
void Integration_oneoverr_GST_double(Vector3D point, GST gst, void* parameters, real* dlp) {
   unsigned int i, j, next;
   Vector3D sphPoint, normalAt, Tvec, pnt;
   Vector3D Tmat[3], X, Y, Z, v1, v2, ac, v2ac, X_l, Y_l, Z_l;
   Vector3D p, pProj, dpdt;
   real t2, t, sphRadius, radius, cosphi, ddu_atan_u, dudt;
   real self_term_tolerance, along_z_axis_tolerance, detJ;
  

   static QuadratureRule DJW_dipole = NULL;
   static QuadratureRule DJW_dipole_self = NULL;
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

   sphPoint = Vector3D_allocate();
   normalAt = Vector3D_allocate();
   Tvec = Vector3D_allocate();
   pnt = Vector3D_allocate();
   Tmat[0] = Vector3D_allocate();
   Tmat[1] = Vector3D_allocate();
   Tmat[2] = Vector3D_allocate();
   X = Vector3D_allocate();
   Y = Vector3D_allocate();
   Z = Vector3D_allocate();
   v1 = Vector3D_allocate();
   v2 = Vector3D_allocate();
   ac = Vector3D_allocate();
   v2ac = Vector3D_allocate();
   X_l = Vector3D_allocate();
   Y_l = Vector3D_allocate();
   Z_l = Vector3D_allocate();
   p  = Vector3D_allocate();
   pProj = Vector3D_allocate();
   dpdt = Vector3D_allocate();

   Tmat[0]->x = 1.0;
   Tmat[1]->y = 1.0;
   Tmat[2]->z = 1.0;
  
   Vector3D_copy(sphPoint, gst->centroid);
   Vector3D_copy(normalAt, gst->normalAtCentroid);

   // translate so that point is at origin
   Vector3D_sub(Tvec, Tvec, point);  
   Vector3D_sub(pnt, point, point);
   Vector3D_sub(sphPoint, sphPoint, point);

   self_term_tolerance = 1e-6;
   along_z_axis_tolerance = 1e-6;
   if (Vector3D_length(sphPoint) < self_term_tolerance){ // ie, selfterm
      Vector3D_copy(Z, normalAt);
      X->z = 1.0;
      Vector3D_addscaled(X, X, -Vector3D_dot(Z,X), Z);
      Vector3D_normalize(X);
      Vector3D_cross(Y, Z, X);
      qr = DJW_dipole_self;
   } else if (sqrt(sphPoint->x*sphPoint->x+sphPoint->y*sphPoint->y)
              > along_z_axis_tolerance) { // 
      Vector3D_copy(Z, sphPoint);
      Vector3D_normalize(Z);
      X->z = 1.0;
      Vector3D_addscaled(X, X, -Vector3D_dot(Z,X), Z);
      Vector3D_normalize(X);
      Vector3D_cross(Y, Z, X);
   } else { // not sure when this should get triggered...
      X->x = 1.0;
      Z->z = (sphPoint->z > 0.0)?1.0:-1.0;
      Vector3D_cross(Y, Z, X);
   } 

   // rotate coordinate system so panel normal is "up"
   Vector3D_transformVecs_inverse(pnt, X, Y, Z, pnt);
   Vector3D_transformVecs_inverse(Tvec, X, Y, Z, Tvec);
   Vector3D_transformVecs_inverse(Tmat[0], X, Y, Z, Tmat[0]);
   Vector3D_transformVecs_inverse(Tmat[1], X, Y, Z, Tmat[1]);
   Vector3D_transformVecs_inverse(Tmat[2], X, Y, Z, Tmat[2]);


   // determine sphere radius
   sphRadius = 2 ;//2 * (Vector3D_distance(gst->center, pnt) + gst->radius);

   for (i = 0; i < 3; i++) {
      next = i + 1;
      if (next == 3)
         next = 0;

      Vector3D_copy(v1, gst->vertices[i]);
      Vector3D_copy(v2, gst->vertices[next]);
      Vector3D_copy(ac, gst->arccenters[i]);
      // transform them into local coordinate system: normal up Z, pnt at origin
      Vector3D_transformVecs(v1, Tmat[0], Tmat[1], Tmat[2], v1);
      Vector3D_transformVecs(v2, Tmat[0], Tmat[1], Tmat[2], v2);
      Vector3D_transformVecs(ac, Tmat[0], Tmat[1], Tmat[2], ac);
      Vector3D_add(v1, v1, Tvec);
      Vector3D_add(v2, v2, Tvec);
      Vector3D_add(ac, ac, Tvec);

      // now find the circle coordinate system
      Vector3D_sub(X_l, v1, ac);
      Vector3D_sub(Y_l, v2, ac);
      Vector3D_copy(v2ac, Y_l);
      radius = Vector3D_length(X_l);
      Vector3D_normalize(X_l);
      Vector3D_addscaled(Y_l, Y_l, -Vector3D_dot(Y_l, X_l), X_l);
      Vector3D_normalize(Y_l);
      Vector3D_cross(Z_l, X_l, Y_l);

      t2 = atan2(Vector3D_dot(Y_l, v2ac), Vector3D_dot(X_l, v2ac));

      // this is the inner loop, optimize it!

      // original loop

      /*
        for (j = 0; j < qr->order; j++) {
        t = t2 * qr->x[j];  // t1 = 0 (x-axis by definition)

        Vector3D_copy(p, ac);
        Vector3D_addscaled(p, p, radius * cos(t), X_l);
        Vector3D_addscaled(p, p, radius * sin(t), Y_l);

        findSpherePoint(pProj, pnt, sphRadius, p); // pnt is origin, remember!!

        // and sphRadius is the PROJECTION SPHERE RADIUS
        Vector3D_copy(dpdt, X_l);
        Vector3D_scale(dpdt, -sin(t) * radius);
        Vector3D_addscaled(dpdt, dpdt, cos(t) * radius, Y_l);

        phi = acos(pProj->z / sphRadius);

        ddu_atan_u = (p->x * p->x) / (p->x*p->x + p->y*p->y);
        dudt = (p->y * dpdt->x - p->x * dpdt->y) / (p->x*p->x);
        detJ = -ddu_atan_u * dudt;

        *dlp +=  t2 * qr->w[j] * sphRadius * sphRadius * detJ * (1 - cos(phi));
        }
      */

      // new super optimized loop

      real rsint, rcost, sphRadius2 = sphRadius * sphRadius;

      for (j = 0; j < qr->order; j++) {
         t = t2 * qr->x[j];  // t1 = 0 (x-axis by definition)

         // collect these often used expensive terms
         rsint = radius * sin(t);
         rcost = radius * cos(t);

         // remove function calls from setting up p
         p->x = ac->x + rcost * X_l->x + rsint * Y_l->x;
         p->y = ac->y + rcost * X_l->y + rsint * Y_l->y;
         p->z = ac->z + rcost * X_l->z + rsint * Y_l->z;
         //         printf("%3.9f  %3.9f  %3.9f\n", p->x, p->y, p->z);
         // remove function calls from setting up dpdt, skip ->z
         dpdt->x = -rsint * X_l->x + rcost * Y_l->x;
         dpdt->y = -rsint * X_l->y + rcost * Y_l->y;

         // compress everything used to determine cos(phi)
         cosphi = p->z / Vector3D_length(p);

         // compress all expressions into detJ
         detJ = -(p->y * dpdt->x - p->x * dpdt->y) / (p->x*p->x + p->y*p->y);

         *dlp +=  t2 * qr->w[j] * sphRadius2 * detJ * (1.0 - cosphi);
      }
    
   }

   *dlp = - *dlp / (sphRadius * sphRadius);  // messed up directions

   // MALTY BOO

   if (Vector3D_equal(point, gst->centroid))
      *dlp += 4.0 * M_PI;  // this gives limit from outside

   Vector3D_free(pnt);
   Vector3D_free(sphPoint);
   Vector3D_free(normalAt);
   Vector3D_free(X);
   Vector3D_free(Y);
   Vector3D_free(Z);

   Vector3D_free(v1);
   Vector3D_free(v2);
   Vector3D_free(ac);
   Vector3D_free(v2ac);
  
   Vector3D_free(Tvec);
   Vector3D_free(Tmat[0]);
   Vector3D_free(Tmat[1]);
   Vector3D_free(Tmat[2]);

   Vector3D_free(X_l);
   Vector3D_free(Y_l);
   Vector3D_free(Z_l);
   Vector3D_free(p);
   Vector3D_free(pProj);
   Vector3D_free(dpdt);
}

// gonna count on maltman to optimize this: probably best route is to add to flatpanel all the stuff
// that's in the current Panel structure
void greenInt_XinWang(real *Integrals, FlatRefPanel flatpanel, Vector3D localpnt, unsigned int order) {
   unsigned int count;
   real *cur;
   real edgeLength[3], diagLength, area; // will have to change
   // edgelength, npanel to be
   // dynamically allocated when
   // verts is allowed to = 4
   unsigned int i, j, ii, jj, k, next, verts;
   Vector3D X, Y, Z, v1, v3, centroid, npanel[3], nside[3], point;
   real x1, x3, y1, y3, xc, yc, zc;
   real off_plane_tolerance = 1e-5;  // SINGLE PRECISION WARNING
   real ct[3], st[3], deltaxi[3], deltaeta[3];
   real avoid_singular_tolerance = 1e-7; // SINGLE PRECISION WARNING
   real zn, znabs, evalDistance;
   Vector3D rho[3];
   unsigned int OK;
   real xmxv[3], ymyv[3];
   real u[3], U[3], v[3], r[3], fe[3];
   Vector3D cross1, cross2, vtp;
   Matrix cc, cc1, x, y;
   real ru, ru2, phi, psi, psix,psiy, phix,phiy, phixy;
   real s1[3],s2[3],s3[3], c1[3],c2[3],c3[3];
   real Qup, Qlow, Q[3], P[3], intr1, intr3, intxir3, intetar3, intxietar3;
   Matrix sto, cto, lintx, linty, PHI, PSI, PHI1, PSI1;
   real lint;

   X = Vector3D_allocate();
   Y = Vector3D_allocate();
   Z = Vector3D_allocate();
   v1 = Vector3D_allocate();
   v3 = Vector3D_allocate();
   centroid = Vector3D_allocate();
   point = Vector3D_allocate();
   cross1 = Vector3D_allocate();
   cross2 = Vector3D_allocate();
   vtp = Vector3D_allocate();
  
   for (i = 0; i < 3; i++) {
      npanel[i] = Vector3D_allocate();
      nside[i] = Vector3D_allocate();
      rho[i] =Vector3D_allocate();
   }
  
   verts = 3;
   for (i = 0; i < 3; i++) {
      next = i + 1;
      if (next == 3)
         next = 0;
      edgeLength[i] = Vector3D_distance(flatpanel->vertices[next], flatpanel->vertices[i]);
   }

   Vector3D_sub(X, flatpanel->vertices[2], flatpanel->vertices[0]);
   diagLength = Vector3D_length(X);
   Vector3D_sub(Y, flatpanel->vertices[1], flatpanel->vertices[0]);
   Vector3D_cross(Z, X, Y);
   area = .5 * Vector3D_length(Z);
   Vector3D_normalize(X);
   Vector3D_normalize(Z);
   Vector3D_cross(Y, Z, X);

   Vector3D_sub(v1, flatpanel->vertices[1], flatpanel->vertices[0]);
   Vector3D_sub(v3, flatpanel->vertices[2], flatpanel->vertices[0]);

   y1 = Vector3D_dot(v1, Y);
   y3 = Vector3D_dot(v3, Y);
   x1 = Vector3D_dot(v1, X);
   x3 = Vector3D_dot(v3, X);

   yc = (y1 + y3) / 3.0;
   xc = (diagLength + ((x1 * y1 - x3 * y3)/(y1 - y3))) / 3.0; // straight from xin wang.. what if y1=y3?

   Vector3D_copy(centroid, flatpanel->vertices[0]);
   Vector3D_addscaled(centroid, centroid, xc, X);
   Vector3D_addscaled(centroid, centroid, yc, Y);

   // put the corners in the newly defined coordinate system
   for (i = 0; i < 3; i++) {
      Vector3D_sub(npanel[i], flatpanel->vertices[i], centroid);
      Vector3D_transformVecs_inverse(npanel[i], X, Y, Z, npanel[i]);
   }

   // check that panel is in x-y plane
   for (i = 0; i < 3; i++) {
      if (realabs(npanel[i]->z) > off_plane_tolerance * diagLength) { 
         perror("Coordinate transform failure in greenint_xinwang!");
         exit(-1);
      }
   }

   // compute the contribution terms for each edge
   for (i = 0; i < 3; i++) {
      next = i + 1;
      if (next == 3)
         next = 0;

      ct[i] = (npanel[next]->x - npanel[i]->x) / edgeLength[i];
      st[i] = (npanel[next]->y - npanel[i]->y) / edgeLength[i];
      Vector3D_sub(nside[i], npanel[next], npanel[i]);
      deltaxi[i] = nside[i]->x;
      deltaeta[i] = nside[i]->y;
   }

   // done with the panel setup, now loop through the evaluation points
   // order == polynomial order

   // for loop starts here if we have multiple eval points
   Vector3D_sub(point, localpnt, centroid);
   Vector3D_transformVecs_inverse(point, X, Y, Z, point);

   if (realabs(point->z) < avoid_singular_tolerance)
      point->z = avoid_singular_tolerance;  // straight from xin wang, i dunno
  
   zn = point->z;
   znabs = realabs(zn);
   evalDistance = Vector3D_length(point);
   point->z = 0; //little cleverness
   for (i = 0; i < 3; i++) {
      Vector3D_sub(rho[i], point, npanel[i]);
   }
   point->z = zn; // restore point to original state

   // once per vertex computation
   OK = 1;
   for (i = 0; i < 3; i++) {
      next = i + 1;
      if (next == 3)
         next = 0;

      u[i] = Vector3D_dot(rho[i],nside[i]) / Vector3D_length(nside[i]);
      Vector3D_copy(cross1, nside[i]);
      Vector3D_scale(cross1, 1.0/Vector3D_length(nside[i]));
      Vector3D_copy(cross2, nside[i]);
      Vector3D_scale(cross2, u[i]/Vector3D_length(nside[i]));
      Vector3D_sub(cross2, cross2, rho[i]);
      Vector3D_cross(vtp, cross1, cross2);
      v[i] = vtp->z;
      U[i] = Vector3D_dot(rho[next], nside[i]) / Vector3D_length(nside[i]);
      xc = point->x - npanel[i]->x;
      yc = point->y - npanel[i]->y;
      zc = point->z - npanel[i]->z;
      xmxv[i] = xc;
      ymyv[i] = yc;
      fe[i] = xc*xc+zc*zc;
      r[i] = sqrt(yc*yc+fe[i]);

      if (r[i] < (1.005 * znabs)) 
         OK = 0;

   }

   // calculate part of line integrals in newman paper 5.2, 5.3, 5.4;
   // fill CC matrix: binomial coefficients
   cc = Matrix_allocate(order + 1, order + 1);
   cc1 = Matrix_allocate(order + 1, order + 1);
   for (i = 0; i < order + 1; i++) {
      for (j = 0; j <= i; j++) {
         cc[i][j] = 1.0;
      }
   } // end "cc = tril(ones(order+1))" line

   for (i = 1; i < order + 1; i++) {
      cc[i][i] = PetscPowRealInt(-1, i+2); // not i+1 because we're using 0 based indexing
      for (j = 1; j < i; j++) { // not i-1 because our test is < as opposed to matlabs <=
         cc[i][j] = cc[j][j] * (realabs(cc[i-1][j-1]) + realabs(cc[i-1][j]));
      }
   } 
  
   for (i = 0; i < order + 1; i++) {
      for (j = 0; j < order + 1; j++) {
         cc1[i][j] = realabs(cc[i][j]);
      }
   }

   x = Matrix_allocate(verts, order + 1);
   y = Matrix_allocate(verts, order + 1);  // verts = 3 here, always (for now... pain later maybe)
  
   for (i = 0; i < verts; i++) {
      next = i + 1;
      if (next == verts)
         next = 0;

      ru = r[i] - u[i];
      ru2 = r[i]*r[i] - u[i]*u[i];
      if (fabs(r[next]-U[i]) < 1e-6)
         x[i][0] = 1e-8;
      else
         x[i][0] = log((r[next]-U[i]) / ru);
      x[i][1] = r[next] - r[i];
      for (j = 3; j <= order + 1; j++) {
         x[i][j-1] = 1.0 / (j-1) * (PetscPowRealInt(-U[i],j-2) * r[next] - PetscPowRealInt(-u[i],j-2)*r[i]) - ((real)j-2)/(j-1) * ru2 * x[i][j-3];
      }
      for (j = 1; j <= order + 1; j++) {
         for (k = 1; k <= j; k++) {
            y[i][j-1] = y[i][j-1] + cc1[j-1][k-1]*x[i][k-1]*PetscPowRealInt(u[i],j-k);
         }
      }
   }

   Matrix_copy(x, y, verts, order + 1);
  

   phi = 0.0;
   psi = 0.0;
   psix = 0.0;
   psiy = 0.0;
   phix = 0.0;
   phiy = 0.0;
   phixy = 0.0;
  
   for (i = 0; i < verts; i++) {
      next = i + 1;
      if (next == verts)
         next = 0;

      s1[i] = deltaeta[i] *(xmxv[i]*xmxv[i]+point->z*point->z) - deltaxi[i]*xmxv[i]*ymyv[i];
      c1[i] = r[i]* point->z * deltaxi[i];
      s2[i] = deltaeta[i] *(xmxv[next]*xmxv[next]+point->z*point->z) - deltaxi[i]*xmxv[next]*ymyv[next];
      c2[i] = r[next] * point->z * deltaxi[i];
      s3[i] = s1[i]*c2[i]-s2[i]*c1[i];
      c3[i] = c1[i]*c2[i]+s1[i]*s2[i];
      phi = phi + atan2(s3[i],c3[i]);
      Qup = r[i]+r[next]+edgeLength[i];
      Qlow = r[i]+r[next]-edgeLength[i];
      Q[i] = log(Qup/Qlow);
      P[i] = .5 * (u[i]*r[i]-U[i]*r[next]+(r[i]*r[i]-u[i]*u[i])*Q[i]);
      psi = psi+(xmxv[i]*st[i]-ymyv[i]*ct[i]) * Q[i];
      psix = psix - P[i]*st[i];
      psiy = psiy + P[i]*ct[i];
      phix = phix + point->z * Q[i] * st[i];
      phiy = phiy - point->z * Q[i] * ct[i];
      phixy = phixy + point->z * ct[i] * (v[i]*Q[i]*st[i]- ct[i]*(r[next]-r[i]));
   }

   psi = psi - point->z * phi;
   psix = psix + point->x * psi;
   psiy = psiy + point->y * psi;
   phix = phix + point->x * phi;
   phiy = phiy + point->y * phi;
   phixy = phixy + point->x * phiy + point->y * phix - point->x * point->y * phi;

   // put into integration

   intr1 = psi;
   intr3 = phi / point->z;
   intxir3 = phix / point->z;
   intetar3 = phiy / point->z;
   intxietar3 = phixy/point->z;

   // use newman's recurrence equation to get more high order phis and psis
   sto = Matrix_allocate(verts, order + 2);
   cto = Matrix_allocate(verts, order + 2);
   for (i = 1; i <= verts; i++) {
      for (j = 1; j <= order + 2; j++) {
         sto[i-1][j-1] = PetscPowRealInt(st[i-1], j-1);
         cto[i-1][j-1] = PetscPowRealInt(ct[i-1], j-1);
      }
   }

   linty = Matrix_allocate(order+1,order+1);
   lintx = Matrix_allocate(order+1,order+1);

   // no way is it worth screwing with indexing in situations like this
   for (ii = 0; ii <= order; ii++) {
      for (jj = 0; jj <= order-ii; jj++) {
         for (i = 0; i <= ii; i++) {
            for (j = 0; j <= jj; j++) {
               for (k = 1; k <= verts; k++) {
                  lint = x[k-1][ii-i+jj-j] * cto[k-1][ii-i] * sto[k-1][jj-j]* PetscPowRealInt((point->x - npanel[k-1]->x),i)
                     *PetscPowRealInt(point->y - npanel[k-1]->y,j) * cc[ii][i]*cc[jj][j];
                  linty[ii][jj] = linty[ii][jj] - lint * st[k-1];
                  lintx[ii][jj] = lintx[ii][jj] + lint * ct[k-1];
               }
            }
         }
      }
   }

   PHI = Matrix_allocate(order + 1, order + 1);
   PSI = Matrix_allocate(order + 1, order + 1);
   PHI[0][0] = intr3;
   PHI[1][0] = intxir3 - point->x * intr3;
   PHI[0][1] = intetar3 - point->y * intr3;
   PHI[1][1] = intxietar3 - point->x*intetar3-point->y*intxir3 + point->x *point->y * intr3;

   PSI[0][0] = psi;
   PSI[1][0] = psix - point->x*psi;
   PSI[0][1] = psiy - point->y*psi;
   if (order > 2)
      PSI[1][1] = 1.0 / 3.0 * (-point->z*point->z * PHI[1][1] + linty[2][1]+lintx[1][2]);

   for (i = 0; i <= order; i++) {
      for (j = 0; j <= order - i; j++) {
         if ((i>1 || j>1) && (i >= j) && (i+j < order))  {
            PHI[i][j] = (i-1)* PSI[i-2][j] - linty[i-1][j];
            PSI[i][j] = 1.0/((real)i+j+1) *(-point->z*point->z*PHI[i][j]+linty[i+1][j]+lintx[i][j+1]);
         } else if ((i>1 || j>1) && (i >= j) && (i+j == order)) {
            PHI[i][j] = (i-1)*PSI[i-2][j] - linty[i-1][j];
         } else if ((i>1 || j>1) && (j > i) && (i + j < order)) {
            PHI[i][j] = (j-1)*PSI[i][j-2]-lintx[i][j-1];
            PSI[i][j] = 1.0/((real)i+j+1) *(-point->z*point->z*PHI[i][j]+linty[i+1][j]+lintx[i][j+1]);
         } else if ((i>1 || j>1) && (j > i) && (i+j == order)) {
            PHI[i][j] = (j-1)*PSI[i][j-2]- lintx[i][j-1];
         }
      }
   }

   PHI1 = Matrix_allocate(order+1, order+1);
   PSI1 = Matrix_allocate(order+1, order+1);

   for (i = 0; i <= order; i++) {
      for (j = 0; j <= order - i; j++) {
         if (i+j != order) {
            for (ii = 0; ii <= i; ii++) {
               real xtimii = PetscPowRealInt(point->x,i-ii);
               for (jj = 0; jj <= j; jj++) {
                  real ytjmjj = PetscPowRealInt(point->y,j-jj);
                  PHI1[i][j] = PHI1[i][j]+cc1[i][ii]*cc1[j][jj]*PHI[ii][jj]*xtimii*ytjmjj;
                  PSI1[i][j] = PSI1[i][j]+cc1[i][ii]*cc1[j][jj]*PSI[ii][jj]*xtimii*ytjmjj;
               }
            }
         } else {
            for (ii = 0; ii <= i; ii++) {
               real xtimii = PetscPowRealInt(point->x,i-ii);
               for (jj = 0; jj <= j; jj++) {
                  PHI1[i][j] = PHI1[i][j]+cc1[i][ii]*cc1[j][jj]*PHI[ii][jj]*xtimii*PetscPowRealInt(point->y,j-jj);
               }
            }
         }
      }
   }

   // skipping the very end of greenint.m because it's not what we care about!! (yay)
  
   // for loop ends here if we have multiple eval points
   count = 0;
   for (i = 0; i < order; i++) {
      cur = &(PSI1[0][i]); // sets up cur as a pointer to the beginning of
      // the ith column, given how the matrix is set
      // up in memory
      Integrals[count] = cur[0];
      count = count + 1;
      for (j = 0; j < i; j++) {
         cur += order;
         Integrals[count] = cur[0];
         count = count + 1;
      }
   }

   Vector3D_free(X);
   Vector3D_free(Y);
   Vector3D_free(Z);
   Vector3D_free(v1);
   Vector3D_free(v3);
   Vector3D_free(centroid);
   for (i = 0; i < 3; i++) {
      Vector3D_free(npanel[i]);
      Vector3D_free(nside[i]);
      Vector3D_free(rho[i]);
   }
   Vector3D_free(point);
   Vector3D_free(cross1);
   Vector3D_free(cross2);
   Vector3D_free(vtp);
   Matrix_free(cc);
   Matrix_free(cc1);
   Matrix_free(x);
   Matrix_free(y);
   Matrix_free(sto);
   Matrix_free(cto);
   Matrix_free(lintx);
   Matrix_free(linty);
   Matrix_free(PHI);
   Matrix_free(PSI);

   // copy things out of PHI1 and into real *integrals

   Matrix_free(PHI1);
   Matrix_free(PSI1);

  
}

