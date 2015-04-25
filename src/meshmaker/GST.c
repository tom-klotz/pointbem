#include "meshmaker.h"
#include "petscblaslapack.h"

void cartToSph(real* rho, real* theta, real* phi, Vector3D pnt) {
   *rho = Vector3D_length(pnt);
   *theta = atan2(pnt->y, pnt->x);
   *phi = 0;
   if (((*rho) * (*rho)) > 0.0)
      *phi = acos(pnt->z/ *rho);
}

void sphToCart(Vector3D pnt, real rho, real theta, real phi) {
   pnt->x = rho * cos(theta) * sin(phi);
   pnt->y = rho * sin(theta) * sin(phi);
   pnt->z = rho * cos(phi);
}

void GST_readfile(unsigned int *numGSTpanels, GST** gstList, FILE* file, unsigned int ignorecav) {
   char line[128];
   unsigned int i = 0;
   int caporpit = 0, cavity = 0, junk;
   unsigned int numLinesPerGST = 9;

   Vector3D sphcenter = Vector3D_allocate();
   real radius;
   Vector3D v1 = Vector3D_allocate();
   Vector3D v2 = Vector3D_allocate();
   Vector3D v3 = Vector3D_allocate();
   Vector3D ac1 = Vector3D_allocate();
   Vector3D ac2 = Vector3D_allocate();
   Vector3D ac3 = Vector3D_allocate();

   while(fgets(line, 128, file))
      i++;
   rewind(file);

   *numGSTpanels = i / numLinesPerGST;

   *gstList = (GST*)calloc(*numGSTpanels, sizeof(GST));
  
   for (i = 0; i < *numGSTpanels; i++) {
      // line 1: sphcenter
      fgets(line, 128, file);
#ifdef REAL_IS_DOUBLE
      sscanf(line, "%lf %lf %lf", &sphcenter->x, &sphcenter->y, &sphcenter->z);
#else
#ifdef REAL_IS_FLOAT
      sscanf(line, "%f %f %f", &sphcenter->x, &sphcenter->y, &sphcenter->z);
#endif
#endif
      // line 2: sphradius
      fgets(line, 128, file);
#ifdef REAL_IS_DOUBLE
      sscanf(line, "%lf %d", &radius, &caporpit);
#else
#ifdef REAL_IS_FLOAT
      sscanf(line, "%f %d", &radius, &caporpit);
#endif
#endif
      // line 3: v1
      fgets(line, 128, file);
#ifdef REAL_IS_DOUBLE
      sscanf(line, "%lf %lf %lf", &v1->x, &v1->y, &v1->z);
#else
#ifdef REAL_IS_FLOAT
      sscanf(line, "%f %f %f", &v1->x, &v1->y, &v1->z);
#endif
#endif
      // line 4: v2
      fgets(line, 128, file);
#ifdef REAL_IS_DOUBLE
      sscanf(line, "%lf %lf %lf", &v2->x, &v2->y, &v2->z);
#else
#ifdef REAL_IS_FLOAT
      sscanf(line, "%f %f %f", &v2->x, &v2->y, &v2->z);
#endif
#endif
      // line 5: v3
      fgets(line, 128, file);
#ifdef REAL_IS_DOUBLE
      sscanf(line, "%lf %lf %lf", &v3->x, &v3->y, &v3->z);
#else
#ifdef REAL_IS_FLOAT
      sscanf(line, "%f %f %f", &v3->x, &v3->y, &v3->z);
#endif
#endif
      // line 6: ac1
      fgets(line, 128, file);
#ifdef REAL_IS_DOUBLE
      sscanf(line, "%lf %lf %lf", &ac1->x, &ac1->y, &ac1->z);
#else
#ifdef REAL_IS_FLOAT
      sscanf(line, "%f %f %f", &ac1->x, &ac1->y, &ac1->z);
#endif
#endif
      // line 7: ac2
      fgets(line, 128, file);
#ifdef REAL_IS_DOUBLE
      sscanf(line, "%lf %lf %lf", &ac2->x, &ac2->y, &ac2->z);
#else
#ifdef REAL_IS_FLOAT
      sscanf(line, "%f %f %f", &ac2->x, &ac2->y, &ac2->z);
#endif
#endif
      // line 8: ac3
      fgets(line, 128, file);
#ifdef REAL_IS_DOUBLE
      sscanf(line, "%lf %lf %lf", &ac3->x, &ac3->y, &ac3->z);
#else
#ifdef REAL_IS_FLOAT
      sscanf(line, "%f %f %f", &ac3->x, &ac3->y, &ac3->z);
#endif
#endif

      // line 9: arc radii, not needed anymore, but need cavity
      fgets(line, 128, file);
      sscanf(line, "%d %d %d", &junk, &cavity, &junk);

      if (ignorecav)
         cavity = 0;

      // actually initialize the GST
      if (cavity)
         (*gstList)[i] = GST_allocate(sphcenter, radius, caporpit, cavity, v2, v1, v3, ac1, ac3, ac2);
      else
         (*gstList)[i] = GST_allocate(sphcenter, radius, caporpit, cavity, v1, v2, v3, ac1, ac2, ac3);
   }

   Vector3D_free(sphcenter);
   Vector3D_free(v1);
   Vector3D_free(v2);
   Vector3D_free(v3);
   Vector3D_free(ac1);
   Vector3D_free(ac2);
   Vector3D_free(ac3);  
}

GST GST_allocate(Vector3D center, real radius, int caporpit, int cavity,
                 Vector3D v1, Vector3D v2, Vector3D v3,
                 Vector3D ac1, Vector3D ac2, Vector3D ac3) {
   GST gst = (GST)calloc(1, sizeof(struct _GST));
   Vector3D v = Vector3D_allocate();
   real junkreal;
   unsigned int i;
   
   real nonmajorarc_tolerance = 1e-6; // SINGLE PRECISION ALERT
   gst->pit = caporpit;
   gst->cavity = cavity;
   gst->center = Vector3D_allocate();
   Vector3D_copy(gst->center, center);
   gst->radius = radius;
   for (i = 0; i < 3; i++) {
      gst->vertices[i] = Vector3D_allocate();
      gst->arccenters[i] = Vector3D_allocate();
      gst->Tmat[i] = Vector3D_allocate();
   }
   gst->Tvec = Vector3D_allocate();
   gst->centroid = Vector3D_allocate();
   gst->normalAtCentroid = Vector3D_allocate();
   Vector3D_copy(gst->vertices[0], v1);
   Vector3D_copy(gst->vertices[1], v2);
   Vector3D_copy(gst->vertices[2], v3);
   Vector3D_copy(gst->arccenters[0], ac1);
   Vector3D_copy(gst->arccenters[1], ac2);
   Vector3D_copy(gst->arccenters[2], ac3);
   
   for (i = 0; i < 3; i++) {
      findSpherePoint(v, gst->center, gst->radius, gst->vertices[i]);
      Vector3D_copy(gst->vertices[i], v);
   }
   
   // i'm skipping the permute in fixGST because we should be able to
   // do all edge conic corrections now
   
   // CUT by jpb 8/28/05: 
   // get panel "centroid"
   /*  Vector3D_copy(v, gst->vertices[0]); */
   /*   Vector3D_addscaled(v, v, 1.0, gst->vertices[1]); */
   /*   Vector3D_addscaled(v, v, 1.0, gst->vertices[2]); */
   /*   Vector3D_scale(v, (real)(1.0 / 3.0)); */
   /*   findSpherePoint(gst->centroid, gst->center, gst->radius, v); */
   
   /*   // get normal at "centroid" */
   /*   Vector3D_sub(gst->normalAtCentroid, gst->centroid, gst->center); */
   /*   Vector3D_normalize(gst->normalAtCentroid); */
   /*   if (gst->pit > 0) */
   /*     Vector3D_scale(gst->normalAtCentroid, -1.0); */
   // end CUT by jpb
   
   // the centroid of the unit triangle is (1/3, 1/3)
   GST_getCentroid(gst);
   // added back in by maltman
   if (gst->pit > 0)
      Vector3D_scale(gst->normalAtCentroid, -1.0);
   if (gst->cavity > 0)
      Vector3D_scale(gst->normalAtCentroid, -1.0);

   // find the reference panel and transforms (parallels
   // find_GST_flat_transform.m)
   initialize_GST_FlatRefPanel(gst);

   for (i = 0; i < 3; i++) {
      if (nonmajorarc_tolerance <
          Vector3D_distance(gst->arccenters[i], gst->center)) {
         // arc is nonmajor
         figure_out_correction_sign(gst, i); // sets gst->factor[i]
         conic_section_identify(gst, i); //allocs and sets gst->conic[i]
      } else {
         // arc is major
         gst->factor[i] = 0;
         gst->conic[i] = NULL;
      }
   }
   Vector3D_free(v);
   
   // this is HORRIBLE we're fixing it to 16 point quad rule 
   gst->numdirectquadpoints = 16;
   gst->directquadpoints = (Vector3D *)calloc(16, sizeof(Vector3D));
   gst->directquadnormals = (Vector3D *)calloc(16, sizeof(Vector3D));
   gst->directquadweights = (real *)calloc(16, sizeof(real));
   {
      real ksiTab[16];
      real etaTab[16];
      real weight[16];
      Vector3D dTmat[3];
      Vector3D dTvec;
      Curve edges[3];
      dTmat[0] = Vector3D_allocate();
      dTmat[1] = Vector3D_allocate();
      dTmat[2] = Vector3D_allocate();
      dTvec = Vector3D_allocate();
      
      GST_getStandardPosition(gst, edges, dTmat, dTvec);
      gen_Stroud_Rule(ksiTab, etaTab, weight);
      for (i = 0; i < 16; i++) {
         gst->directquadpoints[i] = Vector3D_allocate();
         gst->directquadnormals[i] = Vector3D_allocate();
         GST_get_direct_quadPoint(gst->directquadpoints[i],      // x,y,z
                                  gst->directquadnormals[i],     // nx,ny,nz
                                  &(gst->directquadweights[i]),  // detJ
                                  ksiTab[i], etaTab[i], gst, edges,
                                  dTmat, dTvec);
         if (gst->pit > 0)
            Vector3D_scale(gst->directquadnormals[i], -1.0);
         if (gst->cavity > 0)
            Vector3D_scale(gst->directquadnormals[i], -1.0);
         gst->directquadweights[i] *= weight[i];
      }
      Curve_free(edges[0]);
      Curve_free(edges[1]);
      Curve_free(edges[2]);
      Vector3D_free(dTmat[0]);
      Vector3D_free(dTmat[1]);
      Vector3D_free(dTmat[2]);
      Vector3D_free(dTvec);
   }
   
   return gst;
}

void GST_free(GST gst) {
   unsigned int i;
   for (i = 0; i < 3; i++) {
      Vector3D_free(gst->vertices[i]);
      Vector3D_free(gst->arccenters[i]);
      Vector3D_free(gst->Tmat[i]);
      if (gst->factor[i] != 0)
         Conic_free(gst->conic[i]);
   }
   
   free(gst->directquadweights);
   for (i = 0; i < 16; i++) {
      Vector3D_free(gst->directquadpoints[i]);
      Vector3D_free(gst->directquadnormals[i]);
   }
   free(gst->directquadpoints);
   free(gst->directquadnormals);
   
   Vector3D_free(gst->center);
   Vector3D_free(gst->Tvec);
   Vector3D_free(gst->centroid);
   Vector3D_free(gst->normalAtCentroid);
   FlatRefPanel_free(gst->flatpanel);
   
   Matrix_free(gst->vandermondeInverse);
   
   free(gst);
}

void findSpherePoint(Vector3D sphvert, Vector3D center, real radius, Vector3D oldvert) {
   Vector3D tmp = Vector3D_allocate();
   
   Vector3D_sub(tmp, oldvert, center);
   Vector3D_normalize(tmp);
   
   Vector3D_addscaled(sphvert, center, radius, tmp);
   
   Vector3D_free(tmp);
}

/* Why does this function use spherical to cartesian conversions?? See above. */

/*
  void findSpherePoint(Vector3D sphvert, Vector3D center,
  real radius, Vector3D oldvert) {
  real rho, theta, phi;
  Vector3D tmp = Vector3D_allocate();
  
  Vector3D_sub(tmp, oldvert, center);
  cartToSph(&rho, &theta, &phi, tmp);
  rho = radius;
  sphToCart(tmp, rho, theta, phi);
  
  Vector3D_add(sphvert, tmp, center);
  
  Vector3D_free(tmp);
  }
*/

void initialize_GST_FlatRefPanel(GST gst) {
   real rhs, t, maxEdgelength;
   real *ksi;
   real *eta;
   real *weights;
   Vector3D fp_centroid = Vector3D_allocate();
   Vector3D fp_normal = Vector3D_allocate();
   Vector3D X = Vector3D_allocate();
   Vector3D Y = Vector3D_allocate();
   Vector3D Z = Vector3D_allocate();
   Vector3D r0 = Vector3D_allocate();
   Vector3D r1 = Vector3D_allocate();
   Vector3D delta = Vector3D_allocate();
   
   unsigned int numquadpoints, i, next;
   Vector3D fp_vertices[3];
   for (i = 0; i < 3; i++)
      fp_vertices[i] = Vector3D_allocate();
   
   // the below sets gst->tmat = I (since each vec is zero beforehand)
   gst->Tmat[0]->x = 1.0;
   gst->Tmat[1]->y = 1.0;
   gst->Tmat[2]->z = 1.0;
   
   // for LJ calculations, use the other block...
   if (1) {
      Vector3D_copy(fp_centroid, gst->vertices[0]);
      Vector3D_add(fp_centroid, fp_centroid, gst->vertices[1]);
      Vector3D_add(fp_centroid, fp_centroid, gst->vertices[2]);
      Vector3D_scale(fp_centroid, (real) (1.0/3.0));
      findSpherePoint(fp_centroid, gst->center, gst->radius, fp_centroid);
      Vector3D_sub(fp_normal, fp_centroid, gst->center);
      Vector3D_normalize(fp_normal);
      rhs = Vector3D_dot(fp_normal, fp_centroid);
      
      for (i = 0; i < 3; i++) {
         Vector3D_copy(r0, gst->center);
         Vector3D_copy(r1, gst->vertices[i]);
         Vector3D_sub(delta, r1, r0);
         
         t = (rhs - Vector3D_dot(r0, fp_normal)) / Vector3D_dot(fp_normal, delta);
         Vector3D_copy(fp_vertices[i], r0);
         Vector3D_addscaled(fp_vertices[i], fp_vertices[i], t, delta);
      }
   } else {
   }
   
   // translate everything by -fp_centroid
   Vector3D_copy(gst->Tvec, fp_centroid);  // tvec should be zero to
   // start with anyway, before
   // we do this
   Vector3D_scale(gst->Tvec, -1.0);
   for (i = 0; i < 3; i++)
      Vector3D_sub(fp_vertices[i], fp_vertices[i], fp_centroid);
   Vector3D_sub(fp_centroid,fp_centroid,fp_centroid);
   rhs = 0;
   
   // now rotate so normal is inline with z-axis
   Vector3D_sub(X, fp_vertices[1], fp_vertices[0]);
   Vector3D_normalize(X);
   Vector3D_sub(Y, fp_vertices[2], fp_vertices[0]);
   Vector3D_addscaled(Y, Y, -Vector3D_dot(Y,X), X);
   Vector3D_normalize(Y);
   Vector3D_cross(Z, X, Y);
   for (i = 0; i < 3; i++) {
      Vector3D_transformVecs_inverse(gst->Tmat[i], X, Y, Z, gst->Tmat[i]);
      Vector3D_transformVecs_inverse(fp_vertices[i], X, Y, Z, fp_vertices[i]);
   }
   Vector3D_transformVecs_inverse(gst->Tvec, X, Y, Z, gst->Tvec);
   Vector3D_transformVecs_inverse(fp_centroid, X, Y, Z, fp_centroid);
   Vector3D_transformVecs_inverse(fp_normal, X, Y, Z, fp_normal);
   
   
   // now generate quadrature points
   generateQuadraturePoints(fp_vertices, &numquadpoints, &ksi, &eta, &weights);
   
   maxEdgelength = 0;
   for (i = 0; i < 3; i++) {
      next = i + 1;
      if (next==3)
         next = 0;
      if (maxEdgelength <
          Vector3D_distance(fp_vertices[i], fp_vertices[next]))
         maxEdgelength = Vector3D_distance(fp_vertices[i], fp_vertices[next]);
   }
   
   
   gst->flatpanel = FlatRefPanel_allocate(fp_vertices[0], fp_vertices[1],
                                          fp_vertices[2], fp_centroid,
                                          fp_normal, rhs, numquadpoints,
                                          ksi, eta, weights, maxEdgelength);
   
   free(ksi);
   free(eta);
   free(weights);
   Vector3D_free(fp_centroid);
   Vector3D_free(fp_normal);
   Vector3D_free(fp_vertices[0]);
   Vector3D_free(fp_vertices[1]);
   Vector3D_free(fp_vertices[2]);
   Vector3D_free(r0);
   Vector3D_free(r1);
   Vector3D_free(delta);
   Vector3D_free(X);
   Vector3D_free(Y);
   Vector3D_free(Z);
}

FlatRefPanel FlatRefPanel_allocate(Vector3D v1, Vector3D v2, Vector3D v3,
                                   Vector3D centroid, Vector3D normal,
                                   real rhs, unsigned int numquadpoints,
                                   real *ksi, real *eta, real *weights,
                                   real maxEdgelength) {
   FlatRefPanel fp = (FlatRefPanel)calloc(1, sizeof(struct _FlatRefPanel));
   fp->vertices[0] = Vector3D_allocate();
   Vector3D_copy(fp->vertices[0], v1);
   fp->vertices[1] = Vector3D_allocate();
   Vector3D_copy(fp->vertices[1], v2);
   fp->vertices[2] = Vector3D_allocate();
   Vector3D_copy(fp->vertices[2], v3);
   
   fp->centroid = Vector3D_allocate();
   Vector3D_copy(fp->centroid, centroid);
   fp->normal = Vector3D_allocate();
   Vector3D_copy(fp->normal, normal);
   
   fp->rhs = rhs;
   fp->numquadpoints = numquadpoints;
   fp->ksi = (real *)calloc(numquadpoints, sizeof(real));
   memcpy((void*)fp->ksi,(void*)ksi, numquadpoints * sizeof(real));
   fp->eta = (real *)calloc(numquadpoints, sizeof(real));
   memcpy((void*)fp->eta,(void*)eta, numquadpoints * sizeof(real));
   fp->weights = (real *)calloc(numquadpoints, sizeof(real));
   memcpy((void*)fp->weights,(void*)weights, numquadpoints * sizeof(real));
   
   fp->maxEdgelength = maxEdgelength;
   
   return fp;
}

void FlatRefPanel_free(FlatRefPanel fp) {
   Vector3D_free(fp->centroid);
   Vector3D_free(fp->normal);
   Vector3D_free(fp->vertices[0]);
   Vector3D_free(fp->vertices[1]);
   Vector3D_free(fp->vertices[2]);
   free(fp->ksi);
   free(fp->eta);
   free(fp->weights);
   free(fp);

}

void generateQuadraturePoints(Vector3D* vertices, unsigned int *numquadpoints,
                              real **ksi, real **eta, real **weights) {
   real lengthXside, lengthYside, height, curYlength;
   unsigned int i, j;
   static QuadratureRule qr = NULL;
   Vector3D X, Y,curBase,junk,curPoint;
   X = Vector3D_allocate();
   Y = Vector3D_allocate();
   curBase = Vector3D_allocate();
   junk = Vector3D_allocate();
   curPoint = Vector3D_allocate();

#ifdef OMP
#pragma omp critical
#endif
   if (qr == NULL) {
      qr = QuadratureRule_allocate(4);
   }
   
   Vector3D_sub(X, vertices[0], vertices[1]);
   lengthXside = Vector3D_length(X);
   Vector3D_normalize(X);
   
   Vector3D_sub(Y, vertices[2], vertices[0]);
   lengthYside = Vector3D_length(Y);
   Vector3D_normalize(Y);
   
   Vector3D_sub(junk, vertices[0], vertices[1]);
   Vector3D_addscaled(junk, junk, -Vector3D_dot(junk,Y), Y);
   height = Vector3D_length(junk);
   
   *numquadpoints = qr->order * qr->order;
   *ksi = (real *)calloc(*numquadpoints, sizeof(real));
   *eta = (real *)calloc(*numquadpoints, sizeof(real));
   *weights = (real *)calloc(*numquadpoints, sizeof(real));
   for (i = 0; i < qr->order; i++) {
      Vector3D_addscaled(curBase, vertices[1], qr->x[i] * lengthXside, X);
      curYlength = qr->x[i] * lengthYside;
      for (j = 0; j < qr->order; j++) {
         Vector3D_addscaled(curPoint, curBase, qr->x[j] * curYlength, Y);
         (*ksi)[i*qr->order + j] = curPoint->x;
         (*eta)[i*qr->order + j] = curPoint->y;
         (*weights)[i*qr->order + j] = qr->w[i] * height * qr->w[j] * curYlength;
      }
   }
   
   Vector3D_free(X);
   Vector3D_free(Y);
   Vector3D_free(junk);
   Vector3D_free(curBase);
   Vector3D_free(curPoint);
}

void figure_out_correction_sign(GST gst, unsigned int index) {
   unsigned int this, next;
   Vector3D v1, v2, ac, sc, vs1, vs2, dv1v2, intopanel, theta0, dtheta, stn;
   Vector3D Z;
   v1 = Vector3D_allocate();
   v2 = Vector3D_allocate();
   ac = Vector3D_allocate();
   sc = Vector3D_allocate();
   vs1 = Vector3D_allocate();
   vs2 = Vector3D_allocate();
   dv1v2 = Vector3D_allocate();
   intopanel = Vector3D_allocate();
   theta0 = Vector3D_allocate();
   dtheta = Vector3D_allocate();
   stn = Vector3D_allocate();
   Z = Vector3D_allocate();
   Z->z = 1.0;
   this = index; next = this + 1;
   if (next == 3)
      next = 0;
   
   Vector3D_copy(v1, gst->vertices[this]);
   Vector3D_copy(v2, gst->vertices[next]);
   Vector3D_copy(ac, gst->arccenters[this]);
   Vector3D_copy(sc, gst->center);
   
   Vector3D_transformVecs(v1, gst->Tmat[0],gst->Tmat[1],gst->Tmat[2], v1);
   Vector3D_add(v1, v1, gst->Tvec);
   Vector3D_transformVecs(v2, gst->Tmat[0],gst->Tmat[1],gst->Tmat[2], v2);
   Vector3D_add(v2, v2, gst->Tvec);
   Vector3D_transformVecs(ac, gst->Tmat[0],gst->Tmat[1],gst->Tmat[2], ac);
   Vector3D_add(ac, ac, gst->Tvec);
   Vector3D_transformVecs(sc, gst->Tmat[0],gst->Tmat[1],gst->Tmat[2], sc);
   Vector3D_add(sc, sc, gst->Tvec);
   
   Vector3D_sub(dv1v2, v1, v2);
   Vector3D_cross(intopanel, dv1v2, Z);
   Vector3D_normalize(intopanel);
   
   Vector3D_sub(theta0, v1, ac);
   Vector3D_normalize(theta0);
   Vector3D_sub(dtheta, v2, ac);
   Vector3D_addscaled(dtheta, dtheta, -Vector3D_dot(dtheta, theta0), theta0);
   Vector3D_normalize(dtheta);
   
   Vector3D_sub(vs1, v1, sc);
   Vector3D_sub(vs2, v2, sc);
   Vector3D_cross(stn, vs1, vs2);
   Vector3D_normalize(stn);
   
   if (Vector3D_dot(stn, dtheta) > 0)
      gst->factor[index] = -1;
   else if (Vector3D_dot(stn, dtheta) < 0)
      gst->factor[index] = +1;
   else
      gst->factor[index] = 0;

   if (gst->pit)
      gst->factor[index] = -gst->factor[index];
   if (gst->cavity)
      gst->factor[index] = -gst->factor[index];

   Vector3D_free(v1);
   Vector3D_free(v2);
   Vector3D_free(ac);
   Vector3D_free(sc);
   Vector3D_free(dv1v2);
   Vector3D_free(vs1);
   Vector3D_free(vs2);
   Vector3D_free(intopanel);
   Vector3D_free(theta0);
   Vector3D_free(dtheta);
   Vector3D_free(stn);
   Vector3D_free(Z);
}

void conic_section_identify(GST gst, unsigned int index) {
   Vector3D Tmat[3];
   Vector3D Tvec;
   Vector3D vthisac, vnextac, vprev;
   real min_length_tolerance, transform_tolerance;
   Vector3D transformedGSTcenter;
   Vector3D vthis, vnext, ac;
   unsigned int i, j, this, next, prev;
   real rad, theta, v2_1rotx, v2_1roty, factor;
   Vector3D arcCircleNormal, circlePoint;
   real phi;
   real a, b, theta1Int, theta2Int;
   real V[6], A, B, C, D, E, F, invdet;
   real A2, B2 ,C2, D2, E2, F2;
   real A3, B3, C3, D3, E3, F3;
   real newx, newy, u, v, c, s;
   Vector3D *M, V1, V2, V3;
   
   unsigned int numpoints = 100;
   
   this = index;
   next = index + 1;
   if (index == 0)
      prev = 2;
   else
      prev = index - 1;
   if (next == 3)
      next = 0;
   
   vthis = Vector3D_allocate();
   vnext = Vector3D_allocate();
   ac = Vector3D_allocate();
   transformedGSTcenter = Vector3D_allocate();
   vprev = Vector3D_allocate();
   
   for (i = 0; i < 3; i++) {
      Tmat[i] = Vector3D_allocate();
   }
   vthisac = Vector3D_allocate();
   vnextac = Vector3D_allocate();
   arcCircleNormal = Vector3D_allocate();
   circlePoint = Vector3D_allocate();
   M = (Vector3D *)calloc(numpoints, sizeof(Vector3D));
   for (i = 0; i < numpoints; i++)
      M[i] = Vector3D_allocate();
   Tvec = Vector3D_allocate();
   V1 = Vector3D_allocate();
   V2 = Vector3D_allocate();
   V3 = Vector3D_allocate();
   
   Vector3D_copy(vthis, gst->vertices[this]);
   Vector3D_copy(vnext, gst->vertices[next]);
   Vector3D_copy(vprev, gst->vertices[prev]);
   Vector3D_copy(ac, gst->arccenters[this]);
   Vector3D_copy(transformedGSTcenter, gst->center);
   
   Vector3D_transformVecs(vthis, gst->Tmat[0], gst->Tmat[1], gst->Tmat[2], vthis);
   Vector3D_add(vthis, vthis, gst->Tvec);
   Vector3D_transformVecs(vnext, gst->Tmat[0], gst->Tmat[1], gst->Tmat[2], vnext);
   Vector3D_add(vnext, vnext, gst->Tvec);
   Vector3D_transformVecs(vprev, gst->Tmat[0], gst->Tmat[1], gst->Tmat[2], vprev);
   Vector3D_add(vprev, vprev, gst->Tvec);
   Vector3D_transformVecs(ac, gst->Tmat[0], gst->Tmat[1], gst->Tmat[2], ac);
   Vector3D_add(ac, ac, gst->Tvec);
   Vector3D_transformVecs(transformedGSTcenter, gst->Tmat[0], gst->Tmat[1], gst->Tmat[2], transformedGSTcenter);
   Vector3D_add(transformedGSTcenter, transformedGSTcenter, gst->Tvec);
   
   Tmat[0]->x = 1.0;
   Tmat[1]->y = 1.0;
   Tmat[2]->z = 1.0;  // this sets Tmat = eye(3)
   
   
   
   rad = Vector3D_distance(vthis, ac);
   Vector3D_sub(vthisac, vthis, ac);
   Vector3D_normalize(vthisac); // X axis
   Vector3D_sub(vnextac,  vnext, ac);
   Vector3D_addscaled(vnextac, vnextac, -Vector3D_dot(vnextac,vthisac), vthisac);
   Vector3D_normalize(vnextac); // Y axis
   Vector3D_cross(arcCircleNormal, vthisac, vnextac);
   Vector3D_normalize(arcCircleNormal); // Z axis
   
   for (i = 1; i <= numpoints; i++) {
      theta = i * (real) 2.0 * M_PI / numpoints;
      Vector3D_copy(circlePoint, ac);
      Vector3D_addscaled(circlePoint, circlePoint, rad * cos(theta), vthisac);
      Vector3D_addscaled(circlePoint, circlePoint, rad * sin(theta), vnextac);
      projectToPlane(M[i-1], transformedGSTcenter, circlePoint);
   }
   Vector3D_free(transformedGSTcenter);
   GST_get_conic(V, M, numpoints);
   
   A = V[0]; B = V[1]; C = V[2];
   D = V[3]; E = V[4]; F = V[5];
   invdet = (real)1.0 / (4.0 * A *C  - B* B);
   u = invdet * ((2 * C) * -D + -B * -E);  // center=[2*A B; B 2*C] \ [-D; -E]
   v = invdet * (-B * -D + (2 * A)  *-E);
   
   Tvec->x -= u; // Tvec = Tvec - [center; 0]
   Tvec->y -= v;
   
   // translate conic
   A2 = A;
   B2 = B;
   C2 = C;
   D2 = D + 2 * A * u + B * v;
   E2 = E + 2 * C * v + B * u;
   F2 = F + A * u * u + B * u * v + C * v * v + D * u + E * v;
   
   // rotate conic
   theta  = .5 * atan2(B, A-C);
   c = cos(theta);
   s = sin(theta);
   A3 = A2 * c * c + B2 * c * s + C2 * s * s;
   B3 = -2 * A2 * c * s - B2 * s * s + B2 * c * c + 2 * C2 * s * c;
   C3 = A2 * s * s- B2 *c*s + C2 * c*c;
   D3 = D2 * c + E2*s;
   E3 = -D2 *s + E2*c;
   F3 = F2;
   V1->x = c; V1->y = s; V1->z = 0;
   V2->x = -s; V2->y = c; V2->z = 0;
   V3->x = 0; V3->y = 0; V3->z = 1;
   for (i = 0; i < 3; i++) 
      Vector3D_transformVecs_inverse(Tmat[i], V1, V2, V3, Tmat[i]);
   Vector3D_transformVecs_inverse(Tvec, V1, V2, V3, Tvec);
   
   if (((A3 > 0.0) && (F3 > 0.0)) ||
       ((A3 < 0.0) && (F3 < 0.0)) ) {
      A2 = A3; B2 = B3; C2 = C3; D2 = D3; E2 = E3; F2 = F3;
      
      A3 = C2;
      B3 = -B2;
      C3 = A2;
      D3 = -E2;
      E3 = D2;
      F3 = F2;
      
      V1->x = 0; V1->y = 1; V1->z = 0;
      V2->x = -1; V2->y = 0; V2->z = 0;
      V3->x = 0; V3->y = 0; V3->z = 1;
      for (i = 0; i < 3; i++) 
         Vector3D_transformVecs_inverse(Tmat[i], V1, V2, V3, Tmat[i]);
      Vector3D_transformVecs_inverse(Tvec, V1, V2, V3, Tvec);
   }
   
   a = 1.0 / (-A3 / F3);
   b = 1.0 / (C3 / F3);
   Vector3D_copy(vthis, gst->flatpanel->vertices[this]);
   Vector3D_transformVecs(vthis, Tmat[0], Tmat[1], Tmat[2], vthis);
   Vector3D_add(vthis, vthis, Tvec);
   Vector3D_copy(vnext, gst->flatpanel->vertices[next]);
   Vector3D_transformVecs(vnext, Tmat[0], Tmat[1], Tmat[2], vnext);
   Vector3D_add(vnext, vnext, Tvec);
   theta1Int = atan2(vthis->y, vthis->x);
   newx = cos(theta1Int) * vnext->x + sin(theta1Int) * vnext->y;
   newy = -sin(theta1Int) * vnext->x + cos(theta1Int) * vnext->y;
   theta2Int = theta1Int + atan2(newy, newx);
   
   gst->conic[index] = Conic_allocate(a, b, theta1Int, theta2Int, Tmat, Tvec);
   for (i = 0; i < 3; i++) 
      Vector3D_free(Tmat[i]);
   Vector3D_free(Tvec);
   Vector3D_free(V1);
   Vector3D_free(V2);
   Vector3D_free(V3);
   Vector3D_free(vthisac);
   Vector3D_free(vnextac);
   Vector3D_free(vprev);
   Vector3D_free(vthis);
   Vector3D_free(vnext);
   Vector3D_free(ac);
   Vector3D_free(arcCircleNormal);
   Vector3D_free(circlePoint);
   for (i = 0; i < numpoints; i++)
      Vector3D_free(M[i]);
   free(M);
}

// find the point at which the line defined by start-end intersects the X-Y plane
void projectToPlane(Vector3D onplane, Vector3D start, Vector3D end) {
   real close_to_plane_tolerance = 1e-7;
   Vector3D dv = Vector3D_allocate();
   Vector3D_sub(dv, end, start);
   //if ( (fabs(dv->z) < close_to_plane_tolerance) ||
   //  (start->z > close_to_plane_tolerance)) {
   if (fabs(dv->z) < close_to_plane_tolerance) {
      printf("flat line will never hit x-y plane within reasonable distance %e\n", fabs(dv->z));
      //exit(-1);
   }
   Vector3D_addscaled(onplane, start, -start->z/dv->z, dv);
   Vector3D_free(dv);
}
void GST_get_conic(real* Vret, Vector3D *M, unsigned int numpoints) {
   char JOBU = 'A', JOBVT = 'A';
   numpoints = 6;
   unsigned int m = numpoints, n = 6, LDA = numpoints, LDU = numpoints, LDVT = 6;
   unsigned int mindim = PetscMin(6, numpoints);
   unsigned int LWORK = PetscMax(3*PetscMin(m,n)+PetscMax(m,n),5*PetscMin(m,n));
   int INFO;
   Matrix UT = Matrix_allocate(numpoints, numpoints);
   Matrix S = Matrix_allocate(numpoints, 6);
   Matrix V = Matrix_allocate(6, 6);
   Matrix temp = Matrix_allocate(numpoints, 6);
   real D[mindim], WORK[LWORK], tol;
   
   unsigned int i;
   real x, y;
   Matrix A = Matrix_allocate(numpoints,6);
   for (i = 0; i < numpoints; i++) {
      x = M[i]->x;
      y = M[i]->y;
      A[i][0] = x*x;
      A[i][1] = x*y;
      A[i][2] = y*y;
      A[i][3] = x;
      A[i][4] = y;
      A[i][5] = 1.0;  // actually building column major matrix for SVD
   }
   
#ifdef OMP
#pragma omp critical (lapack)
#endif
   LAPACKgesvd_(&JOBU, &JOBVT, &m, &n, A[0], &LDA, D,
                UT[0], &LDU, V[0], &LDVT, WORK, &LWORK, &INFO);
   
   
   if (INFO) {
      printf("s/dgesvd returned error code %u\n", INFO);
      exit(-2);
   }
   
   
   for (i = 0; i < 6; i++)
      Vret[i] = UT[5][i];
   Matrix_free(V);
   Matrix_free(S);
   Matrix_free(temp);
   Matrix_free(UT);
   Matrix_free(A);
}

Conic Conic_allocate(real asq, real bsq, real theta1Int, real theta2Int,
                     Vector3D *Tmat, Vector3D Tvec) {
   unsigned int i;
   Conic c = (Conic)calloc(1, sizeof(struct _Conic));
   
   c->asquared = asq;
   c->bsquared = bsq;
   c->theta1Int = theta1Int;
   c->theta2Int = theta2Int;
   
   for (i = 0; i < 3; i++) {
      c->Tmat[i] = Vector3D_allocate();
      Vector3D_copy(c->Tmat[i], Tmat[i]);
   }
   c->Tvec = Vector3D_allocate();
   Vector3D_copy(c->Tvec, Tvec);
   
   return c;
}

void Conic_free(Conic c) {
   unsigned int i;
   for (i = 0; i < 3; i++) {
      Vector3D_free(c->Tmat[i]);
   }
   Vector3D_free(c->Tvec);
   free(c);
}

real safeacos(real arg) {
   real phi;
   if ((arg > -1.0) && (arg < 1.0))
      phi = acos(arg);
   else { 
      if (arg <= -1.0)
         phi = M_PI;
      else
         phi = 0.0;
   }
   return phi;
}

void GST_get_direct_quadPoint(Vector3D point, Vector3D normal, real *detJ,
                              real ksi, real eta, GST gst,
                              Curve *edges, Vector3D *Tmat, Vector3D Tvec) {
   real phiStart, phiEnd, phi;
   real z_threshold = 1e-6;
   Vector3D centerCurEtaCircle = Vector3D_allocate();
   Vector3D normalCurEtaCircle = Vector3D_allocate();
   real radiusCurEtaCircle;
   
   Vector3D startPoint = Vector3D_allocate();
   Vector3D endPoint = Vector3D_allocate();
   real alphaStartCurve;
   Vector3D dStart_dAlpha = Vector3D_allocate();
   real alphaEndCurve;
   Vector3D dEnd_dAlpha = Vector3D_allocate();
   
   real thetaStart, thetaEnd, theta;
   Vector3D dR_dTheta = Vector3D_allocate();
   Vector3D dR_dPhi = Vector3D_allocate();
   Vector3D  dStart_dEta = Vector3D_allocate();
   Vector3D dEnd_dEta = Vector3D_allocate();
   real dTheta_dKsi, dPhi_dKsi, dPhi_dEta, dX_dEta;
   real dAlphaStart_dEta, uStart, duStart_dEta, dThetaStart_dEta;
   real dAlphaEnd_dEta, uEnd, duEnd_dEta, dThetaEnd_dEta;
   
   real dTheta_dEta;
   Vector3D dR_dKsi = Vector3D_allocate();
   Vector3D dR_dEta = Vector3D_allocate();
   
   Vector3D Normal2 = Vector3D_allocate();
   
   // okay now start doing work
   phiStart = safeacos(edges[0]->v1->x / gst->radius);
   phiEnd   = safeacos(edges[2]->v1->x / gst->radius);
   phi = phiStart + (phiEnd - phiStart) * eta;
   centerCurEtaCircle->x = gst->radius * cos(phi);
   centerCurEtaCircle->y = 0.0;
   centerCurEtaCircle->z = 0.0;
   normalCurEtaCircle->x = 1.0;
   normalCurEtaCircle->y = 0.0;
   normalCurEtaCircle->z = 0.0;
   radiusCurEtaCircle = gst->radius * sin(phi);
   //   printf("%3.9f  %3.9f  %3.9f\n", phiStart, phiEnd, phi);

   GST_getCircleArcIntersection(centerCurEtaCircle, normalCurEtaCircle,
                                radiusCurEtaCircle, edges[2],
                                startPoint, &alphaStartCurve, dStart_dAlpha);
   GST_getCircleArcIntersection(centerCurEtaCircle, normalCurEtaCircle,
                                radiusCurEtaCircle, edges[1],
                                endPoint, &alphaEndCurve, dEnd_dAlpha);
   
   thetaStart = atan2(startPoint->y, startPoint->z);
   thetaEnd = atan2(endPoint->y, endPoint->z);
   theta = thetaStart + ksi / (1. - eta) * (thetaEnd - thetaStart);
   //   printf("%3.9f  %3.9f  %3.9f\n", thetaStart, thetaEnd, theta);
   GST_sphGSTtoCart(gst->radius, -theta, phi, point, normal, dR_dTheta, dR_dPhi);
   Vector3D_transformVecs(point, Tmat[0], Tmat[1], Tmat[2], point);
   Vector3D_add(point, point, Tvec);
   Vector3D_transformVecs(normal, Tmat[0], Tmat[1], Tmat[2], normal);
   
   
   dTheta_dKsi = (thetaEnd - thetaStart) / (1. - eta);
   dPhi_dKsi = 0;
   dPhi_dEta = (phiEnd - phiStart);
   dX_dEta = dR_dPhi->x * dPhi_dEta;
   
   dAlphaStart_dEta = 1. / (dStart_dAlpha->x / dX_dEta);
   Vector3D_copy(dStart_dEta, dStart_dAlpha);
   Vector3D_scale(dStart_dEta, dAlphaStart_dEta);
   duStart_dEta = (startPoint->z * dStart_dEta->y
                   - startPoint->y * dStart_dEta->z) /
      (startPoint->z * startPoint->z);
   if ( startPoint->z < z_threshold) {  // changed so that startPoint->z doesn't crash the whole thing
      dThetaStart_dEta = 0.0;  //  insert uStart = infty below, you get zero
   } else { 
      uStart = startPoint->y / startPoint->z;
      dThetaStart_dEta = (1. / (1. + uStart * uStart)) * duStart_dEta;
   }
   dAlphaEnd_dEta = 1. / (dEnd_dAlpha->x / dX_dEta);
   Vector3D_copy(dEnd_dEta, dEnd_dAlpha);
   Vector3D_scale(dEnd_dEta, dAlphaEnd_dEta);
   duEnd_dEta = (endPoint->z * dEnd_dEta->y
                 - endPoint->y * dEnd_dEta->z) /
      (endPoint->z * endPoint->z);
   if ( endPoint->z < z_threshold) {  // changed so that startPoint->z doesn't crash the whole thing
      dThetaEnd_dEta = 0.0;  //  insert uStart = infty below, you get zero
   } else { 
      uEnd = endPoint->y / endPoint->z;
      dThetaEnd_dEta = (1. / (1. + uEnd * uEnd)) * duEnd_dEta;
   }

   dTheta_dEta = dThetaStart_dEta
      + (ksi / (1. - eta)) * (dThetaEnd_dEta - dThetaStart_dEta)
      + (ksi / ((1. - eta)  * (1. - eta))) * (thetaEnd - thetaStart);

   Vector3D_scale(dR_dEta, 0.0);
   Vector3D_scale(dR_dKsi, 0.0);
   Vector3D_addscaled(dR_dEta, dR_dEta, dTheta_dEta, dR_dTheta);
   Vector3D_addscaled(dR_dEta, dR_dEta, dPhi_dEta, dR_dPhi);
   Vector3D_addscaled(dR_dKsi, dR_dKsi, dTheta_dKsi, dR_dTheta);
   Vector3D_addscaled(dR_dKsi, dR_dKsi, dPhi_dKsi, dR_dPhi);
   
   Vector3D_cross(Normal2, dR_dKsi, dR_dEta);
   *detJ = Vector3D_length(Normal2);
   
   // clean up, go home
   Vector3D_free(dStart_dEta);
   Vector3D_free(dEnd_dEta);
   Vector3D_free(Normal2);
   Vector3D_free(dR_dEta);
   Vector3D_free(dR_dKsi);
   Vector3D_free(dR_dPhi);
   Vector3D_free(dR_dTheta);
   Vector3D_free(startPoint);
   Vector3D_free(endPoint);
   Vector3D_free(dStart_dAlpha);
   Vector3D_free(dEnd_dAlpha);
   Vector3D_free(centerCurEtaCircle);
   Vector3D_free(normalCurEtaCircle);
}

void GST_getStandardPosition(GST gst, Curve *edges, Vector3D *Tmat, Vector3D Tvec) {
   real edgeLengths[3]; // copy vertices, arccenters to local vars
   real maxEdgelength = -1.0;
   Curve curve, curEdge;
   Vector3D vertices[3];
   Vector3D arccenters[3];
   
   unsigned int i, next, maxIndex;
   for (i = 0; i < 3; i++) {
      vertices[i] = Vector3D_allocate();
      arccenters[i] = Vector3D_allocate();
      Vector3D_copy(vertices[i], gst->vertices[i]);
      Vector3D_copy(arccenters[i], gst->arccenters[i]);
   }
   
   maxIndex = 0;
   for (i = 0; i < 3; i++) {
      next = i + 1;
      if (next > 2)
         next = 0;
      curEdge = Curve_allocate(arccenters[i], vertices[i], vertices[next]);
      edgeLengths[i] = curEdge->radius * curEdge->endTheta;
      if (edgeLengths[i] > maxEdgelength) {
         maxIndex = i;
         maxEdgelength = edgeLengths[i];
      }
      Curve_free(curEdge);
   }
   
  
   if (maxIndex == 0) {
      // no need to permute
   } else if (maxIndex == 1) {
      // v2 -> v1 new
      // v3 -> v2 new
      // v1 -> v3 new
      PermuteVectors(vertices[0], vertices[1], vertices[2]);
      PermuteVectors(arccenters[0], arccenters[1], arccenters[2]);
   } else { // permute twice.
      PermuteVectors(vertices[0], vertices[1], vertices[2]);
      PermuteVectors(arccenters[0], arccenters[1], arccenters[2]);
      PermuteVectors(vertices[0], vertices[1], vertices[2]);
      PermuteVectors(arccenters[0], arccenters[1], arccenters[2]);
   }
   
   Vector3D_copy(Tvec, gst->center);
   for (i = 0; i < 3; i++) {
      Vector3D_sub(vertices[i], vertices[i], gst->center);
      Vector3D_sub(arccenters[i], arccenters[i], gst->center);
   }
   
   curve = Curve_allocate(arccenters[0], vertices[0], vertices[1]);
   Curve_getParamPoint(Tmat[2], NULL, curve, 0.5);
   Vector3D_cross(Tmat[0], curve->X, curve->Y);
   Vector3D_scale(Tmat[0], -1.0);
   Vector3D_sub(Tmat[2], Tmat[2], arccenters[0]);
   Vector3D_scale(Tmat[2], 1.0 / Vector3D_length(Tmat[2]));
   Vector3D_cross(Tmat[1], Tmat[2], Tmat[0]);
   
   /*   Vector3D_transformVecs_inverse(Tvec, Tmat[0], Tmat[1], Tmat[2], Tvec); */
   for (i = 0; i < 3; i++) {
      Vector3D_transformVecs_inverse(vertices[i], Tmat[0], Tmat[1], Tmat[2], vertices[i]);
      Vector3D_transformVecs_inverse(arccenters[i], Tmat[0], Tmat[1], Tmat[2], arccenters[i]);
   }
   
   for (i = 0; i < 3; i++) {
      next = i + 1;
      if (next > 2)
         next = 0;
      edges[i] = Curve_allocate(arccenters[i], vertices[i], vertices[next]);
   }
   
   for (i = 0; i < 3; i++) {
      Vector3D_free(vertices[i]);
      Vector3D_free(arccenters[i]);
   }

   Curve_free(curve);
}

void GST_sphGSTtoCart(real radius, real theta, real phi,
                      Vector3D point, Vector3D normal,
                      Vector3D dPoint_dTheta, Vector3D dPoint_dPhi) {
   Vector3D A[3];
   unsigned int i;
   for (i = 0; i < 3; i++)
      A[i] = Vector3D_allocate();
   
   A[0]->z = 1;
   A[1]->y = -1;
   A[2]->x = 1;
   
   point->x = radius * cos(theta) * sin(phi);
   point->y = radius * sin(theta) * sin(phi);
   point->z = radius * cos(phi);
   Vector3D_transformVecs(point, A[0], A[1], A[2], point);
   
   normal->x = cos(theta) * sin(phi);
   normal->y = sin(theta) * sin(phi);
   normal->z = cos(phi);
   Vector3D_transformVecs(normal, A[0], A[1], A[2], normal);
   
   dPoint_dTheta->x = radius * -sin(theta) * sin(phi);
   dPoint_dTheta->y = radius * cos(theta) * sin(phi);
   dPoint_dTheta->z = 0;
   Vector3D_transformVecs(dPoint_dTheta, A[0], A[1], A[2], dPoint_dTheta);
   
   dPoint_dPhi->x = radius * cos(theta) * cos(phi);
   dPoint_dPhi->y = radius * sin(theta) * cos(phi);
   dPoint_dPhi->z = radius * -sin(phi);
   Vector3D_transformVecs(dPoint_dPhi, A[0], A[1], A[2], dPoint_dPhi);
   
   for (i = 0; i < 3; i++)
      Vector3D_free(A[i]);
}

void GST_getCircleArcIntersection(Vector3D center, Vector3D normal, real radius,
                                  Curve curve, Vector3D intPoint,
                                  real *alpha, Vector3D dIntPoint_dAlpha) {
   Vector3D normal1, normal2;
   Vector3D lTvec, lTmat[3];
   Vector3D curveAC, curveV1, curveV2, curveX, curveY;
   Vector3D circleCenter, circleNormal;
   Vector3D points[2];
   Vector3D x0, line_dir, x1;
   lTvec = Vector3D_allocate();
   lTmat[0] = Vector3D_allocate();
   lTmat[1] = Vector3D_allocate();
   lTmat[2] = Vector3D_allocate();
   points[0] = Vector3D_allocate();
   points[1] = Vector3D_allocate();
   curveX = Vector3D_allocate();
   curveY = Vector3D_allocate();
   curveAC = Vector3D_allocate();
   curveV1 = Vector3D_allocate();
   curveV2 = Vector3D_allocate();
   circleCenter = Vector3D_allocate();
   circleNormal = Vector3D_allocate();
   x0 = Vector3D_allocate();
   x1 = Vector3D_allocate();
   line_dir = Vector3D_allocate();
   
   real p1, p2,a11,a12, a21, a22, invdet, line_dirVec[3];
   real dx, dy, dr, D, discrim, signdy, angle2;
   int yesno;
   real thresholdUP = 1e-6;
   unsigned int others[2];
   real maxElem;
   real onplaneThreshold = 1e-6;
   real threshold = 1e-6;
   real x0_vec[3];
   unsigned int i, maxIndex =0;
   others[0] = 0; others[1] = 0;
	a11 = 0.0; a12 = 0.0; a21 = 0.0; a22 = 0.0;
   normal1 = Vector3D_allocate();
   Vector3D_copy(normal1, normal);
   p1 = Vector3D_dot(normal1, center);
   
   normal2 = Vector3D_allocate();
   Vector3D_cross(normal2, curve->X, curve->Y);
   p2 = Vector3D_dot(normal2, curve->ac);
   if (Vector3D_dot(normal1, normal2) > (1 - threshold))
      printf("normals are almost coincident!!\n");
   
   Vector3D_cross(line_dir, normal1, normal2);
   maxElem = -1.0;
   line_dirVec[0] = fabs(line_dir->x);
   line_dirVec[1] = fabs(line_dir->y);
   line_dirVec[2] = fabs(line_dir->z);
   for (i = 0; i < 3; i++) {
      if (line_dirVec[i] > maxElem) {
         maxIndex = i;
         maxElem = line_dirVec[i];
      }
   }
   
   if (maxIndex == 0) {
      a11 = normal1->y; a12 = normal1->z;
      a21 = normal2->y; a22 = normal2->z;
      others[0] = 1; others[1] = 2;
   } else if (maxIndex == 1) {
      a11 = normal1->x; a12 = normal1->z;
      a21 = normal2->x; a22 = normal2->z;
      others[0] = 0; others[1] = 2;
   } else if (maxIndex == 2) {
      a11 = normal1->x; a12 = normal1->y;
      a21 = normal2->x; a22 = normal2->y;
      others[0] = 0; others[1] = 1;
   }
   invdet = 1.0 / (a11 * a22 - a12 * a21);
   // form N matrix, b vector, solve
   x0_vec[maxIndex] = 0;
   x0_vec[others[0]] = invdet * (a22 * p1 - a12 * p2);
   x0_vec[others[1]] = invdet * (-a21 * p1 + a11 * p2);
   x0->x= x0_vec[0];
   x0->y= x0_vec[1];
   x0->z= x0_vec[2];
   
   if ((fabs(p1 - Vector3D_dot(normal1, x0)) > onplaneThreshold) ||
       (fabs(p2 - Vector3D_dot(normal2, x0)) > onplaneThreshold)) {
      printf("failed to find an x0 on the circle-arc line\n");
      printf("%f %f\n", fabs(p1 - Vector3D_dot(normal1, x0)), fabs(p2 - Vector3D_dot(normal2, x0)));
   }
   
   Vector3D_copy(lTvec, curve->ac);
   Vector3D_sub(curveAC, curve->ac, lTvec);
   Vector3D_sub(curveV1, curve->v1, lTvec);
   Vector3D_sub(curveV2, curve->v2, lTvec);
   Vector3D_sub(circleCenter, center, lTvec);
   Vector3D_sub(x0, x0, lTvec);
   
   Vector3D_copy(lTmat[0], curve->X);
   Vector3D_copy(lTmat[1], curve->Y);
   Vector3D_copy(lTmat[2], normal2);
   Vector3D_transformVecs_inverse(curveAC, lTmat[0], lTmat[1], lTmat[2], curveAC);
   Vector3D_transformVecs_inverse(curveV1, lTmat[0], lTmat[1], lTmat[2], curveV1);
   Vector3D_transformVecs_inverse(curveV2, lTmat[0], lTmat[1], lTmat[2], curveV2);
   Vector3D_transformVecs_inverse(curveX, lTmat[0], lTmat[1], lTmat[2], curve->X);
   Vector3D_transformVecs_inverse(curveY, lTmat[0], lTmat[1], lTmat[2], curve->Y);
   Vector3D_transformVecs_inverse(circleCenter, lTmat[0], lTmat[1], lTmat[2], circleCenter);
   Vector3D_transformVecs_inverse(circleNormal, lTmat[0], lTmat[1], lTmat[2], circleNormal);
   Vector3D_transformVecs_inverse(x0, lTmat[0], lTmat[1], lTmat[2], x0);
   Vector3D_transformVecs_inverse(line_dir, lTmat[0], lTmat[1], lTmat[2], line_dir);
   
   Vector3D_add(x1, x0, line_dir);
   dx = x1->x - x0->x;
   dy = x1->y - x0->y;
   dr = sqrt(dx*dx + dy*dy);
   D = x0->x*x1->y - x1->x * x0->y;
   discrim = curve->radius*curve->radius * dr*dr - D *D;
   if (discrim < 0) {
      printf("no intersection found for arc and circle! aborting...\n");
      exit(-1);
   }
   signdy = (dy >= 0)? 1.0: -1.0;
   
   points[0]->x = (D * dy + signdy * dx * sqrt(discrim))/(dr*dr);  
   points[0]->y = (-D *dx + fabs(dy) * sqrt(discrim))/(dr*dr);
   points[1]->x = (D * dy - signdy * dx * sqrt(discrim))/(dr*dr);
   points[1]->y = (-D*dx -fabs(dy) * sqrt(discrim))/(dr*dr);
   yesno = 0;
   for (i=0; i < 2; i++) {
      if (thresholdUP < abs(radius - Vector3D_distance(points[i], circleCenter)))
         continue;
      angle2 = atan2(Vector3D_dot(curveY, points[i]), Vector3D_dot(curveX, points[i]));
      if ((curve->endTheta < 0) & (angle2 < 0) & (angle2 > curve->endTheta))
         yesno = 1;
      if ((curve->endTheta > 0) & (angle2 > 0) & (angle2 < curve->endTheta))
         yesno = 1;
      
      if (yesno > 0) {
         *alpha = angle2 / curve->endTheta;
         Curve_getParamPoint(intPoint, dIntPoint_dAlpha, curve, *alpha);

         Vector3D_free(normal1);
         Vector3D_free(normal2);
         Vector3D_free(lTvec);
         Vector3D_free(lTmat[0]);
         Vector3D_free(lTmat[1]);
         Vector3D_free(lTmat[2]);
         Vector3D_free(curveAC);
         Vector3D_free(curveV1);
         Vector3D_free(curveV2);
         Vector3D_free(curveX);
         Vector3D_free(curveY);
         Vector3D_free(circleCenter);
         Vector3D_free(circleNormal);
         Vector3D_free(points[0]);
         Vector3D_free(points[1]);
         Vector3D_free(x0);
         Vector3D_free(x1);
         Vector3D_free(line_dir);

         return;
      }
      
   }
   
   Vector3D_free(normal1);
   Vector3D_free(normal2);
   Vector3D_free(lTvec);
   Vector3D_free(lTmat[0]);
   Vector3D_free(lTmat[1]);
   Vector3D_free(lTmat[2]);
   Vector3D_free(curveAC);
   Vector3D_free(curveV1);
   Vector3D_free(curveV2);
   Vector3D_free(curveX);
   Vector3D_free(curveY);
   Vector3D_free(circleCenter);
   Vector3D_free(circleNormal);
   Vector3D_free(points[0]);
   Vector3D_free(points[1]);
   Vector3D_free(x0);
   Vector3D_free(x1);
   Vector3D_free(line_dir);
}


Curve Curve_allocate(Vector3D ac, Vector3D v1, Vector3D v2) {
   real dist_threshold = 1e-4;
   Curve c = (Curve)calloc(1, sizeof(struct _Curve));
   Vector3D v2ac = Vector3D_allocate();
   c->ac = Vector3D_allocate();
   c->v1 = Vector3D_allocate();
   c->v2 = Vector3D_allocate();
   c->X = Vector3D_allocate();
   c->Y = Vector3D_allocate();
   Vector3D_copy(c->ac, ac);
   Vector3D_copy(c->v1, v1);
   Vector3D_copy(c->v2, v2);
   c->radius = Vector3D_distance(c->ac, c->v1);
   if (abs(c->radius - Vector3D_distance(c->ac, c->v2)) > dist_threshold) {
      printf("Curve_allocate called with two points not on circular arc!\n");
      exit(-1);
   }
   Vector3D_sub(c->X, c->v1, c->ac);
   Vector3D_normalize(c->X);
   Vector3D_sub(c->Y, c->v2, c->ac);
   Vector3D_addscaled(c->Y, c->Y, - Vector3D_dot(c->Y, c->X), c->X);
   Vector3D_normalize(c->Y);
   Vector3D_sub(v2ac, c->v2, c->ac);
   c->endTheta = atan2(Vector3D_dot(c->Y, v2ac), Vector3D_dot(c->X, v2ac));
   Vector3D_free(v2ac);
   return c;
}

void Curve_free(Curve c) {
   Vector3D_free(c->ac);
   Vector3D_free(c->v1);
   Vector3D_free(c->v2);
   Vector3D_free(c->X);
   Vector3D_free(c->Y);
   free(c);
}

void Curve_getParamPoint(Vector3D p, Vector3D dpdalpha, Curve c, real alpha) {
   alpha = alpha * c->endTheta;
   Vector3D_copy(p, c->ac);
   Vector3D_addscaled(p, p, c->radius * cos(alpha), c->X);
   Vector3D_addscaled(p, p, c->radius * sin(alpha), c->Y);
   if ( dpdalpha != NULL ) {
      Vector3D_scale(dpdalpha, 0.0);
      Vector3D_addscaled(dpdalpha, dpdalpha, c->radius * -sin(alpha), c->X);
      Vector3D_addscaled(dpdalpha, dpdalpha, c->radius * cos(alpha), c->Y);
   }
}

void PermuteVectors(Vector3D v0, Vector3D v1, Vector3D v2) {
   Vector3D junk = Vector3D_allocate();
   Vector3D_copy(junk, v0);
   Vector3D_copy(v0, v1);
   Vector3D_copy(v1, v2);
   Vector3D_copy(v2, junk);
   Vector3D_free(junk);

}

void GST_getCentroid(GST gst) {
   Vector3D dTmat[3], dTvec;
   real junkreal;
   Curve edges[3];
   dTmat[0] = Vector3D_allocate();
   dTmat[1] = Vector3D_allocate();
   dTmat[2] = Vector3D_allocate();
   
   dTvec = Vector3D_allocate();
   
   GST_getStandardPosition(gst, edges, dTmat, dTvec);
   GST_get_direct_quadPoint(gst->centroid, gst->normalAtCentroid, &junkreal,
                            1./3., 1./3., gst, edges, dTmat, dTvec);
   Curve_free(edges[0]);
   Curve_free(edges[1]);
   Curve_free(edges[2]);
   Vector3D_free(dTmat[0]);
   Vector3D_free(dTmat[1]);
   Vector3D_free(dTmat[2]);
   Vector3D_free(dTvec);
}

