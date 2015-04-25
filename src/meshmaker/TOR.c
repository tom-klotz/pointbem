#include "meshmaker.h"

TOR TOR_allocate(real torusRadius, real probeRadius, real startTheta, real endTheta,
                 real startPsi, real endPsi, Vector3D center, Vector3D normal,
                 Vector3D localX, Vector3D localY,
                 unsigned int thetaIndex, unsigned int psiIndex, int cavity,
                 unsigned int numthetadivs, unsigned int numpsidivs, unsigned int torusID) {
   TOR t = (TOR)calloc(1, sizeof(struct _TOR));

   real tmp, junk;
   t->c = torusRadius;
   t->a = probeRadius;

   if (endPsi < startPsi) {
      tmp = startPsi;
      startPsi = endPsi;
      endPsi = tmp;
   }
   t->startPsi = startPsi;
   t->endPsi = endPsi;

   t->startTheta = startTheta;
   t->endTheta = endTheta;
   t->thetaIndex = thetaIndex;
   t->psiIndex = psiIndex;
   t->cavity = cavity;
   t->numthetadivs = numthetadivs;
   t->numpsidivs = numpsidivs;
   t->torusID = torusID;
   t->isLocal = 0;

   t->center = Vector3D_allocate();
   Vector3D_copy(t->center, center);
   t->normal = Vector3D_allocate();
   Vector3D_copy(t->normal, normal);
   t->localX = Vector3D_allocate();
   Vector3D_copy(t->localX, localX);
   t->localY = Vector3D_allocate();
   Vector3D_copy(t->localY, localY);

   t->Tvec = Vector3D_allocate();
   Vector3D_sub(t->Tvec, t->Tvec, t->center);
   Vector3D_transformVecs_inverse(t->Tvec, t->localX, t->localY, t->normal, t->Tvec);
  
   t->Tmat[0] = Vector3D_allocate();
   t->Tmat[1] = Vector3D_allocate();
   t->Tmat[2] = Vector3D_allocate();
   t->Tmat[0]->x = 1.0;
   t->Tmat[1]->y = 1.0;
   t->Tmat[2]->z = 1.0;
   Vector3D_transformVecs_inverse(t->Tmat[0], t->localX, t->localY, t->normal, t->Tmat[0]);
   Vector3D_transformVecs_inverse(t->Tmat[1], t->localX, t->localY, t->normal, t->Tmat[1]);
   Vector3D_transformVecs_inverse(t->Tmat[2], t->localX, t->localY, t->normal, t->Tmat[2]);

   t->centroid = Vector3D_allocate();
   t->normalAtCentroid =Vector3D_allocate();
   torToCart(t->centroid, t->normalAtCentroid, t->c, t->a, (startTheta+endTheta)/2,
             (startPsi + endPsi)/2);

   if (t->cavity)
      Vector3D_scale(t->normalAtCentroid, -1.0);

   //Tmat takes from global to local; so apply Tmat inverse to go FROM local TO global
   // and add Tvec AFTER rotation
   if (t->isLocal == 0) {
      Vector3D_transformVecs_inverse(t->centroid, t->Tmat[0], t->Tmat[1], t->Tmat[2], t->centroid);   
      Vector3D_add(t->centroid, t->centroid, t->center);
      Vector3D_transformVecs_inverse(t->normalAtCentroid, t->Tmat[0], t->Tmat[1], t->Tmat[2], t->normalAtCentroid);
   }

   t->area = (endTheta - startTheta) * (t->a * t->c * (endPsi - startPsi) + t->a * t->a * (sin(endPsi) - sin(startPsi)));
   if (t->area < 0.0)
      t->area = - t->area;

   t->maxEdgelength = t->a * (endPsi - startPsi);
   junk = (t->c + t->a * cos(endPsi)) * (endTheta - startTheta);
   if (junk > t->maxEdgelength)
      t->maxEdgelength = junk;
   junk  =(t->c + t->a * cos(startPsi)) * (endTheta - startTheta);
   if (junk > t->maxEdgelength)
      t->maxEdgelength = junk;

   // horrible directquad
   t->numdirectquadpoints = 16;
   t->directquadpoints = (Vector3D *)calloc(16, sizeof(Vector3D));
   t->directquadnormals = (Vector3D *)calloc(16, sizeof(Vector3D));
   t->directquadweights = (real *)calloc(16, sizeof(real));
   {
      unsigned int count = 0, i, j;
      static QuadratureRule qr = NULL;

#ifdef OMP
#pragma omp critical
#endif
      if (qr == NULL) {
         qr = QuadratureRule_allocate(4);
      }

      for (i = 0; i < qr->order; i++) {
         real curTheta = t->startTheta + qr->x[i] * (t->endTheta - t->startTheta);
         for (j = 0; j < qr->order; j++) {
            t->directquadpoints[count] = Vector3D_allocate();
            t->directquadnormals[count] = Vector3D_allocate();

            real curPsi = t->startPsi + qr->x[j] * (t->endPsi - t->startPsi);
            torToCart(t->directquadpoints[count], t->directquadnormals[count], t->c, t->a, curTheta, curPsi);
            Vector3D_sub(t->directquadpoints[count], t->directquadpoints[count], t->Tvec);
            Vector3D_transformVecs_inverse(t->directquadpoints[count], t->Tmat[0], t->Tmat[1], t->Tmat[2], t->directquadpoints[count]);
            Vector3D_transformVecs_inverse(t->directquadnormals[count], t->Tmat[0], t->Tmat[1], t->Tmat[2], t->directquadnormals[count]);

            if (t->cavity)
               Vector3D_scale(t->directquadnormals[count], -1.0);

            real detJ = t->a * (t->c + t->a * cos(curPsi));
            t->directquadweights[count] = fabs(t->endTheta - t->startTheta) * qr->w[i] * fabs(t->endPsi - t->startPsi) * qr->w[j] * detJ;
            count++;
         }
      }
   }

   return t;
}

void TOR_free(TOR t) {
   unsigned int i;

   Vector3D_free(t->Tmat[0]);
   Vector3D_free(t->Tmat[1]);
   Vector3D_free(t->Tmat[2]);
   Vector3D_free(t->Tvec);
   Vector3D_free(t->centroid);
   Vector3D_free(t->normalAtCentroid);
   Vector3D_free(t->center);
   Vector3D_free(t->normal);
   Vector3D_free(t->localX);
   Vector3D_free(t->localY);

   for (i = 0; i < t->numdirectquadpoints; i++) {
      Vector3D_free(t->directquadpoints[i]);
      Vector3D_free(t->directquadnormals[i]);
   }
   free(t->directquadpoints);
   free(t->directquadnormals);
   free(t->directquadweights);

   free(t);
}

void TOR_readfile(unsigned int *numTORpanels, TOR** torList, FILE* file, unsigned int ignorecav) {
   char line[128];
   unsigned int i = 0;
   unsigned int numLinesPerTOR = 8;

   Vector3D center, normal, localX, localY;
   real torusRadius, probeRadius, startTheta, endTheta, startPsi, endPsi;
   unsigned int thetaIndex, psiIndex, numthetadivs, numpsidivs, torusID;
   int cavity;
   center = Vector3D_allocate();
   normal = Vector3D_allocate();
   localX = Vector3D_allocate();
   localY = Vector3D_allocate();

   while (fgets(line, 128, file))
      i++;
   rewind(file);

   *numTORpanels = i / numLinesPerTOR;

   *torList = (TOR *)calloc(*numTORpanels, sizeof(TOR));
   for (i = 0; i < *numTORpanels; i++) {
      // line 1: torus center
      fgets(line, 128, file);
#ifdef REAL_IS_DOUBLE
      sscanf(line, "%lf %lf %lf", &center->x, &center->y, &center->z);
#else
#ifdef REAL_IS_FLOAT
      sscanf(line, "%f %f %f", &center->x, &center->y, &center->z);
#endif
#endif
      // line 2: normal
      fgets(line, 128, file);
#ifdef REAL_IS_DOUBLE
      sscanf(line, "%lf %lf %lf", &normal->x, &normal->y, &normal->z);
#else
#ifdef REAL_IS_FLOAT
      sscanf(line, "%f %f %f", &normal->x, &normal->y, &normal->z);
#endif
#endif
      // line 3: local X
      fgets(line, 128, file);
#ifdef REAL_IS_DOUBLE
      sscanf(line, "%lf %lf %lf", &localX->x, &localX->y, &localX->z);
#else
#ifdef REAL_IS_FLOAT
      sscanf(line, "%f %f %f", &localX->x, &localX->y, &localX->z);
#endif
#endif
      // line 4: local Y;
      fgets(line, 128, file);
#ifdef REAL_IS_DOUBLE
      sscanf(line, "%lf %lf %lf", &localY->x, &localY->y, &localY->z);
#else
#ifdef REAL_IS_FLOAT
      sscanf(line, "%f %f %f", &localY->x, &localY->y, &localY->z);
#endif
#endif
      // line 5: torusRadius (c), probeRadius (a), startTheta
      fgets(line, 128, file);
#ifdef REAL_IS_DOUBLE
      sscanf(line, "%lf %lf %lf", &torusRadius, &probeRadius, &startTheta);
#else
#ifdef REAL_IS_FLOAT
      sscanf(line, "%f %f %f", &torusRadius, &probeRadius, &startTheta);
#endif
#endif
    
      // line 6: endTheta, startPsi, endPsi;
      fgets(line, 128, file);
#ifdef REAL_IS_DOUBLE
      sscanf(line, "%lf %lf %lf", &endTheta, &startPsi, &endPsi);
#else
#ifdef REAL_IS_FLOAT
      sscanf(line, "%f %f %f", &endTheta, &startPsi, &endPsi);
#endif
#endif
      // line 7: theta index, psi index, junk
      fgets(line, 128, file);
      sscanf(line, "%d %d %d", &thetaIndex, &psiIndex, &cavity);
      // line 8: numtheatdivs, numpsidivs, torusID
      fgets(line, 128, file);
      sscanf(line, "%d %d %d", &numthetadivs, &numpsidivs, &torusID);

      if (ignorecav)
         cavity = 0;

      // actually initialize
      if (cavity)
         (*torList)[i] = TOR_allocate(torusRadius, probeRadius, endTheta, startTheta,
                                      startPsi, endPsi, center, normal, localX, localY,
                                      thetaIndex, psiIndex, cavity, numthetadivs, numpsidivs, torusID);
      else
         (*torList)[i] = TOR_allocate(torusRadius, probeRadius, startTheta, endTheta,
                                      startPsi, endPsi, center, normal, localX, localY,
                                      thetaIndex, psiIndex, cavity, numthetadivs, numpsidivs, torusID);
   }
  
   Vector3D_free(center);
   Vector3D_free(normal);
   Vector3D_free(localX);
   Vector3D_free(localY);
  
}

void torToCart(Vector3D point, Vector3D normal,
               real c, real a, real theta, real psi) {
   real costheta, cospsi, sintheta, sinpsi, cacospsi;

   sintheta = sin(theta);
   costheta = cos(theta);
   sinpsi = sin(psi);
   cospsi = cos(psi);

   cacospsi = c + a * cospsi;

   point->x = costheta * cacospsi;
   point->y = sintheta * cacospsi;
   point->z = a * sinpsi;

   if (normal) {
      struct _Vector3D thetahat;
      struct _Vector3D psihat;

      thetahat.x = -point->y;
      thetahat.y = point->x;
      thetahat.z = 0;

      psihat.x = -a * sinpsi * costheta;
      psihat.y = -a * sinpsi * sintheta;
      psihat.z = a * cospsi;

      Vector3D_cross(normal, &psihat, &thetahat);
      Vector3D_normalize(normal);
   }
}
