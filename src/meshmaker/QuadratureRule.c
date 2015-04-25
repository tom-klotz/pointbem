#include "meshmaker.h"
#include "petscblaslapack.h"

QuadratureRule QuadratureRule_allocate(unsigned int order) {
   QuadratureRule qr = (QuadratureRule)calloc(1, sizeof(struct _QuadratureRule));
   real fact;
   real STP[7][2];
   unsigned int i;
   real ksi, eta;
   real L1, dL1dKsi, dL1dEta;
   real L2, dL2dKsi, dL2dEta;
   real L3, dL3dKsi, dL3dEta;   

   if (order > MAX_QUADRATURE_ORDER) {
      printf("Warning: Requested quadrature order %u exceeds maximum quadrature order %u\n", order, MAX_QUADRATURE_ORDER);
      order = MAX_QUADRATURE_ORDER;
   }

   qr->order = order;
   qr->x = Vector_allocate(order);
   qr->w = Vector_allocate(order);

   QuadratureRule_generatepoints(qr, 0.0, 1.0);

   /* this seems magical, but look at my notebook, ~ 6/18/03, for
      derivation -- JPB */

   STP[0][0] = 1.;    STP[0][1] = 0.;
   STP[1][0] = 1./4;  STP[1][1] = -sqrt(3)/4;
   STP[2][0] = -1./2; STP[2][1] = -sqrt(3)/2;
   STP[3][0] = -1./2; STP[3][1] = 0;
   STP[4][0] = -1./2; STP[4][1] = sqrt(3)/2;
   STP[5][0] = 1./4;  STP[5][1] = sqrt(3)/4;
   STP[6][0] = 0;     STP[6][1] = 0;

   for (i = 0; i < 7; i++) {
      ksi = STP[i][0];  eta = STP[i][1];
      L1 = 1./3 * (1 + 2 * ksi); dL1dKsi = 2./3; dL1dEta = 0;
      L2 = 1./3 * (1 - ksi - sqrt(3) * eta); dL2dKsi = -1./3; dL2dEta = -sqrt(3)/3;
      L3 = 1./3 * (1 - ksi + sqrt(3) * eta); dL3dKsi = -1./3; dL3dEta = sqrt(3)/3;
      qr->dNdKsi[0][i] = 4 * L1 * dL1dKsi - dL1dKsi;
      qr->dNdKsi[1][i] = 4 * (L1 * dL2dKsi + L2 * dL1dKsi);
      qr->dNdKsi[2][i] = 4 * L2 * dL2dKsi - dL2dKsi;
      qr->dNdKsi[3][i] = 4 * (L2 * dL3dKsi + L3 * dL2dKsi);
      qr->dNdKsi[4][i] = 4 * L3 * dL3dKsi - dL3dKsi;
      qr->dNdKsi[5][i] = 4 * (L3 * dL1dKsi + L1 * dL3dKsi);

      qr->dNdEta[0][i] = 4 * L1 * dL1dEta - dL1dEta;
      qr->dNdEta[1][i] = 4 * (L1 * dL2dEta + L2 * dL1dEta);
      qr->dNdEta[2][i] = 4 * L2 * dL2dEta - dL2dEta;
      qr->dNdEta[3][i] = 4 * (L2 * dL3dEta + L3 * dL2dEta);
      qr->dNdEta[4][i] = 4 * L3 * dL3dEta - dL3dEta;
      qr->dNdEta[5][i] = 4 * (L3 * dL1dEta + L1 * dL3dEta);
   }

   fact = sqrt(3)/80;

   qr->SPW[0] = 3 * fact;
   qr->SPW[1] = 8 * fact;
   qr->SPW[2] = 3 * fact;
   qr->SPW[3] = 8 * fact;
   qr->SPW[4] = 3 * fact;
   qr->SPW[5] = 8 * fact;
   qr->SPW[6] = 27 * fact;

   return qr;
}

void QuadratureRule_free(QuadratureRule qr) {
   Vector_free(qr->x);
   Vector_free(qr->w);
}

void QuadratureRule_generatepoints(QuadratureRule qr, real lower, real upper) {
   unsigned int i, N, LWORK, LDA;
   char JOBZ = 'V';
   char UPLO = 'U';
   Vector A, W, WORK; 
   int INFO;

   N = qr->order; LDA = N;
   A = Vector_allocate(N * N);
   W = Vector_allocate(N);
   LWORK = 10 * N - 1;
   WORK = Vector_allocate(LWORK);

   for (i = 1; i < N; i++) {
      A[i*N+(i-1)] = .5 * pow(1 - pow(2.*(double)i,-2), -.5);
   }

#ifdef OMP
#pragma omp critical (lapack)
#endif
   LAPACKsyev_(&JOBZ, &UPLO, &N, A, &LDA, W, WORK, &LWORK, &INFO);

   if (INFO) {
      printf("s/dsyev returned error code %d\n", INFO);
      exit(-2);
   }

   for (i = 0; i < N; i++) {
      qr->x[i] = W[i]/2 + .5;
      qr->w[i] = pow(A[i*N],2);
   }

   Vector_free(WORK);
   Vector_free(W);
   Vector_free(A);
}
