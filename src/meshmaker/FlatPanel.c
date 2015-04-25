#include "meshmaker.h"

/* Constructors and Destructors */

FlatPanel FlatPanel_allocate(Vector3D v1, Vector3D v2, Vector3D v3) {
   unsigned int i, next, prev;
   real alpha;
   FlatPanel panel = (FlatPanel)calloc(1, sizeof(struct _FlatPanel));

   panel->vertex[0] = Vector3D_allocate();
   panel->vertex[1] = Vector3D_allocate();
   panel->vertex[2] = Vector3D_allocate();
   panel->centroid = Vector3D_allocate();
   panel->panelaxis[0] = Vector3D_allocate();
   panel->panelaxis[1] = Vector3D_allocate();
   panel->panelaxis[2] = Vector3D_allocate();
   panel->panelvertex[0] = Vector3D_allocate();
   panel->panelvertex[1] = Vector3D_allocate();
   panel->panelvertex[2] = Vector3D_allocate();
   panel->panelaxisnum[0] = Vector3D_allocate();
   panel->panelaxisnum[1] = Vector3D_allocate();
   panel->panelaxisnum[2] = Vector3D_allocate();
   panel->panelvertexnum[0] = Vector3D_allocate();
   panel->panelvertexnum[1] = Vector3D_allocate();
   panel->panelvertexnum[2] = Vector3D_allocate();
   panel->edges[0] = Vector3D_allocate();
   panel->edges[1] = Vector3D_allocate();
   panel->edges[2] = Vector3D_allocate();
   Vector3D_copy(panel->vertex[0], v1);
   Vector3D_copy(panel->vertex[1], v2);
   Vector3D_copy(panel->vertex[2], v3);      

   Vector3D_add(panel->centroid, panel->vertex[0], panel->vertex[1]);
   Vector3D_add(panel->centroid, panel->centroid, panel->vertex[2]);
   Vector3D_scale(panel->centroid, (real)1 / (real)3);

   panel->edgelength[0] = Vector3D_distance(panel->vertex[1], panel->vertex[0]);
   panel->edgelength[1] = Vector3D_distance(panel->vertex[2], panel->vertex[1]);
   panel->edgelength[2] = Vector3D_distance(panel->vertex[0], panel->vertex[2]);
   
   Vector3D_sub(panel->panelaxis[0], panel->vertex[2], panel->vertex[0]);
   Vector3D_sub(panel->panelaxis[1], panel->vertex[1], panel->vertex[0]);
   Vector3D_cross(panel->panelaxis[2], panel->panelaxis[0], panel->panelaxis[1]);

   panel->area = 0.5 * Vector3D_length(panel->panelaxis[2]);

   Vector3D_normalize(panel->panelaxis[2]);
   Vector3D_normalize(panel->panelaxis[0]);
   Vector3D_cross(panel->panelaxis[1], panel->panelaxis[2], panel->panelaxis[0]);

   // this section creates a panel coordinate system by determining if
   // panel edge 2 (v2-v3) or panel edge 3 (v3-v1) is less aligned
   // with edge 1 (v1-v2)... and as such should be less susceptible to
   // cancellation error

   Vector3D_sub(panel->edges[0], panel->vertex[1], panel->vertex[0]);
   Vector3D_normalize(panel->edges[0]);
   Vector3D_sub(panel->edges[1], panel->vertex[2], panel->vertex[1]);
   Vector3D_normalize(panel->edges[1]);
   Vector3D_sub(panel->edges[2], panel->vertex[0], panel->vertex[2]);
   Vector3D_normalize(panel->edges[2]);

   Vector3D_copy(panel->panelaxisnum[0], panel->edges[0]);
   Vector3D_copy(panel->panelaxisnum[1], panel->edges[1]);
   if (fabs(Vector3D_dot(panel->edges[2], panel->panelaxisnum[0]))
       < fabs(Vector3D_dot(panel->panelaxisnum[1], panel->panelaxisnum[0])) ) {
      Vector3D_copy(panel->panelaxisnum[1], panel->edges[2]);
      Vector3D_scale(panel->panelaxisnum[1], (real)-1.0);
   }

   alpha = Vector3D_dot(panel->panelaxisnum[1], panel->panelaxisnum[0]);
   Vector3D_scale(panel->edges[0], -alpha);
   Vector3D_add(panel->panelaxisnum[1], panel->panelaxisnum[1], panel->edges[0]);
   Vector3D_normalize(panel->panelaxisnum[1]);
   
   Vector3D_cross(panel->panelaxisnum[2], panel->panelaxisnum[0], panel->panelaxisnum[1]);
   
   // transform all panels into panel local coordinates

   for (i = 0; i < 3; i++) {
      Vector3D vmc = Vector3D_allocate();

      Vector3D_sub(vmc, panel->vertex[i], panel->centroid);

      panel->panelvertex[i]->x = panel->panelaxis[0]->x * vmc->x +
         panel->panelaxis[0]->y * vmc->y +
         panel->panelaxis[0]->z * vmc->z;
      panel->panelvertex[i]->y = panel->panelaxis[1]->x * vmc->x +
         panel->panelaxis[1]->y * vmc->y +
         panel->panelaxis[1]->z * vmc->z;
      panel->panelvertex[i]->z = panel->panelaxis[2]->x * vmc->x +
         panel->panelaxis[2]->y * vmc->y +
         panel->panelaxis[2]->z * vmc->z;

      panel->panelvertexnum[i]->x = panel->panelaxisnum[0]->x * vmc->x +
         panel->panelaxisnum[0]->y * vmc->y +
         panel->panelaxisnum[0]->z * vmc->z;
      panel->panelvertexnum[i]->y = panel->panelaxisnum[1]->x * vmc->x +
         panel->panelaxisnum[1]->y * vmc->y +
         panel->panelaxisnum[1]->z * vmc->z;
      panel->panelvertexnum[i]->z = panel->panelaxisnum[2]->x * vmc->x +
         panel->panelaxisnum[2]->y * vmc->y +
         panel->panelaxisnum[2]->z * vmc->z;

      Vector3D_free(vmc);
   }

   // now form all edge related information, so in calcp we can
   // determine whether an observation point has its projection
   // "inside" the panel or not
   for (i=0; i < 3; i++) {  // this is a rotation by 90 degrees so our
                            // edge becomes the equation for a line
      next = i+1;
      if (i==2)
         next = 0;
      if (i==0)
         prev = 2;
      else
         prev = i -1;
        
      Vector3D_sub(panel->edges[i],
                   panel->panelvertexnum[next],
                   panel->panelvertexnum[i]);
      alpha = panel->edges[i]->x;
      panel->edges[i]->x = -panel->edges[i]->y;
      panel->edges[i]->y = alpha;
      panel->edgeRHS[i] = Vector3D_dot(panel->panelvertexnum[i], panel->edges[i]);
      panel->edgeOV[i] = Vector3D_dot(panel->panelvertexnum[prev], panel->edges[i]) - panel->edgeRHS[i];
   } 
   
   panel->contributionC[0] = (panel->panelvertex[1]->x - panel->panelvertex[0]->x) / panel->edgelength[0];
   panel->contributionS[0] = (panel->panelvertex[1]->y - panel->panelvertex[0]->y) / panel->edgelength[0];
   panel->contributionC[1] = (panel->panelvertex[2]->x - panel->panelvertex[1]->x) / panel->edgelength[1];
   panel->contributionS[1] = (panel->panelvertex[2]->y - panel->panelvertex[1]->y) / panel->edgelength[1];
   panel->contributionC[2] = (panel->panelvertex[0]->x - panel->panelvertex[2]->x) / panel->edgelength[2];
   panel->contributionS[2] = (panel->panelvertex[0]->y - panel->panelvertex[2]->y) / panel->edgelength[2];

   panel->max_diag = panel->edgelength[0];
   if (panel->edgelength[1] > panel->max_diag)
      panel->max_diag = panel->edgelength[1];
   if (panel->edgelength[2] > panel->max_diag)
      panel->max_diag = panel->edgelength[2];

   panel->min_diag = panel->edgelength[0];
   if (panel->edgelength[1] < panel->min_diag)
      panel->min_diag = panel->edgelength[1];
   if (panel->edgelength[2] < panel->min_diag)
      panel->min_diag = panel->edgelength[2];
   
   FlatPanel_moments(panel);

   // Compute direct quad points

   // this is HORRIBLE we're fixing it to 16 point quad rule
   panel->numdirectquadpoints = 16;
   panel->directquadpoints = (Vector3D *)calloc(16, sizeof(Vector3D));
   panel->directquadnormals = (Vector3D *)calloc(16, sizeof(Vector3D));
   panel->directquadweights = (real *)calloc(16, sizeof(real));
   {
      real ksiTab[16];
      real etaTab[16];
      real weight[16];
      gen_Stroud_Rule(ksiTab, etaTab, weight);
      for (i = 0; i < 16; i++) {
         panel->directquadpoints[i] = Vector3D_allocate();
         panel->directquadnormals[i] = Vector3D_allocate();

         Vector3D_addscaled(panel->directquadpoints[i], panel->directquadpoints[i], ksiTab[i], v1);
         Vector3D_addscaled(panel->directquadpoints[i], panel->directquadpoints[i], etaTab[i], v2);
         Vector3D_addscaled(panel->directquadpoints[i], panel->directquadpoints[i], 1.0-etaTab[i]-ksiTab[i], v3);

         Vector3D_copy(panel->directquadnormals[i], panel->panelaxis[2]);
         Vector3D_scale(panel->directquadnormals[i], -1.0);
         panel->directquadweights[i] = weight[i] * 2.0 * panel->area;
      }
   }

   return panel;
}

void FlatPanel_free(FlatPanel panel) {
   Vector3D_free(panel->vertex[0]);
   Vector3D_free(panel->vertex[1]);
   Vector3D_free(panel->vertex[2]);
   Vector3D_free(panel->centroid);
   Vector3D_free(panel->panelaxis[0]);
   Vector3D_free(panel->panelaxis[1]);
   Vector3D_free(panel->panelaxis[2]);
   Vector3D_free(panel->panelaxisnum[0]);
   Vector3D_free(panel->panelaxisnum[1]);
   Vector3D_free(panel->panelaxisnum[2]);
   Vector3D_free(panel->panelvertex[0]);
   Vector3D_free(panel->panelvertex[1]);
   Vector3D_free(panel->panelvertex[2]);
   Vector3D_free(panel->panelvertexnum[0]);
   Vector3D_free(panel->panelvertexnum[1]);
   Vector3D_free(panel->panelvertexnum[2]);
   Vector3D_free(panel->edges[0]);
   Vector3D_free(panel->edges[1]);
   Vector3D_free(panel->edges[2]);

   free(panel);
}

/* Compute the panel moments for a panel
   These moments are used to accelerate panel integration
   They are also passed to Multipole_panel to compute the multipole moments
   This routine is taken directly from FASTCAP */

void FlatPanel_moments(FlatPanel panel) {
   unsigned int panelorder = 4;
   unsigned int i, j, nside, N, M, N1, M1, M2, MN1, MN2;
   real dx, dy, dxdy, dydx, SI, *xp, *yp, *xpn, *ypn;
   real *XP[4], *YP[4], **IM;
   real CS[16] = { 0.0, 1.0, 1.0, 1.5, 1.5, 3.75, 1.0, 3.0,
                   1.5, 7.5, 1.5, 1.5, 3.75, 1.5, 7.5, 3.75 };
 
   for(i = 0; i < 4; i++) {
      XP[i] = (real*)calloc(panelorder+3, sizeof(real));
      YP[i] = (real*)calloc(panelorder+3, sizeof(real));
   }

   /* Allocate the euclidean moments matrix, Imn. */
   IM = (real**)calloc(panelorder+1, sizeof(real*));
   for(i = 0; i <= panelorder; i++)
      IM[i] = (real*)calloc(panelorder+1, sizeof(real));

   /* First zero out the Moments matrix. */
   for(i = 0; i <= panelorder; i++)
      for(j = 0; j <= panelorder; j++)
         IM[i][j] = (real)0;

   /* Compute powers of x and y at corner pts. */
   for(i = 0; i < 3; i++) {
      xp = XP[i];
      yp = YP[i];
      xp[1] = panel->panelvertex[i]->x;
      yp[1] = panel->panelvertex[i]->y;
      for(j = 2; j <= panelorder+2; j++) {
         xp[j] = xp[j-1] * xp[1];
         yp[j] = yp[j-1] * yp[1];
      }
   }

   /* First moment, easy, just the panel area. */
   IM[0][0] = panel->area;

   /* By using centroid, (1,0) and (0,1) are zero, so begin with (2,0). */
   for(nside = 0; nside < 3; nside++) {
      xp = XP[nside];
      yp = YP[nside];
      if(nside == 2) {
         xpn = XP[0];
         ypn = YP[0];
      }
      else {
         xpn = XP[nside + 1];
         ypn = YP[nside + 1];
      }
 
      dx = xpn[1] - xp[1];
      dy = ypn[1] - yp[1];
 
      if (fabs(dx) >= fabs(dy)) {
         dydx = dy/dx;
         for(M = 2; M <= panelorder; M++) {
            M1 = M + 1;
            M2 = M + 2;
 
            SI = ((xpn[M1] * ypn[1]) - (xp[M1] * yp[1])) / M1
               + dydx * (xp[M2] - xpn[M2]) / (M1 * M2);
            IM[M][0] += SI;
 
            for(N = 1; N <= M; N++) {
               N1 = N + 1;
               MN1 = M - N + 1;
               SI = (xpn[MN1] * ypn[N1] - xp[MN1] * yp[N1]) / (MN1 * N1)
                  - (dydx * N * SI) / MN1;
               IM[M-N][N] += SI;
            }
         }
      }
      else {
         dxdy = dx/dy;
         for(M = 2; M <= panelorder; M++) {
            M1 = M + 1;
            M2 = M + 2;
 
            SI = (dxdy / (M1 * M2)) * (ypn[M2] - yp[M2]);
            IM[0][M] += SI;
 
            for(N = 1; N <= M; N++) {
               MN1 = M - N + 1;
               MN2 = MN1 + 1;
               SI = dxdy * ((xpn[N] * ypn[MN2] - xp[N] * yp[MN2]) / (MN1 * MN2)
                            - (N * SI / MN1));
               IM[N][M-N] += SI;
            }
         }
      }
   }
 
   /* Now Create the S vector for calcp. */
   for(i = 0, M = 0; M <= 4; M++) {
      for(N = 0; N <= (4 - M); N++) {
         i++;
         panel->moments[i] = IM[M][N] * CS[i];
      }
   }
 
   for(i = 0; i < 4; i++) {
      free(XP[i]);
      free(YP[i]);
   }

   for(i = 0; i <= panelorder; i++)
      free(IM[i]);
   free(IM);
}

/* Operations */

/* Returns the true memory used by a FlatPanel structure */

unsigned int FlatPanel_memory(unsigned int order, unsigned int doublelayer) {
   unsigned int FlatPanelmem =  sizeof(struct _FlatPanel) +
      10 * sizeof(struct _Vector3D) +
      27 * sizeof(real);

   return FlatPanelmem;
}

// quadrature rules taken from http://people.scs.fsu.edu/~burkardt/m_src/stroud/stroud.html
// (just like the rule in FlatIntegration.c)
void gen_arb_Stroud_rule(unsigned int rule, unsigned int *numquadpoints,
								 real **ytab, real **xtab, real **weight)
{
  real a, b, c, d, e, f, g, h, r, s, t, u, v, w;
  real w1, w2, w3, w4, w5, w6;
  switch (rule) {
  case 1:
	 *numquadpoints = 1;
	 *ytab = (real *)calloc(*numquadpoints, sizeof(real));
	 *xtab = (real *)calloc(*numquadpoints, sizeof(real));
	 *weight = (real *)calloc(*numquadpoints, sizeof(real));
	 (*xtab)[0] = 1.0 / 3.0;
	 (*ytab)[0] = 1.0 / 3.0;
	 (*weight)[0] = 1.0;
	 break;
  case 3:
	 *numquadpoints = 3;
	 *ytab = (real *)calloc(*numquadpoints, sizeof(real));
	 *xtab = (real *)calloc(*numquadpoints, sizeof(real));
	 *weight = (real *)calloc(*numquadpoints, sizeof(real));
	 a = 1.0;  b = 3.0; c = 4.0;  d = 6.0;
	 (*xtab)[0] = c/d; (*xtab)[1] = a/d; (*xtab)[2] = a/d;
	 (*ytab)[0] = a/d; (*ytab)[1] = c/d; (*ytab)[2] = a/d;
	 (*weight)[0] = a/b; (*weight)[1] = a/b; (*weight)[2] = a/b;
	 break;
  case 5:
	 *numquadpoints = 4;
	 *ytab = (real *)calloc(*numquadpoints, sizeof(real));
	 *xtab = (real *)calloc(*numquadpoints, sizeof(real));
	 *weight = (real *)calloc(*numquadpoints, sizeof(real));
	 a = 6.0;  b = 10.0;  c = 18.0;  d = 25.0;
	 e = -27.0;  f = 30.0;  g = 48.0;
	 (*xtab)[0] = b/f; (*xtab)[1] = c/f; (*xtab)[2] = a/f; (*xtab)[3] = a/f;
	 (*ytab)[0] = b/f; (*ytab)[1] = a/f; (*ytab)[2] = c/f; (*ytab)[3] = a/f;
	 (*weight)[0] = e/g; (*weight)[1] = d/g; (*weight)[2] = d/g; (*weight)[3] = d/g;
	 break;
  case 8:
	 *numquadpoints = 6;
	 *ytab = (real *)calloc(*numquadpoints, sizeof(real));
	 *xtab = (real *)calloc(*numquadpoints, sizeof(real));
	 *weight = (real *)calloc(*numquadpoints, sizeof(real));
	 a = 0.816847572980459;
    b = 0.091576213509771;
    c = 0.108103018168070;
    d = 0.445948490915965;
    v = 0.109951743655322;
    w = 0.223381589678011;
	 (*xtab)[0] = a; (*xtab)[1] = b; (*xtab)[2] = b;
	 (*xtab)[3] = c; (*xtab)[4] = d; (*xtab)[5] = d;
	 (*ytab)[0] = b; (*ytab)[1] = a; (*ytab)[2] = b;
	 (*ytab)[3] = d; (*ytab)[4] = c; (*ytab)[5] = d;
	 (*weight)[0] = v; (*weight)[1] = v; (*weight)[2] = v;
	 (*weight)[3] = w; (*weight)[4] = w; (*weight)[5] = w;
	 break;
  case 10:
	 *numquadpoints = 7;
	 *ytab = (real *)calloc(*numquadpoints, sizeof(real));
	 *xtab = (real *)calloc(*numquadpoints, sizeof(real));
	 *weight = (real *)calloc(*numquadpoints, sizeof(real));
	 a = 1.0/3.0;
	 b = (9.0 + 2.0 * sqrt(15.0) ) / 21.0;
	 c = (6.0 - sqrt(15.0)) / 21.0;
	 d = (9.0 - 2.0 * sqrt(15.0) ) / 21.0;
	 e = (6.0 + sqrt(15.0)) / 21.0;
	 u = 0.225;
	 v = (155.0 - sqrt(15.0)) / 1200.0;
	 w = (155.0 + sqrt(15.0)) / 1200.0;
	 (*xtab)[0] = a; (*xtab)[1] = b; (*xtab)[2] = c; (*xtab)[3] = c;
	 (*xtab)[4] = d; (*xtab)[5] = e; (*xtab)[6] = e;
	 (*ytab)[0] = a; (*ytab)[1] = c; (*ytab)[2] = b; (*ytab)[3] = c;
	 (*ytab)[4] = e; (*ytab)[5] = d; (*ytab)[6] = e;
	 (*weight)[0] = u; (*weight)[1] = v; (*weight)[2] = v; (*weight)[3] = v;
	 (*weight)[4] = w; (*weight)[5] = w; (*weight)[6] = w;
	 break;
  case 11:
	 *numquadpoints = 9;
	 *ytab = (real *)calloc(*numquadpoints, sizeof(real));
	 *xtab = (real *)calloc(*numquadpoints, sizeof(real));
	 *weight = (real *)calloc(*numquadpoints, sizeof(real));
	 a = 0.124949503233232;
    b = 0.437525248383384;
	 c = 0.797112651860071;
    d = 0.165409927389841;
    e = 0.037477420750088;
    u = 0.205950504760887;
    v = 0.063691414286223;
	 (*xtab)[0] = a; (*xtab)[1] = b; (*xtab)[2] = b;
	 (*xtab)[3] = c; (*xtab)[4] = c; (*xtab)[5] = d;
	 (*xtab)[6] = e; (*xtab)[7] = e; (*xtab)[8] = e;
	 (*ytab)[0] = b; (*ytab)[1] = a; (*ytab)[2] = b;
	 (*ytab)[3] = d; (*ytab)[4] = e; (*ytab)[5] = c;
	 (*ytab)[6] = e; (*ytab)[7] = c; (*ytab)[8] = d;
	 (*weight)[0] = u; (*weight)[1] = u; (*weight)[2] = u;
	 (*weight)[3] = v; (*weight)[4] = v; (*weight)[5] = v;
	 (*weight)[6] = v; (*weight)[7] = v; (*weight)[8] = v;
	 break;
  case 13:
	 *numquadpoints = 13;
	 *ytab = (real *)calloc(*numquadpoints, sizeof(real));
	 *xtab = (real *)calloc(*numquadpoints, sizeof(real));
	 *weight = (real *)calloc(*numquadpoints, sizeof(real));
	 a = 0.479308067841923;
    b = 0.260345966079038;
    c = 0.869739794195568;
    d = 0.065130102902216;
    e = 0.638444188569809;
    f = 0.312865496004875;
    g = 0.048690315425316;
    h = 1.0 / 3.0;
    t = 0.175615257433204;
	 u = 0.053347235608839;
    v = 0.077113760890257;
    w = -0.149570044467670;
	 (*xtab)[0] = a; (*xtab)[1] = b; (*xtab)[2] = b;
	 (*xtab)[3] = c; (*xtab)[4] = d; (*xtab)[5] = d;
	 (*xtab)[6] = e; (*xtab)[7] = e; (*xtab)[8] = f;
	 (*xtab)[9] = f; (*xtab)[10]= g; (*xtab)[11]= g;
	 (*xtab)[12]= h;
	 (*ytab)[0] = b; (*ytab)[1] = a; (*ytab)[2] = b;
	 (*ytab)[3] = d; (*ytab)[4] = c; (*ytab)[5] = d;
	 (*ytab)[6] = f; (*ytab)[7] = g; (*ytab)[8] = e;
	 (*ytab)[9] = g; (*ytab)[10]= e; (*ytab)[11]= f;
	 (*ytab)[12]= h;
	 (*weight)[0] = t; (*weight)[1] = t; (*weight)[2] = t;
	 (*weight)[3] = u; (*weight)[4] = u; (*weight)[5] = u;
	 (*weight)[6] = v; (*weight)[7] = v; (*weight)[8] = v;
	 (*weight)[9] = v; (*weight)[10]= v; (*weight)[11]= v;
	 (*weight)[12]= w;
	 break;
  case 15:
	 break;
  case 16:
	 break;
  case 18:
	 *numquadpoints = 19;
	 *ytab = (real *)calloc(*numquadpoints, sizeof(real));
	 *xtab = (real *)calloc(*numquadpoints, sizeof(real));
	 *weight = (real *)calloc(*numquadpoints, sizeof(real));
    a = 1.0 / 3.0;
    b = 0.02063496160252593;
    c = 0.4896825191987370;
    d = 0.1258208170141290;
    e = 0.4370895914929355;
    f = 0.6235929287619356;
    g = 0.1882035356190322;
    r = 0.9105409732110941;
    s = 0.04472951339445297;
    t = 0.7411985987844980;
    u = 0.03683841205473626;
    v = 0.22196288916076574;
    w1 = 0.09713579628279610;
    w2 = 0.03133470022713983;
    w3 = 0.07782754100477543;
    w4 = 0.07964773892720910;
    w5 = 0.02557767565869810;
	 w6 = 0.04328353937728940;
	 (*xtab)[0] = a; (*xtab)[1] = b; (*xtab)[2] = c; (*xtab)[3] = c;
	 (*xtab)[4] = d; (*xtab)[5] = e; (*xtab)[6] = e; (*xtab)[7] = f;
	 (*xtab)[8] = g; (*xtab)[9] = g; (*xtab)[10]= r; (*xtab)[11]= s;
	 (*xtab)[12]= s; (*xtab)[13]= t; (*xtab)[14]= t; (*xtab)[15]= u;
	 (*xtab)[16]= u; (*xtab)[17]= v; (*xtab)[18]= v;
	 (*ytab)[0] = a; (*ytab)[1] = c; (*ytab)[2] = b; (*ytab)[3] = c;
	 (*ytab)[4] = e; (*ytab)[5] = d; (*ytab)[6] = e; (*ytab)[7] = g;
	 (*ytab)[8] = f; (*ytab)[9] = g; (*ytab)[10]= s; (*ytab)[11]= r;
	 (*ytab)[12]= s; (*ytab)[13]= u; (*ytab)[14]= v; (*ytab)[15]= t;
	 (*ytab)[16]= v; (*ytab)[17]= t; (*ytab)[18]= u;
	 (*weight)[0] = w1; (*weight)[1] = w2; (*weight)[2] = w2; (*weight)[3] = w2;
	 (*weight)[4] = w3; (*weight)[5] = w3; (*weight)[6] = w3; (*weight)[7] = w4;
	 (*weight)[8] = w4; (*weight)[9] = w4; (*weight)[10]= w5; (*weight)[11]= w5;
	 (*weight)[12]= w5; (*weight)[13]= w6; (*weight)[14]= w6; (*weight)[15]= w6;
	 (*weight)[16]= w6; (*weight)[17]= w6; (*weight)[18]= w6;
	 break;
  case 19:
	 break;
  case 20:
	 break;
  default:
	 printf("Rule %d is not defined!\n", rule);
	 exit(-1);
  }
}

// this should eventually replace the ugliness of the 16 point fixed direct quad rule
void FlatPanel_getquadrature(FlatPanel panel, unsigned int rule, unsigned int* numquadpoints,
									  Vector3D** quadpoints, Vector3D** quadnormals, Vector* quadweights) {
  real *ksiTab;
  real *etaTab;
  real *weight;

  gen_arb_Stroud_rule(rule, numquadpoints, &ksiTab, &etaTab, &weight);  

  *quadpoints = (Vector3D *)calloc(*numquadpoints, sizeof(Vector3D));
  *quadnormals = (Vector3D *)calloc(*numquadpoints, sizeof(Vector3D));
  *quadweights = Vector_allocate(*numquadpoints);
  unsigned int i;
  for (i = 0; i < *numquadpoints; i++) { // this was shamelessly ripped from FlatPanel_allocate()
	 (*quadpoints)[i] = Vector3D_allocate();
	 (*quadnormals)[i] = Vector3D_allocate();
	 
	 Vector3D_addscaled((*quadpoints)[i], (*quadpoints)[i], ksiTab[i], panel->vertex[0]);
	 Vector3D_addscaled((*quadpoints)[i], (*quadpoints)[i], etaTab[i], panel->vertex[1]);
	 Vector3D_addscaled((*quadpoints)[i], (*quadpoints)[i], 1.0-etaTab[i]-ksiTab[i], panel->vertex[2]);
	 
	 Vector3D_copy((*quadnormals)[i], panel->panelaxis[2]);
	 //	 Vector3D_scale((*quadnormals)[i], -1.0);
	 (*quadweights)[i] = weight[i] * panel->area;
  }
  free(ksiTab);
  free(etaTab);
  free(weight);
}
