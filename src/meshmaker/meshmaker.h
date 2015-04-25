#include "petsc.h"
#define REAL_IS_DOUBLE
typedef PetscReal real;
typedef PetscReal sreal;
typedef double complex complexsreal;

#include "stdio.h"
#include "errno.h"
#include "unistd.h"
#include "ctype.h"
#include "float.h"

/* Replacing #include "FFTSVDpbeAPI.h" */

/*   FFTSVD.h */
typedef enum { CONSTANT_KERNEL, X_KERNEL, Y_KERNEL, Z_KERNEL, POISSON_KERNEL, HELMHOLTZ_KERNEL, DESINGULARIZED_HELMHOLTZ_KERNEL, LJ_KERNEL, LJ12_KERNEL, LJ6_KERNEL, MONOMIAL_KERNEL, INVERSEPOWER_KERNEL, GHOSH_KERNEL, GRYCUK_KERNEL, LESYNG_KERNEL, VOLUME_KERNEL, XGB_KERNEL, JUFFER_KERNEL, BIBEE_P_KERNEL, BIBEE_CFA_KERNEL} BEMKernelType;
typedef enum { SINGLE_LAYER_INT, DOUBLE_LAYER_INT, SINGLE_AND_DOUBLE_LAYER_INT, NORMDERIV_SINGLE_LAYER_INT, NORMDERIV_DOUBLE_LAYER_INT } BEMLayerType;

/*   FFTSVDpbeAPI.h */
#define DEBYE_CONSTANT 3.047  /* For compatibility with DelPhi */
#define KT_CONVERSION 561.0  /* For compatibility with DelPhi */
#define GRID_SIZE 4
#define MAX_PANELS_PER_FINEST_CUBE 32
#define SVD_ERROR 1e-5

extern real ionexclusionradius;
extern real proberadius;

typedef struct {
   char record[7];
   unsigned int atomnumber;
   char atomname[5];
   char alternatelocation;
   char residuename[4];
   char chain;
   unsigned int residuenumber;
   char residueinsertion;
   real x;
   real y;
   real z;
   real occupancy;
   real temperature;
   unsigned int footnotenumber;

   real radius;
   real charge;
   real potential;
} PDBentry;

typedef struct {
   char atomlabel[7];
   char residuelabel[4];
   real radius;
} SIZentry;

typedef struct {
   char atomlabel[7];
   char residuelabel[4];
   unsigned int residuenumber;
   char chain;
   real charge;
} CRGentry;

/*   Vector.h */
typedef real *Vector;
Vector Vector_allocate(unsigned int length);
void Vector_free(Vector vector);
void Vector_copy(Vector vectordest, Vector vectorsrc, unsigned int length);
void Vector_copypiece(Vector vectordest, unsigned int destStart, Vector vectorsrc, unsigned int srcStart, unsigned int length);
void Vector_copyscaledpiece(Vector vectordest, unsigned int destStart, Vector vectorsrc, unsigned int srcStart, unsigned int length, real scalefactor);
void Vector_zero(Vector vector, unsigned int length);
void Vector_scale(Vector vector, real scale, unsigned int length);
void Vector_subtractscaledvector(Vector vector1, real scale, Vector vector2, unsigned int length);
void Vector_addscaledvector(Vector vector1, real scale, Vector vector2, unsigned int length);
void Vector_addvector(Vector vector1, Vector vector2, unsigned int length);
void Vector_subtractvector(Vector vector1, Vector vector2, unsigned int length);
void Vector_pointwisescale(Vector dest, Vector scalefactors, unsigned int length);
real Vector_dot(Vector vector1, Vector vector2, unsigned int length);
real Vector_norm(Vector vector, unsigned int length);
void Vector_normalize(Vector vector, unsigned int length);
void Vector_writefile(char* filename, Vector vector, unsigned int length);
/*   SVector.h */
typedef sreal *SVector;
/*   ComplexSVector.h */
typedef complexsreal *ComplexSVector;
/*   Matrix.h */
typedef real **Matrix;
Matrix Matrix_allocate(unsigned int rows, unsigned int columns);
void Matrix_free(Matrix matrix);
void Matrix_copy(Matrix matrixdest, Matrix matrixsrc, unsigned int rows, unsigned int columns);
void Matrix_copypiece(Matrix matrixdest, unsigned int destStartRow, unsigned int destStartCol,
                      Matrix matrixsrc, unsigned int srcStartRow, unsigned int srcStartCol,
                      unsigned int numrows, unsigned int numcols);
void Matrix_zero(Matrix matrix, unsigned int rows, unsigned int columns);
void Matrix_identity(Matrix matrix, unsigned int rows, unsigned int columns);
void Matrix_scale(Matrix A, real factor, unsigned int rows, unsigned int columns);
void Matrix_add(Matrix B, Matrix A, Matrix X, unsigned int rows, unsigned int columns);
void Matrix_stats(Matrix matrix, unsigned int rows, unsigned int columns);
real Matrix_norm(Matrix matrix, unsigned int rows, unsigned int columns);
void Matrix_multiplyvector(Vector b, Matrix A, Vector x, unsigned int rows, unsigned int columns);
void Matrix_multiplyvector_transpose(Vector b, Matrix A, Vector x, unsigned int rows, unsigned int columns);
void Matrix_multiplymatrix(Matrix B, Matrix A, Matrix X, unsigned int rowsA, unsigned int columnsA, unsigned int columnsX);
void Matrix_multiplymatrix_transpose(Matrix B, Matrix A, Matrix X, unsigned int rowsA, unsigned int columnsA, unsigned int columnsX);
void Matrix_multiplytranspose_matrix(Matrix B, Matrix A, Matrix X, unsigned int rowsA, unsigned int columnsA, unsigned int rowsX);
void Matrix_rowbasis_pivotedmgs(Matrix* U, unsigned int* rowrank, Matrix matrix, unsigned int rows, unsigned int columns, real epsilon);
void Matrix_columnbasis_pivotedmgs(Matrix* V, unsigned int* columnrank, Matrix matrix, unsigned int rows, unsigned int columns, real epsilon);
void Matrix_transpose(Matrix* matrix, unsigned int rows, unsigned int columns);
void Matrix_diff(Matrix X, Matrix A, Matrix B, unsigned int rows, unsigned int columns);
void Matrix_pseudoinverse_droptol(Matrix XI, Matrix X, unsigned int rows, unsigned int columns, real droptol);
void Matrix_pseudoinverse(Matrix XI, Matrix X, unsigned int rows, unsigned int columns);
void Matrix_lsq_solve(Vector x, Matrix A, Vector b, unsigned int rows, unsigned int columns);
void Matrix_columnbasis_check(Matrix V, unsigned int rows, unsigned int columns);
void Matrix_eigendecomposition(Matrix A, Matrix V, real* d, unsigned int rows, unsigned int columns);
void Matrix_writefile(char* filename, Matrix matrix, unsigned int rows, unsigned int columns);
void Matrix_readfile(Matrix matrix, FILE* file, unsigned int rows, unsigned int columns);
void Matrix_print(Matrix matrix, unsigned int rows, unsigned int columns);
void Matrix_writebinary(FILE* file, Matrix matrix, unsigned int rows, unsigned int columns);
void Matrix_readbinary(Matrix matrix, FILE* file, unsigned int rows, unsigned int columns);
/*   SMatrix.h */
typedef sreal **SMatrix;

/*   Vector3D.h */
typedef struct _Vector3D *Vector3D;
struct _Vector3D {
   real x, y, z;
};

Vector3D Vector3D_allocate();
void Vector3D_free(Vector3D v);
void Vector3D_copy(Vector3D vectordest, Vector3D vectorsrc);
void Vector3D_add(Vector3D sum, Vector3D v1, Vector3D v2);
void Vector3D_sub(Vector3D diff, Vector3D v1, Vector3D v2);
void Vector3D_scale(Vector3D v, real scale);
void Vector3D_cross(Vector3D cross, Vector3D v1, Vector3D v2);
real Vector3D_dot(Vector3D v1, Vector3D v2);
real Vector3D_length(Vector3D v);
void Vector3D_normalize(Vector3D v);
real Vector3D_distance(Vector3D v1, Vector3D v2);
unsigned int Vector3D_equal(Vector3D v1, Vector3D v2);
void Vector3D_transform(Vector3D dest, real *A, Vector3D src);
void Vector3D_addscaled(Vector3D sum, Vector3D v1, real scale, Vector3D v2);
void Vector3D_transformVecs(Vector3D Ax, Vector3D A1, Vector3D A2, Vector3D A3,Vector3D x);
void Vector3D_transformVecs_inverse(Vector3D Ax, Vector3D A1, Vector3D A2, Vector3D A3, Vector3D x);
void Vector3D_print(Vector3D v);

/*   QuadratureRule.h */
#define MAX_QUADRATURE_ORDER 1024
typedef struct _QuadratureRule *QuadratureRule;
struct _QuadratureRule {
  /* basic integration */
  unsigned int order;
  real *x;
  real *w;
                                                                                
  /* 7 point quadrature */
  real dNdKsi[6][7];
  real dNdEta[6][7];
  real SPW[7];
};

QuadratureRule QuadratureRule_allocate(unsigned int order);
void QuadratureRule_free(QuadratureRule qr);
void QuadratureRule_generatepoints(QuadratureRule qr, real lower, real upper);

/*   FlatPanel.h */
typedef struct _FlatPanel *FlatPanel;
struct _FlatPanel {
   Vector3D vertex[3];  /* Vertex coordinates in global frame */
   Vector3D centroid;  /* Centroid coordinates in global frame */
   real edgelength[3];  /* Edge lengths */
   Vector3D panelaxis[3];   /* Panel reference frame */
   Vector3D panelvertex[3];  /* Vertex coordinates in panel reference frame */
   real contributionC[3], contributionS[3];  /* cos/sin contribution terms */

   // until we figure out new formulas for contribution[CS]
   Vector3D panelaxisnum[3];
   Vector3D panelvertexnum[3];
   Vector3D edges[3];   /* equations for edges in panel ref frame JPB 11/20/04*/
   real edgeRHS[3]; /* added by JPB for more stable panel transforms */
   real edgeOV[3]; /* vector added by JPB--tells us whether opp. vert is inside*/

   real area;  /* Panel area */
   real moments[16];
   real max_diag;
   real min_diag;

   unsigned int numdirectquadpoints;
   Vector3D *directquadpoints;
   Vector3D *directquadnormals;
   real *directquadweights;
};

FlatPanel FlatPanel_allocate(Vector3D v1, Vector3D v2, Vector3D v3);
void FlatPanel_free(FlatPanel panel);
void FlatPanel_moments(FlatPanel panel);
unsigned int FlatPanel_memory(unsigned int order, unsigned int doublelayer);
void FlatPanel_getquadrature(FlatPanel panel, unsigned int rule, unsigned int* numquadpoints,
									  Vector3D** quadpoints, Vector3D** quadnormals, Vector* quadweights);

/*   FlatIntegration.h */
void FlatIntegration_oneoverr_grad(Vector3D point, FlatPanel panel, void* parameters, real* slp, real* dlp);
void FlatIntegration_oneoverr_deriv(Vector3D point, FlatPanel panel, void* parameters, real* slp, real* dlp);
void FlatIntegration_oneoverr_deriv_numerical_qual(Vector3D point, FlatPanel panel, void* parameters, real* slp, real* dlp);
void FlatIntegration_oneoverr_deriv_qual(FlatPanel point, FlatPanel panel, void* parameters, real* slp, real* dlp);
void FlatIntegration_oneoverr_deriv_qual_point(Vector3D point, FlatPanel panel, void* parameters, real* slp, real* dlp);

void FlatIntegration_oneoverr(Vector3D point, FlatPanel panel, void* parameters, real* slp, real* dlp);
void FlatIntegration_oneoverr_numerical(Vector3D point, FlatPanel panel, void* parameters, real* slp, real* dlp);
void FlatIntegration_ekroverr_numerical(Vector3D point, FlatPanel panel, void* parameters, real* slp, real* dlp);
void FlatIntegration_ekroverr_deriv_numerical(Vector3D point, FlatPanel panel, Vector3D normal, void* parameters, real* slp, real* dlp);
void FlatIntegration_ekroverr_desingularized(Vector3D point, FlatPanel panel, void* parameters, real* slp, real* dlp);
void FlatIntegration_LJ(Vector3D point, FlatPanel panel, void* parameters, real *slp, real *dlp);
void FlatIntegration_LJ12(Vector3D point, FlatPanel panel, void* parameters, real *slp, real *dlp);
void FlatIntegration_LJ6(Vector3D point, FlatPanel panel, void* parameters, real *slp, real *dlp);
void FlatIntegration_Ghosh(Vector3D point, FlatPanel panel, void* parameters, real *slp, real *dlp);
void FlatIntegration_Lesyng(Vector3D point, FlatPanel panel, void* parameters, real *slp, real *dlp);
void FlatIntegration_Grycuk(Vector3D point, FlatPanel panel, void* parameters, real *slp, real *dlp);
void FlatIntegration_general(Vector3D point, FlatPanel panel, BEMKernelType kernel, BEMLayerType layer, void* parameters, real* integral);

real FlatIntegration_maltquad(Vector3D point, FlatPanel panel, BEMKernelType kerneltype, void* parameters, BEMLayerType layertype);

void gen_Stroud_Rule(real *xtab, real *ytab, real *wtab);

/*   GST.h */
typedef struct _FlatRefPanel *FlatRefPanel;
struct _FlatRefPanel {
  Vector3D centroid;
  Vector3D normal;
  real rhs;
  Vector3D vertices[3];
  unsigned int numquadpoints;
  real *ksi;
  real *eta;
  real *weights;
  real maxEdgelength;
};

typedef struct _Conic *Conic;
struct _Conic {
  real asquared;
  real bsquared;
  real theta1Int;
  real theta2Int;
  Vector3D Tmat[3];
  Vector3D Tvec;
};

typedef struct _Curve *Curve;
struct _Curve {
  Vector3D ac;
  Vector3D v1;
  Vector3D v2;
  Vector3D X; // points from ac -> v1
  Vector3D Y;

  real radius;
  real endTheta;
};

typedef struct _GST *GST;
struct _GST {
  Vector3D center;
  Vector3D vertices[3];
  Vector3D arccenters[3];
  int pit, cavity;
  real radius;
  real area;

  Vector3D Tmat[3];
  Vector3D Tvec;

  Vector3D centroid; // on GST
  Vector3D normalAtCentroid; 
  FlatRefPanel flatpanel;
  int factor[3];  // add or subtract or zero
  Conic conic[3];

  unsigned int numdirectquadpoints;
  Vector3D *directquadpoints;
  Vector3D *directquadnormals;
  real *directquadweights;
  
  Matrix vandermondeInverse;
};

FlatRefPanel FlatRefPanel_allocate(Vector3D v1, Vector3D v2, Vector3D v3,
									  Vector3D centroid, Vector3D normal,
									  real rhs, unsigned int numquadpoints,
									  real *ksi, real *eta, real *weights,
									  real maxEdgelength);
void FlatRefPanel_free(FlatRefPanel fp);
Conic Conic_allocate(real asq, real bsq, real theta1Int, real theta2Int,
							Vector3D *Tmat, Vector3D Tvec);
void Conic_free(Conic c);
Curve Curve_allocate(Vector3D ac, Vector3D v1, Vector3D v2);
void Curve_getParamPoint(Vector3D p, Vector3D dpdalpha, Curve c, real alpha);
void Curve_free(Curve c);
void PermuteVectors(Vector3D v0, Vector3D v1, Vector3D v2);

GST GST_allocate(Vector3D center, real radius, int caporpit, int cavity,
					  Vector3D v1, Vector3D v2, Vector3D v3,
					  Vector3D ac1, Vector3D ac2, Vector3D ac3);
void GST_free(GST gst);

void GST_readfile(unsigned int *numGSTpanels, GST** gstList, FILE* file, unsigned int ignorecav);
void findSpherePoint(Vector3D sphvert, Vector3D center,
							real radius, Vector3D oldvert);
void initialize_GST_FlatRefPanel(GST gst);
void generateQuadraturePoints(Vector3D* vertices, unsigned int *numquadpoints, real **ksi, real **eta, real **weights);
void figure_out_correction_sign(GST gst, unsigned int index);
void conic_section_identify(GST gst, unsigned int index);
void projectToPlane(Vector3D onplane, Vector3D start, Vector3D end);
void GST_get_conic(real* Vret, Vector3D *M, unsigned int numpoints);
void GST_get_direct_quadPoint(Vector3D point, Vector3D normal, real *detJ,
										real ksi, real eta, GST gst,
										Curve *edges, Vector3D *Tmat, Vector3D Tvec);

void GST_getStandardPosition(GST gst, Curve *edges, Vector3D *Tmat, Vector3D Tvec);
void GST_sphGSTtoCart(real radius, real theta, real phi, Vector3D point, Vector3D normal, Vector3D dPoint_dTheta, Vector3D dPoint_dPhi);
void GST_getCircleArcIntersection(Vector3D center, Vector3D normal, real radius,
											 Curve curve, Vector3D intPoint,
											 real *alpha, Vector3D dIntPoint_dAlpha);
void GST_getCentroid(GST gst);

/*   TOR.h */
typedef struct _TOR *TOR;
struct _TOR {
  Vector3D centroid;
  Vector3D normalAtCentroid;
  int cavity;
  real area;
  real maxEdgelength;
  real c, a;
  real startPsi, endPsi;
  real startTheta, endTheta;
  Vector3D center; // torus center
  Vector3D normal; // normal of donut plane
  Vector3D localX, localY;
  Vector3D Tmat[3];
  Vector3D Tvec;
  int isLocal; // should be deprecated soon...
  unsigned int thetaIndex, psiIndex;
  unsigned int numthetadivs, numpsidivs;
  unsigned int torusID;

  unsigned int numdirectquadpoints;
  Vector3D *directquadpoints;
  Vector3D *directquadnormals;
  real *directquadweights;
};

TOR TOR_allocate(real torusRadius, real probeRadius, real startTheta, real endTheta,
					  real startPsi, real endPsi, Vector3D center, Vector3D normal,
					  Vector3D localX, Vector3D localY, unsigned int thetaIndex, unsigned int psiIndex, int cavity,
					  unsigned int numthetadivs, unsigned int numpsidivs, unsigned int torusID);
void TOR_free(TOR t);
void TOR_readfile(unsigned int *numTORpanels, TOR** torList, FILE* file, unsigned int ignorecav);
void torToCart(Vector3D point, Vector3D normal,
					real c, real a, real theta, real psi);

/*   Panel.h */
typedef enum { FLAT_PANEL, GST_PANEL, TOR_PANEL, POINT_PANEL } PanelType;
typedef struct _Panel *Panel;
struct _Panel {
   Vector3D centroid;  /* Centroid coordinates in global frame */
   Vector3D normal;  /* Normal vector at centroid */
   real area;  /* Panel area */
   real maxedgelength;

   unsigned int numdirectquadpoints;
   Vector3D* directquadpoints;
   Vector3D* directquadnormals;
   Vector directquadweights;

#ifdef ENABLE_GALERKIN
   unsigned int numgalerkinquadpoints;
   Vector3D* galerkinquadpoints;
   Vector3D* galerkinquadnormals;
   Vector galerkinquadweights;
#endif

   PanelType type;
   void* realpanel;
};

Panel Panel_allocate();
void Panel_free(Panel);
void Panel_FlatPanel(Panel panel, FlatPanel flatpanel);
void Panel_GST(Panel Panel, GST gst, unsigned int usecom);
void Panel_TOR(Panel panel, TOR tor, unsigned int usecom);
void Panel_Vector3D(Panel panel, Vector3D v3d);
real Panel_memory(Panel panel, BEMLayerType layertype);
real Panel_quadrature(Vector3D point, Panel panel, BEMKernelType kerneltype, void* parameters, BEMLayerType layertype);

/*   GreensFunction.h */
real GreensFunction(Vector3D dest, Vector3D src, BEMKernelType kerneltype, void* parameters);
real GreensFunction_deriv(Vector3D dest, Vector3D src, BEMKernelType kerneltype, void* parameters, Vector3D direction);
real GreensFunction_doublederiv(Vector3D dest, Vector3D src, BEMKernelType kerneltype,
										  void* parameters, Vector3D direction, Vector3D direction2);

real GreensFunction_oneoverr(Vector3D dest, Vector3D src, void* parameters);
real GreensFunction_oneoverr_deriv(Vector3D dest, Vector3D src, void* parameters, Vector3D direction);
real GreensFunction_ekroverr(Vector3D dest, Vector3D src, void* parameters);
real GreensFunction_ekroverr_deriv(Vector3D dest, Vector3D src, void* parameters, Vector3D direction);
real GreensFunction_ekroverr_desingularized(Vector3D dest, Vector3D src, void* parameters);
real GreensFunction_ekroverr_desingularized_deriv(Vector3D dest, Vector3D src, void* parameters, Vector3D direction);
real GreensFunction_LJ(Vector3D dest, Vector3D src, void* parameters);
real GreensFunction_LJ12(Vector3D dest, Vector3D src, void* parameters);
real GreensFunction_LJ6(Vector3D dest, Vector3D src, void* parameters);
real GreensFunction_LJ12_deriv(Vector3D dest, Vector3D src, void* parameters, Vector3D direction);
real GreensFunction_LJ6_deriv(Vector3D dest, Vector3D src, void* parameters, Vector3D direction);
real GreensFunction_monomial(Vector3D dest, Vector3D src, void* parameters);
real GreensFunction_monomial_deriv(Vector3D dest, Vector3D src, void* parameters, Vector3D direction);
real GreensFunction_inversepower(Vector3D dest, Vector3D src, void* parameters);
real GreensFunction_inversepower_deriv(Vector3D dest, Vector3D src, void* parameters, Vector3D direction);

// generalized-born kernels
real GreensFunction_Ghosh(Vector3D dest, Vector3D src, void* parameters);
real GreensFunction_Ghosh_deriv(Vector3D dest, Vector3D src, void* parameters, Vector3D direction);
real GreensFunction_XGB(Vector3D dest, Vector3D src, void* parameters);
real GreensFunction_XGB_deriv(Vector3D dest, Vector3D src, void* parameters, Vector3D direction);
real GreensFunction_Lesyng(Vector3D dest, Vector3D src, void* parameters);
real GreensFunction_Lesyng_deriv(Vector3D dest, Vector3D src, void* parameters, Vector3D direction);
real GreensFunction_Grycuk(Vector3D dest, Vector3D src, void* parameters);
real GreensFunction_Grycuk_deriv(Vector3D dest, Vector3D src, void* parameters, Vector3D direction);

/*   calcpc_GST.h */
void findIntersection(Vector3D point, Vector3D v1, Vector3D v2, real theta);
void Integration_general_GST_single(Vector3D point, GST gst, BEMKernelType kernel, void* parameters, real* slp);
void Integration_general_GST_double(Vector3D point, GST gst, BEMKernelType kernel, void* parameters, real* dlp);
void Integration_oneoverr_GST(Vector3D point, GST gst, void* parameters, real* slp, real *dlp);
void Integration_oneoverr_GST_single(Vector3D point, GST gst, void* parameters, real* slp);
void Integration_oneoverr_GST_double(Vector3D point, GST gst, void* parameters, real* dlp);
void greenInt_XinWang(real *Integrals, FlatRefPanel flatpanel, Vector3D localpnt, unsigned int order);
void Integration_general_GST_single(Vector3D point, GST gst, BEMKernelType kernel, void* parameters, real *slp);
void GST_computePolynomialIntegrals();
void GST_computeVandermondeMatrix();
void GST_computePolynomialValues(); 
real findJacobianDet(Vector3D center, real radius, Vector3D pnt);
void ellipse_int2(real *ellipse_integrals, GST gst, Vector3D point, unsigned int numcoeffs,
						unsigned int* kpvec, unsigned int* epvec,
						BEMKernelType kernel, BEMLayerType layertype, void* parameters);
void conic_section_integrate(real *ellipse_integrals, unsigned int arcnumber,
									  GST gst, Vector3D point, unsigned int numcoeffs,
									  unsigned int* kpvec, unsigned int* epvec,
									  	BEMKernelType kernel, BEMLayerType layertype, void *parameters);

/*   calcpc_TOR.h */
real realabs(real a);
real callGreensFunction(Vector3D obspnt, Vector3D srcpnt, Vector3D normal,
								BEMKernelType kernel, void *parameters);
void Integration_general_TOR_single(Vector3D point, TOR tor, BEMKernelType kernel, void* parameters, real* slp);
void Integration_oneoverr_TOR(Vector3D point, TOR tor, void* parameters, real* slp, real* dlp);
void Integration_oneoverr_TOR_single(Vector3D point, TOR tor, void* parameters, real* slp);
void Integration_oneoverr_TOR_double(Vector3D point, TOR tor, void* parameters, real* dlp);
void Integration_oneoverr_TOR_double_DJW(Vector3D point, TOR tor, void* parameters, real* dlp);
void Integration_oneoverr_TOR_double_local(Vector3D point, TOR tor, void* parameters, real* dlp);
void Integration_general_TOR_single_local(Vector3D point, TOR tor,
														BEMKernelType kernel, void* parameters, real* slp);
void Integration_general_TOR_single_self(Vector3D point, TOR tor,
													  BEMKernelType kernel, void* parameters, real* slp);
void Integration_general_TOR_single_near(Vector3D point, TOR tor,
													  BEMKernelType kernel, void* parameters, real* slp);

extern unsigned int scount;
extern unsigned int ncount;
extern unsigned int fcount;

/*   Integration.h */
real Integration(Vector3D point, Panel panel, BEMKernelType kerneltype, void* parameters, BEMLayerType layertype);

/*   VertFace.h */
typedef struct _VertFace *VertFace;
struct _VertFace {
   unsigned int numvertices, numfaces;
   Vector3D* vertices;
   unsigned int* facesx;
   unsigned int* facesy;
   unsigned int* facesz;
   unsigned int* facesgenus;
} _VertFace;

VertFace VertFace_allocate();
void VertFace_free(VertFace vf);
void VertFace_readvert(VertFace vf, FILE* file);
void VertFace_readface(VertFace vf, FILE* file);
void VertFace_readface_flip(VertFace vf, FILE* file);
void VertFace_fix(VertFace vf, unsigned int accessible);
void VertFace_getpanels(VertFace vf, Panel* panels);

/*   config.h */
#define BINARY_DIR "."
#define TRUE 1

/* meshutil.h */
void error(char* message, ...);
void readPDB(const char* filename, unsigned int* numPDBentries, PDBentry** PDBentries);
void readCRD(const char* filename, unsigned int* numPDBentries, PDBentry** PDBentries);
void readXYZR(const char* filename, unsigned int* numPDBentries, PDBentry** PDBentries);
void readSIZ(const char* filename, unsigned int* numSIZentries, SIZentry** SIZentries);
void readCRG(const char* filename, unsigned int* numCRGentries, CRGentry** CRGentries);
void assignRadiiCharges(PDBentry* PDBentries, unsigned int numPDBentries, SIZentry* SIZentries, unsigned int numSIZentries, CRGentry* CRGentries, unsigned int numCRGentries);

/* Done Replacing #include "FFTSVDpbeAPI.h" */
