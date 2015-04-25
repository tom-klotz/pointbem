#include "meshmaker.h"

/* Constructors and Destructors */

void Vector3D_print(Vector3D v)
{
   printf("v->x = %f,v->y = %f, v->z = %f\n",v->x, v->y, v->z);
}

Vector3D Vector3D_allocate() {
   return (Vector3D)calloc(1, sizeof(struct _Vector3D));
}
 
void Vector3D_free(Vector3D v) {
   free(v);
}

/* Initialization and Copying */
 
void Vector3D_copy(Vector3D vdest, Vector3D vsrc) {
   vdest->x = vsrc->x;
   vdest->y = vsrc->y;
   vdest->z = vsrc->z;
}
  
/* Arithmetic Operations */

void Vector3D_add(Vector3D sum, Vector3D v1, Vector3D v2) {
   sum->x = v1->x + v2->x;
   sum->y = v1->y + v2->y;
   sum->z = v1->z + v2->z;
}

void Vector3D_sub(Vector3D diff, Vector3D v1, Vector3D v2) {
   diff->x = v1->x - v2->x;
   diff->y = v1->y - v2->y;
   diff->z = v1->z - v2->z;
}

void Vector3D_scale(Vector3D v, real scale) {
   v->x *= scale;
   v->y *= scale;
   v->z *= scale;
}

void Vector3D_cross(Vector3D cross, Vector3D v1, Vector3D v2) {
   cross->x = v1->y*v2->z - v1->z*v2->y;
   cross->y = v1->z*v2->x - v1->x*v2->z;
   cross->z = v1->x*v2->y - v1->y*v2->x;
}

real Vector3D_dot(Vector3D v1, Vector3D v2) {
   return (v1->x*v2->x + v1->y*v2->y + v1->z*v2->z);
}

real Vector3D_length(Vector3D v) {
   return sqrt(Vector3D_dot(v, v));
}

void Vector3D_normalize(Vector3D v) {
   Vector3D_scale(v, (real)1 / Vector3D_length(v));
}

real Vector3D_distance(Vector3D v1, Vector3D v2) {
   return sqrt((v1->x - v2->x) * (v1->x - v2->x) +
               (v1->y - v2->y) * (v1->y - v2->y) +
               (v1->z - v2->z) * (v1->z - v2->z));
}

unsigned int Vector3D_equal(Vector3D v1, Vector3D v2) {
   return ((v1->x == v2->x) && (v1->y == v2->y) && (v1->z == v2->z));
}

void Vector3D_transform(Vector3D dest, real *A, Vector3D src) {
   dest->x = A[0] * src->x + A[1] * src->y + A[2] * src->z;
   dest->y = A[3] * src->x + A[4] * src->y + A[5] * src->z;
   dest->z = A[6] * src->x + A[7] * src->y + A[8] * src->z;
}

void Vector3D_addscaled(Vector3D sum, Vector3D v1, real scale, Vector3D v2) {
   sum->x = v1->x + scale * v2->x;
   sum->y = v1->y + scale * v2->y;
   sum->z = v1->z + scale * v2->z;
}

void Vector3D_transformVecs(Vector3D Ax, Vector3D A1, Vector3D A2, Vector3D A3, Vector3D x) {
   real A[9];
   Vector3D dest = Vector3D_allocate();
   A[0] = A1->x; A[1] = A2->x; A[2] = A3->x;
   A[3] = A1->y; A[4] = A2->y; A[5] = A3->y;
   A[6] = A1->z; A[7] = A2->z; A[8] = A3->z;
   Vector3D_transform(dest, A, x);
   Vector3D_copy(Ax, dest);
   Vector3D_free(dest);
}

void Vector3D_transformVecs_inverse(Vector3D Ax, Vector3D A1, Vector3D A2, Vector3D A3, Vector3D x) {
   real A[9];
   Vector3D dest = Vector3D_allocate();
   A[0] = A1->x; A[3] = A2->x; A[6] = A3->x;
   A[1] = A1->y; A[4] = A2->y; A[7] = A3->y;
   A[2] = A1->z; A[5] = A2->z; A[8] = A3->z;
   Vector3D_transform(dest, A, x);
   Vector3D_copy(Ax, dest);
   Vector3D_free(dest);
}
