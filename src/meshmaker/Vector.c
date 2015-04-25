#include "meshmaker.h"

/* Constructors and Destructors */

Vector Vector_allocate(unsigned int length) {
   return (Vector)calloc(length, sizeof(real));
}

void Vector_free(Vector vector) {
   free(vector);
}

/* Initialization and Copying */

void Vector_copy(Vector vectordest, Vector vectorsrc, unsigned int length) {
   unsigned int i;
   
#pragma ivdep  
   for (i = 0; i < length; i++)
      vectordest[i] = vectorsrc[i];
}

void Vector_copypiece(Vector vectordest, unsigned int destStart, Vector vectorsrc, unsigned int srcStart,unsigned int length) {
   unsigned int i;
#pragma ivdep
   for (i = 0; i < length; i++)
      vectordest[destStart + i] = vectorsrc[srcStart + i];
}

void Vector_copyscaledpiece(Vector vectordest, unsigned int destStart, Vector vectorsrc, unsigned int srcStart,unsigned int length, real scalefactor) {
   unsigned int i;
#pragma ivdep
   for (i = 0; i < length; i++)
      vectordest[destStart + i] = scalefactor * vectorsrc[srcStart + i];
}

void Vector_zero(Vector vector, unsigned int length) {
   unsigned int i;
   
#pragma ivdep  
   for (i = 0; i < length; i++)
      vector[i] = (real)0;
}

/* Arithmetic Operations */

void Vector_scale(Vector vector, real scale, unsigned int length) {
   unsigned int i;
   
#pragma ivdep  
   for (i = 0; i < length; i++)
      vector[i] *= scale;
}

void Vector_subtractscaledvector(Vector vector1, real scale, Vector vector2, unsigned int length) {
   unsigned int i;
   
#pragma ivdep  
   for (i = 0; i < length; i++)
      vector1[i] -= scale * vector2[i];
}

void Vector_addscaledvector(Vector vector1, real scale, Vector vector2, unsigned int length) {
   unsigned int i;
   
#pragma ivdep  
   for (i = 0; i < length; i++)
      vector1[i] += scale * vector2[i];
}

void Vector_addvector(Vector vector1, Vector vector2, unsigned int length) {
   unsigned int i;
 
#pragma ivdep  
   for (i = 0; i < length; i++)
      vector1[i] += vector2[i];
}

void Vector_subtractvector(Vector vector1, Vector vector2, unsigned int length) {
   unsigned int i;
 
#pragma ivdep  
   for (i = 0; i < length; i++)
      vector1[i] -= vector2[i];
}

/* Linear Algebra Operations */

real Vector_dot(Vector vector1, Vector vector2, unsigned int length) {
   unsigned int i;
   real dot = (real)0;
   
#pragma ivdep  
   for (i = 0; i < length; i++)
      dot += vector1[i] * vector2[i];
      
   return dot;
}

real Vector_norm(Vector vector, unsigned int length) {
   unsigned int i;
   real norm = (real)0;
   
#pragma ivdep  
   for (i = 0; i < length; i++)
      norm += vector[i] * vector[i];
   
   norm = sqrt(norm);
   
   return norm;
}

void Vector_normalize(Vector vector, unsigned int length) {
   Vector_scale(vector, (real)1 / Vector_norm(vector, length), length);
}

void Vector_writefile(char *filename, Vector vector, unsigned int length) {
  unsigned int i;
  FILE *file = NULL;
  file = fopen(filename, "w");
  for (i = 0; i < length; i++)
	 fprintf(file, "%20.15e\n", vector[i]);
  fclose(file);
}

void Vector_pointwisescale(Vector dest, Vector scalefactors, unsigned int length) {
  unsigned int i;
  for (i = 0; i < length; i++) {
	 dest[i] *= scalefactors[i];
  }
}
