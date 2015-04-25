#include "meshmaker.h"
#include "petscblaslapack.h"

/* Constructors and Destructors */

Matrix Matrix_allocate(unsigned int rows, unsigned int columns) {
   unsigned int i;
   Matrix matrix = NULL;
   
   matrix = (Matrix)calloc(rows, sizeof(Vector));
   matrix[0] = (Vector)calloc(rows * columns, sizeof(real));
   
#pragma ivdep
   for (i = 1; i < rows; i++)
      matrix[i] = matrix[0] + i * columns;
   
   return matrix;
}

void Matrix_free(Matrix matrix) {
   if (matrix != NULL)
      free(matrix[0]);
   free(matrix);
}

/* Initialization and Copying */

// B = A + X;  B can be either A or X (or both to double heh)
void Matrix_add(Matrix B, Matrix A, Matrix X, unsigned int rows, unsigned int columns) {
  unsigned int i;
  Vector pB = B[0];
  Vector pA = A[0];
  Vector pX = X[0];
#pragma ivdep
  for (i = 0; i < rows * columns; i++)
	 pB[i] = pA[i] + pX[i];
}

void Matrix_scale(Matrix A, real factor, unsigned int rows, unsigned int columns) {
  unsigned int i;
  Vector pA = A[0];
#pragma ivdep
  for (i = 0; i < rows * columns; i++)
	 pA[i] *= factor;
}

void Matrix_copy(Matrix matrixdest, Matrix matrixsrc, unsigned int rows, unsigned int columns) {
   unsigned int i;   
   Vector pdest = matrixdest[0];
   Vector psrc = matrixsrc[0];
   
#pragma ivdep
   for (i = 0; i < rows * columns; i++)
      pdest[i] = psrc[i];
}

// this function is horribly inefficient but is only needed in a few
// minor cases, once or twice per run.
void Matrix_copypiece(Matrix matrixdest, unsigned int destStartRow, unsigned int destStartCol,
                      Matrix matrixsrc, unsigned int srcStartRow, unsigned int srcStartCol,
                      unsigned int numrows, unsigned int numcols) {
   unsigned int i, j;
   for (i = 0; i < numrows; i++) {
      for (j = 0; j < numcols; j++) {
         matrixdest[destStartRow+i][destStartCol+j] = matrixsrc[srcStartRow+i][srcStartCol+j];
      }
   }
}

void Matrix_zero(Matrix matrix, unsigned int rows, unsigned int columns) {
   unsigned int i;   
   Vector pmatrix = matrix[0];
   
#pragma ivdep
   for (i = 0; i < rows * columns; i++)
      pmatrix[i] = (real)0;
}

void Matrix_identity(Matrix matrix, unsigned int rows, unsigned int columns) {
   unsigned int i;
   
   #pragma ivdep
   for (i = 0; i < PetscMin(rows, columns); i++)
      matrix[i][i] = (real)1;
}

/* Linear Algebra Operations */

void Matrix_stats(Matrix matrix, unsigned int rows, unsigned int columns) {
   unsigned int i;   
   Vector pmatrix = matrix[0];
   real min = 1e10, max = -1e10;
   
   #pragma ivdep
   for (i = 0; i < rows * columns; i++) {
      if (pmatrix[i] > max) max = pmatrix[i];
      if (pmatrix[i] < min) min = pmatrix[i];
   }

   printf("MIN: %e MAX: %e\n", min, max);
}

/* Matrix Frobenius norm */

real Matrix_norm(Matrix matrix, unsigned int rows, unsigned int columns) {
   unsigned int i;
   real norm = (real)0;
   Matrix temp = Matrix_allocate(columns, columns);
   
   Matrix_multiplymatrix_transpose(temp, matrix, matrix, rows, columns, columns);
   
#pragma ivdep
   for (i = 0; i < columns; i++)
      norm += temp[i][i];
   
   norm = sqrt(norm);
   
   Matrix_free(temp);
   
   return norm;
}

/* b = A*x */
void Matrix_multiplyvector(Vector b, Matrix A, Vector x, unsigned int rows, unsigned int columns) {
   unsigned int i, j;

   for (i = 0; i < rows; i++) {
      Vector Ai = A[i];
#pragma ivdep
      for (j = 0; j < columns; j++)
         b[i] += Ai[j] * x[j];
   }
}

/* b = A'*x */
void Matrix_multiplyvector_transpose(Vector b, Matrix A, Vector x, unsigned int rows, unsigned int columns) {
   unsigned int i, j;

   for (i = 0; i < rows; i++) {
      Vector Ai = A[i];
      real xi = x[i];
#pragma ivdep
      for (j = 0; j < columns; j++)
         b[j] += Ai[j] * xi;
   }
}

/* B = A*X */
void Matrix_multiplymatrix(Matrix B, Matrix A, Matrix X, unsigned int rowsA, unsigned int columnsA, unsigned int columnsX) {
   unsigned int i, j, k;
   
   Matrix_zero(B, rowsA, columnsX);
   
   for (i = 0; i < rowsA; i++) {
      Vector Bi = B[i];
      for (j = 0; j < columnsA; j++) {
         real Aij = A[i][j];
         Vector Xj = X[j];
#pragma ivdep
         for (k = 0; k < columnsX; k++)
            Bi[k] += Aij * Xj[k];
      }
   }
}

/* B = A'*X */
void Matrix_multiplymatrix_transpose(Matrix B, Matrix A, Matrix X, unsigned int rowsA, unsigned int columnsA, unsigned int columnsX) {
   unsigned int i, j, k;
   
   Matrix_zero(B, columnsA, columnsX);
   
   for (i = 0; i < columnsA; i++) {
      Vector Bi = B[i];
      for (j = 0; j < rowsA; j++) {
         real Aji = A[j][i];
         Vector Xj = X[j];
#pragma ivdep
         for (k = 0; k < columnsX; k++)
            Bi[k] += Aji * Xj[k];
      }
   }
}

/* B = A*X' */
void Matrix_multiplytranspose_matrix(Matrix B, Matrix A, Matrix X, unsigned int rowsA, unsigned int columnsA, unsigned int rowsX) {
   unsigned int i, j, k;
   
   Matrix_zero(B, rowsA, rowsX);
   
   for (i = 0; i < rowsA; i++) {
      Vector Bi = B[i];
      for (j = 0; j < columnsA; j++) {
         real Aij = A[i][j];
#pragma ivdep
         for (k = 0; k < rowsX; k++)
            Bi[k] += Aij * X[k][j];
      }
   }
}

void Matrix_rowbasis_pivotedmgs(Matrix* U, unsigned int* rowrank, Matrix matrix, unsigned int rows, unsigned int columns, real epsilon) {
   Matrix transpose = Matrix_allocate(rows, columns);
   
   Matrix_copy(transpose, matrix, rows, columns);
   
   Matrix_transpose(&transpose, rows, columns);
   
   Matrix_columnbasis_pivotedmgs(U, rowrank, transpose, columns, rows, epsilon);
   
   Matrix_transpose(U, columns, *rowrank);
   
   Matrix_free(transpose);
}

void Matrix_columnbasis_pivotedmgs(Matrix* V, unsigned int* columnrank, Matrix matrix, unsigned int rows, unsigned int columns, real epsilon) {
   unsigned int i, j, k, l, reorth;
   Matrix A = Matrix_allocate(rows, columns);
   Matrix Q = Matrix_allocate(rows, columns);
   Matrix R = Matrix_allocate(rows, columns);
   real matrixnorm = Matrix_norm(matrix, rows, columns);

   Matrix_copy(A, matrix, rows, columns);
   
   for (k = 0; k < columns; k++) {
      real currentnorm, maxnorm = (real)0;
      unsigned int largestnormcolumn = 0;

      for (l = k; l < columns; l++) {
         currentnorm = (real)0;
         
         for (i = 0; i < rows; i++)
            currentnorm += A[i][l] * A[i][l];
         currentnorm = sqrt(currentnorm);
         
         if (currentnorm > maxnorm) {
            maxnorm = currentnorm;
            largestnormcolumn = l;
         }
      }
      
      if ((k >= 1) && (maxnorm < (epsilon * matrixnorm))) {
         *columnrank = k;
         break;
      }
      
      for (i = 0; i < rows; i++) {
         real temp = A[i][largestnormcolumn];
         A[i][largestnormcolumn] = A[i][k];
         A[i][k] = temp;
      }

      for (i = 0; i < rows; i++) {
         real temp = R[i][largestnormcolumn];
         R[i][largestnormcolumn] = R[i][k];
         R[i][k] = temp;
      }

      R[k][k] = maxnorm;
      
      for (i = 0; i < rows; i++)
         Q[i][k] = A[i][k] / maxnorm;
         
      if (k < (columns - 1)) {
         for (reorth = 0; reorth < 2; reorth++) {
            Vector Rk = R[k];
      
            for (j = k + 1; j < columns; j++)
               Rk[j] = (real)0;

            for (i = 0; i < rows; i++) {
               real Qik = Q[i][k];
               Vector Ai = A[i];
#pragma ivdep
               for (j = k + 1; j < columns; j++)
                  Rk[j] += Qik * Ai[j];
            }
         
            for (i = 0; i < rows; i++) {
               real Qik = Q[i][k];
               Vector Ai = A[i];
#pragma ivdep
               for (j = k + 1; j < columns; j++)
                  Ai[j] -= Qik * Rk[j];
            }
         }
      }
      
      if (k == (columns - 1))
         *columnrank = columns;
   }

   *V = Matrix_allocate(rows, *columnrank);
   
   for (i = 0; i < rows; i++) {
      Vector Vi = (*V)[i];
      Vector Qi = Q[i];
      unsigned int r = *columnrank;
#pragma ivdep
      for (j = 0; j < r; j++)
         Vi[j] = Qi[j];
   }

   //Matrix_columnbasis_check(*V, rows, *columnrank);
   
   Matrix_free(A);
   Matrix_free(Q);
   Matrix_free(R);
}

void Matrix_columnbasis_check(Matrix V, unsigned int rows, unsigned int columns) {
   unsigned int c1, c2, r;
   real maxselferror = 0.0, maxcrosserror = 0.0, error;

   for (c1 = 0; c1 < columns; c1++)
      for (c2 = c1; c2 < columns; c2++) {
         real dot = 0.0;
         for (r = 0; r < rows; r++)
            dot += V[r][c1] * V[r][c2];
         if (c1 == c2) {
            error = fabs(1.0 - dot);
            if (error > maxselferror)
               maxselferror = error;
         }
         else {
            error = fabs(dot);
            if (error > maxcrosserror)
               maxcrosserror = error;
         }
      }

   printf("MAXSELFERROR: %g MAXCROSSERROR: %g\n", maxselferror, maxcrosserror);
}

void Matrix_transpose(Matrix* matrix, unsigned int rows, unsigned int columns) {
   unsigned int i, j;

   Matrix transpose = Matrix_allocate(columns, rows);
   
   for (i = 0; i < rows; i++)
      for (j = 0; j < columns; j++)
         transpose[j][i] = (*matrix)[i][j];
         
   Matrix_free(*matrix);
   
   *matrix = transpose;
}

void Matrix_diff(Matrix X, Matrix A, Matrix B, unsigned int rows, unsigned int columns) {
   unsigned int i, j;
   
   for (i = 0; i < rows; i++)
      for (j = 0; j < columns; j++)
         X[i][j] = A[i][j] - B[i][j];
}

void Matrix_pseudoinverse_droptol(Matrix XI, Matrix X, unsigned int rows, unsigned int columns, real droptol) {
   char JOBU = 'A', JOBVT = 'A';
   unsigned int m = rows, n = columns, LDA = rows, LDU = rows, LDVT = columns;
   unsigned int mindim = PetscMin(rows, columns);
   unsigned int LWORK = PetscMax(3*PetscMin(m,n)+PetscMax(m,n),5*PetscMin(m,n));
   int INFO;
   Matrix Xcolumnmajor = Matrix_allocate(columns, rows);
   Matrix UT = Matrix_allocate(rows, rows);
   Matrix S = Matrix_allocate(columns, rows);
   Matrix V = Matrix_allocate(columns, columns);
   Matrix temp = Matrix_allocate(columns, rows);
   real D[mindim], WORK[LWORK], tol;
   unsigned int i, j;

   for (i = 0; i < rows; i++)
      for (j = 0; j < columns; j++)
         Xcolumnmajor[j][i] = X[i][j];

#ifdef OMP
#pragma omp critical (lapack)
#endif
   LAPACKgesvd_(&JOBU, &JOBVT, &m, &n, Xcolumnmajor[0], &LDA, D,
                UT[0], &LDU, V[0], &LDVT, WORK, &LWORK, &INFO);

   if (INFO) {
      printf("s/dgesvd returned error code %d\n", INFO);
      exit(-2);
   }

   Matrix_free(Xcolumnmajor);

   tol = D[0] * droptol;

   for (i = 0; i < mindim; i++)
      if (D[i] > tol)
         S[i][i] = (real)1 / D[i];

   Matrix_multiplymatrix(temp, V, S, columns, columns, rows);

   Matrix_free(V);
   Matrix_free(S);

   Matrix_multiplymatrix(XI, temp, UT, columns, rows, rows);

   Matrix_free(temp);
   Matrix_free(UT);
}

void Matrix_pseudoinverse(Matrix XI, Matrix X, unsigned int rows, unsigned int columns) {
   char JOBU = 'A', JOBVT = 'A';
   unsigned int m = rows, n = columns, LDA = rows, LDU = rows, LDVT = columns;
   unsigned int mindim = PetscMin(rows, columns);
   unsigned int LWORK = PetscMax(3*PetscMin(m,n)+PetscMax(m,n),5*PetscMin(m,n));
   int INFO;
   Matrix Xcolumnmajor = Matrix_allocate(columns, rows);
   Matrix UT = Matrix_allocate(rows, rows);
   Matrix S = Matrix_allocate(columns, rows);
   Matrix V = Matrix_allocate(columns, columns);
   Matrix temp = Matrix_allocate(columns, rows);
   real D[mindim], WORK[LWORK], tol;
   unsigned int i, j;

   for (i = 0; i < rows; i++)
      for (j = 0; j < columns; j++)
         Xcolumnmajor[j][i] = X[i][j];

#ifdef OMP
#pragma omp critical (lapack)
#endif
   LAPACKgesvd_(&JOBU, &JOBVT, &m, &n, Xcolumnmajor[0], &LDA, D,
                UT[0], &LDU, V[0], &LDVT, WORK, &LWORK, &INFO);

   if (INFO) {
      printf("s/dgesvd returned error code %d\n", INFO);
      exit(-2);
   }

   Matrix_free(Xcolumnmajor);

   /* These tolerances are the default for MATLAB's pinv routine */

#ifdef REAL_IS_DOUBLE
   tol = PetscMax(rows, columns) * D[0] * DBL_EPSILON;
#else
#ifdef REAL_IS_FLOAT
   tol = PetscMax(rows, columns) * D[0] * FLT_EPSILON;
#endif
#endif

   for (i = 0; i < mindim; i++)
      if (D[i] > tol)
         S[i][i] = (real)1 / D[i];

   Matrix_multiplymatrix(temp, V, S, columns, columns, rows);

   Matrix_free(V);
   Matrix_free(S);

   Matrix_multiplymatrix(XI, temp, UT, columns, rows, rows);

   Matrix_free(temp);
   Matrix_free(UT);
}

#if 0
void Matrix_lsq_solve(Vector x, Matrix A, Vector b, unsigned int rows, unsigned int columns) {
   char TRANS = 'N';
   unsigned int M = rows, N = columns, NRHS = 1, LDA = rows, LDB = PetscMax(rows,columns);
   unsigned int LWORK = PetscMax(1, M*N + PetscMax(M*N, NRHS));
   int INFO;
   real WORK[LWORK];
   Matrix Acolumnmajor = Matrix_allocate(columns, rows);
   Vector bcopy = Vector_allocate(PetscMax(rows, columns));
   unsigned int i, j;

   for (i = 0; i < rows; i++)
      for (j = 0; j < columns; j++)
         Acolumnmajor[j][i] = A[i][j];

   Vector_copy(bcopy, b, rows);

#ifdef REAL_IS_DOUBLE
#ifdef _AIX
   dgels(&TRANS, &M, &N, &NRHS, Acolumnmajor[0], &LDA, bcopy, &LDB, WORK, &LWORK, &INFO);
#else
   dgels_(&TRANS, &M, &N, &NRHS, Acolumnmajor[0], &LDA, bcopy, &LDB, WORK, &LWORK, &INFO);
#endif
#else
#ifdef REAL_IS_FLOAT
#endif
#ifdef _AIX
   sgels(&TRANS, &M, &N, &NRHS, Acolumnmajor[0], &LDA, bcopy, &LDB, WORK, &LWORK, &INFO);
#else
   sgels_(&TRANS, &M, &N, &NRHS, Acolumnmajor[0], &LDA, bcopy, &LDB, WORK, &LWORK, &INFO);
#endif
#endif

   if (INFO) {
      printf("s/dgels returned error code %d\n", INFO);
      exit(-2);
   }

   Vector_copy(x, bcopy, columns);

   Vector_free(bcopy);
   Matrix_free(Acolumnmajor);
}
#endif

/* I/O */

void Matrix_readfile(Matrix matrix, FILE* file, unsigned int rows, unsigned int columns) {
  unsigned int i,j;
  for (i = 0; i < rows; i++) {
	 for (j = 0; j < columns; j++) {
		fscanf(file, "%lf ",&(matrix[i][j]));
	 }
  }
}

void Matrix_writefile(char* filename, Matrix matrix, unsigned int rows, unsigned int columns) {
   unsigned int i, j;
   FILE* file = NULL;
   
   file = fopen(filename, "w");
   
   for (i = 0; i < rows; i++) {
      for (j = 0; j < columns; j++)
         fprintf(file,"%20.15e ", matrix[i][j]);
      fprintf(file,"\n");
   }
   
   fclose(file);
}

void Matrix_print(Matrix matrix, unsigned int rows, unsigned int columns) {
   unsigned int i, j;
   
   for (i = 0; i < rows; i++) {
      for (j = 0; j < columns; j++)
         printf("%20.15e ", matrix[i][j]);
      printf("\n");
   }
}

void Matrix_writebinary(FILE* file, Matrix matrix, unsigned int rows, unsigned int columns) {
   unsigned int count = fwrite(matrix[0], sizeof(real), rows*columns, file);
   if (count != rows*columns)
      printf("Short Write!\n");
}

void Matrix_readbinary(Matrix matrix, FILE* file, unsigned int rows, unsigned int columns) {
/*
   char filename[1024], buf[64];
   sprintf(buf, "/proc/self/fd/%d", fileno(file));
   unsigned int len = readlink(buf, &filename[0], sizeof(filename));
   filename[len] = '\0';
   //printf("READING %u %u %u %s\n", ftell(file), rows, columns, filename);
*/
   if (feof(file))
      printf("End of File!\n");
   if (ferror(file))
      printf("File Error %d!\n", ferror(file));
   unsigned int count = fread(matrix[0], sizeof(real), rows*columns, file);
   if (count != rows*columns)
      printf("Short Read!\n");
}

#if 0
void Matrix_eigendecomposition(Matrix A, Matrix V_mat, Vector d, unsigned int rows, unsigned int columns) {
  unsigned int i,j;
  assert(rows==columns);
  Matrix T = Matrix_allocate(rows,columns);

  Matrix_copy(T, A, rows, columns);
/*   Matrix_transpose(&T, rows, columns); */
/*   Matrix_add(T, A, T, rows, columns); */
/*   Matrix_scale(T, 0.5, rows, columns); */
  // now T = 0.5 * ( A + A' )

  char JOBZ = 'V', RANGE = 'A', UPLO = 'U';
  unsigned int N = rows;

  unsigned int LDA = N, LDZ = N;
  real VL,VU;
  unsigned int IL, IU;
  real ABSTOL = -1.;
  unsigned int M;
  Matrix S = Matrix_allocate(N, N);
  Matrix V = Matrix_allocate(N,N);
  Matrix Vt = Matrix_allocate(N,N);
  Matrix temp = Matrix_allocate(N,N);
  Vector D = Vector_allocate(N);
  unsigned int *ISUPPZ = (unsigned int *)malloc(N * sizeof(unsigned int));
  real *WORK = (real *)malloc(26 * N * sizeof(real));
  unsigned int LWORK = 26 * N;
  unsigned int *IWORK = (unsigned int *)malloc(10 * N * sizeof(unsigned int));
  unsigned int LIWORK = 10 * N;
  int INFO;

#ifdef REAL_IS_DOUBLE
#ifdef OMP
#pragma omp critical (lapack)
#endif
  printf("JOBZ = %c\nRANGE = %c\nUPLO =  %c\nN = %d\nT = %d\nLDA = %d\nVL = %f\nVU = %f\n", JOBZ, RANGE, UPLO, N, T[0], LDA, VL, VU);
  dsyevr_(&JOBZ, &RANGE, &UPLO, &N,
			 T[0], &LDA, &VL, &VU,
			 &IL, &IU, &ABSTOL,
			 &M, D, Vt[0], &LDZ,
			 ISUPPZ, WORK, &LWORK,
			 IWORK, &LIWORK, &INFO);  // V is returned as V'!!!
#else
#ifdef REAL_IS_FLOAT
#ifdef OMP
#pragma omp critical (lapack)
#endif
  ssyevr_(&JOBZ, &RANGE, &UPLO, &N,
			 T[0], &LDA, &VL, &VU,
			 &IL, &IU, &ABSTOL,
			 &M, D, Vt[0], &LDZ,
			 ISUPPZ, WORK, &LWORK,
			 IWORK, &LIWORK, &INFO);  // V is returned as V'!!!
#endif
#endif
  Matrix_copy(V, Vt, N, N);
  Matrix_transpose(&V, N, N);
  Matrix_copy(V_mat, V, N, N);
  Matrix_free(T);
}
#endif
