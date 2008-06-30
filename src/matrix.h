#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <math.h>
#include <time.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
#include <R_ext/Linpack.h>
#include <R_ext/Applic.h>
#include <R_ext/Utils.h>
#include <R_ext/Random.h>
#include <R.h>

#define ME(matrix,row,col) (((matrix)->entries)[(col) * ((matrix)->nr) + (row)])
#define VE(vector,index) (((vector)->entries)[(index)])

#define oops(s) {error((s));}
#define malloc_safe(s,t) if(((s) = malloc((t))) == NULL) { oops("error: malloc() "); }
#define calloc_safe(s,t,T) if(((s) = calloc((t),(T))) == NULL) { oops("error: calloc() "); }
#define max(a,b) ( ((a) > (b)) ? (a) : (b) )
#define min(a,b) ( ((a) > (b)) ? (b) : (a) )
#define malloc_mat(NR, NC, M) { calloc_safe((M),((size_t) 1),sizeof(matrix)); ((M)->nr) = (NR); ((M)->nc) = (NC); ((M)->entries) = calloc(((size_t) (NR)*(NC)) , sizeof(double));}
#define malloc_vec(L, V) { calloc_safe((V),((size_t) 1),sizeof(vector)); ((V)->length) = (L); ((V)->entries) = calloc(((size_t) (L)), sizeof(double));}


typedef struct{
  int nr;
  int nc;
  double *entries;
} matrix;

typedef struct{
  int length;
  double *entries;
} vector;

typedef struct{
  double timec;
  int callc;
} counter;

/* void malloc_mat(int *nrow, int *ncol, matrix *M); */

void free_mat(matrix *M);

/* void malloc_vec(int *length, vector *V); */

void free_vec(vector *V);

int nrow_matrix(matrix *M);

int ncol_matrix(matrix *M);

int length_vector(vector *v);

void print_a_matrix(matrix *M);

extern void F77_SUB(dpotri)(const char* uplo, const int* n,
		 double* a, const int* lda, int* info);

extern void F77_SUB(dpotrf)(const char* uplo, const int* n,
		 double* a, const int* lda, int* info);

extern void F77_SUB(dgemm)(const char *transa, const char *transb, const int *m, \
		const int *n, const int *k, const double *alpha,\
		const double *a, const int *lda,\
		const double *b, const int *ldb,\
		const double *beta, double *c, const int *ldc);

extern void F77_SUB(dgemv)(const char *trans, const int *m, const int *n,
		const double *alpha, const double *a, const int *lda,
		const double *x, const int *incx, const double *beta,
		double *y, const int *incy);

extern void F77_SUB(dgetrf)(const int* m, const int* n, double* a, const int* lda,
                 int* ipiv, int* info);

extern void F77_SUB(dgetri)(const int* n, double* a, const int* lda,
                 int* ipiv, double* work, const int* lwork, int* info);

extern void F77_SUB(dqrdc2)(double *x, int *ldx, int *n, int *p,
                      double *tol, int *rank,
                      double *qraux, int *pivot, double *work);

extern void F77_SUB(dtrco)(double*, int*, int*, double*, double*, int*);

void MtM(matrix *M, matrix *A);

void invertSPD(matrix *A, matrix *AI);

void Mv(matrix *M, vector *v1, vector *v2);

void vM(matrix *M, vector *v1, vector *v2);

vector *vec_star(vector *v1, vector *v2, vector *v3);
  
double vec_sum(vector *v);

double vec_min(vector *v, int *imin);
  
void mat_zeros(matrix *M);
  
void vec_zeros(vector *v);
 
void print_mat(matrix *M);

void print_vec(vector *v);
 
vector *extract_row(matrix *M, int row_to_get, vector *v);

void replace_row(matrix *M, int row_to_set, vector *v);

void vec_add(vector *v1, vector *v2, vector *v3);

vector *scl_vec_mult(double scalar, vector *v1, vector *v2);

matrix *scl_mat_mult(double scalar, matrix *m1, matrix *m2);

matrix *mat_copy(matrix *m1, matrix *m2);

vector *vec_copy(vector *v1, vector *v2);

void mat_subsec(matrix *m1, int rowStart, int colStart,
		       int rowStop, int colStop, matrix *m2);

matrix *mat_transp(matrix *m1, matrix *m2);

void vec_subtr(vector *v1, vector *v2, vector *v3);

void mat_subtr(matrix *m1, matrix *m2, matrix *m3);

void mat_add(matrix *m1, matrix *m2, matrix *m3);

void vec_add_mult(vector *v1, vector *v2, double s, vector *v3);

void MtA(matrix *M, matrix *A, matrix *Mout);

void MAt(matrix *M, matrix *A, matrix *Mout);

void invert(matrix *A, matrix *AI);

void MxA(matrix *M, matrix *A, matrix *Mout);

void R_CheckUserInterrupt(void);

void print_clock(clock_t *intime,int i);

void update_clock(clock_t *intime, counter *C);

void zcntr(counter *C);

void print_counter(int i, counter *C);

void head_matrix(matrix *M);

void head_vector(vector *V);

void identity_matrix(matrix *M);

void malloc_mats(int nrow, int ncol, ...);

void malloc_vecs(int length, ...);

void free_mats(matrix **M, ...);

void free_vecs(vector **V, ...);

vector *vec_ones(vector *v);

void replace_col(matrix *M, int col_to_set, vector *v);

vector *extract_col(matrix *M, int col_to_get, vector *v);
