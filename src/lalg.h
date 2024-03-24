#ifndef BVSR_LALG
#define BVSR_LALG

#include "global.h"
#include "ftime.h"
#include "generic.h"

#define ZERO 1e-12

double* allocate1D (int);
float* allocate1D_float (int);
float** allocate1D_float_pointer (int);
double* allocate1D (int, int);
double** allocate2D (int);
void free1D (double*);
void free1D (float*);
void free1D (float**);
void free2D (double**);
void copy1D (double*, const double*, int);
void copy1D (double*, const float*, int);
void copy1D (float*, const float*, int);
void copy1D (float**, float**, int);
void mat_vec_mul (double*, double*, int, int, double*, int);
void mat_vec_mul (float*, double*, int, int, double*, int);
void mat_vec_mul (float**, double*, int, int, double*, int);
void mat_vec_mul (float**, float*, int, int, double*, int);
double vec_vec_mul (double*, double*, int);
double vec_vec_mul (float*, double*, int);
double vec_vec_mul (float*, float*, int);
double L2_norm (double*, int);
void center_rows(double*, int , int);
void center_rows(float*, int , int);
double center_array(double*, int);
void array_del (const std::vector<int>&, int, double*, double*);
void matrix_del (const std::vector<int>&, int, int, double*, double* );
void matrix_del (const std::vector<int>&, int, int, float*, float* );
double array_sst(double*, int);

#endif

