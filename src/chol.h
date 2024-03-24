#ifndef BVSR_CHOL
#define BVSR_CHOL

#include "xy.h"
#include "generic.h"
#include "global.h"


void forward_sub (int, int, double*, double*, double);
void Chol_add (int, int, double*, double*, double*);
void Chol_del (std::vector<int>, int, double*, double*);	
void Chol_init (int, double*, double*);
double Chol_decomp_A (int, double*, double*, double*);
double Chol_decomp_A (int, double*, double*, double);
void Chol_decomp_g_cov(void);
void Chol_solve (int, double*, double*, double*);

#endif

