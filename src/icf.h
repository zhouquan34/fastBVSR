#ifndef BVSR_ICF
#define BVSR_ICF

#include "lalg.h"


/* 
   We assume the first c columns of X represent confounding covariates.
   In total X contains p columns. 
   The diagonal regularization is a vector with first c elements = 0. 

   ICF(R, z, p, g_n_cov, d, b);
    R is p times p. Solve: (R^t R + d^2 I)b = z

*/

int ICF(const double*, const double*, int, int, const double, double*);

#endif

