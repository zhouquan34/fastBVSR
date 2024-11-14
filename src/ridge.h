#ifndef BVSR_RIDGE
#define BVSR_RIDGE

#include "global.h"
#include "generic.h"
#include "icf.h"
#include "lalg.h"
#include "chol.h"

/*
 calc_globals();
    if g_n_cov > 0, regress out confounding covariates.
    compute g_yy, g_xty. 

 calc_xtx();
    precompute part of xtx.

 filter_snp(sigma); 
    filter identical snps if filtering is turned on.
    compute single SNP bayes factor and order BF. 
    compute g_single_id, g_single_ord.  

 order: calc_globals => filter_snp => calc_xtx. 
 
 calc_ridge(Ra, R, z, p, d, b, XX);
     solve (RtR + d^2 I) b = z or Ab = z. 
     Ra is the Chol factor of A; if it exists, use it. 
     if p is small, use exact algorithm and return log_determinant of A.
     otherwise, use ICF. if ICF fails to converge in 10 iterations, do exact.
     Note that if RA is calculated, it is also initialized.  


*/


double calc_ridge(double*&, double*, double*, int, double*, double*, double*);
double calc_ridge(double*&, double*, double*, int, double, double*, double*);
void calc_globals(void);
void filter_snp(const double);
void calc_xtx(void);

#endif

