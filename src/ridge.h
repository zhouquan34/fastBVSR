#ifndef BVSR_RIDGE
#define BVSR_RIDGE

#include "global.h"
#include "generic.h"
#include "icf.h"
#include "lalg.h"
#include "chol.h"

double calc_ridge(double*&, double*, double*, int, double*, double*, double*);
double calc_ridge(double*&, double*, double*, int, double, double*, double*);
void calc_globals(void);
void filter_snp(const double);
void calc_xtx(void);

#endif

