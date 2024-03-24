#ifndef BVSR_READIN
#define BVSR_READIN


#include "global.h"
#include "generic.h"
#include "lalg.h"

void read_meang(std::string);
void read_plink(std::string);
void read_matrix(std::string);
void read_pheno(std::string);
void read_cov(std::string);
void read_last_bvsr(std::vector<int>&, double&, int);

#endif

