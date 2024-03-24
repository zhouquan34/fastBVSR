#ifndef BVSR_XY
#define BVSR_XY

#include "global.h"
#include "lalg.h"
#include "generic.h"

void X_init(const std::vector<int>&, double*);
void X_init(const std::vector<int>&, float**);
void X_add(const std::vector<int>&, int, double*, double*);
void X_add(const std::vector<int>&, int, float**, float**);
void X_del(const std::vector<int>&, int, double*, double*);
void X_del(const std::vector<int>&, int, float**, float**);
void Xy_init(double*, int, double*);
void Xy_init(float**, int, double*);
void Xy_add(const std::vector<int>&, int, double*, double*);
void Xy_del(const std::vector<int>&, int, double*, double*);
void XX_init(double*, int, double*);
void XX_init(float*, int, double*);
void XX_init(const std::vector<int>&, double*, int, double*);
void XX_init(const std::vector<int>&, float**, int, double*);
void XX_add(const std::vector<int>&, int, double*, double*, double*);
void XX_add(const std::vector<int>&, const std::vector<int>&, int, double*, double*, double*);
void XX_add(const std::vector<int>&, const std::vector<int>&, int, float**, double*, double*);
void XX_del(const std::vector<int>&, int, double*, double*);

#endif
