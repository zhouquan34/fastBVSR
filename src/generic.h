#ifndef BVSR_GEN
#define BVSR_GEN

#include "global.h"

void gsl_init(void);
void gsl_exit(void);
bool compare_pair(std::pair<int, double>, std::pair<int, double>);
double log_sum_log(double, double);
double vector_mean(std::vector<double>&);
void get_row_col (const std::string, int&, int&);
void print_array(const double*, int);
void print_array(const float*, int);
void print_matrix(const double*, int, int);
void print_matrix(const float*, int, int);
void write_log_cmds(int, char**);
double sample_beta(double, double);
double calc_mean_scale(double);

template <typename T>
void vector_print(std::vector<T>& v){
	typename std::vector<T>::size_type i;
	for (i=0; i<v.size(); i++){
		std::cout << v[i] << " ";
	}
	std::cout << std::endl;
}

#endif

