#ifndef BVSR_GLOBAL
#define BVSR_GLOBAL

#include <iostream>
#include <fstream>
#include <cstring>
#include <vector>
#include <map>
#include <sys/time.h>
#include <set>
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <algorithm>
#include <sstream>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_cdf.h>

extern int g_n_sub;
extern int g_n_cov;
extern int g_n_snp;
extern std::vector<int> g_n_grp;
extern int g_group;
extern int g_missing;
extern float* g_data;
extern float* g_cov;
extern double* g_pheno;
extern double* g_xty;
extern double g_yy;
extern double g_ypwy;
extern double* g_cov_R;
extern int g_no_filter;
extern double g_corr_cut;
extern std::vector<std::string> g_snp_names;
extern std::vector<double> g_snp_var;
extern std::vector<double> g_snp_mean;
extern std::vector<int> g_snp_map;
extern std::set<int> g_dup_keep;
extern std::map<int, int> g_dup_map;
extern std::map<int, int> g_dup_count;
extern std::map<int, double> g_dup_bf;
extern int g_precomp;
extern double* g_xtx;
extern std::vector<double> g_dd_moments;

extern std::string g_out_prefix;
extern std::ofstream LOG;
extern std::ofstream PATH;
extern std::ofstream MOD;
extern std::ofstream OUT;
extern std::ofstream YHAT;
extern std::map<std::string, int> g_paras;

extern int g_output_time;
extern std::map<std::string, double> g_time;
extern struct timeval g_t_beg;
extern struct timeval g_t_end;

extern double g_diag_c;
extern double g_icf_abs_tol;
extern int g_chol_call;
extern int g_icf_min_p;
extern std::map<int, int> g_icf_iter;

extern std::string g_last_bvsr;
extern double g_default_sigma;
extern int g_start_size;
extern int g_rb_thin;
extern int g_max_jump;
extern double g_prop_uniform;
extern double g_geometric;
extern double g_prop_add;
extern double g_long_range;
extern double g_long_exchange;
extern double g_jump_h2;
extern double g_jump_rho;
extern double g_prop_h2_uniform;

extern double g_log_min_h;
extern double g_log_max_h;
extern double g_min_h;
extern double g_max_h;
extern double g_min_pi;
extern double g_max_pi;
extern double g_min_pi_2;
extern double g_max_pi_2;
extern int g_pi_range_set;
extern int g_pi_range_set_2;
extern int g_log_uniform_h;
extern int g_pi_prior_alpha;
extern int g_pi_prior_beta;
extern int g_mcmc_warm;
extern int g_max_model_size;
extern int g_min_model_size;
extern int g_expect_size_a;
extern int g_pi_alpha_2;
extern int g_pi_beta_2;
extern int g_pseudo_n;
extern int g_pset;
extern int g_pmin;
extern int g_pmax;

extern int g_mcmc_i;
extern int g_mcmc_init_try;
extern int g_sbf_perm;
extern int g_sbf;
extern int g_rb_count;
extern int g_mcmc_iter;
extern int g_last_iter;
extern int g_mcmc_warm;
extern int g_mcmc_stay;
extern int g_continue;
extern int g_exact_bf;
extern std::vector<double> g_single_bf;
extern std::vector<double> g_single_h2;
extern std::vector<int> g_pip;
extern std::vector<double> g_beta;
extern std::vector<double> g_rb_pip;
extern std::vector<double> g_rb_beta;
extern std::vector<int> g_single_id;
extern std::vector<int> g_single_ord;
extern std::vector<int> g_propose;
extern std::vector<int> g_accept;
extern std::map<int, int> g_model_size;
extern std::map<double, int> g_h2;
extern std::map<double, int> g_pve;

extern int g_output_yhat;

extern const gsl_rng_type* gsl_type;
extern gsl_rng* gsl_r;
extern int gsl_seed;

#endif

