#include "global.h"

using namespace std;

/*******  data *******/
int g_n_sub = 0;
int g_n_cov = 0;
int g_n_snp = 0;
vector<int> g_n_grp;
int g_group = 1;
int g_missing = 0;
float* g_data = NULL;
float* g_cov = NULL;
double* g_pheno = NULL;
double* g_xty = NULL;
double* g_cov_R = NULL;
double g_ypwy = 0; // no use
double g_yy = 0;
int g_no_filter = 1;
double g_corr_cut = 1;
vector<string> g_snp_names;
vector<double> g_snp_var;
vector<double> g_snp_mean;
vector<int> g_snp_map;
set<int> g_dup_keep;
map<int, int> g_dup_map;
map<int, int> g_dup_count;
map<int, double> g_dup_bf;
int g_precomp = 5000;
double* g_xtx = NULL;
vector<double> g_dd_moments;

/******** IO *********/
string g_out_prefix = "bvsr_test";
ofstream LOG;
ofstream PATH;
ofstream MOD;
ofstream YHAT;
ofstream OUT;
map<string, int> g_paras;
int g_output_yhat = 0;

/******** time *********/
int g_output_time = 0;
struct timeval g_t_beg;
struct timeval g_t_end;
map<string, double> g_time;


/******** ridge *******/
double g_diag_c = 1e-2;
double g_icf_abs_tol = 1e-6;
int g_icf_min_p = 30;
map<int, int> g_icf_iter;
int g_chol_call = 0;

/******* mcmc settings *******/
string g_last_bvsr;
double g_default_sigma = 0.2;
int g_start_size = 0;
int g_rb_thin = 1000;
int g_max_jump = 5;
double g_prop_uniform = 0.3;
double g_geometric = 0.002;
double g_prop_add = 0.5; // CANNOT BE CHANGED
double g_long_range = 0.1;
double g_long_exchange = 0.3;
double g_jump_h2 = 0.1;
double g_jump_rho = 0.1;
double g_prop_h2_uniform = 0.8;
double g_min_pi = 0;
double g_max_pi = 1;
double g_min_pi_2 = 0;
double g_max_pi_2 = 1;
int g_pi_range_set = 0;
int g_pi_range_set_2 = 0;
double g_min_h = 1e-6;
double g_max_h = 0.999;
double g_log_min_h = -10.0;
double g_log_max_h = 0.0;
int g_log_uniform_h = 0;
int g_pi_prior_alpha = 0;
int g_pi_prior_beta = 1;
int g_max_model_size = 4000;
int g_min_model_size = 1;
int g_pi_alpha_2 = 1;
int g_pi_beta_2 = 1;
int g_expect_size_a = 0;
int g_pseudo_n = 0; 
int g_pmax = -1;
int g_pmin = -1;
int g_pset = 0;

/******* mcmc ********/
int g_mcmc_i = 0;
int g_mcmc_init_try = -1;
int g_sbf_perm = 10;
int g_sbf = 0;
int g_rb_count = 0;
int g_mcmc_iter = 1000;
int g_last_iter = 0;
int g_mcmc_warm = 0;
int g_mcmc_stay = 0;
int g_continue = 0;
int g_exact_bf = 0;
vector<double> g_single_bf;
vector<double> g_single_h2;
vector<int> g_pip;
vector<double> g_beta;
vector<double> g_rb_pip;
vector<double> g_rb_beta;
vector<int> g_single_ord;
vector<int> g_single_id;
vector<int> g_propose;
vector<int> g_accept;
map<int, int> g_model_size;
map<double, int> g_h2;
map<double, int> g_pve;

/****** gsl *******/
const gsl_rng_type* gsl_type;
gsl_rng* gsl_r;
int gsl_seed;



