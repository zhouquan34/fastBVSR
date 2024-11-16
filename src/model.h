#ifndef BVSR_MODEL
#define BVSR_MODEL

#define FAIL_ALPHA -1e15

#include "chol.h"
#include "ridge.h"
#include "generic.h"
#include "global.h"
#include "lalg.h"

extern std::vector<double> g_log_score;

class model{
	private:
		double sigma;
		double pi;
		double sse;
		double h2; 
		double pve;
		double br2; 
		double tau;
		double like;
		double var_sy;
		double sbf;
		double gamma;
		double alpha;
		int msize;
		std::vector<int> subset;
		std::vector<int> selected;
	public:
		double* R;
		double* RA;
		float** X;
		double* XX;
		double* Xy;
		double* Beta;
		double* sY;
		double* sYhat;
		double* hatY;
		model();
		~model();
		void allocate(int);
		void initialize(std::vector<int>&, double=0.2);
		void copy(class model*);
		int model_size(void);
		void update_marginal(int);
		std::string para_str(void);
		std::string model_str(void);
		void calc_sigma(void);
		void calc_rss(void);
		double calc_rss(double*);
		void sample_pve(void);
		double sample_h2(void);
		void sample_pi(void);
		double add_one_snp(std::vector<int>&);
		double add_snps(int, model*);
		double delete_one_snp(std::vector<int>&);
		double delete_snps(int, model*);
		void sample_y (void);
		int compare_models(double, model*);
		int compare_models_sbf(double, model*);
		void rao_blackwell(void);
		void sample_tau(void);
		void calc_sbf(int);
		double calc_alpha (void);
		void print_pointers(void);
};

#endif

