#ifndef BVSR_MODEL_2
#define BVSR_MODEL_2

#include "chol.h"
#include "ridge.h"
#include "generic.h"
#include "global.h"
#include "lalg.h"
#include "ftime.h"
#include "model.h"

//extern std::map<double, int> g_h2;
extern int g_size_g1;
extern int g_size_g2;
extern std::map<double, int> g_rho;

class model_two_comp{
	private:
		int size_a;  
		int size_b;
		double sigma_a;
		double sigma_b;
		double pi_a;  
		double pi_b;
		double sse;
		double h2; 
		double rho;
		double pve;
		double tau;
		double like;
		double var_sy;
		double sbf;
		int msize;
		double* diags;
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
		model_two_comp();
		~model_two_comp();
		void allocate(int);
		void initialize(std::vector<int>&, double=0.2);
		void copy(class model_two_comp*);
		int model_size(void);
		void update_marginal(int);
		std::string para_str(void);
		std::string model_str(void);
		void calc_sigmas(void);
		void calc_rss(void);
		double calc_rss(double*);
		void sample_pve(void);
		double sample_h2(void);
		void sample_pi(void);
		double add_one_snp(std::vector<int>&);
		double add_snps(int, model_two_comp*);
		double delete_one_snp(std::vector<int>&);
		double delete_snps(int, model_two_comp*);
		void sample_y (void);
		int compare_models(double, model_two_comp*);
		int compare_models_sbf(double, model_two_comp*);
		void rao_blackwell(void);
		void sample_tau(void);
		void calc_sbf(int);
};

#endif

