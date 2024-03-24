#include "bvsr.h"

using namespace std;

void exit_bvsr(void){
	gsl_exit();
	free1D(g_data);
	free1D(g_cov);
	free1D(g_pheno);
	free1D(g_xty);
	free1D(g_cov_R);
	free1D(g_xtx);
	return;
}

void allocate_all(float**& X, double*& XX, double*& R, double*& z, int p){
	X = allocate1D_float_pointer(p);
	XX = allocate1D(p * (p+1) / 2);
	R = allocate1D(p * p);
	z = allocate1D(p);
}

void free_all(float** X, double* XX, double* R, double* z){
	free1D(X);
	free1D(XX);
	free1D(R);
	free1D(z);
}

void tests(){
	cout.precision(15);
	float x (1.23456789f);
	double dx = (double) x;
//	cout << x << "\t" << dx << endl;
	vector<double> tt;
	tt.push_back(1.0/0.0); tt.push_back(11);
	if (isinf(tt[0])){cout << "Inf" << endl;}
	else{cout << "nonINF" << endl;}

	float** X0; float** X1; float** X2;
	double* XX0; double* XX1; double* XX2;
	double* R0; double* R1; double* R2;
	double* z0; double* z1; double* z2;
	int p0 = 4 + g_n_cov; int p1 = 6 + g_n_cov; int p2 = 4 + g_n_cov;
	allocate_all(X0, XX0, R0, z0, p0);
	allocate_all(X1, XX1, R1, z1, p1);
	allocate_all(X2, XX2, R2, z2, p2);
	
	vector<int> s1; 
	s1.push_back(0); s1.push_back(1); 
	s1.push_back(6); s1.push_back(9);
	vector<int> add; add.push_back(2); add.push_back(8);
	vector<int> del; del.push_back(1); del.push_back(3);
	vector<int> s2 = s1; s2.push_back(2); s2.push_back(8);

	X_init(s1, X0);
	X_add(add, p0, X0, X1);
	X_del(del, p1, X1, X2);

	Xy_init(X0, p0, z0);
	Xy_add(add, p0, z0, z1);
	Xy_del(del, p1, z1, z2);
	
	XX_init(s1, X0, p0, XX0);
	XX_add(s2, add, p0, X1, XX0, XX1);
	XX_del(del, p1, XX1, XX2);
	
	Chol_init(p0, XX0, R0);
	Chol_add(2, p0, XX1, R0, R1);
	Chol_del(del, p1, R1, R2);

	double* diags = allocate1D(p2); 
	double d1 = sqrt(4.0 - g_diag_c);
	double d2 = sqrt(25.0 - g_diag_c);
	for (int i=0; i<g_n_cov; i++){diags[i]=0;}
	for (int i=g_n_cov; i<p2-2; i++){diags[i]=d1;}
	for (int i=p2-2; i<p2; i++){diags[i]=d2;}
	double* b = allocate1D(p2);
	double* RA = NULL; 
	calc_ridge(RA, R2, z2, p2, diags, b, XX2);
	print_array(b, p2);
//	calc_ridge(RA, R2, z2, p2, sqrt(4.0 - g_diag_c), b, XX2);
//	print_array(b, p2);
//	print_matrix(RA, p2 * p2, p2);
//	print_array(XX2, p2 * (p2 + 1)/2);
//	print_matrix(X2, p2 * g_n_sub, g_n_sub);
//	print_array(z2, p2);
	free1D(diags);
	free1D(RA);
	free1D(b);

	free_all(X0, XX0, R0, z0);
	free_all(X1, XX1, R1, z1);
	free_all(X2, XX2, R2, z2);
	return;
}

int main(int argc, char** argv){
	cout.precision(15);
/*	
	cout << gsl_cdf_beta_P(1e-9, 200, 200000) << endl;
	cout << gsl_cdf_beta_P(1e-8, 200, 200000) << endl;
	cout << gsl_cdf_beta_P(1e-7, 200, 200000) << endl;
	cout << gsl_cdf_beta_P(1e-6, 200, 200000) << endl;
	cout << gsl_cdf_beta_P(1e-5, 200, 200000) << endl;
	cout << gsl_cdf_beta_P(1e-4, 200, 200000) << endl;
	cout << gsl_cdf_beta_P(1e-3, 200, 200000) << endl;
*/
	read_args(argc, argv);
	gsl_init();
	open_outputs(g_out_prefix);
	write_log_cmds(argc, argv);
	calc_globals();
	filter_snp(0.2);
	calc_xtx();

//	tests()	;	
	clock_t begin_t = clock();
	if (g_group == 1){
		mcmc();
	}else if (g_group == 2){
		mcmc_two_comp();
	}
	clock_t end_t = clock();
	double used_t =( (double) end_t - begin_t)/CLOCKS_PER_SEC;
	write_time(used_t);
	
	if (g_output_time == 1){
		output_time();
	}
	
	close_outputs();
	exit_bvsr();
	return 0;
}


