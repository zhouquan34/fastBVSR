#include "model.h"

using namespace std;

vector<double> g_log_score; 
int g_upper_add = 2;
int g_lower_add = 0;
int g_upper_del = 1;
int g_lower_del = 0;

void gsl_matrix_print (const gsl_matrix* m){
	for (int i=0; i<m->size1; i++){
		for (int j=0; j<m->size2; j++){
			cout << gsl_matrix_get(m, i, j) << "\t";	
		}	
		cout << endl;
	}	
	return;
}


// int g_fix_sigma = 1;

model::model(void){
	sse = 0.0;
	sigma = 0.0;
	pi = 0.0;
	h2 = 0.0;
	pve = 0.0;
	pve_rb = 0.0;
	tau = 0.0;
	var_sy = 0.0;
	like = 0.0;
	sbf = 0.0;
	br2 = 0.0;

	R = NULL;
	RA = NULL;
	XX = NULL;
	X = NULL;
	Xy = NULL;
	Beta = NULL;
	hatY = NULL;
	sYhat = NULL;
	sY = NULL;
}

model::~model(){
	free1D(R);	
	free1D(RA);	
	free1D(XX);	
	free1D(Xy);
	free1D(Beta);
	free1D(sY);
	free1D(sYhat);
	free1D(hatY);
	free1D_pointer(X);
}


int model::model_size(void){return subset.size();}

string model::model_str(void){
	stringstream s;
	s << subset[0];
	for (int i=1; i<subset.size(); i++){
		s << "," << subset[i];
	}
	if (g_exact_bf == 1){
		double l = like + gsl_sf_lngamma(subset.size()) + gsl_sf_lngamma(g_pseudo_n + 1 - subset.size()) ;
		s << "\t" << sigma << "\t" << l;
	}
	return s.str();
}

string model::para_str(void){
	stringstream s;
	s << subset.size() << "\t" << h2 << "\t"  << pve << "\t" << pi << "\t" << sigma << "\t" << br2 << "\t" << pve_rb;
	return s.str();
}

void model::update_marginal(int k){
//	record_time();
	for (int i=0; i<subset.size(); i++){
		g_pip[subset[i]] += k;
		g_beta[subset[i]] += k * Beta[i + g_n_cov];
	}
	g_model_size[subset.size()] += k;
	g_h2[h2] += k;
	g_pve[pve] += k;
//	record_time("UpdateMarg");
	return;
}

void model::allocate(int p){
// RA is allocated in calc_ridge
//	record_time();
	X = allocate1D_float_pointer(p);
	XX = allocate1D(p * (p+1) / 2);
	R = allocate1D(p * p);
	Xy = allocate1D(p);
	Beta = allocate1D(p);
	sY = allocate1D(g_n_sub);
	sYhat = allocate1D(g_n_sub); // this must be initiated to zero
	hatY = allocate1D(g_n_sub);
//	record_time("Allocate");
}

void model::initialize(vector<int>& start, double h2_start){
	subset = start;
	h2 = h2_start; 
	msize = subset.size() + g_n_cov;
	allocate(msize);	
	X_init(subset, X);
	XX_init(subset, X, msize, XX);
	Xy_init(X, msize, Xy);
	Chol_init(msize, XX, R);
	
	selected.clear();
	selected.assign(g_n_snp,0);
	for (int i=0; i<subset.size(); i++){
		selected[subset[i]] = 1;
	}
	
	calc_sigma();
	calc_rss();
	sample_pi();	
	sample_y();
	sample_pve();
	
	return;
}

void model::copy(model* m){  // sbf pve tau var_sy not to be copied
//	record_time();
	sigma = m->sigma;
	pi = m->pi;
	h2 = m->h2;
	sse = m->sse;
	like = m->like;
	msize = m->msize;
	subset = m->subset;
	selected = m->selected;
//	record_time("CopyModel");
}

void model::calc_sigma(void){
	double ss = 0.0;
	for (int i=0; i<subset.size(); i++){
		ss += g_snp_var.at(subset[i]);
	}
	sigma = sqrt(h2 / (1-h2) / ss);
	return;
}

double model::calc_alpha(void){
	return 1.0 - g_gamma * h2; 
}

void model::print_pointers(void){
	cerr << "RA address" << RA << endl;
	cerr << "R address" << R << endl;
	cerr << "Xy address" << Xy << endl;
	cerr << "Beta address" << Beta << endl;
	cerr << "XX address" << XX << endl;	
}


void model::calc_rss(void){
	double log_det = calc_ridge(RA, R, Xy, msize, sqrt(1.0/(sigma*sigma*calc_alpha()) - g_diag_c), Beta, XX);	
	sse = g_yy - vec_vec_mul(Xy, Beta, msize); 
	like = -0.5 * g_n_sub * calc_alpha() * log(sse);
	if (g_exact_bf == 1){
		like += ( (-0.5 * log_det)  - subset.size() * log(sigma) );
	}
	return;
}

double model::calc_rss(double* tilde_y){
	double* z1 = allocate1D(msize);
	record_time();
	mat_vec_mul(X, tilde_y, msize, g_n_sub, z1, 0);
	record_time("X-y");

	double tyy = vec_vec_mul(tilde_y, tilde_y, g_n_sub);
	double* beta1 = allocate1D(msize);
	calc_ridge(RA, R, z1, msize, sqrt(1.0/(sigma*sigma*calc_alpha()) - g_diag_c), beta1, XX);
	double rss = tyy - vec_vec_mul(z1, beta1, msize);
	free1D(z1);
	free1D(beta1);
	return -0.5 * g_n_sub * calc_alpha() * log(rss);
}

void model::sample_pve(void){
	record_time();	
	mat_vec_mul(X, Beta, msize, g_n_sub, hatY, 1);
	record_time("X-beta");
	
	double v_pred = array_sst(hatY, g_n_sub);	
	double* res = allocate1D(g_n_sub);
	for (int i=0; i<g_n_sub; i++){res[i] = g_pheno[i] - hatY[i];}
	double v_res = array_sst(res, g_n_sub);
	br2 = v_pred/(v_pred + v_res);
	free1D(res);

	sample_tau();
	double sum_yy = 0; 
	for (int i=0; i<g_n_sub; i++){
		double y = hatY[i] + gsl_ran_gaussian(gsl_r, sqrt(((double) subset.size())/(g_n_sub * tau)) );
		sum_yy += y * y;
	}
	double v1 = sum_yy * tau / (double) g_n_sub;
	pve = v1 / (1.0 + v1);
//	record_time("SamplePVE");
	return;
}

double model::sample_h2(void){
	double r = 0.0;
	if (g_log_uniform_h == 0){
		h2 = h2 + 2.0 * g_jump_h2 * gsl_rng_uniform(gsl_r) - g_jump_h2; 
		if (h2 < g_min_h){ h2 = 2.0 * g_min_h - h2; }
		if (h2 > g_max_h){ h2 = 2.0 * g_max_h - h2; }
	}else{
		if (gsl_rng_uniform(gsl_r) < g_prop_h2_uniform){
			r = log(1.0/h2);
			h2 = h2 + 2.0 * g_jump_h2 * gsl_rng_uniform(gsl_r) - g_jump_h2; 
			if (h2 < g_min_h){ h2 = 2.0 * g_min_h - h2; }
			if (h2 > g_max_h){ h2 = 2.0 * g_max_h - h2; }
			r = log(1.0/h2) - r;
		}else{
			double log_h2 = log(h2) + 0.3 * gsl_rng_uniform(gsl_r) - 0.15 ; 
			if (log_h2 > g_log_max_h) {log_h2 = 2.0*g_log_max_h - log_h2; }
			if (log_h2 < g_log_min_h) {log_h2 = 2.0*g_log_min_h - log_h2;}
			h2 = exp(log_h2);
		}
	}
	return r;
}


void model::sample_pi(void){
	record_time();
	pi = sample_beta_distr(subset.size(), g_pseudo_n - subset.size() + 1.0);
	record_time("SamplePI");
}

double add_weight (int k){
	double q = g_log_score[k];
	double lp = log((double) g_n_snp);
	if (q > g_upper_add * lp ){ q = g_upper_add * lp;}
	if (q < g_lower_add * lp ){ q = g_lower_add * lp;}
	return q; 
}

double del_weight (int k){
	double q = -g_log_score[k]; 
	double lp = log((double) g_n_snp);
	if (q > g_upper_del * lp ){ q = g_upper_del * lp;}
	if (q < g_lower_del * lp ){ q = g_lower_del * lp;}
	return q; 
}


double model::add_one_snp(vector<int>& add){
	vector<int> ids; 
	vector<double> pp; 	
	for (int i=0; i<g_n_snp; i++){
		if (selected[i] == 0){
			ids.push_back(i);
			pp.push_back(add_weight(i));
		}
	}
	double q_old2new; 
	int k = log_weighted_sample(pp, q_old2new);
	int pick = ids[k]; 
	selected[pick] = 1;
	subset.push_back(pick);
	add.push_back(pick);
	msize ++;
	
	double s = -1e10;
	double qk = 0;
	
	for (int j=0; j<subset.size(); j++){
		double q = del_weight(subset[j]);
		s = log_sum_log(s, q); 
		if (j == subset.size() - 1){qk = q;}
	}
	double q_new2old = qk - s;
	return q_new2old - q_old2new;
}

double model::add_snps (int add_size, model* last_model){
	if (subset.size() + add_size > g_max_model_size){return FAIL_ALPHA;} 
	vector<int> add;
	int ss1 = subset.size();	
	double alpha = 0.0;
	allocate(last_model->msize + add_size);
	for (int i=0; i<add_size; i++){
		alpha += add_one_snp(add); // this is the proposal ratio of gamma
	}
	X_add(add, last_model->msize, last_model->X, X);
	XX_add(subset, add, last_model->msize, X, last_model->XX, XX);
	Xy_add(add, last_model->msize, last_model->Xy, Xy);
	Chol_add(add_size, last_model->msize, XX, last_model->R, R);
	alpha += (gsl_sf_lnpoch(ss1 + g_pi_prior_alpha, add_size) - gsl_sf_lnpoch(g_pseudo_n + g_pi_prior_beta - ss1 - add_size, add_size));
	if (g_pi_range_set == 1){
		int a1 = ss1 + add_size + g_pi_prior_alpha;
		int b1 = g_pseudo_n - ss1 - add_size + g_pi_prior_beta;
		int a2 = ss1 + g_pi_prior_alpha;
		int b2 = g_pseudo_n - ss1 + g_pi_prior_beta;
		alpha += log( (gsl_cdf_beta_P(g_max_pi, a1, b1) - gsl_cdf_beta_P(g_min_pi, a1, b1)) / (gsl_cdf_beta_P(g_max_pi, a2, b2) - gsl_cdf_beta_P(g_min_pi, a2, b2)) ); 
	}

//	gsl_rng_uniform(gsl_r);
	if (g_mcmc_i % 50000 == 0){
		Chol_decomp_A(msize, XX, R, 0.0);
	}
	return alpha;
}


double model::delete_one_snp(vector<int>& del){
	vector<double> pp; 	
	for (int j=0; j<subset.size(); j++){
		pp.push_back(del_weight(subset[j]));
	}
	double q_old2new; 
	int pick_r = log_weighted_sample(pp, q_old2new);
	int pick = subset[pick_r]; 	
	selected[pick] = 0;
	subset.erase(subset.begin() + pick_r);
	msize --;
	
	set<int> tmp;
	int del_r = pick_r;
	int max_del = del_r + del.size();
	for (int i=0; i<del.size(); i++){
		int m = del.at(i);
		if (m <= del_r){
			del_r ++;
			while (tmp.find(del_r) != tmp.end()){
				del_r ++;
			}
		}
		if (m > del_r && m <= max_del){
			tmp.insert(m);
		}
	}	
	del.push_back(del_r);


	double s = -1e10;
	double qk = 0;
	for (int i=0; i<g_n_snp; i++){
		if (selected[i] == 0){
			double q = add_weight(i);
			s = log_sum_log(s, q);
			if (i == pick){ qk = q; }
		}
	}
	double q_new2old = qk - s; 
	return q_new2old - q_old2new;
}

double model::delete_snps (int del_size, model* last_model){	
	if (subset.size() < (del_size + g_min_model_size) ){return FAIL_ALPHA;}	
	double alpha = 0.0;
	int ss1 = subset.size();
	allocate(last_model->msize - del_size);
	vector<int> del;
	for (int i=0; i<del_size; i++){
		alpha += delete_one_snp(del); // this is the proposal ratio of gamma
	}	
	sort(del.begin(), del.end()); // important. 
	X_del(del, last_model->msize, last_model->X, X);
	XX_del(del, last_model->msize, last_model->XX, XX);
	Xy_del(del, last_model->msize, last_model->Xy, Xy);
	Chol_del(del, last_model->msize, last_model->R, R);	
	alpha += ( -gsl_sf_lnpoch(ss1 + g_pi_prior_alpha - del_size, del_size) + gsl_sf_lnpoch(g_pseudo_n + g_pi_prior_beta - ss1, del_size) );
	if (g_pi_range_set == 1){
		int a1 = ss1 - del_size + g_pi_prior_alpha;
		int b1 = g_pseudo_n - ss1 + del_size + g_pi_prior_beta;
		int a2 = ss1 + g_pi_prior_alpha;
		int b2 = g_pseudo_n - ss1 + g_pi_prior_beta;
		alpha += log( (gsl_cdf_beta_P(g_max_pi, a1, b1) - gsl_cdf_beta_P(g_min_pi, a1, b1)) / (gsl_cdf_beta_P(g_max_pi, a2, b2) - gsl_cdf_beta_P(g_min_pi, a2, b2)) ); 
	}	
	return alpha;
}


void model::sample_y (void){
	// when we sample Y which is used for calculating BF, we don't need to sample tau. This is because tau is cancelled out in BF. So assume tau = 1.	
	double* sbeta = allocate1D(subset.size()); // beta for cov is zero	
	for (int i=0; i<subset.size(); i++){ sbeta[i] = gsl_ran_gaussian(gsl_r, sigma * sqrt(calc_alpha()));}
	
	record_time();	
	mat_vec_mul(&X[g_n_cov], sbeta, subset.size(), g_n_sub, sYhat, 1);
	record_time("X-beta");
	
	for (int i=0; i<g_n_sub; i++){ sY[i] = sYhat[i] +  gsl_ran_gaussian(gsl_r, 1.0);}		
	var_sy = center_array(sY, g_n_sub);
	free1D(sbeta);
	return;
}

int model::compare_models(double alpha, model* last_model){
//	cerr << "C1a" << endl;
	calc_sigma();
//	cerr << "C1b" << endl;
	calc_rss();
//	cerr << "C1c" << endl;
	double mh = alpha;
	mh += (like - last_model->like);
//	cerr << "C1d" << endl;
	if (g_exact_bf == 0){
//		cerr << "C2" << endl;
		sample_y();
		double f_new = calc_rss(sY);
		double f_old = last_model->calc_rss(sY); // These two calculations could be faster since they are almost the same. Might save 2% time.
		mh += (f_old - f_new);   //cout << alpha << "\t" << f_old << "\t" << f_new << "\t" << mh << endl;
//		cerr << "C3" << endl;
	}

	double r = gsl_rng_uniform(gsl_r);
	if (exp(mh) > r){
		return 1; // accept
	}else{
		return 0; // reject
	}
}


void model::sample_tau(void){
	tau = gsl_ran_gamma(gsl_r, g_n_sub * calc_alpha()/2, 2.0/sse);	
	return;
}


void model::rao_blackwell (void){
	record_time();
//	sample_tau();
	
	// assuming g_n_cov = 0; need revision later.  
	double* Beta_rb = allocate1D(g_n_snp);
	double lambda = pi/(1-pi);
	double ss = h2 /( (1-h2) * sigma * sigma);
	double* tilde_y = allocate1D(g_n_sub);
	double* residual = allocate1D(g_n_sub);
	copy1D(residual, g_pheno, g_n_sub);
	for (int i=0; i<g_n_sub; i++){ residual[i] -= hatY[i];}	

	
	double correct = 0.0; 
	for (int i=0; i<g_n_snp; i++){
		double sigma_1, sigma_0;
		double v = g_snp_var.at(i);
		double xx = v * g_n_sub;
		if (selected[i] == 1){ // purely guess 
			sigma_1 = sigma * sqrt(calc_alpha());		
			sigma_0 = sqrt(h2 / (1-h2) / (ss - v)  * calc_alpha() );
		}else{
			sigma_0 = sigma * sqrt(calc_alpha());		
			sigma_1 = sqrt(h2 / (1-h2) / (ss + v)  * calc_alpha() );
		}
		double beta_j = 0.0;
		double sum_b2 = 0.0;
		double k1 = subset.size();
		for (int j=0; j<subset.size(); j++){
			if (subset[j] == i){
				beta_j = Beta[j + g_n_cov]; // In fact this Beta should be sampled. But sampling from multivariate normal is too time-consuming.
				k1 -- ;
				continue;
			}
			double b = Beta[j + g_n_cov];
			sum_b2 += (b*b);
		}	
		double d = ( - k1 * log(sigma_1 / sigma_0) - 0.5 * tau * sum_b2 * ( 1.0/sigma_1/sigma_1 - 1.0/sigma_0/sigma_0 )  );
		double c;  // log scale
		// ss = v; only one SNP;
		if ( (sigma_0 != sigma_0) || isinf(sigma_0)){c = 1e10; }
		else{c = log(lambda) + d; }

		float* xj = &g_data[i * g_n_sub];
		double xj_ty = 0.0;
		if (selected[i] == 1){
			copy1D(tilde_y, xj, g_n_sub);
			for (int l=0; l<g_n_sub; l++){
				tilde_y[l] *= beta_j;
				tilde_y[l] += residual[l];
			}
			xj_ty = vec_vec_mul(xj, tilde_y, g_n_sub);
		}else{
			xj_ty = vec_vec_mul(xj, residual, g_n_sub);
		}
	
		double beta_rb = xj_ty/ (xx + 1/sigma_1/sigma_1);
	
		c = c - log(sigma_1) - 0.5 * log(xx + 1.0/sigma_1/sigma_1) + tau * xj_ty * beta_rb / 2; 
		g_log_score[i] = c; 
		double cp = 1.0 / (1.0 + exp(-c));
		g_rb_pip[i] += cp;
		g_rb_beta[i] += beta_rb * cp;
		Beta_rb[i] = beta_rb * cp; 
		
		double vxx = 1 / (xx + 1/sigma_1/sigma_1) / tau;
		double var_beta_rb = cp * (beta_rb * beta_rb + vxx) - cp * cp * beta_rb * beta_rb;
		correct += var_beta_rb * xx; 
	}
	free1D(tilde_y);
	free1D(residual);
	g_rb_count ++;
	
	// in sst: mean should always be zero
	double* yhat = allocate1D(g_n_sub);
	mat_vec_mul(g_data_mat, Beta_rb, g_n_snp, g_n_sub, yhat, 1);
	pve_rb = (correct + array_sst(yhat, g_n_sub)) / g_yy; 
	free1D(yhat);		
	free1D(Beta_rb);

	record_time("RaoBlackwell");
	return;
} 


