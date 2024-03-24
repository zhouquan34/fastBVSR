#include "model2c.h"

using namespace std;

int g_size_g1 = 0;
int g_size_g2 = 0;
map<double, int> g_rho;

model_two_comp::model_two_comp(void){
	sse = 0.0;
	sigma_a = 0.0;
	sigma_b = 0.0;
	size_a = 0;
	size_b = 0;
	pi_a = 0.0;
	pi_b = 0.0;
	h2 = 0.0;
	pve = 0.0;
	rho = 0.0;
	tau = 0.0;
	var_sy = 0.0;
	like = 0.0;
	sbf = 0.0;
	diags = NULL;
	R = NULL;
	RA = NULL;
	XX = NULL;
	X = NULL;
	Xy = NULL;
	Beta = NULL;
	sYhat = NULL;
	sY = NULL;
	hatY = NULL;
}

model_two_comp::~model_two_comp(){
	free1D(diags);
	free1D(R);	
	free1D(RA);	
	free1D(XX);	
	free1D(X);
	free1D(Xy);
	free1D(Beta);
	free1D(sY);
	free1D(sYhat);
	free1D(hatY);
}


int model_two_comp::model_size(void){return subset.size();}

string model_two_comp::model_str(void){
	stringstream s;
	s << subset[0];
	for (int i=1; i<subset.size(); i++){
		s << ", " << subset[i];
	}	
	return s.str();
}

string model_two_comp::para_str(void){
	stringstream s;
	s << size_a << "\t" << size_b << "\t" << h2 << "\t"  << rho <<  "\t" << pve << "\t"  << pi_a << "\t" << pi_b << "\t" << sigma_a << "\t" << sigma_b;
	return s.str();
}

void model_two_comp::update_marginal(int k){
	for (int i=0; i<subset.size(); i++){
		g_pip[subset[i]] += k;
		g_beta[subset[i]] += k * Beta[i + g_n_cov];
	}
	g_model_size[subset.size()] += k;
	g_h2[h2] += k;
	g_rho[rho] += k;
	g_pve[pve] += k;
	return;
}

void model_two_comp::allocate(int p){
	X = allocate1D_float_pointer(p);
	XX = allocate1D(p * (p+1) / 2);
	R = allocate1D(p * p);
	Xy = allocate1D(p);
	Beta = allocate1D(p);
	sY = allocate1D(g_n_sub);
	sYhat = allocate1D(g_n_sub); 
	hatY = allocate1D(g_n_sub);
}

void model_two_comp::initialize(vector<int>& start, double h2_start){
	subset = start;
	h2 = h2_start; 
//	rho = rho_start;
	msize = subset.size() + g_n_cov;
	
	allocate(msize);	
	X_init(subset, X);
	XX_init(subset, X, msize, XX);
	Xy_init(X, msize, Xy);
	Chol_init(msize, XX, R);
	
	selected.clear();
	selected.assign(g_n_snp,0);
	size_a = 0; 
	size_b = 0;
	for (int i=0; i<subset.size(); i++){
		selected[subset[i]] = 1;
		if (subset[i] < g_size_g1){
			size_a ++;
		}else{
			size_b ++;
		}
	}
	
	calc_sigmas();
	calc_rss();
	sample_pi();	
	sample_y();
	sample_pve();
	
	if (g_sbf == 1){
		calc_sbf(g_sbf_perm);
	}
	return;
}

void model_two_comp::copy(model_two_comp* m){  // sbf pve tau var_sy not to be copied
	size_a = m->size_a;
	size_b = m->size_b;
	sigma_a = m->sigma_a;
	sigma_b = m->sigma_b;
	rho = m->rho;
	pi_a = m->pi_a;
	pi_b = m->pi_b;
	h2 = m->h2;
	sse = m->sse;
	like = m->like;
	msize = m->msize;
	subset = m->subset;
	selected = m->selected;
}

void model_two_comp::calc_sigmas(void){
	double ha = 0.0, hb = 0.0; 	
	double ssa = 0.0, ssb = 0.0;
	for (int i=0; i<subset.size(); i++){
		if (subset[i] < g_size_g1){
			ssa += g_snp_var.at(subset[i]);
			ha += g_single_h2.at(subset[i]);
		}else{	
			ssb += g_snp_var.at(subset[i]);
			hb += g_single_h2.at(subset[i]);
		}
	}	
	
	rho = ha / (ha + hb);
	double hh = h2 / (1-h2);
	ha = hh * rho;
	hb = hh - ha;
		
	sigma_a = sqrt(ha / ssa);
	sigma_b = sqrt(hb / ssb);

	diags = allocate1D(msize);
	double d1 = sqrt(1.0/sigma_a/sigma_a - g_diag_c);
	double d2 = sqrt(1.0/sigma_b/sigma_b - g_diag_c);
	for (int i=0; i<g_n_cov; i++){diags[i] = 0;}
	for (int i=0; i<subset.size(); i++){
		if (subset[i] < g_size_g1){
			diags[i + g_n_cov] = d1;
		}else{	
			diags[i + g_n_cov] = d2;
		}
	}	
	
	return;
}

void model_two_comp::calc_rss(void){
	calc_ridge(RA, R, Xy, msize, diags, Beta, XX);	
	sse = g_yy - vec_vec_mul(Xy, Beta, msize); 
	like = -0.5 * g_n_sub * log(sse);
	return;
}

double model_two_comp::calc_rss(double* tilde_y){
	double* z1 = allocate1D(msize);
	record_time();
	mat_vec_mul(X, tilde_y, msize, g_n_sub, z1, 0);
	record_time("X-y");

	double tyy = vec_vec_mul(tilde_y, tilde_y, g_n_sub);
	double* beta1 = allocate1D(msize);
	calc_ridge(RA, R, z1, msize, diags, beta1, XX);
	double rss = tyy - vec_vec_mul(z1, beta1, msize);
	free1D(z1);
	free1D(beta1);
	return -0.5 * g_n_sub * log(rss);
}

void model_two_comp::sample_pve(void){
	record_time();	
	mat_vec_mul(X, Beta, msize, g_n_sub, hatY, 1);
	record_time("X-beta");
	
	sample_tau();
	double sum_yy = 0; 
	for (int i=0; i<g_n_sub; i++){
		double y = hatY[i] + gsl_ran_gaussian(gsl_r, sqrt(((double) subset.size())/(g_n_sub * tau)) );
		sum_yy += y * y;
	}
	double v1 = sum_yy * tau / (double) g_n_sub;
	pve = v1 / (1.0 + v1);
//	pve = array_sst(hatY, g_n_sub) / g_yy;
	return;
}

double model_two_comp::sample_h2(void){
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
	
	// sample rho
/*	rho = rho + 2.0 * g_jump_rho * gsl_rng_uniform(gsl_r) - g_jump_rho; 
	if (rho < 0){ rho = -rho; }
	if (rho > 1){ rho = 2.0 - rho; }
*/	
	return r;
}


void model_two_comp::sample_pi(void){
//	cout << size_a << "\t" << size_b << "\t" << g_size_g1 << "\t" << g_size_g2 << endl;
	pi_a = sample_beta(size_a, g_size_g1 - size_a + 1.0);
	pi_b = sample_beta(size_b, g_size_g2 - size_b + 1.0);
/*	while (1){
		double p = gsl_ran_beta(gsl_r, size_a, g_size_g1 - size_a + 1.0);
		if (p > g_min_pi && p < g_max_pi){	
			pi_a = p;
			break;
		}
	}
	while (1){
		double p = gsl_ran_beta(gsl_r, size_b, g_size_g2 - size_b + 1.0);
		if (p > g_min_pi && p < g_max_pi){	
			pi_b = p;
			break;
		}
	}
*/
//	cout << pi_b << endl;
}


double model_two_comp::add_one_snp(vector<int>& add){
	double r = gsl_rng_uniform(gsl_r);
	int pick = -1;
	int pick_rank = -1;
	if (r < g_prop_uniform){ // uniform
		pick = gsl_rng_uniform_int(gsl_r, g_n_snp);
		while (selected[pick] == 1){ pick = gsl_rng_uniform_int(gsl_r, g_n_snp); }
		pick_rank = g_single_ord[pick];
		int k = pick_rank;
		for (int i=0; i<subset.size(); i++){
			if ( g_single_ord[subset[i]] <= k){ 
				pick_rank --; 
			}
		}		
	}else{ // geometric  Pr(K = k) = 0.02 * 0.98^k (truncated)
		double r2 = log(gsl_rng_uniform(gsl_r))/log(1 - g_geometric);
		int k = (int) r2;
		while (k >= g_n_snp - subset.size() ){
			r2 = log(gsl_rng_uniform(gsl_r))/log(1 - g_geometric);
			k = (int) r2;
		}
		pick_rank = k;
		// perhaps we can also save the order of the current SNPs in the model. But I think it won't be more efficient (we need copy).
		set<int> tmp;
		int max_k = k + subset.size();
		for (int i=0; i<subset.size(); i++){
			int m = g_single_ord[subset[i]];
			if ( m <= k){ 
				k ++; 
				while ( tmp.find(k) != tmp.end() ){
					k++;
				}
			}
			if ( m > k && m <= max_k ){
				tmp.insert(m);
			}
		}		
		pick = g_single_id[k];
	}	
	selected[pick] = 1;
	subset.push_back(pick);
	add.push_back(pick);
	msize ++;
	if (pick < g_size_g1){size_a ++;}
	else {size_b ++;}

	double sum_geom = 1 - pow(1-g_geometric, g_n_snp - subset.size() + 1); 
	double pa = log(g_prop_uniform) - log(g_n_snp - subset.size() + 1);
	double pb = log(1-g_prop_uniform) - log(sum_geom) + log(g_geometric) + pick_rank * log(1-g_geometric);
	double q_new2old = -log(subset.size());
	double q_old2new = log_sum_log(pa, pb); 
	
	return q_new2old - q_old2new;
}

double model_two_comp::add_snps (int add_size, model_two_comp* last_model){
	if (subset.size() + add_size > g_max_model_size){return FAIL_ALPHA;} 
	vector<int> add;
	int ssa = size_a;	
	int ssb = size_b;
	double alpha = 0.0;
	allocate(last_model->msize + add_size);
	for (int i=0; i<add_size; i++){
		alpha += add_one_snp(add); // this is the proposal ratio of gamma
	}
	X_add(add, last_model->msize, last_model->X, X);
	XX_add(subset, add, last_model->msize, X, last_model->XX, XX);
	Xy_add(add, last_model->msize, last_model->Xy, Xy);
	Chol_add(add_size, last_model->msize, XX, last_model->R, R);
	
//	cout << g_mcmc_i << "\tAdd\t" << ssa << "\t" << ssb << "\t" << size_a << "\t" << size_b << endl;
	int add_a = size_a - ssa;
	int add_b = size_b - ssb;
	if (add_a > 0){
		alpha += (gsl_sf_lnpoch(ssa + g_pi_alpha_2, add_a) - gsl_sf_lnpoch(g_size_g1 + g_pi_beta_2 - ssa - add_a, add_a));
	}
	if (g_pi_range_set_2 == 1){
		int a1 = ssa + add_a + g_pi_alpha_2;
		int b1 = g_size_g1 - ssa - add_a + g_pi_beta_2;
		int a2 = ssa + g_pi_alpha_2;
		int b2 = g_size_g1 - ssa + g_pi_beta_2;
		alpha += log( (gsl_cdf_beta_P(g_max_pi, a1, b1) - gsl_cdf_beta_P(g_min_pi, a1, b1)) / (gsl_cdf_beta_P(g_max_pi, a2, b2) - gsl_cdf_beta_P(g_min_pi, a2, b2)) ); 
	}		
	if (add_b > 0){
		alpha += (gsl_sf_lnpoch(ssb + g_pi_prior_alpha, add_b) - gsl_sf_lnpoch(g_size_g2 + g_pi_prior_beta - ssb - add_b, add_b));
	}
	if (g_pi_range_set == 1){
		int a1 = ssb + add_b + g_pi_prior_alpha;
		int b1 = g_n_snp - ssb - add_b + g_pi_prior_beta;
		int a2 = ssb + g_pi_prior_alpha;
		int b2 = g_n_snp - ssb + g_pi_prior_beta;
		alpha += log( (gsl_cdf_beta_P(g_max_pi, a1, b1) - gsl_cdf_beta_P(g_min_pi, a1, b1)) / (gsl_cdf_beta_P(g_max_pi, a2, b2) - gsl_cdf_beta_P(g_min_pi, a2, b2)) ); 
	}		
	return alpha;
}


double model_two_comp::delete_one_snp(vector<int>& del){
	int pick_r = gsl_rng_uniform_int(gsl_r, subset.size());
	int pick = subset[pick_r];			
	int pick_rank = g_single_ord[pick];
	int k = pick_rank;
	for (int i=0; i<subset.size(); i++){
		if ( g_single_ord[subset[i]] < k){   // skip m == k 
			pick_rank --; 
		}
	}

	msize --;
	if (subset[pick_r] < g_size_g1){size_a --;}
	else{size_b --;}
	selected[pick] = 0;
	subset.erase(subset.begin() + pick_r);

	double sum_geom = 1 - pow(1-g_geometric, g_n_snp - subset.size() ); 
	double pa = log(g_prop_uniform) - log(g_n_snp - subset.size());
	double pb = log(1-g_prop_uniform) - log(sum_geom) + log(g_geometric) + pick_rank * log(1-g_geometric);
	double q_old2new = -log(subset.size() + 1);
	double q_new2old = log_sum_log(pa, pb); 

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

	return q_new2old - q_old2new;
}

double model_two_comp::delete_snps (int del_size, model_two_comp* last_model){	
	if (subset.size() < (del_size + 1) ){return FAIL_ALPHA;}	
	double alpha = 0.0;
	int ssa = size_a; 
	int ssb = size_b;
	allocate(last_model->msize - del_size);
	vector<int> del;
	for (int i=0; i<del_size; i++){
		alpha += delete_one_snp(del); // this is the proposal ratio of gamma
	}
	if (size_a == 0 || size_b == 0){return FAIL_ALPHA;}

	sort(del.begin(), del.end()); // important. 
	X_del(del, last_model->msize, last_model->X, X);
	XX_del(del, last_model->msize, last_model->XX, XX);
	Xy_del(del, last_model->msize, last_model->Xy, Xy);
	Chol_del(del, last_model->msize, last_model->R, R);	
	
//	cout << g_mcmc_i << "\tDel\t" << ssa << "\t" << ssb << "\t" << size_a << "\t" << size_b << endl;
	int del_a = ssa - size_a; 
	int del_b = ssb - size_b;
	if (del_a > 0){
		alpha += ( -gsl_sf_lnpoch(ssa + g_pi_alpha_2 - del_a, del_a) + gsl_sf_lnpoch(g_size_g1 + g_pi_beta_2 - ssa, del_a) );
	}
	if (g_pi_range_set_2 == 1){
		int a1 = ssa - del_a + g_pi_alpha_2;
		int b1 = g_size_g1 - ssa + del_a + g_pi_beta_2;
		int a2 = ssa + g_pi_alpha_2;
		int b2 = g_size_g1 - ssa + g_pi_beta_2;
		alpha += log( (gsl_cdf_beta_P(g_max_pi, a1, b1) - gsl_cdf_beta_P(g_min_pi, a1, b1)) / (gsl_cdf_beta_P(g_max_pi, a2, b2) - gsl_cdf_beta_P(g_min_pi, a2, b2)) ); 
	}		
	if (del_b > 0){
		alpha += ( -gsl_sf_lnpoch(ssb + g_pi_prior_alpha - del_b, del_b) + gsl_sf_lnpoch(g_size_g2 + g_pi_prior_beta - ssb, del_b) );
	}
	if (g_pi_range_set == 1){
		int a1 = ssb - del_b + g_pi_prior_alpha;
		int b1 = g_n_snp - ssb + del_b + g_pi_prior_beta;
		int a2 = ssb + g_pi_prior_alpha;
		int b2 = g_n_snp - ssb + g_pi_prior_beta;
		alpha += log( (gsl_cdf_beta_P(g_max_pi, a1, b1) - gsl_cdf_beta_P(g_min_pi, a1, b1)) / (gsl_cdf_beta_P(g_max_pi, a2, b2) - gsl_cdf_beta_P(g_min_pi, a2, b2)) ); 
	}		

	return alpha;
}


void model_two_comp::sample_y (void){
	double* sbeta = allocate1D(subset.size()); // beta for cov is zero	
	for (int i=0; i<subset.size(); i++){ 
		if (subset[i] < g_size_g1){
			sbeta[i] = gsl_ran_gaussian(gsl_r, sigma_a);
		}else{
			sbeta[i] = gsl_ran_gaussian(gsl_r, sigma_b);
		}
	}
	
	record_time();	
	mat_vec_mul(&X[g_n_cov], sbeta, subset.size(), g_n_sub, sYhat, 1);
	record_time("X-beta");
	
	for (int i=0; i<g_n_sub; i++){ sY[i] = sYhat[i] +  gsl_ran_gaussian(gsl_r, 1.0);}		
	var_sy = center_array(sY, g_n_sub);
	free1D(sbeta);
	return;
}

void model_two_comp::calc_sbf(int n){
	// permute n times; assume y is already the residual after regressing out W.
	double* tilde_y = allocate1D(g_n_sub);
	copy1D(tilde_y, g_pheno, g_n_sub);
	
	double sum = 0.0;
	for (int i=0; i<n; i++){
		gsl_ran_shuffle(gsl_r, tilde_y, g_n_sub, sizeof(double));
		double* z1 = allocate1D(msize);
		mat_vec_mul(X, tilde_y, msize, g_n_sub, z1, 0);
		double* beta1 = allocate1D(msize);
		calc_ridge(RA, R, z1, msize, diags, beta1, XX);
		double rss = g_yy - vec_vec_mul(z1, beta1, msize);
		sum += (-0.5 * g_n_sub * log(rss));
		free1D(z1);
		free1D(beta1);
	}
	
	free1D(tilde_y);
	sbf = like - sum / (double) n;
	return;
}

int model_two_comp::compare_models(double alpha, model_two_comp* last_model){
	calc_sigmas();
	calc_rss();
	double mh = alpha;
	mh += (like - last_model->like);
	sample_y();
	double f_new = calc_rss(sY);
	double f_old = last_model->calc_rss(sY); // These two calculations could be faster since they are almost the same. Might save 2% time.
	mh += (f_old - f_new);

	double r = gsl_rng_uniform(gsl_r);
	if (exp(mh) > r){
		return 1; // accept
	}else{
		return 0; // reject
	}
}

int model_two_comp::compare_models_sbf(double alpha, model_two_comp* last_model){
	calc_sigmas();
	calc_rss();
	calc_sbf(g_sbf_perm);

	double mh = alpha + sbf - last_model->sbf;
	double r = gsl_rng_uniform(gsl_r);
	if (exp(mh) > r){
		return 1; // accept
	}else{
		return 0; // reject
	}
}


void model_two_comp::sample_tau(void){
	tau = gsl_ran_gamma(gsl_r, g_n_sub/2, 2.0/sse);	
	return;
}


void model_two_comp::rao_blackwell (void){
	record_time();
//	sample_tau();

	double lambda_a = pi_a/(1-pi_a);
	double lambda_b = pi_b/(1-pi_b);
	double hh = h2 / (1-h2);
	double ha = hh * rho;
	double hb = hh - ha;
	double ss_a = ha /( sigma_a * sigma_a);
	double ss_b = hb /( sigma_b * sigma_b);

	double* tilde_y = allocate1D(g_n_sub);
	double* residual = allocate1D(g_n_sub);
	copy1D(residual, g_pheno, g_n_sub);
	for (int i=0; i<g_n_sub; i++){ residual[i] -= hatY[i];}

	for (int i=0; i<g_n_snp; i++){
		double sigma, hh, ss, lambda;
		if (i < g_size_g1){
			sigma = sigma_a;
			hh = ha;
			ss = ss_a;
			lambda = lambda_a;
		}else{
			sigma = sigma_b;
			hh = hb;
			ss = ss_b;
			lambda = lambda_b;
		}
	
		double sigma_1, sigma_0;
		double v = g_snp_var.at(i);
		double xx = v * g_n_sub;
		if (selected[i] == 1){
			sigma_1 = sigma;		
			sigma_0 = sqrt(hh / (1-hh) / (ss - v) );
		}else{
			sigma_0 = sigma;		
			sigma_1 = sqrt(hh / (1-hh) / (ss + v) );
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
		double c; 
		// ss = v; only one SNP; 
		if (sigma_0 != sigma_0) {  c = 1.0 / 0.0;}
		else{	c = lambda * exp(d);}

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
		c = c / sigma_1 / sqrt(xx + 1/sigma_1/sigma_1) * exp(tau * xj_ty / 2 );
		if (isinf(c) == 1){ // check inf
			g_rb_pip[i] += 1.0;
			g_rb_beta[i] += beta_rb;
		}else{
			g_rb_pip[i] += c / (1+c);
			g_rb_beta[i] += beta_rb * c / (1+c);
		}
	}


	free1D(tilde_y);
	free1D(residual);
	g_rb_count ++;
	
	record_time("RaoBlackwell");
	return;
} 






