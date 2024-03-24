#include "model.h"

using namespace std;
// int g_fix_sigma = 1;

model::model(void){
	sse = 0.0;
	sigma = 0.0;
	pi = 0.0;
	h2 = 0.0;
	pve = 0.0;
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
	free1D(X);
	free1D(Xy);
	free1D(Beta);
	free1D(sY);
	free1D(sYhat);
	free1D(hatY);
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
	s << subset.size() << "\t" << h2 << "\t"  << pve << "\t" << pi << "\t" << sigma << "\t" << br2;
//	s << "\t" << sbf;
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
	
	if (g_sbf == 1){
		calc_sbf(g_sbf_perm);
	}
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
/*	if (g_fix_sigma == 1){
		sigma = 0.2; 
		double s1 = ss * sigma * sigma;
		h2 = s1 / (1.0 + s1);
		return;
	}else{
*/	
	sigma = sqrt(h2 / (1-h2) / ss);
	return;
}

void model::calc_rss(void){
	double log_det = calc_ridge(RA, R, Xy, msize, sqrt(1.0/sigma/sigma - g_diag_c), Beta, XX);	
//	record_time();
	sse = g_yy - vec_vec_mul(Xy, Beta, msize); 
	like = -0.5 * g_n_sub * log(sse);
	if (g_exact_bf == 1){
//		if (log_det == 0){cout << "Wrong!" << endl;}
		like += ( (-0.5 * log_det)  - subset.size() * log(sigma) );
//		cout << like << "\t" << sigma << "\t" << log_det << "\t"; vector_print(subset);
	}
//	record_time("CalcRSS");
	return;
}

double model::calc_rss(double* tilde_y){
	double* z1 = allocate1D(msize);
	record_time();
	mat_vec_mul(X, tilde_y, msize, g_n_sub, z1, 0);
	record_time("X-y");

	double tyy = vec_vec_mul(tilde_y, tilde_y, g_n_sub);
	double* beta1 = allocate1D(msize);
	calc_ridge(RA, R, z1, msize, sqrt(1.0/sigma/sigma - g_diag_c), beta1, XX);
	double rss = tyy - vec_vec_mul(z1, beta1, msize);
	free1D(z1);
	free1D(beta1);
	return -0.5 * g_n_sub * log(rss);
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
	pi = sample_beta(subset.size(), g_pseudo_n - subset.size() + 1.0);
/*
	while (1){
		double p = gsl_ran_beta(gsl_r, subset.size(), g_n_snp - subset.size() + 1.0);
		if (p > g_min_pi && p < g_max_pi){	
			pi = p;
			break;
		}
	}
*/
	record_time("SamplePI");
}


double model::add_one_snp(vector<int>& add){
//	record_time();
//	cout << "Add: " ; vector_print(subset);
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
	
	double sum_geom = 1 - pow(1-g_geometric, g_n_snp - subset.size() + 1); 
	double pa = log(g_prop_uniform) - log(g_n_snp - subset.size() + 1);
	double pb = log(1-g_prop_uniform) - log(sum_geom) + log(g_geometric) + pick_rank * log(1-g_geometric);
	double q_new2old = -log(subset.size());
	double q_old2new = log_sum_log(pa, pb); 
	
//	cout << pick << "\t" << sum_geom << "\t" << pa << "\t" << pb << "\t" << q_new2old << "\t" << q_old2new << endl;
//	record_time("AddSNP");
	return q_new2old - q_old2new;
}

void gsl_matrix_print (const gsl_matrix* m){
	for (int i=0; i<m->size1; i++){
		for (int j=0; j<m->size2; j++){
			cout << gsl_matrix_get(m, i, j) << "\t";	
		}	
		cout << endl;
	}	
	return;
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
/*	if (subset.size() == 3){
		cout << "\n" << msize << "\n";
		print_array(XX, msize * (msize + 1)/2);
		print_matrix(R, msize * msize, msize);
	}
*/
/*	if (g_mcmc_i > 10000){
		cout << "Iter " << g_mcmc_i << endl;
		print_matrix(R, msize * msize, msize);	
		double* R_tmp_1 = allocate1D(msize * msize);
		Chol_decomp_A(msize, XX, R_tmp_1, 0.0);
		cout << "Method 2" << endl;
		print_matrix(R_tmp_1, msize*msize, msize);
		gsl_matrix* R_tmp_2 = gsl_matrix_alloc(msize, msize);
		int i=0, j=0;
		for (int k=0; k<msize*(msize+1)/2; k++){
			gsl_matrix_set(R_tmp_2, i, j,  XX[k]);
			j ++ ;
			if (j > i){
				j = 0;
				i ++;
			}
		} 
		gsl_linalg_cholesky_decomp1(R_tmp_2);
		cout << "Method 3" << endl;
		gsl_matrix_print(R_tmp_2);
	}
*/	

//	gsl_rng_uniform(gsl_r);
	if (g_mcmc_i % 50000 == 0){
	//	cout << "Reset Chol" << endl;
		Chol_decomp_A(msize, XX, R, 0.0);
	}
	return alpha;
}


double model::delete_one_snp(vector<int>& del){
//	record_time();
//	cout << "Del: " ; vector_print(subset);
	int pick_r = gsl_rng_uniform_int(gsl_r, subset.size());
	int pick = subset[pick_r];			
	int pick_rank = g_single_ord[pick];
	int k = pick_rank;
	for (int i=0; i<subset.size(); i++){
		if ( g_single_ord[subset[i]] < k){   // skip m == k 
			pick_rank --; 
		}
	}

	selected[pick] = 0;
	subset.erase(subset.begin() + pick_r);
	msize --;

	double sum_geom = 1 - pow(1-g_geometric, g_n_snp - subset.size() ); 
	double pa = log(g_prop_uniform) - log(g_n_snp - subset.size());
	double pb = log(1-g_prop_uniform) - log(sum_geom) + log(g_geometric) + pick_rank * log(1-g_geometric);
	double q_old2new = -log(subset.size() + 1);
	double q_new2old = log_sum_log(pa, pb); 
//	cout << pick << "\t" << sum_geom << "\t" << pa << "\t" << pb << "\t" << q_new2old << "\t" << q_old2new << endl;
	
	set<int> tmp;
//	int del_r = pick_r + g_n_cov;
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

//	record_time("AddSNP");
//	record_time("DelSNP");
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
	for (int i=0; i<subset.size(); i++){ sbeta[i] = gsl_ran_gaussian(gsl_r, sigma);}
	
	record_time();	
	mat_vec_mul(&X[g_n_cov], sbeta, subset.size(), g_n_sub, sYhat, 1);
	record_time("X-beta");
	
	for (int i=0; i<g_n_sub; i++){ sY[i] = sYhat[i] +  gsl_ran_gaussian(gsl_r, 1.0);}		
	var_sy = center_array(sY, g_n_sub);
	free1D(sbeta);
	return;
}

void model::calc_sbf(int n){
	// permute n times; assume y is already the residual after regressing out W.
	double* tilde_y = allocate1D(g_n_sub);
	copy1D(tilde_y, g_pheno, g_n_sub);
	
	double sum = 0.0;
	for (int i=0; i<n; i++){
		gsl_ran_shuffle(gsl_r, tilde_y, g_n_sub, sizeof(double));
		double* z1 = allocate1D(msize);
		mat_vec_mul(X, tilde_y, msize, g_n_sub, z1, 0);
		double* beta1 = allocate1D(msize);
		calc_ridge(RA, R, z1, msize, sqrt(1.0/sigma/sigma - g_diag_c), beta1, XX);
		double rss = g_yy - vec_vec_mul(z1, beta1, msize);
		sum += (-0.5 * g_n_sub * log(rss));
		free1D(z1);
		free1D(beta1);
	}
	
	free1D(tilde_y);
	double sbf_null = sum / (double) n;
	sbf = like - sbf_null;
	double mean_scale = calc_mean_scale(sigma);
//	cout << sbf << "\t" << sbf_null << "\t" << sigma << "\t" << mean_scale << "\t";
	sbf = sbf - subset.size() * log(mean_scale);
	return;
}

int model::compare_models(double alpha, model* last_model){
	calc_sigma();
	calc_rss();
	double mh = alpha;
	mh += (like - last_model->like);
	if (g_exact_bf == 0){
		sample_y();
		double f_new = calc_rss(sY);
		double f_old = last_model->calc_rss(sY); // These two calculations could be faster since they are almost the same. Might save 2% time.
		mh += (f_old - f_new);   //cout << alpha << "\t" << f_old << "\t" << f_new << "\t" << mh << endl;
	}

	double r = gsl_rng_uniform(gsl_r);
	if (exp(mh) > r){
		return 1; // accept
	}else{
		return 0; // reject
	}
}

int model::compare_models_sbf(double alpha, model* last_model){
	calc_sigma();
	calc_rss();
	calc_sbf(g_sbf_perm);

	double mh = alpha + sbf - last_model->sbf;
	double r = gsl_rng_uniform(gsl_r);
	if (exp(mh) > r){
//		cout << "\t" << last_model->sbf << "\t" << 1 << endl;
		return 1; // accept
	}else{
//		cout << "\t" << last_model->sbf << "\t" << 0 << endl;
		return 0; // reject
	}
}


void model::sample_tau(void){
	tau = gsl_ran_gamma(gsl_r, g_n_sub/2, 2.0/sse);	
	return;
}


void model::rao_blackwell (void){
	record_time();
//	sample_tau();

	double lambda = pi/(1-pi);
	double ss = h2 /( (1-h2) * sigma * sigma);
	double* tilde_y = allocate1D(g_n_sub);
	double* residual = allocate1D(g_n_sub);
	copy1D(residual, g_pheno, g_n_sub);
	for (int i=0; i<g_n_sub; i++){ residual[i] -= hatY[i];}

	for (int i=0; i<g_n_snp; i++){
		double sigma_1, sigma_0;
		double v = g_snp_var.at(i);
		double xx = v * g_n_sub;
		if (selected[i] == 1){
			sigma_1 = sigma;		
			sigma_0 = sqrt(h2 / (1-h2) / (ss - v) );
		}else{
			sigma_0 = sigma;		
			sigma_1 = sqrt(h2 / (1-h2) / (ss + v) );
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
		if ( (sigma_0 != sigma_0) || isinf(sigma_0)  ){  c = 1.0 / 0.0;}
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
		if (g_sbf == 0){
			c = c / sigma_1 / sqrt(xx + 1/sigma_1/sigma_1) * exp(tau * xj_ty * beta_rb / 2 );
		}else{
			c = c * exp (-0.5 * xx / (xx + 1/sigma_1/sigma_1) ) * exp(tau * xj_ty * beta_rb/ 2 );
			double mean_scale = calc_mean_scale(sigma_1);
			c = c / mean_scale;
		}
//		if (c != c){
//			cout << c << "\t" << sigma_0 << "\t" << sigma_1 << "\t" << lambda << "\t" << d << "\t" << xj_ty << endl; 
//		}
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
//	free1D(hatY);
	g_rb_count ++;
	
	record_time("RaoBlackwell");
	return;
} 






