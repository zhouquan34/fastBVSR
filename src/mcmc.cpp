#include "mcmc.h"

using namespace std;


void proposal_type (int& add, int& del, int& prop_type){
	prop_type = 0;
	if (gsl_rng_uniform(gsl_r) < g_long_range){ 
		prop_type = 2;
		int csize = gsl_rng_uniform_int(gsl_r, g_max_jump - 1) + 2; 
		add = (int) (csize * gsl_rng_uniform(gsl_r) + 0.5);
		del = csize - add;
	}else{
		if (gsl_rng_uniform(gsl_r) > g_prop_add){prop_type = 1;}
	}
}

void mcmc_sample(model*& last_model){
	int add = 0, del = 0;
	int prop_type = 0;
	proposal_type(add, del, prop_type);
//	cout << add << " " << del << " " << prop_type << endl;
	double alpha = 0.0;
	model* proposed_model = new model();
	if (prop_type == 0){ // add one 
		proposed_model->copy(last_model);
		alpha = proposed_model->add_snps(1, last_model);	
	}else if (prop_type == 1){ // delete one 
		proposed_model->copy(last_model);	
		alpha = proposed_model->delete_snps(1, last_model);
	}else{ // long-range
		if (add == 0){
			proposed_model->copy(last_model);
			alpha = proposed_model->delete_snps(del, last_model);	
		}else if(del == 0){
			proposed_model->copy(last_model);
			alpha = proposed_model->add_snps(add, last_model);	
		}else{
			model* tmp_model = new model();
			tmp_model -> copy(last_model);	
			double alpha_1 = tmp_model->add_snps(add, last_model);
			if (alpha_1 == FAIL_ALPHA){
				alpha = FAIL_ALPHA;
			}else{
				proposed_model -> copy(tmp_model);
				double alpha_2 = proposed_model->delete_snps(del, tmp_model);
				if (alpha_2 == FAIL_ALPHA){
					alpha = FAIL_ALPHA;
				}else{
					alpha += (alpha_1 + alpha_2);
				}
			}
			delete(tmp_model);
		}
	}


	int acc = 0; 
	if (alpha > FAIL_ALPHA){
		double ah = proposed_model->sample_h2();
		if (g_sbf == 0){
			acc = proposed_model->compare_models(alpha + ah, last_model);
		}else{
			acc = proposed_model->compare_models_sbf(alpha + ah, last_model);
		}
	}

//	string model = proposed_model->model_str();
//	cout << g_mcmc_i << "\t" << csize << "\t" << prop_type << "\tM = " << model <<  endl;

	if (acc == 1){	
		if (g_mcmc_i > 0 ){
			string paras = last_model->para_str();
			string model = last_model->model_str();
			PATH << g_mcmc_i - 1 << "\t" << g_mcmc_stay << "\t"  << paras << endl;
			MOD << model << endl;
			if (g_output_yhat == 1){
				for (int i=0; i<g_n_sub; i++){ YHAT << last_model->hatY[i] << " ";}
				YHAT << endl; 
			}
		}
		proposed_model->sample_pi(); // pi can be sampled at last. It doesn't enter the calculation of MH
		proposed_model->sample_pve();
		if (g_mcmc_i > g_mcmc_warm) {
			int s = g_mcmc_i - g_mcmc_warm - 1; 
			if (s > g_mcmc_stay) {s = g_mcmc_stay;}
			last_model -> update_marginal(s);
		}
		delete(last_model);
		last_model = proposed_model;
		g_mcmc_stay = 1;
	}else{
		delete(proposed_model);
		g_mcmc_stay ++ ;
	}
	
	
	if (g_mcmc_i>g_mcmc_warm){
	//	if (add + del > 1 && prop_type < 2){prop_type += 3;}
		g_propose[prop_type] ++;
		if (acc == 1) {g_accept[prop_type] ++;}
	}	
	
//	cout << "Leaving" << endl;
	return;
}


void mcmc_start(model*& m){
	vector<int> s;
	if (g_start_size > 0){
		for (int i=0; i<g_start_size; i++){
			s.push_back(g_single_id[i]);
		}
		m->initialize(s);
		LOG << "Starting model size = " << s.size() << endl;
	}else{
		s.push_back(g_single_id[0]);
		m->initialize(s);
		double g1 = g_prop_uniform;
		double g2 = g_geometric; 
		g_prop_uniform = 0.001;
		g_geometric = 0.02;
		g_prop_add = 1;
		if (g_mcmc_init_try == -1){
			g_mcmc_init_try = 1000;
			if (g_sbf == 1){g_mcmc_init_try = 200;}
		}
		for (int i = 1; i<= g_mcmc_init_try; i++){
			mcmc_sample(m);
			if (m->model_size()>=g_n_snp){break;}
		}
		g_prop_uniform = g1;
		g_geometric = g2;
		g_prop_add =0.5;	
		LOG << "Starting model size = " << m->model_size() << endl;	
	}
	g_chol_call=0;
	return;
}

void calc_sing_distr(void){
	record_time();
	int* snps = new int[g_n_snp];
	int test_p = 100;
	if (test_p > g_n_snp) {test_p = g_n_snp / 2;}
	for (int i=0; i<g_n_snp; i++){snps[i] = i;}
	
	gsl_matrix* X = gsl_matrix_alloc(g_n_sub, test_p);
	gsl_matrix* M = gsl_matrix_alloc(test_p, test_p);
	gsl_matrix* V = gsl_matrix_alloc(test_p, test_p);
	gsl_vector* S = gsl_vector_alloc(test_p);
	gsl_vector* work = gsl_vector_alloc(test_p);
	gsl_matrix_set_zero(X);
	
	vector<double> dd;
	int nsim = 1000;
	for (int sim = 0; sim < nsim; sim++){	
//		cout << sim << endl;
		int* sel_snp = new int[test_p];	
		gsl_ran_choose(gsl_r, sel_snp, test_p, snps, g_n_snp, sizeof(int));
		gsl_matrix_set_zero(M);
		gsl_matrix_set_zero(V);
		gsl_vector_set_zero(S);
		gsl_vector_set_zero(work);
		for (int j=0; j<test_p; j++){
			int s = sel_snp[j];
			for (int i=0; i<g_n_sub; i++){
				gsl_matrix_set(X, i, j, g_data[s * g_n_sub + i]);
			}
		}		
		gsl_linalg_SV_decomp_mod(X, M, V, S, work);		
		for (int k=0; k<test_p; k++){
			dd.push_back(gsl_vector_get(S, k) * gsl_vector_get(S, k));
		}
	}
	
	double dd_s1 = 0.0, dd_s2 = 0.0;
	for (int i=0; i<dd.size(); i++){
		dd_s1 += dd[i];
		dd_s2 += (dd[i] * dd[i]);
	}
	double nl  = (double) nsim * test_p; 
	double dd_mean = dd_s1 / nl;
	g_dd_moments.assign(4, 0);
	g_dd_moments[0] = dd_mean;
	for (int i=0; i<dd.size(); i++){
		for (int j=1; j<=3; j++){
			g_dd_moments[j] += pow(dd[i] - dd_mean, j + 1);
		}
	}		
	for (int j=1; j<=3; j++){
		g_dd_moments[j] /= nl; 
	}

//	cerr << "Mean dd = " << g_dd_mean << "\tVar dd = " << g_dd_var << "\n";
	gsl_matrix_free(X);
	gsl_matrix_free(M);
	gsl_matrix_free(V);
	gsl_vector_free(S);
	gsl_vector_free(work);	
	record_time("EstSingularDistr");
	return;
}

void mcmc_init(model*& model){
	if (g_max_model_size > (g_n_snp-1)) {g_max_model_size = g_n_snp - 1;}
	if (g_min_model_size < 1){g_min_model_size = 1;}
	if (g_pseudo_n == 0) {g_pseudo_n = g_n_snp;}
	if (g_exact_bf == 1) {g_icf_min_p = g_n_snp + 1;}
	g_log_max_h = log(g_max_h);
	g_log_min_h = log(g_min_h);

	g_pip.clear(); g_pip.assign(g_n_snp, 0);
	g_rb_pip.clear(); g_rb_pip.assign(g_n_snp, 0.0);
	g_beta.clear(); g_beta.assign(g_n_snp, 0.0);
	g_rb_beta.clear(); g_rb_beta.assign(g_n_snp, 0.0);	
	g_propose.assign(5,0);
	g_accept.assign(5,0);

	if (g_sbf == 1){
		calc_sing_distr();	
	}

	if (g_continue == 0){
		mcmc_start(model);
	}else{
	//	cout << "g_last_iter = " << g_last_iter << endl;
		vector<int> start_subset; double start_h2;
		read_last_bvsr(start_subset, start_h2, g_last_iter);
		model->initialize(start_subset, start_h2);	
	}
	return;
}

int calc_refresh(void){
	int u = 0;
	if (g_continue == 0){
		u = (g_mcmc_warm + g_mcmc_iter) / 50;
	}else{
		u = g_mcmc_iter / 50;
	}
	if (u <= 0){u = 1;}
	return u;
}

void calc_start_end(int& s, int& e){
	if (g_continue == 0){
		s = 1;
		e = g_mcmc_iter + g_mcmc_warm;
	}else{
		s = g_last_iter + g_mcmc_warm + 1;
		e = g_mcmc_iter + g_mcmc_warm + g_last_iter;
	}
}

void mcmc (void){
	cerr << "Initializing MCMC" << endl;
	model* my_model = new model();
	mcmc_init(my_model);
	int print_unit = calc_refresh();
	g_mcmc_stay = 0;
	cerr << "MCMC sampling" << endl;
	PATH << "No.\tLife\tModelSize\th2\tPVE\tpi\tsigma\tBayesianR2" << endl;
	MOD << "#Iteration index matched with the PATH file" << endl;
//	if (g_output_yhat == 1){YHAT << "#Iteration index matched with the PATH file" << endl;}
	int start_iter = 0, total_iter = 0;
	calc_start_end(start_iter, total_iter);
	for (g_mcmc_i = start_iter; g_mcmc_i <= total_iter; g_mcmc_i ++){
		if (g_mcmc_i % print_unit == 0){cerr << '='; fflush(stdout); }
		mcmc_sample(my_model);
		if (g_mcmc_i > g_mcmc_warm){
			if (g_mcmc_i % g_rb_thin == 0){ my_model->rao_blackwell();}
		}
	}
	string paras = my_model->para_str();
	string model = my_model->model_str();
	PATH << g_mcmc_i - 1  << "\t" << g_mcmc_stay  << "\t" << paras << endl;	
	MOD << model << endl;
	if (g_output_yhat == 1){
		for (int i=0; i<g_n_sub; i++){ YHAT << my_model->hatY[i] << " ";}
		YHAT << endl; 
	}

	my_model -> update_marginal(g_mcmc_stay);
	delete(my_model);	
	mcmc_output();	
	return;
}	



// two component
void mcmc_sample(model_two_comp*& last_model){
	int add = 0, del = 0;
	int prop_type = 0;
	proposal_type(add, del, prop_type);
//	cout << g_mcmc_i << "\t" << prop_type << "\t" << csize << "\n";
	double alpha = 0.0;
	model_two_comp* proposed_model = new model_two_comp();
	if (prop_type == 0){ // add one or more
		proposed_model->copy(last_model);
		alpha = proposed_model->add_snps(1, last_model);	
	}else if (prop_type == 1){ // delete one or more
		proposed_model->copy(last_model);	
		alpha = proposed_model->delete_snps(1, last_model);
	}else{ // exchange
		model_two_comp* tmp_model = new model_two_comp();
		tmp_model -> copy(last_model);	
		double alpha_1 = tmp_model->add_snps(add, last_model);
		if (alpha_1 == FAIL_ALPHA){
			alpha = FAIL_ALPHA;
		}else{
			proposed_model -> copy(tmp_model);
			double alpha_2 = proposed_model->delete_snps(del, tmp_model);
			if (alpha_2 == FAIL_ALPHA){
				alpha = FAIL_ALPHA;
			}else{
				alpha += (alpha_1 + alpha_2);
			}
		}
		delete(tmp_model);
	}
//	cout << "Done" << endl;
		
	int acc = 0; 
	if (alpha > FAIL_ALPHA){
		double ah = proposed_model->sample_h2();
		if (g_sbf == 0){
			acc = proposed_model->compare_models(alpha + ah, last_model);
		}else{
			acc = proposed_model->compare_models_sbf(alpha + ah, last_model);
		}
	}
//	cout << g_mcmc_i << endl;
	if (acc == 1){	
		if (g_mcmc_i > 0 ){
			string paras = last_model->para_str();
			string model = last_model->model_str();
			PATH << g_mcmc_i - 1 << "\t" << g_mcmc_stay << "\t"  << paras << endl;
			MOD << model << endl;
		}
		proposed_model->sample_pi(); // pi can be sampled at last. It doesn't enter the calculation of MH
		proposed_model->sample_pve();
		if (g_mcmc_i > g_mcmc_warm) {
			int s = g_mcmc_i - g_mcmc_warm - 1; 
			if (s > g_mcmc_stay) {s = g_mcmc_stay;}
			last_model -> update_marginal(s);
		}
		delete(last_model);
		last_model = proposed_model;
		g_mcmc_stay = 1;
	}else{
		delete(proposed_model);
		g_mcmc_stay ++ ;
	}
	
	
	if (g_mcmc_i>g_mcmc_warm){
//		if (csize > 1 && prop_type < 2){prop_type += 3;}
		g_propose[prop_type] ++;
		if (acc == 1) {g_accept[prop_type] ++;}
	}	
	
	return;
}

int first_snp(int g){
	for (int i=0; i<g_n_snp; i++){
//		cout << i << "\t" << g_single_id[i] << endl;
		if (g_single_id[i] < g_size_g1){
			if (g == 0){return g_single_id[i];}
		}else{
			if (g == 1){return g_single_id[i];}
		}
	}
	return -1;
}

void mcmc_start(model_two_comp*& m){
	vector<int> s;
//	cout << "P1" << endl;
	if (g_start_size > 0){
		int sa = 0, sb = 0;
		for (int i=0; i<g_start_size; i++){
			s.push_back(g_single_id[i]);
			if (g_single_id[i] < g_size_g1){sa ++;}
			else {sb ++;}
		}
		if (sa == 0){s.push_back(first_snp(0));}
		if (sb == 0){s.push_back(first_snp(1));}
		m->initialize(s);
		LOG << "Starting model size = " << s.size() << endl;
	}else{
//		cout << "P2" << endl;
		s.push_back(first_snp(0));
		s.push_back(first_snp(1));
		m->initialize(s);
		double g1 = g_prop_uniform;
		double g2 = g_geometric; 
		g_prop_uniform = 0.001;
		g_geometric = 0.02;
		g_prop_add = 1;
		for (int i = 1; i<= 1000; i++){
	//		cout << i << endl;
			mcmc_sample(m);
			if (m->model_size()>=g_n_snp){break;}
		}
		g_prop_uniform = g1;
		g_geometric = g2;
		g_prop_add = 0.5;	
		LOG << "Starting model size = " << m->model_size() << endl;	
	}
	g_chol_call=0;
	return;
}

void mcmc_init(model_two_comp*& model){
	g_size_g1 = g_n_grp[0];
	g_size_g2 = g_n_snp - g_size_g1;
	
	if (g_expect_size_a < g_size_g1 && g_expect_size_a > 0){
		g_pi_alpha_2 = g_expect_size_a; 
		g_pi_beta_2 = g_size_g1 - g_expect_size_a;
	}
	
	if (g_max_model_size > (g_n_snp - 1) ) {g_max_model_size = g_n_snp - 1;}
	g_log_max_h = log(g_max_h);
	g_log_min_h = log(g_min_h);

	g_pip.clear(); g_pip.assign(g_n_snp, 0);
	g_rb_pip.clear(); g_rb_pip.assign(g_n_snp, 0.0);
	g_beta.clear(); g_beta.assign(g_n_snp, 0.0);
	g_rb_beta.clear(); g_rb_beta.assign(g_n_snp, 0.0);	
	g_propose.assign(5,0);
	g_accept.assign(5,0);

	if (g_continue == 0){
		mcmc_start(model);
	}else{
		vector<int> start_subset; double start_h2;
		read_last_bvsr(start_subset, start_h2, g_last_iter);
		model->initialize(start_subset, start_h2);	
	}
	return;
}

void mcmc_two_comp (void){
	cerr << "Initializing MCMC" << endl;
	model_two_comp* my_model = new model_two_comp();
	mcmc_init(my_model);
	int print_unit = calc_refresh();
	g_mcmc_stay = 0;
	cerr << "MCMC sampling" << endl;
	PATH << "No.\tLife\tSizeA\tSizeB\th2\trho\tPVE\tpi_A\tpi_B\tsigma_A\tsigma_B" << endl;
	MOD << "#Iteration index matched with the PATH file" << endl;
	if (g_output_yhat == 1){YHAT << "#Iteration index matched with the PATH file" << endl;}
	int start_iter = 0, total_iter = 0;
	calc_start_end(start_iter, total_iter);
	for (g_mcmc_i = start_iter; g_mcmc_i <= total_iter; g_mcmc_i ++){
		if (g_mcmc_i % print_unit == 0){cerr << '='; fflush(stdout); }
		mcmc_sample(my_model);
		if (g_mcmc_i > g_mcmc_warm){
			if (g_mcmc_i % g_rb_thin == 0){ my_model->rao_blackwell();}
		}
	}
	string paras = my_model->para_str();
	string model = my_model->model_str();
	PATH << g_mcmc_i - 1  << "\t" << g_mcmc_stay  << "\t" << paras << endl;	
	MOD << model << endl;
	my_model -> update_marginal(g_mcmc_stay);
	delete(my_model);	
	mcmc_output();	
	return;
}	



