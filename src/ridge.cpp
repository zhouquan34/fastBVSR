#include "ridge.h"

using namespace std;


double calc_ridge(double*& RA, double* R, double* z, int p, double* diags, double* b, double* XX){
	if (RA != NULL){
		Chol_solve(p, RA, z, b);		
		return 0;
	}else if (p <= g_icf_min_p){
		RA = allocate1D(p * p);
		double log_det = Chol_decomp_A(p, XX, RA, diags);	
		Chol_solve(p, RA, z, b);
		return log_det;
	}else{	
		record_time();
		int iter = ICF(R, z, p, g_n_cov, diags, b);
		record_time("ICF");

		if (iter == 0){
			g_chol_call ++ ;
			RA = allocate1D(p * p);
			Chol_decomp_A(p, XX, RA, diags);	
			Chol_solve(p, RA, z, b);		
		}else{
			g_icf_iter[iter] ++ ;
		}
		return 0;
	}
}

double calc_ridge(double*& RA, double* R, double* z, int p, double diag, double* b, double* XX){
	if (RA != NULL){
		Chol_solve(p, RA, z, b);		
		return 0;
	}else if (p <= g_icf_min_p){
		RA = allocate1D(p * p);
		double log_det = Chol_decomp_A(p, XX, RA, diag);	
		Chol_solve(p, RA, z, b);	
		return log_det;
	}else{	
		record_time();
		int iter = ICF(R, z, p, g_n_cov, diag, b);
		record_time("ICF");

		if (iter == 0){
			g_chol_call ++ ;
			RA = allocate1D(p * p);
			Chol_decomp_A(p, XX, RA, diag);	
			Chol_solve(p, RA, z, b);		
		}else{
			g_icf_iter[iter] ++ ;
		}
		return 0;
	}
//	print_array(RA, p * p);
//	print_array(z, p);
//	print_array(b, p);
}



void calc_globals (void){
	
	record_time();
	if (g_n_cov > 0){	
		Chol_decomp_g_cov();	
		
		double* wy = allocate1D(g_n_cov);
		double* cb = allocate1D(g_n_cov);
		mat_vec_mul(g_cov, g_pheno, g_n_cov, g_n_sub, wy, 0);
		Chol_solve(g_n_cov, g_cov_R, wy, cb);
		double* hwy = allocate1D(g_n_sub);
		mat_vec_mul(g_cov, cb, g_n_cov, g_n_sub, hwy, 1);
		for (int i=0; i<g_n_sub; i++){
			g_pheno[i] -= hwy[i];
		}
		free1D(hwy);
		free1D(wy);
		free1D(cb);
	
	}
	
	g_yy = vec_vec_mul(g_pheno, g_pheno, g_n_sub);
	g_xty = allocate1D(g_n_snp);
	mat_vec_mul(g_data, g_pheno, g_n_snp, g_n_sub, g_xty, 0);		
	
	record_time("CalcGlobal");
	return;
}

void calc_single_bf(const double sb, vector< pair<int,double> >& bfs){
	record_time();		
	for (int i=0; i<g_n_snp; i++){
		double xy = g_xty[i];
		double xx = g_n_sub * g_snp_var.at(i);
		double yhy = xy * xy / (xx + 1.0/sb/sb);	
		double log_bf = -0.5 * g_n_sub * log(1 - yhy/g_yy) - log(sb) - 0.5 * log(xx + 1/sb/sb); 
		g_single_bf.push_back(log_bf);
		g_single_h2.push_back(yhy);
		pair<int, double> sing = make_pair(i,  log_bf);
		bfs.push_back(sing);
	}
	sort(bfs.begin(), bfs.end(), compare_pair);
	record_time("SingleBF");
	return;
}

void sort_single_bf_2(void){
	record_time();
	vector< pair<int, double> > single_prob; 
	double shift = ( (double) g_n_grp[1] - g_n_grp[0] ) / (double) g_n_grp[0];
	shift = log(shift);
	for (int i=0; i<g_n_grp[0]; i++){
		pair<int, double> sing = make_pair(i, g_single_bf[i] + shift);
		single_prob.push_back(sing);
	}	
	for (int i=g_n_grp[0]; i<g_n_snp; i++){
		pair<int, double> sing = make_pair(i, g_single_bf[i]);
		single_prob.push_back(sing);
	}	

	sort(single_prob.begin(), single_prob.end(), compare_pair);
	g_single_id.clear();
	g_single_ord.clear();
	g_single_ord.assign(g_n_snp, 0);
	for (int i=0; i<g_n_snp; i++){
		g_single_id.push_back(single_prob[i].first);	
		g_single_ord[single_prob[i].first] = i;
	}
	record_time("SingleBF");
	return;
}

void sort_single_bf(void){
	if (g_group == 2){
		sort_single_bf_2();
		return;
	}

	record_time();
	vector< pair<int, double> > single_prob; 
	for (int i=0; i<g_n_snp; i++){
		pair<int, double> sing = make_pair(i, g_single_bf[i]);
		single_prob.push_back(sing);
	}	
	sort(single_prob.begin(), single_prob.end(), compare_pair);
	g_single_id.clear();
	g_single_ord.clear();
	g_single_ord.assign(g_n_snp, 0);
	for (int i=0; i<g_n_snp; i++){
		g_single_id.push_back(single_prob[i].first);	
		g_single_ord[single_prob[i].first] = i;
	}
	record_time("SingleBF");
	return;
}

int if_same_snps (int a, int b){
	float* ai = &g_data[a * g_n_sub];
	float* bi = &g_data[b * g_n_sub];
	int sign = 1;
	if (abs(ai[0] + bi[0]) < ZERO ){sign = -1;}
	
	for (int i=0; i<g_n_sub; i++){
		if (abs(ai[i] - bi[i] * sign) > ZERO){
			return 0;
		}
	}
	return 1;
}

void filter_snp(const double sb){
	vector<int> identical(g_n_snp, 0);
	vector< pair<int, double> > bfs;
	calc_single_bf(sb, bfs);
/*
	double tmp = 0.0;
	for (int i=0; i<g_n_sub; i++){	
		tmp = gsl_rng_uniform(gsl_r); 
		tmp = gsl_rng_uniform(gsl_r); 
	}
*/
	
	g_snp_map.assign(g_n_snp, -1);	
	if (g_no_filter == 0) {
		cerr << "Removing identical SNPs" << endl;	
		record_time();
		vector<int> to_remove;
		set<int> rsnps;
		int nc = 0;
		for (int i=0; i<g_n_snp; i++){
			int a = bfs[i].first;
			if (rsnps.find(a) != rsnps.end()){ continue;}
			for (int j=i+1; j<g_n_snp; j++){
				if (bfs[i].second - bfs[j].second > ZERO){ break;}
				else{
					int b = bfs[j].first;
					if (abs (g_snp_var[a] - g_snp_var[b]) > ZERO ){ continue; }
					if (rsnps.find(b) == rsnps.end()){
						nc ++;
						if (if_same_snps(a,b) == 1 ){
							g_dup_map[b] = a; 
							rsnps.insert(b);
							to_remove.push_back(b);
						}
					}
				}	
			}	
		}
//		cerr << "Compared " << nc << " pairs" << endl;

		if (to_remove.size()>0){
			string file = g_out_prefix;
			file.append(".remove.txt");
			ofstream REM(file.c_str()); 
			REM << "Keep\tRemove" << endl;	
			LOG << "Remove " << to_remove.size() << " identical SNPs" << endl;
			
		//	map< string, vector<string> > group;
			map< int, vector<int> > group;
			for (map<int, int>::iterator it=g_dup_map.begin(); it!=g_dup_map.end(); it++ ){
				g_dup_bf[it->first] = g_single_bf[it->first];
			//	group[g_snp_names[it->second]].push_back(g_snp_names[it->first]);	
				group[it->second].push_back(it->first);	
			}
			
		//	map< string, vector<string> >::iterator ssit;
	//		vector<int> dup_keep;
			map< int, vector<int> >::iterator ssit;
			for (ssit = group.begin(); ssit != group.end(); ssit ++){
				g_dup_keep.insert(ssit->first);
				REM << g_snp_names[ssit->first] << "\t";
				int nr = ssit->second.size();
				g_dup_count[ssit->first] = nr + 1; 
				for (int j=0; j < nr-1; j++){
					g_dup_count[ssit->second.at(j)] = nr + 1; 
					REM << g_snp_names[ssit->second.at(j)] << ", ";
				}
				g_dup_count[ssit->second.at(nr - 1)] = nr + 1; 
				REM << g_snp_names[ssit->second.at(nr-1)] << endl;
			}		
			REM.close();

			// recompute g_data, g_data_t_pheno
			sort(to_remove.begin(), to_remove.end());
			int new_n_snp = g_n_snp - to_remove.size();
			float* new_data = allocate1D_float(new_n_snp * g_n_sub);	
			double* new_xty = allocate1D(new_n_snp);
			matrix_del(to_remove, g_n_snp, g_n_sub, g_data, new_data);
			array_del(to_remove, g_n_snp, g_xty, new_xty); // no need to recompute
			free1D(g_data);
			free1D(g_xty);
			g_data = new_data;
			g_xty = new_xty;
		
			// recompute g_n_snp, g_snp_var, g_single_bf; g_snp_mean won't be used
			g_n_snp = new_n_snp;
			for (int i=(to_remove.size()-1); i>=0; i--){
				identical[to_remove[i]] = 1;
				g_snp_var.erase(g_snp_var.begin() + to_remove[i]);
				g_snp_mean.erase(g_snp_mean.begin() + to_remove[i]);
				g_single_bf.erase(g_single_bf.begin() + to_remove[i]);
				g_single_h2.erase(g_single_h2.begin() + to_remove[i]);
			}
		}
		record_time("FilterSNP");
	}

	int si = 0;
	for (int i=0; i<identical.size(); i++){
		if (identical[i] == 0){
			g_snp_map[i] = si;
			si ++;
		}
	}	
	sort_single_bf();

	return;
}

// i(row) >= j(col)
void calc_xtx (void){
	record_time();
	cerr << "Precalculation of part of the gram matrix" << endl;
	if (g_precomp > g_n_snp){g_precomp = g_n_snp;}
	int p = g_precomp + g_n_cov;
	g_xtx = allocate1D( p * p , 0 );
	
	for (int i = 0; i < g_n_cov; i++){
		int r = i * p;
		for (int j = 0; j <= i; j++){
			g_xtx[r + j] = vec_vec_mul(&g_cov[i * g_n_sub], &g_cov[j * g_n_sub], g_n_sub);
		}
	}

	if (g_n_cov > 0){
		for (int i = g_n_cov; i < p; i++){
			int r = i * p;
			for (int j = 0; j < g_n_cov; j++){
				g_xtx[r + j] = vec_vec_mul(&g_data[g_single_id[i - g_n_cov] * g_n_sub], &g_cov[j * g_n_sub], g_n_sub);
			}
		}
	}

	for (int i = g_n_cov; i < p; i++){
		int r = i * p;
		for (int j = g_n_cov; j < i; j++){
			g_xtx[r + j] = vec_vec_mul(&g_data[g_single_id[i - g_n_cov] * g_n_sub], &g_data[g_single_id[j - g_n_cov] * g_n_sub], g_n_sub);
		}
	}

	for (int i=0; i<p; i++){
		for (int j=i+1; j<p; j++){
			g_xtx[i * p + j] = g_xtx[j * p + i];
		}
	}	
	record_time("PrecalcXX");
//	cout.precision(15);
//	print_matrix(g_cov, g_n_sub * 2, g_n_sub);
//	print_matrix(g_data, g_n_sub * g_n_snp, g_n_sub);
//	print_matrix(g_xtx, p*p, p);
}





