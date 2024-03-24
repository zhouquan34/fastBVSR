#include "generic.h"

using namespace std;

int g_trunc_beta_iter = 10;

bool compare_pair(pair<int, double> x, pair<int, double> y){
	return x.second>y.second;
}

double log_sum_log (double a, double b){
	if (a < b){
		double c = a;
		a = b;
		b = c;
	}	// so a > b;
	return a + log(1.0 + exp(b-a));
}

double vector_mean (vector<double>& v){
	double s = 0.0;
	int n = v.size();
	for (int i=0; i<v.size(); i++){
		s += v.at(i);	
	}
	return s / (double) n;
}


void gsl_init (void){
	if (g_pset == 1){
		if (g_pmin > -1){
			g_min_pi = (double) g_pmin / (double) g_n_snp;	
			g_pi_range_set = 1;
		}
		if (g_pmax > -1){
			g_max_pi = (double) g_pmax / (double) g_n_snp;	
			g_pi_range_set = 1;
		}
		cerr << "Setting pi-max = " << g_max_pi << " and pi-min = " << g_min_pi << endl;
	}else{
	//	cerr << "pi-max = " << g_max_pi << " and pi-min = " << g_min_pi << endl;
	}

	
	gsl_rng_env_setup();
	gsl_type = gsl_rng_default;
	gsl_r = gsl_rng_alloc(gsl_type);
	gsl_rng_set(gsl_r, gsl_seed);
	return;
}

void gsl_exit (void){
	gsl_rng_free(gsl_r);
	return;
}


void get_row_col (const string file, int& nrow, int& ncol ){
	ifstream DAT(file.c_str() ); 
	string s;
	int row = 0;
	while (getline(DAT, s))	{
		replace(s.begin(), s.end(), ',', ' ');
		row ++ ;
		if (s.empty()){row --;}
		if (row == 2){
			stringstream ss(s);
			int c = 0; string tmp;
			while (ss >> tmp){
				c++;
			}
			ncol = c;
		}
	}
	nrow = row;
	DAT.close();
	return;
}

void print_array(const float* v, int p){
	for (int i = 0; i<p; i++){
		cout << v[i] << " ";
	}
	cout << endl;
}


void print_array(const double* v, int p){
	for (int i = 0; i<p; i++){
		cout << v[i] << " ";
	}
	cout << endl;
}

void print_matrix(const float* v, int p, int k){
	for (int i = 0; i < p; i += k){
		for (int j=0; j<k; j++){
			cout << v[i + j] << " ";
		}
		cout << endl;
	}
}

void print_matrix(const double* v, int p, int k){
	for (int i = 0; i < p; i += k){
		for (int j=0; j<k; j++){
			cout << v[i + j] << " ";
		}
		cout << endl;
	}
}


void write_log_cmds(int argc, char** argv){
	LOG << "# Command: "; 
	for (int i=0; i<argc; i++){ LOG << argv[i] << " ";} 
	LOG  << endl;
	if (g_continue == 1){
		cerr << "# Continue MCMC sampling of dataset " << g_last_bvsr << endl;
		cerr << "# Please note that only '-b', '-s', '-o', '-r' are enabled." << endl;
		cerr << "# " << g_out_prefix << ".beta.txt will be computed using all runs; all the other information is only for this run." << endl;	
		LOG << "\n# This is the continuation of run " << g_last_bvsr << endl;
		LOG << "# " << g_out_prefix << ".beta.txt will be computed using all runs; all the other information is only for this run." << endl;	
		LOG << "# Try 'fastBVSR-post n_burn_in path_file_1 path_file_2 ...' to recompute the posterior for multiple runs." << endl;
	}
	LOG << endl;
	LOG << "N_subject = " << g_n_sub  << "; N_predicator = " << g_n_snp << endl;
}

// currently we consider only g_min_pi = 0; 
double sample_beta(double alpha, double beta){
	if (g_pi_range_set == 0){
		return gsl_ran_beta(gsl_r, alpha, beta);
	}else{
		if (g_min_pi > 1e-12){
			return alpha/(alpha + beta);
		}
		double x = 0.0;
		double ymax = pow(1.0 - g_min_pi, beta - 1);
		double y = gsl_rng_uniform(gsl_r) *  ymax;
		double a = pow(g_min_pi, alpha);
		double ex1 = 1.0/(beta - 1.0);
		double ex2 = 1.0/(double) alpha;
//		cout << alpha << "\t" << beta << "\t" << g_min_pi << "\t" << g_max_pi << endl;
		for (int i=0; i<g_trunc_beta_iter; i++){
			double b = g_max_pi;
			double b1 = 1.0 - pow(y, ex1);
			if (b1 < b){ b = b1;}
			double r = gsl_rng_uniform(gsl_r);
			// if g_min_pi > 0; probably we need to check (a/b)^alpha == 0
			// double c = pow(b, alpha) - a;
			// x = pow(r * c + a, ex2);
			x = pow(r, ex2) * b;
			y = gsl_rng_uniform(gsl_r) * pow(1.0 - x, beta - 1.0);
//			cout << i << "\t" << x << endl;
		}
		return x;
	}
}

double calc_mean_scale (double sb){
	double d2 = g_dd_moments[0];
	double sb2 = sb * sb;
	double sb4 = sb2 * sb2;
	double sb8 = sb4 * sb4;
	double d4 = d2 * d2;
	double d8 = d4 * d4;
	double sb16 = sb8 * sb8;

	double lambda = d2 /(d2 + 1.0/sb2);
	double expl = exp(-0.5 * lambda);
	double f0 = expl * sqrt(1 + sb2 * d2);
	double f2 = 0.25 * expl * (2 * sb4 - d4 * sb8) /pow(d2 * sb2 + 1.0, 3.5);
	double f3 = 0.125 * expl * (-16 * sb4 * sb2 - 18 * d2 * sb8 + 3 * d4 * d2 * sb8 * sb4) / pow(d2 * sb2 + 1.0, 5.5 );
	double f4 = expl * ( 156 * sb8 + 320 * d2 * sb8 * sb2 + 180 * d4 * sb8 * sb4 - 15 * d8 * sb16) / pow(d2 * sb2 + 1.0, 7.5 ) / 16.0;

	return f0 + g_dd_moments[1] * f2 / 2.0 + g_dd_moments[2] * f3 / 6.0 + g_dd_moments[3] * f4 / 24.0;
}




