#include "output.h"

using namespace std;

void open_outputs (string prefix){
	string log_file = prefix;
	log_file.append(".log.txt");
	string out_file = prefix;
	out_file.append(".beta.txt");
	string path_file = prefix;
	path_file.append(".path.txt");
	string model_file = prefix;
	model_file.append(".model.txt");	
	
	LOG.open(log_file.c_str(), ofstream::out);
	OUT.open(out_file.c_str(), ofstream::out);
	PATH.open(path_file.c_str(), ofstream::out);
	MOD.open(model_file.c_str(), ofstream::out);
	LOG.precision(4);
	
	if (g_output_yhat == 1){
		string y_file = prefix;
		y_file.append(".yhat.txt");	
		YHAT.open(y_file.c_str(), ofstream::out);		
	}

	return;
}

void close_outputs(void){
	LOG.close();
	OUT.close();
	PATH.close();
	MOD.close();

	if (g_output_yhat == 1){ YHAT.close(); }

	return;
}

void post_quantile(map<double, int>& v){
	map<double, int>::iterator it;
	int n = 0;
	double ss = 0.0;
	for (it=v.begin(); it!=v.end(); it++){
		n += it->second;
		ss += (double) it->second * it->first;
	}
	if (n > 0){
		double m = ss / n;
		vector<double> qs(5,0);
		double quant[] = {0.025, 0.05, 0.5, 0.95, 0.975};
		double k = 0.0;
		int qi = 0;
		for (it=v.begin(); it!=v.end(); it++){
			k += it->second;
			while (k/n > quant[qi] - ZERO){ 
				qs[qi] = it->first; 
				qi ++;
				if (qi == 5) {break;}
			}
			if (qi == 5) {break;}
		}
		LOG << "Mean = " << m << ";  Median = " << qs[2] << endl;
		LOG << "90% credible interval = (" << qs[1] <<  ", " << qs[3] << ")" << endl;
		LOG << "95% credible interval = (" << qs[0] <<  ", " << qs[4] << ")" << endl;	
	}
	return;
}


void post_quantile(map<int, int>& v){
	map<int, int>::iterator it;
	int n = 0;
	double ss = 0.0;
	for (it=v.begin(); it!=v.end(); it++){
		n += it->second;
		ss += (double) it->second * it->first;
	}
	if (n > 0){
		double m = ss / n;
		vector<double> qs(5,0);
		double quant[] = {0.025, 0.05, 0.5, 0.95, 0.975};
		double k = 0.0;
		int qi = 0;
		for (it=v.begin(); it!=v.end(); it++){
			k += it->second;
			while (k/n > quant[qi]){ 
				qs[qi] = it->first; 
				qi ++;
				if (qi == 5) {break;}
			}
			if (qi == 5) {break;}
		}
		LOG << "Mean = " << m << ";  Median = " << qs[2] << endl;
		LOG << "90% credible interval = (" << qs[1] <<  ", " << qs[3] << ")" << endl;
		LOG << "95% credible interval = (" << qs[0] <<  ", " << qs[4] << ")" << endl;	
	}
	return;
}


void mcmc_output(void){
	record_time();
	cerr << "\nOutput results" << endl;
	cerr << "Try 'fastBVSR-trace " << g_out_prefix << "' to draw the traceplots" << endl;
	if (g_continue == 1){
		cerr << "Try 'fastBVSR-post n_burn_in path_file_1 path_file_2 ...' to recompute the posterior for multiple runs" << endl;
	}

	vector<double> acr (5, 0);
	for (int i=0; i<5; i++){
		if (g_propose[i] == 0){acr[i] = 0;}
		else{
			acr[i] =( (double) g_accept[i])/g_propose[i]; 
		}
	}
	
	LOG << "\nModel size:" << endl;
	post_quantile(g_model_size);

	LOG << "\nHeritability:" << endl;
	post_quantile(g_h2);
	
	LOG << "\nPVE:" << endl;
	post_quantile(g_pve);
	LOG << endl;


	LOG << "Acceptance rates:" << endl;
	LOG << "Add = " << acr[0] << ";  Delete = " << acr[1] << "; Long-range = " << acr[2] << endl;
//	LOG << "Long-range = " << acr[2] << ", " << acr[3] << ", " << acr[4] << endl;
	LOG << endl;

	// output icf summary
	map<int, int>::iterator it;
	int icf_n = 0;
	for (it=g_icf_iter.begin(); it!=g_icf_iter.end(); it++){ icf_n += it->second; }
	int icf_n0 = icf_n + g_chol_call;
	if (icf_n0 == 0){LOG << "ICF was not called" << endl;}
	else{
		double fail_rate = (double) g_chol_call / (double) icf_n0;
		LOG << "ICF calls = " << icf_n0 << "; Cholesky calling rate = " << fail_rate << endl;
		LOG << "ICF iterations:" << endl;
		post_quantile(g_icf_iter);
	}
	LOG << endl;

	OUT << "SNP\tlog10(BF)\tPIP\tRB-PIP\tBeta\tRB-Beta\n";
	if (g_rb_count == 0) {g_rb_count = 1;}
	int total_iter = g_mcmc_iter  + g_last_iter; 
	int total_rb = g_rb_count + g_last_iter / g_rb_thin; 
	cout << "total_iter = " << total_iter << endl;
	cout << "total_rb = " << total_rb << endl;
	for (int i=0; i<g_snp_names.size(); i++){
		int smap = g_snp_map[i];
		if (smap >= 0){
			double pip  = (double) g_pip[smap] / total_iter;
			double beta = g_beta[smap] / total_iter;
			double rb_pip = g_rb_pip[smap] / total_rb;
			double rb_beta = g_rb_beta[smap] / total_rb;
			string snp = g_snp_names[i];
			if (g_dup_keep.find(i) != g_dup_keep.end()){
				double nc = (double) g_dup_count[i];
				OUT << snp << "\t" << g_single_bf[smap]/log(10.0) << "\t" << pip/nc << "\t" << rb_pip/nc << "\t" << beta/nc << "\t" << rb_beta/nc << endl;
			}else{
				OUT << snp << "\t" << g_single_bf[smap]/log(10.0) << "\t" << pip << "\t" << rb_pip << "\t" << beta << "\t" << rb_beta << endl;
			}
		}else{
			int k = g_snp_map[g_dup_map[i]];
			double nc = (double) g_dup_count[i];
			double pip  = (double) g_pip[k] / total_iter;
			double beta = g_beta[k] / total_iter;
			double rb_pip = g_rb_pip[k] / total_rb;
			double rb_beta = g_rb_beta[k] / total_rb;
			OUT <<  g_snp_names[i] << "\t" << g_dup_bf[i]/log(10.0) << "\t" << pip/nc << "\t" << rb_pip/nc << "\t" << beta/nc << "\t" << rb_beta/nc << endl;
		}
	}
	record_time("Output");
	return;
}





