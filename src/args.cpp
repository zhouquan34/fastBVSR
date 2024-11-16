#include "args.h"

using namespace std;

void print_help(void){
	cout << "Example command: " << endl;
	cout << "\t\t./fastBVSR -g test.mg.txt -p test.ph -o bvsr -s 10000" << endl;
	cout << endl;
	cout << "Arguments:" << endl;
	cout << "-g: mean genotype file (same format as piMASS input)" << endl;
	cout << "-a: PLINK 1.9 A-transpose file" << endl;
	cout << "-m: design matrix file" << endl;
	cout << "-p: phenotype file (same format as piMASS input)" << endl;
	cout << "-o: output prefix" << endl;
	cout << "-w: number of burn-in iterations" << endl;
	cout << "-s: number of MCMC iterations of this run" << endl;
	cout << "-R: do Rao-Blackwellization every R iterations" << endl; 
	cout << "-r: random seed" << endl;
	exit(EXIT_FAILURE);
	return;
}


void print_manual(void){
	cout << "All arguments:" << endl;
	cout << "--mg:        input mean genotype file" << endl;
	cout << "--plink:     input plink A-transpose genotype file" << endl;
	cout << "--mat:       input matrix genotype file" << endl;
	cout << "--pheno:     input phenotype file" << endl;
	cout << "--out:       output prefix" << endl;
	cout << "--cov:       covariate file" << endl;
	cout << "--last:      output prefix of the dataset to continue" << endl;
	cout << "--mcmc:      MCMC iteration" << endl;
	cout << "--burn:      burn-in iteration" << endl;
	cout << "--rb:        Rao-Blackwell frequency" << endl;
	cout << "--gamma:        set gamma used in fractional likelihood" << endl;
	cout << "--start:     (-nstart) starting model size" << endl;
	cout << "--lunif-h2:  use uniform prior on log-h2" << endl;
	cout << "--time:      output time usage details" << endl;
	cout << "--long:      long-range proposal proportion" << endl;
	cout << "--long-ex:   exchange long-range proportion" << endl;
	cout << "--geom:      geometric rate of add proposal" << endl;
	cout << "--jump:      max jump size" << endl;
	cout << "--h2-walk:   max h2 walk distance" << endl;
	cout << "--h2-unif:   proportion of random walk on h2" << endl;
	cout << "--add-unif:  proportion of uniform add" << endl;
	cout << "--icf-err:   ICF max error" << endl;
	cout << "--cor:       correlation cutoff (may cause problems with -b)" << endl;
	cout << "--pi-alpha:  prior alpha of pi" << endl;
	cout << "--pi-beta:   prior beta of pi" << endl;
	cout << "--pi-min:    min pi" << endl;
	cout << "--pi-max:    max pi" << endl;
	cout << "--h2-min:    (-hmin) min h2" << endl;
	cout << "--h2-max:    (-hmax) max h2" << endl;
	cout << "--filter:    remove duplicate SNPs" << endl; // automatically true for matrix data
	cout << "--prec:      dim of XtX to be precalculated" << endl;
	cout << "--max-size:  (-smax) max model size allowed" << endl;
	cout << "--min-size:  (-smin) min model size allowed" << endl;
	cout << "--exact:     exact calculation of BF" << endl;
	cout << "--seed:      GSL seed" << endl;
	cout << "--help:      print help" << endl;
	cout << "--HELP:      print full help" << endl;
	cout << "--n-snp:     the original total number of SNPs" << endl;
	cout << "--yhat:       output fitted values" << endl;
	exit(EXIT_FAILURE);
	return;
}


void init_paras(void){
	g_paras["-g"] = PARA_MG;
	g_paras["--mg"] = PARA_MG;
	g_paras["-a"] = PARA_PLINK;
	g_paras["--plink"] = PARA_PLINK;
	g_paras["-m"] = PARA_MAT;
	g_paras["--mat"] = PARA_MAT;
	g_paras["-p"] = PARA_PHENO;
	g_paras["--pheno"] = PARA_PHENO;
	g_paras["-o"] = PARA_OUT;
	g_paras["--out"] = PARA_OUT;
	g_paras["-c"] = PARA_COV;
	g_paras["--cov"] = PARA_COV;
	g_paras["-s"] = PARA_ITER;
	g_paras["--mcmc"] = PARA_ITER;
	g_paras["-w"] = PARA_BURN;
	g_paras["--burn"] = PARA_BURN;
	g_paras["-R"] = PARA_RB;
	g_paras["--rb"] = PARA_RB;
	g_paras["-k"] = PARA_START;
	g_paras["--start"] = PARA_START;
	g_paras["-t"] = PARA_TIME;
	g_paras["--time"] = PARA_TIME;
	g_paras["-b"] = PARA_BVSR;
	g_paras["--last"] = PARA_BVSR;
	g_paras["-z"] = PARA_SEED;
	g_paras["-r"] = PARA_SEED;
	g_paras["--seed"] = PARA_SEED;
	g_paras["-h"] = PARA_MAN;
	g_paras["--help"] = PARA_MAN;	
	g_paras["--long"] = PARA_LONG;
	g_paras["--long-ex"] = PARA_LEX;
	g_paras["--geom"] = PARA_GEOM;
	g_paras["--icf-err"] = PARA_ICFE;
	g_paras["--cor"] = PARA_LD;
	g_paras["--pi-alpha"] = PARA_PI_A;
	g_paras["--pi-beta"] = PARA_PI_B;
	g_paras["--pi-min"] = PARA_MINPI; 
	g_paras["--pi-max"] = PARA_MAXPI; 
	g_paras["--h2-min"] = PARA_MINH2;
	g_paras["--h2-max"] = PARA_MAXH2;
	g_paras["--gamma"] = PARA_GAMMA;
	g_paras["--jump"] = PARA_MAXJ;
	g_paras["--lunif-h2"] = PARA_LOGH;
	g_paras["--h2-unif"] = PARA_UH2;
	g_paras["--add-unif"] = PARA_UADD;
	g_paras["--h2-walk"] = PARA_HJUMP;
	g_paras["--filter"] = PARA_NOF;
	g_paras["--prec"] = PARA_PREC;
	g_paras["--max-size"] = PARA_MAXS;
	g_paras["--min-size"] = PARA_MINS;
	g_paras["--n-snp"] = PARA_NSNP;
	g_paras["--n-try"] = PARA_NTRY;
	g_paras["--exact"] = PARA_EXACT;
	g_paras["--HELP"] = PARA_MAN;
	g_paras["--yhat"] = PARA_YHAT;
	// PIMASS parameters //
	g_paras["-num"] = PARA_RB1;
	g_paras["-hmin"] = PARA_MINH2;
	g_paras["-hmax"] = PARA_MAXH2;
	g_paras["-smax"] = PARA_MAXS;
	g_paras["-smin"] = PARA_MINS;
	g_paras["-pmax"] = PARA_PMAX;
	g_paras["-pmin"] = PARA_PMIN;
	g_paras["-nstart"] = PARA_START;
	
	// PIMASS parameteres not defined //
	g_paras["-pos"] = PARA_NOT;
	g_paras["-cc"] = PARA_NOT;
	g_paras["-v"] = PARA_NOT2;
	g_paras["-exclude-maf"] = PARA_NOT2;
	g_paras["-exclude-nopos"] = PARA_NOT2;
	g_paras["-silence"] = PARA_NOT2;
	return;
}

void print_continue_help(void){
	cerr << "The continuation mode only accepts the following arguments: \n" << endl;
	cerr << "-b: prefix of the MCMC to continue" << endl;
	cerr << "-o: output prefix" << endl;
	cerr << "-s: number of MCMC iterations of this run" << endl;
	cerr << "-z: random seed" << endl;
	cerr << "\nAll the other arguments are forced to take their values in the first run" << endl;
	exit(EXIT_FAILURE);
}

int set_args_part_one (string para, char* val){
	if (g_paras.find(para) == g_paras.end()){ return 0;}
	switch(g_paras[para]){
		case PARA_BVSR:
			g_last_bvsr = val; return 1;
		case PARA_SEED:
			cerr << "WARN: -R is for Rao-Blackwellization and -r is for random seed" << endl;
			gsl_seed = atoi(val); return 1;
		case PARA_OUT:
			g_out_prefix = val; return 1;
		case PARA_ITER:
			g_mcmc_iter = atoi(val); return 1;	
		default:
			return 0;
	}
	return 0;
}

int set_args_part_two (string para, char* val){
	if (g_paras.find(para) == g_paras.end()){ return 0;}
	switch(g_paras[para]){
		case PARA_MG:
			read_meang(val); return 1;
		case PARA_PLINK:
			read_plink(val); return 1;
		case PARA_MAT:
			read_matrix(val); return 1;
		case PARA_PHENO:
			read_pheno(val); return 1;
		case PARA_COV:
			read_cov(val); return 1;
		case PARA_ITER:
			g_mcmc_iter = atoi(val); return 1;  // may be g_last_iter
		case PARA_BURN:
			g_mcmc_warm = atoi(val); return 1;
		case PARA_RB:
			g_rb_thin = atoi(val); return 1;
		case PARA_RB1:
			g_rb_thin = atoi(val) * 10; return 1;
		case PARA_START:
			g_start_size = atoi(val); return 1; 
		case PARA_TIME:
			g_output_time = 1; return 1;
		case PARA_LOGH: // turn on
			g_log_uniform_h = 1; return 1;
		case PARA_LONG:	
			g_long_range = atof(val); return 1;
		case PARA_LEX:
			g_long_exchange = atof(val); return 1;
		case PARA_ICFE:
			g_icf_abs_tol = atof(val); return 1; 
		case PARA_LD:
			g_corr_cut = atof(val); return 1;
		case PARA_GEOM:
			g_geometric = atof(val); return 1;
		case PARA_MINH2:
			g_min_h = atof(val); return 1;
		case PARA_MAXH2:
			g_max_h = atof(val); return 1;
		case PARA_MINPI:
			g_min_pi = atof(val); g_pi_range_set = 1; return 1;
		case PARA_MAXPI:
			g_max_pi = atof(val); g_pi_range_set = 1; return 1;
		case PARA_PMAX:
			g_pmax = atof(val); g_pset = 1; return 1;
		case PARA_PMIN:
			g_pmin = atof(val); g_pset = 1; return 1;		
		case PARA_GAMMA:
			g_gamma = atof(val); return 1;
		case PARA_MAXJ:
			g_max_jump = atoi(val); return 1;
		case PARA_HJUMP:
			g_jump_h2 = atof(val); return 1;
		case PARA_UADD:
			g_prop_uniform = atof(val); return 1;
		case PARA_UH2:
			g_prop_h2_uniform = atof(val); return 1;
		case PARA_PI_A:
			g_pi_prior_alpha = atoi(val); return 1;
		case PARA_PI_B:
			g_pi_prior_beta = atoi(val); return 1;
		case PARA_PREC:
			g_precomp = atoi(val); return 1;
		case PARA_MAXS:
			g_max_model_size = atoi(val); return 1;
		case PARA_MINS:
			g_min_model_size = atoi(val); return 1;
		case PARA_NOF:
			g_no_filter = 0; return 1;  // originally it was to set g_no_filter = 1
		case PARA_YHAT:
			g_output_yhat = 1; return 1; 
		case PARA_NSNP:
			g_pseudo_n = atoi(val); return 1;
		case PARA_NTRY:
			g_mcmc_init_try = atoi(val); return 1;
		case PARA_EXACT:
			g_exact_bf = 1; return 1;
		case PARA_HELP:
			print_help(); return 1;
		case PARA_MAN:
			print_manual(); return 1;	
		case PARA_NOT:
			cerr << "NOTE: -pos, -cc are not enabled in fastBVSR" << endl; return 1;
		case PARA_NOT2:
			cerr << "NOTE: -v, -silence, -exclude-maf, -exclude-nopos are not enabled in fastBVSR"<< endl; return 1;
		default:
			return 0;
	}
	return 0;
}


int scan_continue(int argc, char** argv){
	for (int i=1; i<argc; i++){
		if (argv[i][0] == '-'){
			string para = argv[i];
			if (g_paras.find(para) != g_paras.end()){
				if (g_paras[para] == PARA_BVSR){
					return 1;
				}	
			}
		}
	}
	return 0;
}

void read_continue_args (int argc, char** argv){
	for (int i=1; i<argc; i++){
		if (argv[i][0] == '-'){
			string para = argv[i];
			char* val = NULL;
			if (i < argc - 1) { val = argv[i+1]; }
			if (set_args_part_one(para, val) == 0){
				print_help();
			}
		}
	}
	return;
}

void read_last_args (void){
	string log_file = g_last_bvsr;
	log_file.append(".log.txt");		
	ifstream flog(log_file.c_str());
	if (flog.fail()){
		cerr << "Log file of last run does not exist!" << endl;
		exit(EXIT_FAILURE);
	}

	string s;
	int mcmc_iter_copy = g_mcmc_iter; // already set
//	int gsl_seed_copy  = gsl_seed;
	g_mcmc_iter = 1000; // default value
	while (getline(flog, s)){
		stringstream ss(s); 
		vector<string> args;	
		char pound; ss >> pound;  
		if (pound != '#'){break;} // scan all the previous commands
		if (s.empty()){break;}
		string tmp; ss >> tmp;
		while (ss >> tmp){ args.push_back(tmp); }
		LOG << "# Previous: "; 
		for (int i=0; i<args.size(); i++){ LOG << args[i] << " ";} 
		LOG << endl;
		
		for (int i=1; i<args.size(); i++){
			if (args[i].at(0) != '-'){ continue; }
			char* val = NULL;
			if (i < args.size() - 1) { val = const_cast<char*> (args[i+1].c_str()); }
			set_args_part_two(args[i], val);
		}
		g_last_iter += g_mcmc_iter; 
	}
	flog.close();
	g_mcmc_iter = mcmc_iter_copy; // recover
//	gsl_seed = gsl_seed_copy;
	return;
}


void read_args (int argc, char** argv){
	if (argc <= 1){ print_help(); }
	init_paras();
	if (scan_continue(argc, argv) == 1){
		g_continue = 1;
		read_continue_args(argc, argv);
		read_last_args();
		return;
	}

	for (int i=1; i<argc; i++){
		if (argv[i][0] == '-'){
			string para = argv[i];
			char* val = NULL;
			if (i < argc - 1) { val = argv[i+1]; }
			int defined = 0; 
			defined += set_args_part_one(para, val);
			defined += set_args_part_two(para, val);
			if (defined == 0){		
				print_help();
			}
		}
	}
	
	return;
}



