#include "ftime.h"

using namespace std;

void write_time(double used_t){
	LOG << "Time used by MCMC = " ;
	if (used_t > 2000){LOG << used_t/60.0 << " m" << endl;}
	else{LOG << used_t << " s" << endl;}
	return;
}


void record_time (string x){
	if (x == ""){ 
		gettimeofday(&g_t_beg, NULL);
		return;
	}else{
		gettimeofday(&g_t_end, NULL);
		double dt = (double) (g_t_end.tv_usec - g_t_beg.tv_usec) / 1000.0 + (g_t_end.tv_sec - g_t_beg.tv_sec) * 1000.0;	
		if (g_time.find(x) != g_time.end()){
			g_time[x] += dt / 1000.0 ;
		}else{
			g_time[x] = dt/1000.0;	
		}
		return;
	}
}


void output_time(void){
	LOG << "\nTime usage details (in seconds):" << endl;
	map<double, string> tmp;
	map<string, double>::iterator it;
	for (it=g_time.begin(); it!=g_time.end(); it++){
		tmp[it->second] = it->first;
	}
	map<double, string>::iterator it2;
	for (it2=tmp.begin(); it2!=tmp.end(); it2++){
		LOG << it2->second << "\t" << it2->first << endl;
	}
	return;
}





