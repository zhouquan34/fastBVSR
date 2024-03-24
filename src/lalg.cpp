#include "lalg.h"

using namespace std;

float* allocate1D_float (int dim){
	float* m = (float*) malloc((size_t) (dim * sizeof(float)));
	return m;
}

float** allocate1D_float_pointer (int dim){
	float** m = (float**) malloc((size_t) (dim * sizeof(float*)));
	return m;
}

double* allocate1D (int dim){
	double* m = (double*) malloc((size_t) (dim * sizeof(double)));
//	memset(m, 0, dim * sizeof(double)); // NOT necessary. Forbid for max speed. 
	return m;
}

double* allocate1D (int dim, int v){
	double* m = (double*) malloc((size_t) (dim * sizeof(double)));
	memset(m, v, dim * sizeof(double));
	return m;
}

double** allocate2D (int dim){
	double** m;
	m = (double**) malloc((size_t) (2 * sizeof(double*))); // every row is double*
	m[0] = (double*) malloc((size_t) (2 * dim * sizeof(double)));
	memset(m[0], 0, 2 * dim * sizeof(double));
	m[1] = m[0] + dim;
	return m;
}

void free1D(float** m){
	if (m == NULL){return;}
	free(m);
	m = NULL;
	return;
}

void free1D(float* m){
	if (m == NULL){return;}
	free(m);
	m = NULL;
	return;
}

void free1D(double* m){
	if (m == NULL){return;}
	free(m);
	m = NULL;
	return;
}

void free2D (double** m){
	if (m == NULL) {return;}
	free(m[0]);
	free(m);
	m = NULL;
	return;
}

void copy1D (float* m1, const float* m2, int p){
	for (int i=0; i<p; i++){
		m1[i] = m2[i];
	}
	return;
}

void copy1D (float** m1, float** m2, int p){
	for (int i=0; i<p; i++){
		m1[i] = m2[i];
	}
	return;
}

void copy1D (double* m1, const float* m2, int p){
	for (int i=0; i<p; i++){
		m1[i] = m2[i];
	}
	return;
}

void copy1D (double* m1, const double* m2, int p){
	for (int i=0; i<p; i++){
		m1[i] = m2[i];
	}
	return;
}

// no need to set z zero
void mat_times_vec_opt (double* m, double* v, int p, int n, double* z){	
	double* zpos = &z[0];
	double* mp1 = &m[0]; double *mp2, *mp3, *mp4;
	if (p>3){
		mp2 = &m[n];
		mp3 = &m[2*n];
		mp4 = &m[3*n];
	}

	int k = p/4;
	for (int i=0; i<k; i++){
		double z1 = 0, z2 = 0, z3 = 0, z4 = 0;
		double* vpos = &v[0];
		for (int j=0; j<n; j++){
			z1 += (*mp1) * (*vpos);
			z2 += (*mp2) * (*vpos);
			z3 += (*mp3) * (*vpos);
			z4 += (*mp4) * (*vpos);
			vpos ++;
			mp1 ++;
			mp2 ++;
			mp3 ++;
			mp4 ++;
		}
		*zpos = z1;  zpos ++;
		*zpos = z2;  zpos ++;
		*zpos = z3;  zpos ++;
		*zpos = z4;  zpos ++;
		
		mp1 += 3 * n;
		mp2 += 3 * n;
		mp3 += 3 * n;
		mp4 += 3 * n;
	}
	
	for (int i=k*4; i<p; i++){
		double z0 = 0;
		double* vpos = &v[0];
		for (int j=0; j<n; j++){
			z0 += (*mp1) * (*vpos);
			vpos ++;
			mp1 ++;
		}
		*zpos = z0; zpos ++;
	}

}

void mat_t_times_vec_opt (double* m, double* v, int p, int n, double* z){
	double* vpos = &v[0];
	double* mp1 = &m[0]; double *mp2, *mp3, *mp4;
	if (p>3){
		mp2 = &m[n];
		mp3 = &m[2*n];
		mp4 = &m[3*n];
	}
	int k = p/4;
	for (int i=0; i<k; i++){
		double* zpos = &z[0];
		double v1 = *vpos; vpos ++; 
		double v2 = *vpos; vpos ++; 
		double v3 = *vpos; vpos ++; 
		double v4 = *vpos; vpos ++; 

		for (int j=0; j<n; j++){
			*zpos += (*mp1) * v1;
			*zpos += (*mp2) * v2;
			*zpos += (*mp3) * v3;
			*zpos += (*mp4) * v4;
			zpos ++;
			mp1 ++;
			mp2 ++;
			mp3 ++;
			mp4 ++;
		}
		
		mp1 += 3 * n;
		mp2 += 3 * n;
		mp3 += 3 * n;
		mp4 += 3 * n;
	}
	
	for (int i=k*4; i<p; i++){
		double* zpos = &z[0];
		double v0 = *vpos; vpos ++;
		for (int j=0; j<n; j++){
			*zpos += (*mp1) * v0;
			zpos ++;
			mp1 ++;
		}
	}

	return;
}


// M (p x n) V (n x 1) = Z
void mat_times_vec (double* m, double* v, int p, int n, double* z){	
	int i = 0, j = 0, r = 0;
	int np = n*p;
	double s = 0.0;
	for (i = 0; i < np; i++){
		s += m[i] * v[j];
		j ++; 
		if ( j == n){
			z[r] = s;
			r ++;
			s = 0.0;
			j = 0;
		} 
	}
	return;
}

// t(M (p x n)) V (p x 1) = Z
void mat_t_times_vec (double* m, double* v, int p, int n, double* z){
	int i = 0, j = 0, r = 0;
	for (i = 0; i < n; i++){z[i] = 0;}
	int np = n*p;
	for (i = 0; i < np; i++){
		z[j] += m[i] * v[r];
		j ++; 
		if ( j == n){
			r ++;
			j = 0;
		} 
	}
	return;
}

void mat_vec_mul (double* m, double* v, int p, int n, double* z, int trans){
	if (trans == 0){
		mat_times_vec_opt (m, v, p, n, z);
	}else{
		memset(z, 0, sizeof(double) * n);
		mat_t_times_vec_opt (m, v, p, n, z);
	}
}


double L2_norm (double* v, int p){
	double s = 0.0;
	for (int i=0; i<p; i++){
		s += v[i] * v[i];
	}
	return sqrt(s);
}


double vec_vec_mul (double* v1, double* v2, int p){
	double* vp1 = &v1[0]; double* vp2 = &v2[0];
	double s = 0.0;
	for (int i=0; i<p; i++){
		s += (*vp1) * (*vp2);
		vp1 ++;
		vp2 ++;
	}
	return s;
}

double vec_vec_mul (float* v1, double* v2, int p){
	float* vp1 = &v1[0]; double* vp2 = &v2[0];
	double s = 0.0;
	for (int i=0; i<p; i++){
		s += (*vp1) * (*vp2);  // precision
		vp1 ++;
		vp2 ++;
	}
	return s;
}


double vec_vec_mul (float* v1, float* v2, int p){
//	cout.precision(15); cout << "V-V" << endl;
	float* vp1 = &v1[0]; float* vp2 = &v2[0];
	double s = 0.0;
	for (int i=0; i<p; i++){
//		cout << s << " ";
		s += (*vp1) * (*vp2); // precision 
		vp1 ++;
		vp2 ++;
	}
//	cout << endl;
	return s;
}

void center_row(double* m, int n, int calc_var){
	double sum = 0.0;
	for (int i=0; i<n; i++){
		sum += m[i];
	}
	double mean = sum / n;
	if (calc_var == 1){
		double var = 0.0;
		for (int i=0; i<n; i++){
			m[i] = m[i] - mean;
			var += m[i] * m[i];
		}
		var = var/n;
		g_snp_mean.push_back(mean);
		g_snp_var.push_back(var);
	}else{
		for (int i=0; i<n; i++){
			m[i] = m[i] - mean;
		}
	}
	return;
}

// X is p x n
void center_rows(double* m, int p, int n){
	record_time();
	for (int i=0; i<p; i++){
		center_row(&m[i * n], n, 1);
	}
	record_time("Centering");
	return;
}	

void center_row(float* m, int n, int calc_var){
	double sum = 0.0;
	for (int i=0; i<n; i++){
		sum += m[i];
	}
	double mean = sum / n;
	if (calc_var == 1){
		double var = 0.0;
		for (int i=0; i<n; i++){
			m[i] = m[i] - mean;
			var += m[i] * m[i];
		}
		var = var/n;
		g_snp_mean.push_back(mean);
		g_snp_var.push_back(var);
	}else{
		for (int i=0; i<n; i++){
			m[i] = m[i] - mean;
		}
	}
	return;
}

// X is p x n
void center_rows(float* m, int p, int n){
	record_time();
	for (int i=0; i<p; i++){
		center_row(&m[i * n], n, 1);
	}
	record_time("Centering");
	return;
}	


double center_array(double* v, int p){
	double ss1 = 0.0;
	double ss2 = 0.0;
	for (int i=0; i<p; i++){
		ss1 += v[i];
		ss2 += v[i] * v[i];
	}
	double sst =  ss2 - ss1 * ss1 / (double) p;
	double mean = ss1 / p;
	for (int i=0; i<p; i++){
		v[i] -= mean;	
	}
	return sst / p;
}


void array_del (const vector<int>& del, int p_last, double* v_last, double* v){
	int s = del.size();
	int c = 0;
	int d = 0;
	for (int i=0; i< s; i++){
		int r = del.at(i);
		if (r - c > 0){
			copy1D(&v[d], &v_last[c], r - c);
		}
		d = d + r - c;
		c = r + 1;		
	}
	if (c < p_last){	
		copy1D(&v[d], &v_last[c], p_last - c);
	}
	return;
}

void matrix_del (const vector<int>& del, int p_last, int n_last, double* m_last, double* m){
	int s = del.size();
	int c = 0;
	int d = 0;
	for (int i=0; i< s; i++){
		int r = del.at(i);
		if (r - c > 0){
			copy1D(&m[d * n_last], &m_last[c * n_last], (r-c) * n_last);
		}
		d = d + r - c;
		c = r + 1;		
	}
	if (c < p_last){	
		copy1D(&m[d * n_last], &m_last[c * n_last], (p_last - c) * n_last);
	}
	return;
}

void matrix_del (const vector<int>& del, int p_last, int n_last, float* m_last, float* m){
	int s = del.size();
	int c = 0;
	int d = 0;
	for (int i=0; i< s; i++){
		int r = del.at(i);
		if (r - c > 0){
			copy1D(&m[d * n_last], &m_last[c * n_last], (r-c) * n_last);
		}
		d = d + r - c;
		c = r + 1;		
	}
	if (c < p_last){	
		copy1D(&m[d * n_last], &m_last[c * n_last], (p_last - c) * n_last);
	}
	return;
}


double array_sst(double* v, int p){
	double ss1 = 0.0;
	double ss2 = 0.0;
	for (int i=0; i<p; i++){
		ss1 += v[i];
		ss2 += v[i] * v[i];
	}
	return ss2 - ss1 * ss1 / (double) p;
}



void mat_times_vec_opt (float* m, double* v, int p, int n, double* z){	
	double* zpos = &z[0];
	float* mp1 = &m[0]; float *mp2, *mp3, *mp4;
	if (p>3){
		mp2 = &m[n];
		mp3 = &m[2*n];
		mp4 = &m[3*n];
	}

	int k = p/4;
	for (int i=0; i<k; i++){
		double z1 = 0, z2 = 0, z3 = 0, z4 = 0;
		double* vpos = &v[0];
		for (int j=0; j<n; j++){
			z1 += (*mp1) * (*vpos);
			z2 += (*mp2) * (*vpos);
			z3 += (*mp3) * (*vpos);
			z4 += (*mp4) * (*vpos);
			vpos ++;
			mp1 ++;
			mp2 ++;
			mp3 ++;
			mp4 ++;
		}
		*zpos = z1;  zpos ++;
		*zpos = z2;  zpos ++;
		*zpos = z3;  zpos ++;
		*zpos = z4;  zpos ++;
		
		mp1 += 3 * n;
		mp2 += 3 * n;
		mp3 += 3 * n;
		mp4 += 3 * n;
	}
	
	for (int i=k*4; i<p; i++){
		double z0 = 0;
		double* vpos = &v[0];
		for (int j=0; j<n; j++){
			z0 += (*mp1) * (*vpos);
			vpos ++;
			mp1 ++;
		}
		*zpos = z0; zpos ++;
	}

}

void mat_t_times_vec_opt (float* m, double* v, int p, int n, double* z){
	double* vpos = &v[0];
	float* mp1 = &m[0]; float *mp2, *mp3, *mp4;
	if (p>3){
		mp2 = &m[n];
		mp3 = &m[2*n];
		mp4 = &m[3*n];
	}
	int k = p/4;
	for (int i=0; i<k; i++){
		double* zpos = &z[0];
		double v1 = *vpos; vpos ++; 
		double v2 = *vpos; vpos ++; 
		double v3 = *vpos; vpos ++; 
		double v4 = *vpos; vpos ++; 

		for (int j=0; j<n; j++){
			*zpos += (*mp1) * v1;
			*zpos += (*mp2) * v2;
			*zpos += (*mp3) * v3;
			*zpos += (*mp4) * v4;
			zpos ++;
			mp1 ++;
			mp2 ++;
			mp3 ++;
			mp4 ++;
		}
		
		mp1 += 3 * n;
		mp2 += 3 * n;
		mp3 += 3 * n;
		mp4 += 3 * n;
	}
	
	for (int i=k*4; i<p; i++){
		double* zpos = &z[0];
		double v0 = *vpos; vpos ++;
		for (int j=0; j<n; j++){
			*zpos += (*mp1) * v0;
			zpos ++;
			mp1 ++;
		}
	}

	return;
}

void mat_vec_mul (float* m, double* v, int p, int n, double* z, int trans){
	if (trans == 0){
		mat_times_vec_opt (m, v, p, n, z);
	}else{
		memset(z, 0, sizeof(double) * n);
		mat_t_times_vec_opt (m, v, p, n, z);
	}
}

// float** double* 
void mat_times_vec_opt (float** m, double* v, int p, int n, double* z){	
	double* zpos = &z[0];
	int k = p/4;
	for (int i=0; i<k; i++){
		float* mp1 = m[i*4];
		float* mp2 = m[i*4 + 1];
		float* mp3 = m[i*4 + 2];
		float* mp4 = m[i*4 + 3];	
		double z1 = 0, z2 = 0, z3 = 0, z4 = 0;
		double* vpos = &v[0];
		for (int j=0; j<n; j++){
			z1 += (*mp1) * (*vpos);
			z2 += (*mp2) * (*vpos);
			z3 += (*mp3) * (*vpos);
			z4 += (*mp4) * (*vpos);
			vpos ++;
			mp1 ++;
			mp2 ++;
			mp3 ++;
			mp4 ++;
		}
		*zpos = z1;  zpos ++;
		*zpos = z2;  zpos ++;
		*zpos = z3;  zpos ++;
		*zpos = z4;  zpos ++;
	}
	
	for (int i=k*4; i<p; i++){
		double z0 = 0;
		double* vpos = &v[0];
		float* mp1 = m[i];
		for (int j=0; j<n; j++){
			z0 += (*mp1) * (*vpos);
			vpos ++;
			mp1 ++;
		}
		*zpos = z0; zpos ++;
	}

}

void mat_t_times_vec_opt (float** m, double* v, int p, int n, double* z){
	double* vpos = &v[0];
	int k = p/4;
	for (int i=0; i<k; i++){
		double* zpos = &z[0];
		double v1 = *vpos; vpos ++; 
		double v2 = *vpos; vpos ++; 
		double v3 = *vpos; vpos ++; 
		double v4 = *vpos; vpos ++; 
		float* mp1 = m[i*4];
		float* mp2 = m[i*4 + 1];
		float* mp3 = m[i*4 + 2];
		float* mp4 = m[i*4 + 3];		
		for (int j=0; j<n; j++){
			*zpos += (*mp1) * v1;
			*zpos += (*mp2) * v2;
			*zpos += (*mp3) * v3;
			*zpos += (*mp4) * v4;
			zpos ++;
			mp1 ++;
			mp2 ++;
			mp3 ++;
			mp4 ++;
		}
	}
	
	for (int i=k*4; i<p; i++){
		double* zpos = &z[0];
		double v0 = *vpos; vpos ++;
		float* mp1 = m[i];
		for (int j=0; j<n; j++){
			*zpos += (*mp1) * v0;
			zpos ++;
			mp1 ++;
		}
	}

	return;
}

void mat_vec_mul (float** m, double* v, int p, int n, double* z, int trans){
//	print_array(v, 3);
//	cout << "Mat-vec\t" << p << "\t" << n << "\t" << trans << endl;
	if (trans == 0){
		mat_times_vec_opt (m, v, p, n, z);
	}else{
		memset(z, 0, sizeof(double) * n);
		mat_t_times_vec_opt (m, v, p, n, z);
	}
//	cout << "Mat-vec-done" << endl;
}

// float** float* 
void mat_times_vec_opt (float** m, float* v, int p, int n, double* z){	
	double* zpos = &z[0];
	int k = p/4;
	for (int i=0; i<k; i++){
		float* mp1 = m[i*4];
		float* mp2 = m[i*4 + 1];
		float* mp3 = m[i*4 + 2];
		float* mp4 = m[i*4 + 3];	
		double z1 = 0, z2 = 0, z3 = 0, z4 = 0; // float?
		float* vpos = &v[0];
		for (int j=0; j<n; j++){
			z1 += (*mp1) * (*vpos);
			z2 += (*mp2) * (*vpos);
			z3 += (*mp3) * (*vpos);
			z4 += (*mp4) * (*vpos);
			vpos ++;
			mp1 ++;
			mp2 ++;
			mp3 ++;
			mp4 ++;
		}
		*zpos = z1;  zpos ++;
		*zpos = z2;  zpos ++;
		*zpos = z3;  zpos ++;
		*zpos = z4;  zpos ++;
	}
	
	for (int i=k*4; i<p; i++){
		float z0 = 0;
		float* vpos = &v[0];
		float* mp1 = m[i];
		for (int j=0; j<n; j++){
			z0 += (*mp1) * (*vpos);
			vpos ++;
			mp1 ++;
		}
		*zpos = z0; zpos ++;
	}

}

void mat_t_times_vec_opt (float** m, float* v, int p, int n, double* z){
	float* vpos = &v[0];
	int k = p/4;
	for (int i=0; i<k; i++){
		double* zpos = &z[0];
		float v1 = *vpos; vpos ++; 
		float v2 = *vpos; vpos ++; 
		float v3 = *vpos; vpos ++; 
		float v4 = *vpos; vpos ++; 
		float* mp1 = m[i*4];
		float* mp2 = m[i*4 + 1];
		float* mp3 = m[i*4 + 2];
		float* mp4 = m[i*4 + 3];		
		for (int j=0; j<n; j++){
			*zpos += (*mp1) * v1;
			*zpos += (*mp2) * v2;
			*zpos += (*mp3) * v3;
			*zpos += (*mp4) * v4;
			zpos ++;
			mp1 ++;
			mp2 ++;
			mp3 ++;
			mp4 ++;
		}
	}
	
	for (int i=k*4; i<p; i++){
		double* zpos = &z[0];
		float v0 = *vpos; vpos ++;
		float* mp1 = m[i];
		for (int j=0; j<n; j++){
			*zpos += (*mp1) * v0;
			zpos ++;
			mp1 ++;
		}
	}

	return;
}

void mat_vec_mul (float** m, float* v, int p, int n, double* z, int trans){
//	cout << "Mat-vec" << endl;
	if (trans == 0){
		mat_times_vec_opt (m, v, p, n, z);
	}else{
		memset(z, 0, sizeof(double) * n);
		mat_t_times_vec_opt (m, v, p, n, z);
	}
//	cout << "Mat-vec-done" << endl;
}






