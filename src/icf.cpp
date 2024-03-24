/* we assume the first c columns of X represent confounding covariates.
   In total X contains p columns. 
   The diagonal regularization is a vector with first c elements = 0. */
#include "icf.h"

using namespace std;

// R^t K x - K R x  = b
void D_times_vec(const double* R, const double* x, double* b, const double* diags, const int p, const int c){
	double* x1 = allocate1D(p);
	double* b1 = allocate1D(p, 0);
	int a=0; 
	for (a=0; a<p; a++){ 
		x1[a] = x[a] * diags[a];
		b[a] = 0;
	}
	// R[k] = (matrix)R[i, j]
	int i=c, j=c+1, k=c*p+c+1;
	
	while (i < (p-1) ){
		b1[i] -= x[j] * R[k]; 
		b[j] += x1[i] * R[k];  
		k ++;
		j ++;
		if (j == p ){
			i ++ ;
			j = i + 1;
			k += (i + 1);
		}	
	}

	for (a=0; a<p; a++){ 
		b[a] += b1[a] * diags[a]; 
	}

	free1D(b1);
	free1D(x1);
	return;
}

// (a + bi) (x + yi) = c + di;  x = (ac + bd)/(a2 + b2)   y = (ad - bc)/(a2 + b2).
void R_solve(const double* R, double** b, const int trans, const int p, const int c, const double* diags){
	int i, j;
	if (trans == 0){ // no conj trans	
		for (i = p - 1; i >= c; i--){
			double re = b[0][i], im = b[1][i];
			for (j = p - 1; j > i; j --){
				re -= b[0][j] * R[i * p + j];				
				im -= b[1][j] * R[i * p + j];
			}
			double m = diags[i] * diags[i] + R[i * p + i] * R[i * p + i];
			b[0][i] = ( R[i * p + i] * re + im * diags[i] )/m;
			b[1][i] = ( R[i * p + i] * im - re * diags[i] )/m;
		}
		for (i = c - 1; i >= 0; i--){
			double re = b[0][i], im = b[1][i];
			for (j = p - 1; j > i; j --){
				re -= b[0][j] * R[i * p + j];				
				im -= b[1][j] * R[i * p + j];
			}
			b[0][i] = re / R[i * p + i];
			b[1][i] = im / R[i * p + i];
		}
	}else{
		for (i = 0; i < c; i++){
			double re = b[0][i], im = b[1][i];
			for (j = 0; j < i; j ++){
				re -= b[0][j] * R[j * p + i];				
				im -= b[1][j] * R[j * p + i];
			}
			b[0][i] = re / R[i * p + i];
			b[1][i] = im / R[i * p + i];
		}	
		for (i = c; i < p; i++){
			double re = b[0][i], im = b[1][i];
			for (j = 0; j < i; j ++){
				re -= b[0][j] * R[j * p + i];				
				im -= b[1][j] * R[j * p + i];
			}
			double m = diags[i]*diags[i] + R[i * p + i] * R[i * p + i];
			b[0][i] = ( R[i * p + i] * re - im * diags[i] )/m;
			b[1][i] = ( R[i * p + i] * im + re * diags[i] )/m;
		}	
	}	
	return;
}


int whether_to_stop(const double* a, const double* b, int p){
	for (int i=0; i<p; i++){
		if (fabs(a[i] - b[i]) > g_icf_abs_tol){
			return 0;
		}
	}
	return 1;
}

double average_vec(double* v1, double* v2, double w1, double w2, int p){
	double err = 0.0;
	for (int i=0; i<p; i++){
		double e = v1[i] - v2[i]; 
		err += e*e;		
		v1[i] = v1[i] * w1 + v2[i] * w2;
	}
	return w1 * sqrt(err);
}

int ICF
(const double* R, const double* z, int p, int c, const double* diags, double* b){
	double w = 1.0, cw = 0.0;
	int i=0, j=0;

	double* beta_last = allocate1D(p);
	double** beta = allocate2D(p);	
	double error = 0.0, error_last = 0.0;
	
	int cf_iter = 1;
	copy1D(beta[0], z, p);
	R_solve(R, beta, 1, p, c, diags);
	R_solve(R, beta, 0, p, c, diags);
	error = L2_norm(beta[0], p);
	
	while (1){
		cf_iter ++; 
		error_last = error;
		copy1D(beta_last, beta[0], p);  
		D_times_vec(R, beta[0], beta[1], diags, p, c); 
		copy1D(beta[0], z, p); 
		R_solve(R, beta, 1, p, c, diags); 
		R_solve(R, beta, 0, p, c, diags); 
		error = average_vec(beta[0], beta_last, w, cw, p);
		w = 2*w / (1 + w + error / error_last);
		cw = 1.0 - w;
		
		int stop = whether_to_stop(beta[0], beta_last, p);
		if (stop == 1){break;}	
		if (cf_iter == 10){
			if (error/error_last > 0.8){
				free1D(beta_last);
				free2D(beta);
				return 0;	
			}
		}
	}

	copy1D(b, beta[0], p);
	free1D(beta_last);
	free2D(beta);
	return cf_iter;
}



void D_times_vec(const double* R, const double* x, double* b, const double diag, const int p, const int c){
	int a=0; 
	for (a=0; a<p; a++){ b[a] = 0;}
	int i=c, j=c+1, k=c*p+c+1;
	
	while (i < (p-1) ){
		b[i] -= x[j] * R[k]; 
		b[j] += x[i] * R[k];  
		k ++;
		j ++;
		if (j == p ){
			i ++ ;
			j = i + 1;
			k += (i + 1);
		}	
	}

	for (a=c; a<p; a++){ b[a] *= diag;}
	return;
}

void R_solve(const double* R, double** b, const int trans, const int p, const int c, const double diag){
	int i, j;
	if (trans == 0){ // no conj trans	
		for (i = p - 1; i >= c; i--){
			double re = b[0][i], im = b[1][i];
			for (j = p - 1; j > i; j --){
				re -= b[0][j] * R[i * p + j];				
				im -= b[1][j] * R[i * p + j];
			}
			double m = diag * diag + R[i * p + i] * R[i * p + i];
			b[0][i] = ( R[i * p + i] * re + im * diag )/m;
			b[1][i] = ( R[i * p + i] * im - re * diag )/m;
		}
		for (i = c - 1; i >= 0; i--){
			double re = b[0][i], im = b[1][i];
			for (j = p - 1; j > i; j --){
				re -= b[0][j] * R[i * p + j];				
				im -= b[1][j] * R[i * p + j];
			}
			b[0][i] = re / R[i * p + i];
			b[1][i] = im / R[i * p + i];
		}
	}else{
		for (i = 0; i < c; i++){
			double re = b[0][i], im = b[1][i];
			for (j = 0; j < i; j ++){
				re -= b[0][j] * R[j * p + i];				
				im -= b[1][j] * R[j * p + i];
			}
			b[0][i] = re / R[i * p + i];
			b[1][i] = im / R[i * p + i];
		}	
		for (i = c; i < p; i++){
			double re = b[0][i], im = b[1][i];
			for (j = 0; j < i; j ++){
				re -= b[0][j] * R[j * p + i];				
				im -= b[1][j] * R[j * p + i];
			}
			double m = diag*diag + R[i * p + i] * R[i * p + i];
			b[0][i] = ( R[i * p + i] * re - im * diag )/m;
			b[1][i] = ( R[i * p + i] * im + re * diag )/m;
		}	
	}	
	return;
}


int ICF
(const double* R, const double* z, int p, int c, const double diag, double* b){
	double w = 1.0, cw = 0.0;
	int i=0, j=0;

	double* beta_last = allocate1D(p);
	double** beta = allocate2D(p);	
	double error = 0.0, error_last = 0.0;
	
	int cf_iter = 1;
	copy1D(beta[0], z, p);
	R_solve(R, beta, 1, p, c, diag);
	R_solve(R, beta, 0, p, c, diag);
	error = L2_norm(beta[0], p);
	
	while (1){
		cf_iter ++; 
		error_last = error;
		copy1D(beta_last, beta[0], p);  
		D_times_vec(R, beta[0], beta[1], diag, p, c); 
		copy1D(beta[0], z, p); 
		R_solve(R, beta, 1, p, c, diag); 
		R_solve(R, beta, 0, p, c, diag); 
		error = average_vec(beta[0], beta_last, w, cw, p);
		w = 2*w / (1 + w + error / error_last);
		cw = 1.0 - w;
		
		int stop = whether_to_stop(beta[0], beta_last, p);
		if (stop == 1){break;}	
		if (cf_iter == 10){
			if (error/error_last > 0.8){
				free1D(beta_last);
				free2D(beta);
				return 0;	
			}
		}
	}

	copy1D(b, beta[0], p);
	free1D(beta_last);
	free2D(beta);
	return cf_iter;
}






