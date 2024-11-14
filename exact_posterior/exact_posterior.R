## calculate the exact posterior distribution
args = commandArgs(TRUE)
folder = args[1]

# prior for pi is Beta(alpha, beta)
pi_min = exp(-10)
pi_max = 1
pi_alpha = 0
pi_beta = 1

dat = read.table("toy.txt", header=F, row.names=1)
X = t(as.matrix(dat[,-(1:2)]))
p = ncol(X)
n = nrow(X)
y = scan("toy.ph", quiet=T)
X = scale(X, scale=F, center=T) # centering but no scaling 
y = y - mean(y) # centering 
yy = as.numeric(t(y) %*% y)
g_var = apply(X, 2, var) * (n-1)/n # MLE estimates for covariate variances 

# I do not consider the empty model 
# assume uniform prior on h2, and integrate out pi ~ Beta(pi_alpha, pi_beta) restricted to [pi_min, pi_max]
prior = function(gamma, h2){
	k = length(gamma)
	a1 = lbeta(pi_alpha + k, pi_beta + p - k) 
	a2 = log(pbeta(pi_max, pi_alpha + k, pi_beta + p - k) - pbeta(pi_min, pi_alpha + k, pi_beta + p -k))
	return(a1 + a2)
}

# log-posterior of a model gamma and herit h2
post = function(gamma, h2){
	Z = X[,gamma] 
	ss = sum(g_var[gamma])
	sigma = sqrt(h2/(1-h2)/ss) # sigma is obtained from h2 and ss 
	k = length(gamma)
	A = t(Z) %*% Z + diag(sigma^(-2), k)
	RSS = as.numeric(yy - t(y) %*% Z %*% solve(A,  t(Z) %*% y))
	log_post = -n/2 * log(RSS/yy) - k * log(sigma) - 0.5 * as.numeric(determinant(A, log=TRUE)$mod) + prior(gamma, h2)
	return(log_post)
}

# for each model, integrate out h2  by simple Riemann sum 
h2_seq = seq(0.01, 0.99, by = 0.0005)
n_h2 = length(h2_seq)
NM = 2^p - 1
models = character(NM)
posts = numeric(NM)
models_num = matrix(0, nrow=NM, ncol=p)
for (j in 1:NM){
	binary = as.logical(intToBits(j)[1:p])
	gamma = which(binary)
	models_num[j, ] = as.numeric(binary)
	lpost = numeric(n_h2)
	for (hi in 1:n_h2){
		lpost[hi] = post(gamma, h2_seq[hi])
	}
	lp_max = max(lpost)
	post_gamma = log(sum(exp(lpost - lp_max))) + lp_max - log(length(n_h2))
	models[j] = paste(as.character(gamma), collapse=",")
	posts[j] = post_gamma
}

cat("Marginal PIP\n")
pp = exp(posts)/sum(exp(posts))
for (j in 1:p){
	sel = which(models_num[,j] == 1)
	cat(j, sum(pp[sel]), "\n")
}

cat("\nModel Posterior\n")
for (i in 1:length(models)){
	cat("model = ", models[i], "\tposterior = ", pp[i], "\n", sep="")
}


