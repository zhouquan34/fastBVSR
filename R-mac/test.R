source("bvsr.R")
set.seed(1)
n = 100
p = 20
X <- matrix(rnorm(n * p), nrow=n)
y <- X %*% c(rexp(5), rep(0, p-5)) + rnorm(n)
X <- scale(X, scale=FALSE, center = TRUE)
y <- y - mean(y)
my.bvsr = bvsr(X, y, iter=10000, burn=1000,  model=TRUE, yhat=TRUE, suppress=FALSE, save=0)

### Check the dimensions of the returned objects ## 
sum(my.bvsr$freq)
length(my.bvsr$freq)
dim(my.bvsr$path)
length(my.bvsr$model)
dim(my.bvsr$yhat)

### Here are two ways to compute the posterior mean model size ##
weighted.mean(my.bvsr$path[,1], my.bvsr$freq)
weighted.mean(as.numeric(lapply(my.bvsr$model, length)), my.bvsr$freq)

### So this is the posterior mean of sigma ## 
weighted.mean(my.bvsr$path[,2], my.bvsr$freq)
### and this is the posterior mean of Bayesian R-squared ##
weighted.mean(my.bvsr$path[,3], my.bvsr$freq)

## Actually, 'yhat' can be recoverd using 'model' and 'path'. Here is an example. ## 
x0 = as.matrix(X[,my.bvsr$model[[1]]])
s = my.bvsr$path[1,2]
y1 = x0 %*% solve(t(x0) %*% x0 + diag(s^(-2), ncol(x0)), t(x0)) %*% y
head(cbind(y1, my.bvsr$yhat[1,]))
 
## Suppose we have a new dataset Z. To do prediction, one may try the following code. ## 
Z = matrix(rnorm(10 * p), nrow=10)
y.new = Z %*% my.bvsr$beta[, 4]
## One can also use beta[, 3]. The difference should be negligible in most cases. ##

