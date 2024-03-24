########### Input ############
# Please put the executable 'fastBVSR' in the R working directory 
# If there is a permission problem, try 'chmod 755 fastBVSR' in the terminal
# 'iter' is the number of MCMC sampling iterations
# 'burn' is the number of MCMC burn-in iterations
# 'save = 1' if you want to keep all the BVSR output files
# 'file' is the prefix for all the BVSR output files
# 'seed' is the initial random seed used by BVSR 
# 'suppress = TRUE/FALSE' is used to suppress the messages of BVSR  
# 'yhat = TRUE/FALSE' means whether to return the fitted values of y 
# 'model = TRUE/FALSE' means whether to return the sampled models

########### Output ############
# When iter or p is large, the objects 'yhat' and 'model' could be extremely large and use up the memory. 
# A list is returned which consists of up to 5 objects.
# 'beta' is a matrix with p rows. They are the posterior mean estimates. 
# 'freq' is a vector that stores the frequencies of each sampled model in MCMC. The sum should be equal to 'iter'. We use k to denote its length, which may change each run. 
# 'path' is a matrix with k rows. It gives the model size, sigma, and Bayesian R-squared of each sampled model. 
# 'yhat' is a matrix with k rows and n columns. The i-th row gives the fitted values of Y using the i-th sampled model. 
# 'model' is a list with k elements. The i-th element gives the IDs of the selected covariates for the i-th sampled model. 
# Please try the test.R file, which explains how the posterior means of the hyperparameters can be calculated. 

bvsr <- function(X, y, iter=10000, burn=1000, yhat=FALSE, model=FALSE, file='', save=0, seed=-1, suppress=TRUE){
	if (file == ''){
		file = paste("try", Sys.Date(), as.integer(Sys.time()), sep="-" ) 
	}	
	if (seed == -1){
		seed = floor(as.numeric(Sys.time())*100)
	}
	mat.file=paste(file, "mat", sep=".")
	ph.file=paste(file, "ph", sep=".") 
	write.table(X, mat.file, row.names=F, col.names=T, quote=F)
	write.table(y, ph.file, row.names=F, col.names=F, quote=F)
	
	cmd <- paste( "./fastBVSR -m",  mat.file, "-p",  ph.file, "-r", seed, "-o", file, "-w", format(burn, scientific=FALSE), "-s", format(iter, scientific=F), sep=" ")
	if (yhat == TRUE){
		cmd = paste(cmd, "--yhat", sep = " ")
	}
	system(cmd, ignore.stdout=TRUE, ignore.stderr=suppress, intern=TRUE)
	
	beta = as.matrix(read.table(paste(file, ".beta.txt", sep=""), header=TRUE, row.names=1))
	beta = beta[, -1]
	colnames(beta) = c('post-inclusion-prob', 'post-inclusion-prob-RB', 'beta', 'beta-RB')
	
	path = as.matrix(read.table(paste(file, ".path.txt", sep=""), header=TRUE))
	iter = path[,1]
	freq = path[,2]
	start = 1 
	for (i in 1:length(iter)){
		if (iter[i] <= burn){freq[i] = 0}
		if (iter[i] > burn){freq[i] = iter[i] - burn; start = i; break;}
	}
	freq = freq[start:length(freq)]	
	path = path[start:nrow(path), c(3,7,8)]
	
	res = list("beta" = beta, "freq" = freq, "path" = path)

	if (model == TRUE){
		mod.list <- list()
		mod.file <- paste(file, ".model.txt", sep="")
		con  <- file(mod.file, open = "r")
		ln = 0 
		while (length(oneLine <- readLines(con, n = 1, warn = FALSE)) > 0) {
  			ln = ln + 1
			if (ln <= start){next}
			ids <- strsplit(oneLine, ",")
			ids <- as.numeric(ids[[1]])
			ids = sort(ids + 1)
			mod.list[[ln - start]] = ids
		} 
		close(con)			
		res[["model"]] = mod.list
	}
	
	if (yhat == TRUE){
		ys = as.matrix(read.table(paste(file, ".yhat.txt", sep=""), header=FALSE))
		ys = ys[start:nrow(ys), ]	
		res[["yhat"]] = ys 
	}

	if (save == 0){
		system(paste("rm ", file, "*", sep=""))
	}

	return(res)
}


