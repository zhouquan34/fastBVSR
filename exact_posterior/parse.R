data <- readLines('tmp')
n <- length(data)
count <- list()
total <- 0
for (i in 1:n){
	row <- data[i]
	tmp <- unlist(strsplit(row, "\t"))
	set <- tmp[2]
	life <- as.numeric(tmp[1])
	sorted <- sort(as.numeric(unlist(strsplit(set, ","))))
	sorted <- sorted + 1 
	model <- paste(sorted, collapse=",")
	if (is.null(count[[model]])){
		count[[model]] = 0
	}
	count[[model]] = count[[model]] + life
	total = total + life
}

for (m in names(count)){
	cat(m, count[[m]]/total, "\n", sep="\t")
}

