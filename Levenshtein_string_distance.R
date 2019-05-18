StringDistance <- function(s,t) {
	m <- nchar(s)
	n <- nchar(t)
	d <- matrix(0,nrow=m+1,ncol=n+1)

	for (i in 1:m) {
		d[i+1,1] <- i
	}
	for (j in 1:n) {
		d[1,j+1] <- j
	}

	for (j in 1:n) {
		for (i in 1:m) {
			cost <- 0
			if (substr(s,i,i) == substr(t,j,j)) {
			} else {
				cost <- 1
			}
			d[i+1,j+1] <- min(
				d[i,j+1] + 1, # deletion
				d[i+1,j] + 1, # insertion
				d[i,j] + cost # substitution
			)
		}
	}
	return(d[m+1,n+1])
}

# res <- StringDistance("sitting","kitten")