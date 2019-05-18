# Custom.t.test is a helper function to use
# the apply function to calculate the t.test across
# columns in a large dataframe
# Only send dataframes without string columns so apply can safely convert
# to matrix
Custom.t.test <- function(df,cols1,cols2,na.value=NA) {
  if ( all(df[cols1] == df[cols2]) ) {
    return(1)
  } 
  # In some rare cases, all values in cols1 are equal to each other and
  # all values in cols2 are equal to each other
  # In this case there is no sd, so p.Value can't be calculated
    obj<-try(t.test(df[cols1],df[cols2],var.equal=TRUE), silent=TRUE)
  if (is(obj, "try-error")) return(na.value) else return(obj$p.value)
#    return(t.test(df[cols1],df[cols2],var.equal=TRUE)$p.value)
}
# apply(df[,2:ncol(df)], 1, Custom.t.test, col.wt, col.02)



# RemoveZeroRows removes all rows from dataframe df in which 
# all columns of vector column have 0 values
RemoveZeroRows <- function(df,columns) {
  l <- length(columns)
  df$RemoveZeroRows <- rep(FALSE,length(df[,1]))
  for (i in 1:l) {
    df$RemoveZeroRows[ df[,columns[i]] > 0] <- TRUE
  }
  df <- df[(df$RemoveZeroRows),]
  df$RemoveZeroRows <- NULL
  return(df)
}

# Variance returns a dataframe containing the means and variance within events
Variance <- function(v.df) {
  means   <- rowMeans(v.df)
  std_dev <- apply(v.df,1,sd,na.rm=TRUE)
  variance <- data.frame(
  MEANS=means,
  fluctuation= log(std_dev/means)/log(2)
  )
  return(variance)
}

# DiffExpressed returns dataframes of genes that have different expression
DiffExpressed <- function(df.temp,col.base,col.intr) {
  df.temp$ind <- ( df.temp[,col.intr] / df.temp[,col.base] )
  x2  <- df.temp[ ( (df.temp$ind  >= 2) & (df.temp$ind  < 5) ),]
  x5  <- df.temp[ ( (df.temp$ind  >= 5) & (df.temp$ind  < 10) ),]
  x10 <- df.temp[ (  df.temp$ind  >= 10) ,]
  m <- list(x2,x5,x10)
  return(m) 
}

# MatchExpressionWithGeneList returns dataframes containing genes that
# are present in both Dataframe and list
# e.g Matrix: List of dataframes with differential expressed genes
# list: GO class genelist
MatchExpressionWithGeneList <- function(matrix,list,identifier="gene_short_name") {
  m <-list()
  for ( i in 1:length(matrix) ) {
    df.temp <- data.frame(matrix[i])
    hits <- df.temp[(df.temp[,identifier] %in% list), ]
    m[i] <- list(V=hits)
  }
  return(m)
}

# MakeSubsetLists first calculates fold change between two conditions and then returns a table 
# for each subset of combined with genelists
MakeSubsetLists <- function(df, col1, col2, list.dir, title=paste(col2,col1,sep=".vs."), list.pattern=".txt") {
  dst <- file.path("overlaps/Tables",title) # relative to working directory
  dir.create(dst)
  # Calculate fold change and print it to files
  m1 <- DiffExpressed(df,col1,col2)
  for (i in 1:length(m1)) {
    write.table(data.frame(m1[i]), paste(dst,"/", paste(title, i, sep="_"), ".xls", sep=""), sep="\t",quot=F,row.names=F )
  }
  
  # Analyze what kind of genes have increased expression after stimulation
  # Print the lists for each of those genes
  lists <- dir(path = list.dir, pattern = list.pattern)

  for (i in 1:length(lists)) {
    file <- file.path(list.dir,lists[i])
    list <- unique( scan(file, what="character") )
    m <- MatchExpressionWithGeneList(m1,list)
    for (j in 1:length(m)) {
      write.table(data.frame(m[j]), paste(dst,"/",lists[i],"_",j,".xls",sep=""), sep="\t",quot=F,row.names=F )
    }
  }
  return(title)
}


Classify <- function(df,
  columns,
  input,
  fold=1.5,
  threshold=5

  ) {

  flag = 0
  # FLAG:
  # 1: present at 12h
  # 2: present at 8h
  # 3: present at 4h
  # 4: prsent at 0h

  # present at 0h
  if ( (as.numeric(df[columns[1]]) >= threshold) &
       (as.numeric(df[columns[1]]) > as.numeric(df[input[1]]) * fold)
   ) {
    flag <- flag + 1*2**0 # %% 2 = 1
  }

  # present at 4h
  if ( (as.numeric(df[columns[2]]) >= threshold) &
       (as.numeric(df[columns[2]]) > as.numeric(df[input[2]]) * fold)
   ) {
    flag <- flag + 1*2**1  # %% 4 = 1
  }

  # present at 8h
  if ( (as.numeric(df[columns[3]]) >= threshold) &
       (as.numeric(df[columns[3]]) > as.numeric(df[input[3]]) * fold)
   ) {
    flag <- flag + 1*2**2
  }

  # present at 4h
  if ( (as.numeric(df[columns[4]]) >= threshold) &
       (as.numeric(df[columns[4]]) > as.numeric(df[input[4]]) * fold)
   ) {
    flag <- flag + 1*2**3
  }

  return(flag)
}

checkflag <- function(column,position) {
  return(ifelse(bitwShiftR(column,position) %% 2 == 1, TRUE,FALSE))
}

checkflags <- function(column, positions,results) {
  if ( length(positions) != length(results) ) {
    return(FALSE)
  }
  return( ifelse(
    (
    bitwShiftR(column,positions[1]) %% 2 == results[1] &
    bitwShiftR(column,positions[2]) %% 2 == results[2] &
    bitwShiftR(column,positions[3]) %% 2 == results[3] &
    bitwShiftR(column,positions[4]) %% 2 == results[4] 
    ),TRUE,FALSE) )

}


# df.n$class <- apply(df.n,1,Classify,cols.brd3,cols.input)
# sub3 <- subset(df.n, ifelse(bitwShiftR(class,3) %% 2 == 1, TRUE,FALSE))
# sub3.2 <- subset(df.n, checkflag(class,3))

# Creates a data frame with sliding bins. Receives df with two columns: x & y
# x will be used for sliding, 
# y will be used to calculate mean, extremes and 95% confidence interval values for each bin
slidingBins <- function(df,binsize=200,step=40,column="y",sorting="x") {
  new.df <- df[order(df[,sorting]),]
  x_vector <- c()
  y_vector <- c()
  y_top_vector <- c()
  y_low_vector <- c()
  y_top_extreme <- c()
  y_low_extreme <- c()
  for (i in seq(1,nrow(df)-binsize,by=step)) {
    x <- new.df[i:(i+binsize),sorting]
    y <- new.df[i:(i+binsize),column]
    stat <- t.test(y)$conf.int
    x_vector <- append(x_vector,mean(x))
    y_vector <- append(y_vector,mean(y))
    y_top_vector <- append(y_top_vector,ifelse(is.na(stat[2]), mean(y), stat[2]))
    y_low_vector <- append(y_low_vector,ifelse(is.na(stat[1]), mean(y), stat[1]))
    y_top_extreme <- append(y_top_extreme,max(y))
    y_low_extreme <- append(y_low_extreme,min(y))
  }
  return(data.frame(
    x=x_vector,
    y=y_vector,
    y_conf_min=y_low_vector,
    y_conf_max=y_top_vector,
    y_min=y_low_extreme,
    y_max=y_top_extreme
  ))
}

# Like slidingBins(), but instead of returning a 95% confidence interval 
# form t.test, it returns the 95% highest and lowest value in the set
slidingBins.2 <- function(df,binsize=200,step=40,column="y",sorting="x") {
  new.df <- df[order(df[,sorting]),]
  x_vector <- c()
  y_vector <- c()
  y_top_vector <- c()
  y_low_vector <- c()
  y_top_extreme <- c()
  y_low_extreme <- c()
  for (i in seq(1,nrow(df)-binsize,by=step)) {
    x <- new.df[i:(i+binsize),sorting]
    y <- new.df[i:(i+binsize),column]
    stat <- t.test(y)$conf.int
    x_vector <- append(x_vector,mean(x))
    y_vector <- append(y_vector,mean(y))
    y_top_vector <- append(y_top_vector,quantile(y,0.95))
    y_low_vector <- append(y_low_vector,quantile(y,0.05))
    y_top_extreme <- append(y_top_extreme,max(y))
    y_low_extreme <- append(y_low_extreme,min(y))
  }
  return(data.frame(
    x=x_vector,
    y=y_vector,
    y_conf_min=y_low_vector,
    y_conf_max=y_top_vector,
    y_min=y_low_extreme,
    y_max=y_top_extreme
  ))
}
scale0To10 <- function(df,columns) {
  max_value <- max(df[,columns])
  df[,columns] <- df[,columns]*10/max_value
  return(df)
}
