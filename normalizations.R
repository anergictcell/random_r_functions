# Used for QuantileNormalization
library('preprocessCore')


# QuantileNormalize normalizes a set of columns within a dataframe 
# to each other for comparison
# based on B. M. Bolstad, et al., Bioinformatics, 2003
# kept the old version for compatibility reasons, in case I need to 
# repeat old analysis
# Basically everything before 9. May 2016 was done with the old one
# All Ezh1,2 Data
# THe initial MMSet Bcell data (RNAseq)
QuantileNormalizeOLD <- function(df, columns) {
  # temp is a temporary dataframe that will contain the sorted values
  # from each column in df
  # initiates temp with the first column from df and 
  # renames the column to give it the same name as in df
  temp <- data.frame(col1=sort( df[,columns[1]]) ) 
  names(temp)[1] <- columns[1]

  # copies each sorted column from df to temp and
  # renames the column to give it the same name as in df
  for ( i in 2:length(columns) ) {
    temp$coltemp <- sort( df[,columns[i]]) 
    names(temp)[i] <- columns[i]
  }
  # calculates the means of each row (this value will be used for normalization)
  means <- rowMeans(temp)
  
  # temporarily orders each column in df by value and then overwrites
  # the values with the normalized means
  for ( i in 1:length(columns) ) {
    df[order(df[,columns[i]]), ][,columns[i]] <- means
  }
  return(df)
}

# Uses preProcessCore library QuantileNormalize
QuantileNormalize2 <- function(df, columns) {
  df.n <- data.frame(normalize.quantiles(as.matrix(df[,columns])))
  names(df.n) <- columns
  df.n$symbol <- df$symbol
  return(df.n[,c("symbol",columns)])
}

SingleQuantileNormalize <- function(df, base, adjust) {
  df$indexnormalize <- row(df)[,1]        # index the rows according to current order so it can be returned in the same order
  ref <- sort(df[, base])                 # generate a reference list with all the values from the reference column
  sorted.df <- df[order(df[, adjust]), ]  # order the dataframe by the column to be normalized
  sorted.df[, adjust] <- ref              # change values in adjust columns to the reference values
  df <- sorted.df[ order(sorted.df$indexnormalize),]  # bring order of dataframe back to original order
  df$indexnormalize <- NULL                           # remove temporary index
  return(df)
}

InvariantNormalize <- function(
  df,         # data frame  
  base,       # column name ==> compare/base = fold induction 
  compare,    # column name ==> compare/base = fold induction 
  adjust,     # column name for data that should be adjusted
  fold=1.5) { # defines range for values to be considered)

  temp <- df[ ( (df[,compare]/df[,base]) < fold & (df[,compare]/df[,base]) > (1/fold) ), ]
  return(temp)
}

Cumulative <- function(df,x_axis) {
  # Receives a dataframe and a name for the column of interest
  # Returns a dataframe with two columns:
  # x => data values sorted ascending
  # y => index of the row devided by numer of rows, 
  #      equals the cumulative fraction of values equal or smaller than current
  t   <- data.frame(x = sort(df[,x_axis]))
  t$y <- row(t)/length(t[,1])
  return(t)
}

# SORT : sorts values in a list
# ORDER : returns the indexes of a sorted list

# DetectionLimit returns the MDL (LOQ) for a dataset
# LOQ describes the smallest concentration of a measurand that can be 
# reliably measured by an analytical procedure
# - Analytic accuracy is calculated by dividing the mean measurand of events by 
# their standard deviation.
# - LOQ is defined as the lowest value in a set of values whose accuray passes
# a defined limit (fold) in a significant manner (percentage)
# Skoog, Douglas A. & Leary, James J. Principles of Instrumental Analysis, 
#  Fourth Edition, Saunders College Publishing, 1992
# Parameters:
# df => Dataframe that contains datasets in columns and each row represents one event
# columns => Vecor of names of the columns that are to be considered
# sort => Boolean value determines if each column should be sorted independently
# percentage => Percent of events that have to have low variability
# fold => Determines threshold for ( mean > fold x SD) 
DetectionLimit <- function(df,columns=FALSE,sort=FALSE,percentage=70,fold=1.5) {
  if (!columns) {
    columns <- names(df)
  }
  if (sort) {
    # each column of the dataframe is sorted independently
     temp <- data.frame(col1=sort(df[,columns[1]]) )
     names(temp)[1] <- columns[1]
     # copies each sorted column from df to temp and
     # renames the column to give it the same name as in df
     for ( i in 2:length(columns) ) {
       temp$coltemp <- sort(df[,columns[i]])
       names(temp)[i] <- columns[i]
    }
  } else {
    # the orders of the columns stays constant
    temp <- df[columns]
  }
  # create working dataframe dist that contains describes the significance for each mean
  means   <- rowMeans(temp)
  std_dev <- apply(temp,1,sd,na.rm=TRUE)
  dist <- data.frame(MEANS=means,SIG=(means > std_dev * fold) )
  dist <- dist[order(dist$MEANS),]
  # Make maximum 100 tests per list
  steps <- 10^floor( log10( length(dist[,1]) ) )/100
  return( CalculateThreshold(dist, percentage, steps=steps ) )
}

# CalculateThreshold is a helper function for DetectionLimit
CalculateThreshold <-function(dist,percentage,steps) {
  l <- length(dist[,1])
  # defines return value as NA, in case no significant results is found
  threshold <- NA
  for ( i in 1:(l/steps) ) {
    n = (i-1)*steps + 1   # 1=>1, 2=>101(2), 3=>201(3)...
    hi <- dist[(n):l,]
    l.hi <- length(hi[,1])
    perc <- ( length( hi[,1][hi$SIG] )/l.hi ) * 100
    if (perc > percentage) {
      threshold <- dist$MEANS[n]
      break
    }
  }
  if ( (steps == 1) | (n <= steps) ) {
    return( threshold )
  } else {
    return( CalculateThreshold(dist[(n-steps):l,] ,percentage, steps=(steps/10) ) )
  }
}

# RemoveLowest removes the rows with the lowest values
# It assumes an equal distribution between all columns passed, and only uses
# the first column to identify a threshold
RemoveLowest <- function(df,cols,percent) {
  # Calculating the threshold value.
  # 100-"percent" of rows have higher values than min
  l <- nrow(df)
  min <- df[order(df[,cols[1]]),cols[1]][round(percent/100*l)]
  # return all rows in which at least one of the given columns has a higher value
  # than the threshold "min"
  return( df[apply(df[,cols],1,max) > min,] )
}
# Usage: 
# df.high <- RemoveLowest(df,c("Column1","Column2","Column3"),5)

