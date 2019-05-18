library('hexbin')
FacsPlot <- function(
	black=NULL,
	gray=NULL,
	red=NULL,
	green=NULL,
	blue=NULL,
	orange=NULL,
	log="xy",
	xlim=NULL,
	ylim=NULL,
	cex=0.3,
	pch=19,
	saturation=0.4,
	main=""
	) {

	if (!is.null(black)) {
		names(black) <- rep(NA,ncol(black))
	}
	if (!is.null(gray)) {
		names(gray) <- rep(NA,ncol(gray))
	}
	if (!is.null(red)) {
		names(red) <- rep(NA,ncol(red))
	}
	if (!is.null(green)) {
		names(green) <- rep(NA,ncol(green))
	}
	if (!is.null(blue)) {
		names(blue) <- rep(NA,ncol(blue))
	}
	if (!is.null(orange)) {
		names(orange) <- rep(NA,ncol(orange))
	}

	names <- c("black","gray","red","green","blue","orange")
	cols <- data.frame(
		c(0.1,0.1,0.1),
		c(0.4,0.4,0.4),
		c(0.8,0.2,0.2),
		c(0.2,0.8,0.2),
		c(0.2,0.2,0.8),
		c(0.9,0.4,0.1)
		)
	all <- rbind(black,gray,red,green,blue,orange)
	
	if (is.null(xlim)) {
		xlim <- c(min(all[,1]), ceiling(max(all[,1])))
		if (grepl("x",log)) {
			xlim[1] <- xlim[1]/10
		} else {
			xlim[1] <- floor(xlim[1])
		}
	}
	if (is.null(ylim)) {
		ylim <- c(min(all[,2]), ceiling(max(all[,2])))
		if (grepl("y",log)) {
			ylim[1] <- ylim[1]/10
		} else {
			ylim[1] <- floor(ylim[1])
		}
	}
	
	plot(NA,NA,xlim=xlim,ylim=ylim,log=log,main=main)

	if (length(saturation) == 1) {
		saturation <- rep(saturation,6)
	}
	if (length(cex) == 1) {
		cex <- rep(cex,6)
	}

	for (i in 1:6) {
		points(
			eval(as.name(names[i]))[,1],
			eval(as.name(names[i]))[,2],
			pch=pch,cex=cex[i],
			col=(rgb(cols[1,i],cols[2,i],cols[3,i],saturation[i]))
		)
	}
}

Compensate <- function(v1,v2,percent) {
	return(v1-(v2*(percent/100)))
}

DefineGates <- function(df,xcol,ycol,a,b,c,d) {
	# rect:
	# x > [1]
	# y > [2]
	# x <= [3]
	# y <= [4] 
	df[,"gatecol"] <- 0
	df[ (
		df[,xcol] > a[1] & 
		df[,ycol] > a[2] & 
		df[,xcol] <= a[3] & 
		df[,ycol] <= a[4]
		)	,"gatecol"] <- 1
	df[ (
		df[,xcol] > b[1] & 
		df[,ycol] > b[2] & 
		df[,xcol] <= b[3] & 
		df[,ycol] <= b[4]
		)	,"gatecol"] <- 2
	df[ (
		df[,xcol] > c[1] & 
		df[,ycol] > c[2] & 
		df[,xcol] <= c[3] & 
		df[,ycol] <= c[4]
		)	,"gatecol"] <- 3
	df[ (
		df[,xcol] > d[1] & 
		df[,ycol] > d[2] & 
		df[,xcol] <= d[3] & 
		df[,ycol] <= d[4]
		)	,"gatecol"] <- 4

	return(df[,"gatecol"])
}



#	df.n$brd2gates <- DefineGates(df.n,
#	xcol=cols.brd2.c[1],
#	ycol=cols.pol2.c[1],
#	c(1.4,1.5,8,8),
#	c(3,8,8,16),
#	c(5,16,500,500),
#	c(8,1.5,500,16)
#	)

# Cumulative returns a dataframe that can be used to plot Cumulative Plots
# It receives a dataframe and the name of the column to be plotted
Cumulative <- function(df,x_axis) {
  t   <- data.frame(x = sort(df[,x_axis]))
  t$y <- row(t)/length(t[,1])
  return(t)
}

# /**
#  * @param  {data.frame} df data frame that contains at least one column with pValues and one with fold difference (log2)
#  * @param  {string} pValue name of the column that contains the pValues
#  * @param  {string} diff name of the column that contains the fold difference values
#  * @param  {string} highlightRegex optional string to do regex search and highlight certain datapoints on the plot
#  * @param  {Boolean} labels optional Flag wether to write the names of the highlighted datapoints next to it 
# 										(in this case, df must contain one column labeled "symbol")
#  * @return {list}
# */
VolcanoPlot <- function(df,pValue,diff,highlightRegex=FALSE,labels=FALSE) {
	# data frame with all genes that don't significantly change
	df.no <- df[
		(
			df[,diff] < 1 &
			df[,diff] > (-1)
		) |
		df[,pValue] > 0.05
	,]

	# Increasing gene
	df.up <- df[df[,diff] >= 1 & df[,pValue] <= 0.05 ,]
	df.dn <- df[df[,diff] <= (-1) & df[,pValue] <= 0.05 ,]

	if (nrow(df.no)+nrow(df.up)+nrow(df.dn) != nrow(df)) {
		stop("Numbers of unchanged, increased and decreased don't add up")
	}

	# ensure that X axis is centered
	# by setting the max and min based on 
	# significantly up and downregulated genes only
	# allow unsignificantly changed genes to not be included
	maxValue <- max(abs(min(df.up[,diff],df.dn[,diff])),max(df.up[,diff],df.dn[,diff]))

	plot(NA,NA,xlim=c(maxValue*-1,maxValue),ylim=c(1,min(df.up[,pValue],df.dn[,pValue])),
		pch=19,cex=0.3,log="y",axes=F) 
	points(df.no[,diff],df.no[,pValue],pch=19,cex=0.3,col=rgb(0,0,0,0.3))
	points(df.up[,diff],df.up[,pValue],pch=19,cex=0.4,col=rgb(1,0,0,0.6))
	points(df.dn[,diff],df.dn[,pValue],pch=19,cex=0.4,col=rgb(0,1,0,0.6))
	axis(2)
	axis(2,at=c(1,0.05))
	axis(1)
	box()

	if (is.character(highlightRegex)) {
		# Add a second plot with highlighted points
		plot(NA,NA,xlim=c(maxValue*-1,maxValue),ylim=c(1,min(df.up[,pValue],df.dn[,pValue])),
			pch=19,cex=0.3,log="y",axes=F) 
		points(df.no[,diff],df.no[,pValue],pch=19,cex=0.3,col=rgb(0,0,0,0.3))
		points(df.up[,diff],df.up[,pValue],pch=19,cex=0.4,col=rgb(1,0,0,0.6))
		points(df.dn[,diff],df.dn[,pValue],pch=19,cex=0.4,col=rgb(0,1,0,0.6))

		df.hl.up <- df.up[grepl(highlightRegex,df.up$symbol),]
		df.hl.dn <- df.dn[grepl(highlightRegex,df.dn$symbol),]
		
		points(df.hl.up[,diff],df.hl.up[,pValue],pch=19,cex=0.7,col=rgb(1,0,0,1))
		points(df.hl.dn[,diff],df.hl.dn[,pValue],pch=19,cex=0.7,col=rgb(0,1,0,1))
		if (labels) {
			text(df.hl.up[,diff],df.hl.up[,pValue],df.hl.up$symbol)	
		}

		axis(2)
		axis(2,at=c(1,0.05))
		axis(1)
		box()

	}
	# Returning a list with the three dataframes without printing
	temp <- list(df.up,df.no,df.dn)
}

# This function does NOT convert data into sliding window
# Use slidingBins() from >>analysis_function.R<< to convert data
PlotSlidingData <- function(data,filename="temp",main,xlab,ylab) {
	pdf(paste(filename,".pdf",sep="")) 
		plot(NA, xlab=xlab, ylab=ylab,main=main, xlim=c(min(data$x),max(data$x)),ylim=c(min(data$y_conf_min),max(data$y_conf_max)))
		polygon(
			c(data$x,rev(data$x)),
			c(data$y_conf_min,rev(data$y_conf_max)),
			col=rgb(0.3,0.3,0.3,0.4),border=NA)
		lines(y~x,data)
		plot(NA, xlab=xlab, ylab=ylab,main=main,xlim=c(min(data$x),max(data$x)),ylim=c(min(data$y_min),max(data$y_max)))
		polygon(
			c(data$x,rev(data$x)),
			c(data$y_max,rev(data$y_min)),
			col=rgb(0.9,0.9,0.9,1))
		polygon(
			c(data$x,rev(data$x)),
			c(data$y_conf_min,rev(data$y_conf_max)),
			col=rgb(0.3,0.3,0.3,0.4),border=NA)
		lines(y~x,data)
	dev.off()
}

# /**
#  * [description]
#  * @param  {data.frame}  	df.plot [containg the data to plot. $x $y mandatory]
#  * @param  {list}    			breaks  [Fraction of datapoints for each color, cumulative]
#  * @param  {Number}    		bins    [Detail of binning function]
#  * @param  {String}    		main    [Title of plot]
#  * @param  {...}			 		...     [arguments passed to plot()]
#  * @return {[type]}            [description]
#  */
FlowJoPlot <- function(
  df.plot,
  breaks=c(0.5,0.7,0.9,1),
  bins=600,
  main="Jonas' FlowJoPlot",
  ...
  ) {
  colors <- colorRampPalette(c("darkblue","green","yellow", "red"))(length(breaks))
  bin<-hexbin(df.plot$x, df.plot$y, xbins=bins) 
  actual_breaks <- quantile(bin@count,breaks,names=FALSE)
  
  mapColors <- function(n){
    for (i in 1:length(actual_breaks)){
      if (n <= actual_breaks[i]){
        return(colors[i])
      }
    }
  }

  plot(
    bin@xcm,
    bin@ycm,
    pch=".",
    col=sapply(bin@count,mapColors),
    main=main,
    ...
  )
  return(list(
    breaks=actual_breaks,
    points=bin@ncells
  ))
}


# 
# /**
#  * Prints a profile plot for the TSS
#  * @param  {[data.frame/string]} filename  data.frame or path/to/TSS.profile/file
#  * @param  {[vector]} colors    color gradient for the plot
#  * @param  {[vector]} colValues values binding the color gradient
#  * @param  {[function]} transform function applied to value for coloring
#  * @return {[data.frame]}           profile data
#  */
TSSProfilePlot <- function(
  filename,
  colors,
  colValues,
  transform=function(x){return(x)},
  legend.pos="right"
  ) {

  if (class(filename) == "character") {
    df <- read.table(filename,
      header=T, sep="\t",  as.is=T)
  } else if (class(df) == "data.frame") {
    df <- filename
  } else {
    stop("Invalid data.frame or name of file")
  }
  if ("sort" %in% names(df)) {
    profileColumns <- names(df)[6:(ncol(df)-1)]
  } else {
    profileColumns <- names(df)[6:ncol(df)]
    df$sort <- rowMeans(df[,profileColumns])
  }


  df.x <- melt(df[,c("symbol",profileColumns)],id="symbol")

  (p <- ggplot(df.x,aes(variable,symbol)) +
    geom_raster(aes(fill=transform(value))) +
    # geom_raster(aes(fill=value)) +
    scale_y_discrete( limits=(df[order(df$sort),"symbol"]) ) +
    scale_fill_gradientn("% max",
      colours = colors,
      values = colValues) +
    theme(axis.title.y=element_blank(),
      axis.text.y=element_blank(),
      axis.ticks.y=element_blank(),
      axis.title.x=element_blank(),
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank(),
      legend.position=legend.pos
    )
  )
  print(p)
  return(df)
}
