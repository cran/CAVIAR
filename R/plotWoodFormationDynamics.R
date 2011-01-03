
plotWoodFormationDynamics <- function(data, title=" ", x.axis.labels="DOY", wrap=FALSE) {

 # *******************************************************************************************************
 # plotWoodFormationDynamics() function definition
 # -----------------------------------------------
 #
 # Drawing function for plotting intra-annual dynamics of wood formation
 # 
 # Arguments:
 # 		- data: data.frame with imposed column names, typically output from averageRF() function
 #		- title (optional): plot tittle
 #		- x.axis.labels: a character for selecting between day of year (DOY, default) or Sample number (Sample)
 #		- wrap (optional): draws an enveloppe around the main curve if TRUE
 #
 # Output:
 #		- plot
 #
 # Version: 1.1-1
 # Started: 10 Jully 2009
 # Last modifications: 4 November 2010
 # Author: Cyrille RATHGEBER - INRA Nancy
 #
 # *******************************************************************************************************
	
	
	DF <- data
	
	# Months characteristic
	MonthId <- c("J", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D")
	FirstDay <- c(1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335)
	MiddleDay <- c(16, 45, 75, 105, 136, 166, 197, 228, 258, 289, 319, 350)
	
	# Setting the plot region
	ymax <- max(DF$nM, na.rm=TRUE)
	plot(DF$DY, DF$nM, type="n", xlim=c(min(DF$DY), max(DF$DY)), ylim=c(0, ymax), ann=FALSE, axes=FALSE)
	
	# Drawing enveloppes around the main curves
	if (wrap == TRUE) {
		X <- c(DF$DY, DF$DY[rev(order(DF$DY))])
		
		YC <- c(DF$nC.inf, DF$nC.sup[rev(order(DF$DY))])		
		polygon(X, YC, density=NA, col=rgb(0, 1, 0, 0.25))
		
		YE <- c(DF$nE.inf, DF$nE.sup[rev(order(DF$DY))])
		polygon(X, YE, density=NA, col=rgb(0, 0, 1, 0.25))
	
		YL <- c(DF$nL.inf, DF$nL.sup[rev(order(DF$DY))])
		polygon(X, YL, density=NA, col=rgb(1, 0, 0, 0.25))
		
		YM <- c(DF$nM.inf, DF$nM.sup[rev(order(DF$DY))])
		polygon(X, YM, density=NA, col=rgb(0.37, 0.07, 0.55, 0.25))
	}
	
	# Drawing main curves
	points(DF$DY, DF$nC, type="o", pch=17, col=rgb(0, 1, 0, 1))
	points(DF$DY, DF$nE, type="o", pch=15, col=rgb(0, 0, 1, 1))
	points(DF$DY, DF$nL, type="o", pch=18, col=rgb(1, 0, 0, 1))
	points(DF$DY, DF$nM, type="o", pch=19, col=rgb(0.37, 0.07, 0.55, 1))
	
	# Customising axes 1 in function of the require option
	if(x.axis.labels=="Sample") {
		axis(1, at=DF$DY, labels=DF$Sample)
		mtext("Sample number", side=1, line=2.5)
	}
	else {
		axis(3, at=FirstDay, labels=FALSE)
		axis(3, at=MiddleDay, tick=FALSE, labels=MonthId)
		axis(1)
		mtext("Day of year", side=1, line=2.5)
	}
	
	# Customising axes 2
	axis(2)
	mtext("Number of cells", side=2, line=2.5)

	# Writting plot title
	mtext(title, side=3, line=2, adj=0.0, cex=1.25)
	
	# Legend
	xmin <- min(DF$DY)
	r <- ymax / 15
	points(xmin, ymax, pch=17, col="green")
	text(xmin + 1, ymax, "Dividing", pos=4)
	points(xmin, ymax - r, pch=15, col="blue")
	text(xmin + 1, ymax - r, "Enlarging", pos=4)
	points(xmin, ymax - 2*r, pch=18, col="red")
	text(xmin + 1, ymax - 2*r, "Thickening", pos=4)
	points(xmin, ymax - 3*r, pch=19, col="purple")
	text(xmin + 1, ymax - 3*r, "Mature", pos=4)
	
	box()
	
} # End plotWoodFormationDynamics function