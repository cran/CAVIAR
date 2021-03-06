plotWoodFormationCalendar <- function(data, main=" ", sub=" ", xlim=NULL, plotype=4) {

 # *********************************************************************************************
 # plotWoodFormationCalendar() function definition
 # -----------------------------------------------
 #
 # Drawing function for plotting wood formation calendar
 # 
 # Arguments:
 # 		- data: data.frame with imposed column names, typically output from computeCriticalDates()
 #			function
  #		- main (optional): plot tittle
 #		- sub (optional): plot subtittle
 #		- xlim (optional): vector indicating x-axis min and max values
 #		- plotype (optional): 	type 1 --> individual critical dates plot
 #								type 2 --> goup critical dates plot
 #								type 3 --> goup critical duration plot
 #								type 4 --> goup critical dates & duration plot
 #
 # Output:
 #		- plot
 #		- dataframe of the median and median absolute deviation
 # 
 # Version 4.0
 # 		1. plot individual critical dates
 #		2. plot group critical dates
 #		3. plot group phases durations
 #		4. Compute median and inter-quartile range
 #			4.1. Modifying presentation of maature phase and missing values (15/12/2011)
 #			4.2. Changing title name argument for main (01/01/2012)
 # Started: 7 Juillet 2009
 # Last modifications: 01 January 2010
 # Author: Cyrille RATHGEBER - INRA Nancy
 #
 # *********************************************************************************************
	
	DF <- data
	DFa  <- data
	DFl <- is.na(DF)	

	# Months characteristics
	MonthId <- c("J", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D")
	FirstDay <- c(1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335)
	MiddleDay <- c(16, 45, 75, 105, 136, 166, 197, 228, 258, 289, 319, 350)
	
	# x-axis range
	if (is.null(xlim) == TRUE) xlim <- c(min(DF$bE, na.rm=TRUE), max(DF$cL, na.rm=TRUE))

	
	# Replacing NA by min or max
	# ---------------------------
	DFa$bE <- ifelse(is.na(DF$bE), 0, DF$bE)
	DFa$bL <- ifelse(is.na(DF$bL), 0, DF$bL)
	DFa$bM <- ifelse(is.na(DF$bM), 0, DF$bM)
	DFa$cE <- ifelse(is.na(DF$cE), 400, DF$cE)
	DFa$cL <- ifelse(is.na(DF$cL), 400, DF$cL)
	
	
	# Associated functions
	# --------------------
	
	# Definition of plotPlayMark function
	plotPlayMarks <- function(x, y, sx=5, dir, col) {
	
		sy  <- 0.5 # parameter for adjusting the size of the mark in y
	
		if (dir == "left")	{
			X <- c(x - sx, x + sx, x + sx)
			Y <- c(y, y + sy, y - sy)
			polygon(X, Y, density=NA, col=col)
		}

		if (dir == "right")	{
			X <- c(x + sx, x - sx, x - sx)
			Y <- c(y, y + sy, y - sy)
			polygon(X, Y, density=NA, col=col)
		}
	} # End plotPlayMarks function
	
	# Definition of plotDateMarks function
	plotDateMarks <- function(F, L, y, c) {
		
		bw  <- 0.25 # parameter for horizontal bar width
		lw <- 1 # parameter for horizontal bar line width
		tw <- 2 # parameter for vertical delimiters width
		
		X <- c(F[2], F[3], F[4], F[3])
		Y <- c(y, y + bw, y, y - bw)
		lines(c(F[1], F[5]), c(y, y), col=c, lwd=tw)
		polygon(X, Y, density=NA, col=c, border= c, lwd=lw)

		X <- c(L[2], L[3], L[4], L[3])
		Y <- c(y, y + bw, y, y - bw)
		lines(c(L[1], L[5]), c(y, y), col=c, lwd=tw)
		polygon(X, Y, density=NA, col=c, border=c, lwd=lw)
	} # End plotDateMarks function
	
	# Definition of plotDurationBars function
	plotDurationsBars <- function(D, y, c) {
	
		bw  <- 0.25 # parameter for horizontal bar width
		lw <- 2 # parameter for horizontal bar line width
	
		X <- c(0, D[3], D[3], 0)
		Y <- c(y + bw, y + bw, y - bw, y - bw)
		polygon(X, Y, density=NA, col=c, border=c, lwd=lw)
		lines(c(D[1], D[3]), c(y, y), col="white", lwd=lw)
		lines(c(D[2], D[2]), c(y - 0.25*bw, y + 0.25*bw), col="white", lwd=lw)
		lines(c(D[3], D[5]), c(y, y), col=c, lwd=lw)
		lines(c(D[4], D[4]), c(y - 0.25*bw, y + 0.25*bw), col=c, lwd=lw)
	} # End of plotDurationBars function
	
	
	# Plotting indicidual critical dates if option 1 selected
	# -------------------------------------------------------
	
	if (plotype==1) {
		
		# Setting the plot region
		ymax <- 25
		ymin <- 1
		YTicks <- c(1, 2, 3, 4, 6, 7, 8, 9, 11, 12, 13, 14, 16, 17, 18, 19, 21, 22, 23, 24)
		YLabs <- rep(c("M", "L", "E", "D"), times=5)

		plot(0, 0, xlim=xlim, ylim=c(ymin, ymax), type="n", ann=FALSE, axes=FALSE)


		# Plotting the phases for the 5 trees
		for (i in 1:5) {
		
			lw <- 3 # parameter for line width
		
			TreePosition <- 5 * (i - 1)
		
			# Dividing Phase (C)
			y <- TreePosition + 4
			if (length(DFa$bD) > 0 && length(DFa$cD) > 0) {
				plotPlayMarks(DFa$bD[i], y, DFa$bD.sd[i], dir="left", col="green")
				plotPlayMarks(DFa$cD[i], y, DFa$cD.sd[i], dir="right", col="green")
				lines(c(DFa$bD[i], DFa$cD[i]), c(y, y), col="green", lwd=lw)
			}

	
			# Expansion Phase (E)
			y <- TreePosition + 3
			plotPlayMarks(DFa$bE[i], y, DFa$bE.sd[i], dir="left", col="blue")
			plotPlayMarks(DFa$cE[i], y, DFa$cE.sd[i], dir="right", col="blue")
			if ( DFl[i, "bE"] == TRUE | DFl[i, "cE"] == TRUE) {
				X <- seq(DFa$bE[i]+10, DFa$cE[i]-10, 10)
				Y <- rep(y, length(X))
				points(X, Y, pch=3, col="blue")
			}
			else {
				lines(c(DFa$bE[i], DFa$cE[i]), c(y, y), col="blue", lwd=lw)
			}
		
			# Maturing Phase (L)
			y <- TreePosition + 2
			plotPlayMarks(DFa$bL[i], y, DFa$bL.sd[i], dir="left", col="red")
			plotPlayMarks(DFa$cL[i], y, DFa$cL.sd[i], dir="right", col="red")
			if ( DFl[i, "bL"] == TRUE | DFl[i, "cL"] == TRUE) {
				X <- seq(DFa$bL[i]+10, DFa$cL[i]-10, 10)
				Y <- rep(y, length(X))
				points(X, Y, pch=3, col="red")
			}
			else {
				lines(c(DFa$bL[i], DFa$cL[i]), c(y, y), col="red", lwd=lw)
			}

			# Mature phase (M)
			y <- TreePosition + 1
			if (DFl[i, "bM"] == TRUE) {
				X <- seq(DFa$bM[i]+10, 400, 10)
				Y <- rep(y, length(X))
				points(X, Y, pch=3, col="brown")
			}
			else {	
				plotPlayMarks(DFa$bM[i], y, DFa$bM.sd[i], dir="left", col="brown")	
			}
			
		}

		# Customising axes 1
		axis(side=1, at=FirstDay, labels=FALSE)
		axis(side=1, at=MiddleDay, tick=FALSE, labels=MonthId, cex.axis=1)
		mtext("Month", side=1, line=2.5, cex=1)

		# Customising axes 2
		axis(side=2, at=YTicks, labels=YLabs, cex.axis=0.75)
		for (i in 1:5) {
			mtext(paste("Tree", DF$Tree[i]), side=2, line=2, at=(2 + 5*(i-1)), cex=1)
		}

		# Drawing separation lines between trees
		abline(h=c(5, 10, 15, 20), lty=3)

		# Drawing separation between month
		for (i in 1:12) {
			abline(v=FirstDay[i], lty=3)
		}

		# Writting additionnal labels on the plot

		mtext(paste(main), side=3, line=2.5, adj=0, cex=1.25)
		mtext(paste(paste(sub)), side=3, line=1.5, adj=0, cex=1)
		mtext(paste("Onset"), side=3, line=0.5, adj=0.05, cex=1)
		mtext(paste("Cessation"), side=3, line=0.5, adj=0.95, cex=1)

		box()
	}
	
	
	# Plotting critical dates if option selected
	# ------------------------------------------
	
	if (plotype==2 | plotype==4) {
	
		# Setting the plot region
		YLabsCal <- c("M", "L", "E")
		YTicks <- c(1, 2, 3)	
		ymin <- 0
		ymax <- 4

		plot(0, 0, xlim=xlim, ylim=c(ymin, ymax), type="n", ann=FALSE, axes=FALSE)

		# Drawing separation between month
		for (i in 1:12) {
			abline(v=FirstDay[i], lty=3)
		}

		# Expension phase
		y <- ymax - 1
		QF <- quantile(DF$bE, probs=seq(0, 1, 0.25), na.rm = TRUE)
		QL <- quantile(DF$cE, probs=seq(0, 1, 0.25), na.rm = TRUE)
		plotDateMarks(QF, QL, y, "blue")

		# Maturing phase
		y <- y - 1
		QF <- quantile(DF$bL, probs=seq(0, 1, 0.25), na.rm = TRUE)
		QL <- quantile(DF$cL, probs=seq(0, 1, 0.25), na.rm = TRUE)
		plotDateMarks(QF, QL, y, "red")

		# Mature phase
		y <- y - 1
		QF <- quantile(DF$bM, probs=seq(0, 1, 0.25), na.rm = TRUE)
		plotDateMarks(QF, rep(NA, times=5), y, "purple")
	
		# Customising axes 1
		axis(side=1, at=FirstDay, labels=FALSE)
		axis(side=1, at=MiddleDay, tick=FALSE, labels=MonthId, cex.axis=1)
		mtext("Month", side=1, line=2.5, cex=1)

		# Customising axes 2
		axis(side=2, at=YTicks, labels=YLabsCal, cex.axis=1)
		mtext("Phases", side=2, line=2, adj=0.5, cex=1)
	
		# Writting additionnal labels on the plot
		mtext(paste(main), side=3, line=2.5, adj=0, cex=1.25)
		mtext(paste(paste(sub)), side=3, line=1, adj=0, cex=1)
		mtext(paste("Onset"), side=3, line=-1.5, adj=0.05, cex=1)
		mtext(paste("Cessation"), side=3, line=-1.5, adj=0.95, cex=1)

		box()
	}
	
	
	# Plotting durations if option selected
	# -------------------------------------
	
	if (plotype==3 | plotype==4) {
				
		# Setting the plot region
		YLabsDur <- c("X", "L", "E")
		YTicks <- c(1, 2, 3)
		ymax <- 4
	

		plot(0, 0, xlim=xlim, ylim=c(0, ymax), type="n", ann=FALSE, axes=FALSE)

		# Drawing vertical lines
		abline(v=c(0, 50, 100, 150, 200, 250, 300, 350), lty=3)

		# Expension phase
		y <- ymax - 1
		Q <- quantile(DF$dE, probs=seq(0, 1, 0.25), na.rm = TRUE)
		plotDurationsBars(Q, y, "blue")

		# Maturing phase
		y <- y - 1
		Q <- quantile(DF$dL, probs=seq(0, 1, 0.25), na.rm = TRUE)
		plotDurationsBars(Q, y, "red")

		# Wood formation period
		y <- y - 1
		Q <- quantile(DF$dX, probs=seq(0, 1, 0.25), na.rm = TRUE)
		plotDurationsBars(Q, y, "brown")
	
		# Customising axes 1
		axis(side=1, cex.axis=1)
		mtext("Number of days", side=1, line=2.5, cex=1)

		# Customising axes 2
		axis(side=2, at=YTicks, labels=YLabsDur, cex.axis=1)
		mtext("Phases", side=2, line=2, adj=0.5, cex=1)

		# Writting additionnal labels
		mtext(main, side=3, line=2.5, adj=0, cex=1.25)
		mtext(sub, side=3, line=1, adj=0, cex=1)

		box()	
	}
	
	# Computing the median and the median absolute deviation
	RDF <- DF[, c("bE", "bL", "bM", "cE", "cL", "dE", "dL", "dX")]
	Median <- apply(RDF, 2, median)
	Q1 <- apply(RDF, 2, quantile, probs=0.25, na.rm = TRUE)
	Q3 <- apply(RDF, 2, quantile, probs=0.75, na.rm = TRUE)
	IQR <- Q3 - Q1
	MAD <- IQR/2
	
	# Writting the results in a list
	ODF <- rbind(Median, Q1)
	ODF <- rbind(ODF, Q3)
	ODF <- rbind(ODF, IQR)
	ODF <- rbind(ODF, MAD)	

	# Return the output list
	return(ODF)

} # End plotWoodCalendar function
# *********************************************************************************************