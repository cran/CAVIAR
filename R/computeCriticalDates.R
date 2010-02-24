# **************************************************************************************
# Author: Cyrille RATHGEBER - LERFoB UMR1092 - INRA Nancy
# Date: 21 Octobre 2009, 4 Novembre 2009, 29 January 2010, 1, 5 & 11 February 2010
# Purpose: Extracting wood formation phenology critical dates
# Version:
#		1. Linear interpolation
#		2. Binary transformation of input data
#		3. Logistic regression
#		4. Confidence interval computation
# **************************************************************************************


computeCriticalDates <- function(data, plot=TRUE) {
 # *************************************************************************************
 # computeCriticalDates() function definition
 # ------------------------------------------
 #
 # Data processing function for computing wood formation critical dates and durations
 #
 # Arguments:
 # 		data: data frame with imposed column names (see exemple)
 #		plot: TRUE (or FALSE) if output plot (pdf) is desired (or not)
 #
 # Output:
 #		dataframe of critical dates and duration with associated 95% CI
 #		optional pdf file with the plots
 #
 # *************************************************************************************

	# General settings for the plots
	# ------------------------------
	pdfName <- paste("Critical dates computation report - v4.0 - ", Sys.Date(), ".pdf", sep="")
	pdf(file=pdfName)

	# Transforming the raw data frame into a binary data frame
	# --------------------------------------------------------
	DF  <- data
	DF$Tree <- as.factor(DF$Tree)

	# Detecting beginning and ending of xylem development phases (E, L, M)
	# --------------------------------------------------------------------
	# bE: beginning of enlargement phase
	# bL: beginning of thickening and lignification phase
	# bM: beginning of mature phase
	# cE: ending of enlargement phase
	# cL: ending of thickening and lignification phase

	T <- numeric()
	bE <- numeric()
	bE.sd <- numeric()
	cE <- numeric()
	cE.sd <- numeric()
	bL <- numeric()
	bL.sd <- numeric()
	cL <- numeric()
	cL.sd <- numeric()
	bM <- numeric()
	bM.sd <- numeric()

	i <- 0

	for(t in levels(DF$Tree)) {

		# Preparing data
		# --------------
		i <- i + 1
		T[i] <- t
		TDF <- DF[DF$Tree==t, ]

		# Binarization of nE, NL & nM
		BDF <- TDF
		BDF$E <- ifelse(BDF$nE > 0.5, 1, 0)
		BDF$L <- ifelse(BDF$nL > 0.5, 1, 0)
		BDF$M <- ifelse(BDF$nM > 0.5, 1, 0)

		# Computing the mean for the tree t
		BDF$DY <- as.factor(BDF$DY)
		MDF <- aggregate(BDF[, c("E", "L", "M")],
			by=list(BDF$Tree, BDF$DY), FUN=mean, na.rm=TRUE)
		names(MDF) <- c("Tree", "DY", "E", "L", "M")

		# Rounding means
		MDF[, c("E", "L", "M")] <- round(MDF[, c("E", "L", "M")], digits=2)

		# Ordering mean data frame
		MDF <- MDF[order(MDF$Tree, MDF$DY), ]
		MDF$DY <- as.numeric(levels(MDF$DY))[MDF$DY]
		BDF$DY <- as.numeric(levels(BDF$DY))[BDF$DY]

		# Enlarging phase
		# ---------------
		# Delimiting the plain active phase for E
		datE <- MDF$DY[MDF$E > 0.99]
		datE1 <- datE[1]
		datE2 <- datE[length(datE)]

		# spliting dataset
		EDF.b <- BDF[BDF$DY <= datE2, c("DY", "E")]
		EDF.c <- BDF[BDF$DY >= datE1, c("DY", "E")]

		# looking for the beginning using a logistic regression
		DY <- EDF.b$DY
		E <- EDF.b$E
		
		model.bE <- glm(E ~ DY, binomial)
		bE[i] <- as.numeric(round(-model.bE$coef[1] / model.bE$coef[2], digits=0))
		bE1 <- as.numeric((log(97.5/2.5) - model.bE$coef[1]) / model.bE$coef[2])
		bE2 <- as.numeric((log(2.5/97.5) - model.bE$coef[1]) / model.bE$coef[2])
		bE.sd[i] <- round(abs(bE1 - bE2)/2, digit=1)


		# plot the results for the beginning if option selected
		if(plot==TRUE) {
			plot(DY, E, type="n", ann=FALSE)
			abline(h=c(0.025, 0.5, 0.975), col="grey", lty=2)
			xp <- c(bE1, bE2, bE2, bE1)
			yp <- c(0.025, 0.025, 0.975, 0.975)
			polygon(xp, yp, density=NA, col="yellow", border="yellow")
			x <- seq(min(DY), max(DY), 1)
			y <- predict(model.bE, list(DY=x), type="response")
			lines(x, y, col="blue")
			points(DY, E, type="p", pch=1, col="blue")
			points(bE[i], 0.5, pch=15, col="red")

			# Writting additionnal labels on the plot
			title(paste("Tree", t, " - bE"))
			mtext("Day of Year", side=1, line=2)
			mtext("Probability", side=2, line=2)
			mtext(paste("bE = ", bE, " +/- ", bE.sd, "days"), side=3, line=0.5, adj=0)
		}

		# looking for the end using a logistic regression
		DY <- EDF.c$DY
		E <- EDF.c$E

		model.cE <- glm(E ~ DY, binomial)
		cE[i] <- as.numeric(round(-model.cE$coef[1] / model.cE$coef[2], digits=0))
		cE1 <- as.numeric((log(97.5/2.5) - model.cE$coef[1]) / model.cE$coef[2])
		cE2 <- as.numeric((log(2.5/97.5) - model.cE$coef[1]) / model.cE$coef[2])
		cE.sd[i] <- round(abs(cE1 - cE2)/2, digit=1)

		# plot the results for the end if option selected
		if(plot==TRUE) {
			plot(DY, E, type="n", ann=FALSE)
			abline(h=c(0.025, 0.5, 0.975), col="grey", lty=2)
			xp <- c(cE1, cE2, cE2, cE1)
			yp <- c(0.025, 0.025, 0.975, 0.975)
			polygon(xp, yp, density=NA, col="yellow", border="yellow")
			x <- seq(min(DY), max(DY), 1)
			y <- predict(model.cE, list(DY=x), type="response")
			lines(x, y, col="blue")
			points(DY, E, type="p", pch=1, col="blue")
			points(cE[i], 0.5, pch=15, col="red")

			# Writting additionnal labels on the plot
			title(paste("Tree", t, " - cE"))
			mtext("Day of Year", side=1, line=2)
			mtext("Radial file ratio", side=2, line=2)
			mtext(paste("cE = ", cE, " +/- ", cE.sd, "days"), side=3, line=0.5, adj=0)
		}


		# Maturing phase
		# --------------
		# Delimiting the plain active phase for L
		datL <- MDF$DY[MDF$L > 0.99]
		datL1 <- datL[1]
		datL2 <- datL[length(datL)]

		# spliting dataset
		LDF.b <- BDF[BDF$DY <= datL2, c("DY", "L")]
		LDF.c <- BDF[BDF$DY >= datL1, c("DY", "L")]

		# looking for the beginning using a logistic regression
		DY <- LDF.b$DY
		L <- LDF.b$L
		
		model.bL <- glm(L ~ DY, binomial)
		bL[i] <- as.numeric(round((log(1) - model.bL$coef[1]) / model.bL$coef[2], digits=0))
		bL1 <- as.numeric((log(97.5/2.5) - model.bL$coef[1]) / model.bL$coef[2])
		bL2 <- as.numeric((log(2.5/97.5) - model.bL$coef[1]) / model.bL$coef[2])
		bL.sd[i] <- round(abs(bL1 - bL2)/2, digit=1)

		# plot the results for the beginning if option selected
		if(plot==TRUE) {
			plot(DY, L, type="n", ann=FALSE)
			abline(h=c(0.025, 0.5, 0.975), col="grey", lty=2)
			xp <- c(bL1, bL2, bL2, bL1)
			yp <- c(0.025, 0.025, 0.975, 0.975)
			polygon(xp, yp, density=NA, col="yellow", border="yellow")
			x <- seq(min(DY), max(DY), 1)
			y <- predict(model.bL, list(DY=x), type="response")
			lines(x, y, col="blue")
			points(DY, L, type="p", pch=1, col="blue")
			points(bL[i], 0.5, pch=15, col="red")

			# Writting additionnal labels on the plot
			title(paste("Tree", t, " - bL"))
			mtext("Day of Year", side=1, line=2)
			mtext("Radial file ratio", side=2, line=2)
			mtext(paste("bL = ", bL, " +/- ", bL.sd, "days"), side=3, line=0.5, adj=0)
		}


		# looking for the end using a logistic regression
		DY <- LDF.c$DY
		L <- LDF.c$L

		model.cL <- glm(L ~ DY, binomial)
		cL[i] <- as.numeric(round((log(1) - model.cL$coef[1]) / model.cL$coef[2], digits=0))
		cL1 <- as.numeric((log(97.5/2.5) - model.cL$coef[1]) / model.cL$coef[2])
		cL2 <- as.numeric((log(2.5/97.5) - model.cL$coef[1]) / model.cL$coef[2])
		cL.sd[i] <- round(abs(cL1 - cL2)/2, digit=1)

		# plot the results for the end if option selected
		if(plot==TRUE) {
			plot(DY, L, type="n", ann=FALSE)
			abline(h=c(0.025, 0.5, 0.975), col="grey", lty=2)
			xp <- c(cL1, cL2, cL2, cL1)
			yp <- c(0.025, 0.025, 0.975, 0.975)
			polygon(xp, yp, density=NA, col="yellow", border="yellow")
			x <- seq(min(DY), max(DY), 1)
			y <- predict(model.cL, list(DY=x), type="response")
			lines(x, y, col="blue")
			points(DY, L, type="p", pch=1, col="blue")
			points(cL[i], 0.5, pch=15, col="red")

			# Writting additionnal labels on the plot
			title(paste("Tree", t, " - cL"))
			mtext("Day of Year", side=1, line=2)
			mtext("Radial file ratio", side=2, line=2)
			mtext(paste("cL = ", cL, " +/- ", cL.sd, "days"), side=3, line=0.5, adj=0)
		}


		# Mature phase
		# ------------
		# spliting dataset
		MDF.b <- BDF[, c("DY", "M")]

		# looking for the beginning using a logistic regression
		DY <- MDF.b$DY
		M <- MDF.b$M

		model.bM <- glm(M ~ DY, binomial)
		bM[i] <- as.numeric(round((log(1) - model.bM$coef[1]) / model.bM$coef[2], digits=0))
		bM1 <- as.numeric((log(97.5/2.5) - model.bM$coef[1]) / model.bM$coef[2])
		bM2 <- as.numeric((log(2.5/97.5) - model.bM$coef[1]) / model.bM$coef[2])
		bM.sd[i] <- round(abs(bM1 - bM2)/2, digit=1)

		# plot the results for the beginning if option selected
		if(plot==TRUE) {
			plot(DY, M, type="n", ann=FALSE)
			abline(h=c(0.025, 0.5, 0.975), col="grey", lty=2)
			xp <- c(bM1, bM2, bM2, bM1)
			yp <- c(0.025, 0.025, 0.975, 0.975)
			polygon(xp, yp, density=NA, col="yellow", border="yellow")
			x <- seq(min(DY), max(DY), 1)
			y <- predict(model.bM, list(DY=x), type="response")
			lines(x, y, col="blue")
			points(DY, M, type="p", pch=1, col="blue")
			points(bM[i], 0.5, pch=15, col="red")

			# Writting additionnal labels on the plot
			title(paste("Tree", t, " - bM"))
			mtext("Day of Year", side=1, line=2)
			mtext("Radial file ratio", side=2, line=2)
			mtext(paste("bM = ", bM, " +/- ", bM.sd, "days"), side=3, line=0.5, adj=0)
		}
	}

	dev.off()


	# Computing durations of critical periods and associated standard deviation
	# -------------------------------------------------------------------------
	# dE: total duration of enlarging phase
	# dL: total duration of maturing phase
	# dX: total duration of xylem differentiation

	ODF <- data.frame(T, bE, bE.sd, bL, bL.sd, bM, bM.sd, cE, cE.sd, cL, cL.sd)
	names(ODF) <- c("Tree", "bE", "bE.sd", "bL", "bL.sd", "bM", "bM.sd", "cE", "cE.sd", "cL", "cL.sd")
	ODF$dE <- ODF$cE - ODF$bE
	ODF$dE.sd <- round(sqrt(ODF$cE.sd^2 + ODF$bE.sd^2), digit=1)
	ODF$dL <- ODF$cL - ODF$bL
	ODF$dL.sd <- round(sqrt(ODF$cL.sd^2 + ODF$bL.sd^2), digit=1)
	ODF$dX <- ODF$cL - ODF$bE
	ODF$dX.sd <- round(sqrt(ODF$cL.sd^2 + ODF$bE.sd^2), digit=1)

	# Returning the results
	# ---------------------
	return(ODF)

} # end function computeCriticalDates
