# **************************************************************************************
# FUNCTION computeCriticalDates
# -----------------------------
#
# Author: Cyrille RATHGEBER - LERFoB UMR1092 - INRA Nancy
# Purpose: Computing wood formation phenology critical dates
# Version:
#		1. Linear interpolation
#		2. Binary transformation of input data
#		3. Logistic regression
#		4. Confidence interval computation
#			4.1. Values in plain active phase forced to 1
#				4.1-1. Bug concerning plot subtitle corrected
#				4.1-2. Replacing aberrant value by NA
#				4.1-3. Bug concerning tree number handling corrected (12/12/2011)
#				4.1-4. Missing data handling improvement (16, 19, 21/12/2011)
#
# Development started: 21 October 2009
# Last modifications: 21 December 2012
#
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
 #
 # *************************************************************************************

	# General settings for the plots
	# ------------------------------
	pdfName <- paste("Critical dates computation report - ", Sys.Date(), ".pdf", sep="")
	pdf(file=pdfName)

	# Transforming the raw data frame into a binary data frame
	# --------------------------------------------------------
	DF  <- data
	DF$Tree <- factor(DF$Tree, exclude=NULL)

	# Detecting beginning and ending of xylem development phases (E, L, M)
	# --------------------------------------------------------------------
	# bE: beginning of enlargement phase
	# bL: beginning of thickening and lignification phase
	# bM: beginning of mature phase
	# cE: ending of enlargement phase
	# cL: ending of thickening and lignification phase

	Tree <- character()
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
		Tree[i] <- t
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
		# Checking if ther is any values for E
		if (length(na.omit(MDF$E)) == 0){
			bE[i] <-  NA
			bE.sd[i] <- NA
			cE[i] <-  NA
			cE.sd[i] <- NA	
			
		}
		else {
			# Delimiting the plain active phase for E		
			datE <- MDF$DY[MDF$E > 0.99]
			datE.net <- na.omit(datE)		
			datE1 <- datE.net[1]
			datE2 <- datE.net[length(datE.net)]
		
			# spliting dataset
			EDF.b <- BDF[BDF$DY <= datE2, c("DY", "E")]
			EDF.c <- BDF[BDF$DY >= datE1, c("DY", "E")]
		
			# Forcing plain active phase values to 1
			# --> in order to ensure that begining and end are correctly computed!
			EDF.b[EDF.b$DY >= datE1, "E"] <- 1
			EDF.c[EDF.c$DY <= datE2, "E"] <- 1

			# looking for the beginning using a logistic regression
			DY <- EDF.b$DY
			E <- EDF.b$E
			EDF.bm <- aggregate(EDF.b[, c("E")], by=list(EDF.b$DY), FUN=mean, na.rm=TRUE)
			names(EDF.bm) <- c("DY", "E.m")
			MV <- EDF.bm[is.na(EDF.bm$E.m)==TRUE, ]
			if (nrow(MV) >=1) MV$E <- 0.5
		
			model.bE <- glm(E ~ DY, binomial)
			bE[i] <- as.numeric(round(-model.bE$coef[1] / model.bE$coef[2], digits=0))
			
			# Replacing NA by (O, 1) couple to compute SD more realistacally, if needed
			if (bE[i] >= min(MV$DY) & bE[i] <= max(MV$DY)) {
				DY <- c(DY, rep(MV$DY, 2))
				E <- c(E, rep(0, length(MV$DY)), rep(1, length(MV$DY)))
				model.bE <- glm(E ~ DY, binomial)
				bE[i] <- as.numeric(round(-model.bE$coef[1] / model.bE$coef[2], digits=0))
			}
			
			# Estimating the standard error
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
				points(EDF.b$DY, EDF.b$E, type="p", pch=1, col="blue")	
				points(bE[i], 0.5, pch=16, col="red")
				points(MV$DY, MV$E, type="p", pch = 4, col="black")

				# Writting additionnal labels on the plot
				title(paste("Tree", t, " - bE"))
				mtext(paste("bE = ", bE[i], " +/- ", bE.sd[i], "days"), side=3, line=0.5, adj=0)
				mtext("Day of Year", side=1, line=2)
				mtext("Probability", side=2, line=2)
			} # end if

			# looking for the end using a logistic regression
			DY <- EDF.c$DY
			E <- EDF.c$E
			EDF.cm <- aggregate(EDF.c[, c("E")], by=list(EDF.c$DY), FUN=mean, na.rm=TRUE)
			names(EDF.cm) <- c("DY", "E.m")
			MV <- EDF.cm[is.na(EDF.cm$E.m)==TRUE, ]
			if (nrow(MV) >=1) MV$E <- 0.5

			model.cE <- glm(E ~ DY, binomial)
			cE[i] <- as.numeric(round(-model.cE$coef[1] / model.cE$coef[2], digits=0))
			
			# Replacing NA by (O, 1) couple to compute SD more realistacally, if needed
			if (cE[i] >= min(MV$DY) & cE[i] <= max(MV$DY)) {
				DY <- c(DY, rep(MV$DY, 2))
				E <- c(E, rep(0,length(MV$DY)), rep(1,length(MV$DY)))
				model.cE <- glm(E ~ DY, binomial)
				bE[i] <- as.numeric(round(-model.cE$coef[1] / model.cE$coef[2], digits=0))
			}
			
			# Estimating the standard error			
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
				points(EDF.c$DY, EDF.c$E, type="p", pch=1, col="blue")
				points(cE[i], 0.5, pch=15, col="red")
				points(MV$DY, MV$E, type="p", pch = 4, col="black")

				# Writting additionnal labels on the plot
				title(paste("Tree", t, " - cE"))
				mtext(paste("cE = ", cE[i], " +/- ", cE.sd[i], "days"), side=3, line=0.5, adj=0)
				mtext("Day of Year", side=1, line=2)
				mtext("Radial file ratio", side=2, line=2)
			} # end if
		} # end else 

		# Maturing phase
		# --------------
		# Checking if ther is any values for L
		if (length(na.omit(MDF$L)) == 0){
			bL[i] <-  NA
			bL.sd[i] <- NA
			cL[i] <-  NA
			cL.sd[i] <- NA
		}
		else {
			# Delimiting the plain active phase for L
			datL <- MDF$DY[MDF$L > 0.99]
			datL.net <- na.omit(datL)
			datL1 <- datL.net[1]
			datL1 <- datL[1]
			datL2 <- datL.net[length(datL.net)]

			# spliting dataset
			LDF.b <- BDF[BDF$DY <= datL2, c("DY", "L")]
			LDF.c <- BDF[BDF$DY >= datL1, c("DY", "L")]
		
			# Forcing plain active phase values to 1
			# --> in order to ensure that begining and end are correctly computed!
			LDF.b[LDF.b$DY >= datL1, "L"] <- 1
			LDF.c[LDF.c$DY <= datL2, "L"] <- 1

			# looking for the beginning using a logistic regression
			DY <- LDF.b$DY
			L <- LDF.b$L
			LDF.bm <- aggregate(LDF.b[, c("L")], by=list(LDF.b$DY), FUN=mean, na.rm=TRUE)
			names(LDF.bm) <- c("DY", "L.m")
			MV <- LDF.bm[is.na(LDF.bm$L)==TRUE, ]
			if (nrow(MV) >=1) MV$L <- 0.5
		
			model.bL <- glm(L ~ DY, binomial)
			bL[i] <- as.numeric(round((log(1) - model.bL$coef[1]) / model.bL$coef[2], digits=0))
			
			# Replacing NA by (O, 1) couple to compute SD more realistacally, if needed
			if (bL[i] >= min(MV$DY) & bL[i] <= max(MV$DY)) {
				DY <- c(DY, rep(MV$DY, 2))
				L <- c(L, rep(0, length(MV$DY)), rep(1, length(MV$DY)))
				model.bL <- glm(L ~ DY, binomial)
				bL[i] <- as.numeric(round(-model.bL$coef[1] / model.bL$coef[2], digits=0))
			}
			
			# Estimating the standard error
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
				points(MV$DY, MV$L, type="p", pch = 4, col="black")

				# Writting additionnal labels on the plot
				title(paste("Tree", t, " - bL"))
				mtext(paste("bL = ", bL[i], " +/- ", bL.sd[i], "days"), side=3, line=0.5, adj=0)
				mtext("Day of Year", side=1, line=2)
				mtext("Radial file ratio", side=2, line=2)
			} # end if


			# looking for the end using a logistic regression
			DY <- LDF.c$DY
			L <- LDF.c$L
			LDF.cm <- aggregate(LDF.c[, c("L")], by=list(LDF.c$DY), FUN=mean, na.rm=TRUE)
			names(LDF.cm) <- c("DY", "L.m")
			MV <- LDF.cm[is.na(LDF.cm$L)==TRUE, ]
			if (nrow(MV) >=1) MV$L <- 0.5

			model.cL <- glm(L ~ DY, binomial)
			cL[i] <- as.numeric(round((log(1) - model.cL$coef[1]) / model.cL$coef[2], digits=0))
			
			# Replacing NA by (O, 1) couple to compute SD more realistacally, if needed
			if (cL[i] >= min(MV$DY) & cL[i] <= max(MV$DY)) {
				DY <- c(DY, rep(MV$DY, 2))
				L <- c(L, rep(0, length(MV$DY)), rep(1, length(MV$DY)))
				model.cL <- glm(L ~ DY, binomial)
				cL[i] <- as.numeric(round(-model.cL$coef[1] / model.cL$coef[2], digits=0))
			}
			
			# Estimating the standard error
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
				points(MV$DY, MV$L, type="p", pch = 4, col="black")

				# Writting additionnal labels on the plot
				title(paste("Tree", t, " - cL"))
				mtext(paste("cL = ", cL[i], " +/- ", cL.sd[i], "days"), side=3, line=0.5, adj=0)
				mtext("Day of Year", side=1, line=2)
			mtext("Radial file ratio", side=2, line=2)
			}# end if
		} # end else


		# Mature phase
		# ------------
		# Checking if ther is any values for M
		if (length(na.omit(MDF$M)) == 0){
			bM[i] <-  NA
			bM.sd[i] <- NA
		}
		
		else {
			# spliting dataset
			MDF.b <- BDF[, c("DY", "M")]

			# looking for the beginning using a logistic regression
			DY <- MDF.b$DY
			M <- MDF.b$M
			MDF.bm <- aggregate(MDF.b[, c("M")], by=list(MDF.b$DY), FUN=mean, na.rm=TRUE)
			names(MDF.bm) <- c("DY", "M.m")
			MV <- MDF.bm[is.na(MDF.bm$M)==TRUE, ]
			if (nrow(MV) >=1) MV$M <- 0.5

			model.bM <- glm(M ~ DY, binomial)
			bM[i] <- as.numeric(round((log(1) - model.bM$coef[1]) / model.bM$coef[2], digits=0))
			
			# Replacing NA by (O, 1) couple to compute SD more realistacally, if needed
			if (bM[i] >= min(MV$DY) & bM[i] <= max(MV$DY)) {
				DY <- c(DY, rep(MV$DY, 2))
				M <- c(M, rep(0, length(MV$DY)), rep(1, length(MV$DY)))
				model.bM <- glm(M ~ DY, binomial)
				bM[i] <- as.numeric(round(-model.bM$coef[1] / model.bM$coef[2], digits=0))
			}
			
			# Estimating the standard error
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
				points(MV$DY, MV$M, type="p", pch = 4, col="black")

				# Writting additionnal labels on the plot
				title(paste("Tree", t, " - bM"))
				mtext(paste("bM = ", bM[i], " +/- ", bM.sd[i], "days"), side=3, line=0.5, adj=0)
				mtext("Day of Year", side=1, line=2)
				mtext("Radial file ratio", side=2, line=2)
			} # end if
		} # end else	
	}

	dev.off()


	# Computing durations of critical periods and associated standard deviation
	# -------------------------------------------------------------------------
	# dE: total duration of enlarging phase
	# dL: total duration of maturing phase
	# dX: total duration of xylem differentiation

	RDF <- data.frame(bE, bE.sd, bL, bL.sd, bM, bM.sd, cE, cE.sd, cL, cL.sd)
	names(RDF) <- c("bE", "bE.sd", "bL", "bL.sd", "bM", "bM.sd", "cE", "cE.sd", "cL", "cL.sd")
	RDF$dE <- RDF$cE - RDF$bE
	RDF$dE.sd <- round(sqrt(RDF$cE.sd^2 + RDF$bE.sd^2), digit=1)
	RDF$dL <- RDF$cL - RDF$bL
	RDF$dL.sd <- round(sqrt(RDF$cL.sd^2 + RDF$bL.sd^2), digit=1)
	RDF$dX <- RDF$cL - RDF$bE
	RDF$dX.sd <- round(sqrt(RDF$cL.sd^2 + RDF$bE.sd^2), digit=1)

	# Removing aberrant values
	# ------------------------
	abr.rm <- function(x) {
		ifelse(x >= 0 & x <= 365, x, NA)
	}
	RDF.net <- as.data.frame(apply(RDF, 2, abr.rm))
	
	# Returning the results
	# ---------------------
	ODF <-  data.frame(Tree, RDF.net)

	return(ODF)

} # end function computeCriticalDates
# **************************************************************************************