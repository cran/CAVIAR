# **************************************************************************************
# FUNCTION computeRingCellNb
# ---------------------------
#
# Author: Cyrille RATHGEBER - LERFoB UMR1092 - INRA Nancy
# Purpose: Computing characteristic ring cell numbers
# Versions:
#		1.0. Compute RCN, ICN, FCN
#		1.1. Bug on RCN computation corrected
#		1.2. adding se and mae
#
# Started: 25 March 2010
# Last modified: 3 December 2010
#
# **************************************************************************************

computeRingCellNb <- function(data, stat="median") {
 # *************************************************************************************
 # computeRingCellNb() function definition
 # ---------------------------------------
 #
 # Data processing function for computing final total ring cell number (RCN)
 #
 # Arguments:
 # 		data: data frame with imposed column names (see exemple)
 #		stat: "mean" or "median" (default)
 #
 # Output:
 #		dataframe of RCN, ICN & FCN (median with median absolute deviation or mean with sd)
 #
 #
 # *************************************************************************************

	# Manipulating the raw data frame
	# -------------------------------
	DF  <- data
	DF$Tree <- as.factor(DF$Tree)
	
	# data frame declaration
	RCN <- data.frame()
	ICN <- data.frame()
	FCN <- data.frame()
	
	# Selecting the data needed for the calculation of RCN, ICN and FCN
	# -----------------------------------------------------------------
	# RCN
	RCN1  <- DF[is.na(DF$nE) == FALSE & DF$nE < 0.1
		& is.na(DF$nM) == FALSE & DF$nM >= 1, c("Tree", "DY", "RF", "nL", "nM")]
		
	RCN1$nT <- RCN1$nL + RCN1$nM
	RCN1$c <- 1
	
	# ICN
	ICN1  <- DF[is.na(DF$nC) == FALSE & is.na(DF$nE) == FALSE & DF$nE < 0.1
		& is.na(DF$nM) == FALSE & DF$nM < 0.1, c("Tree", "DY", "RF", "nC")]
	ICN1$c <- 1
	
	# FCN
	FCN1  <- DF[is.na(DF$nC) == FALSE & is.na(DF$nE) == FALSE & DF$nE < 0.1
		& is.na(DF$nM) == FALSE & DF$nM > 1, c("Tree", "DY", "RF", "nC")]
	FCN1$c <- 1	
	
	
	# Option "mean"
	# -------------
	if (stat=="mean") {
		
		# Computing the mean tree-ring cell number (RCN) and its variation for each tree
		# ------------------------------------------------------------------------------
		RCN2 <- aggregate(RCN1[, "nT"], by=list(Tree= RCN1$Tree), mean, na.rm=TRUE)
		names(RCN2) <- c("Tree", "RCN.mean")
		RCN2$RCN.mean <- round(RCN2$RCN.mean, digits=1)
	
		RCN3 <- aggregate(RCN1[, "nT"], by=list(Tree= RCN1$Tree), sd, na.rm=TRUE)
		names(RCN3) <- c("Tree", "RCN.sd")
		
		RCN4 <- aggregate(RCN1[, "c"], by=list(Tree= RCN1$Tree), sum, na.rm=TRUE)
		names(RCN4) <- c("Tree", "RCN.nb")
	
		RCNa <- merge(RCN2, RCN3, by="Tree")
		RCN <- merge(RCNa, RCN4, by="Tree")
		
		RCN$RCN.se <- RCN$RCN.sd /sqrt(RCN$RCN.nb)
		
		RCN$RCN.sd <- round(RCN$RCN.sd, digits=1)
		RCN$RCN.se <- round(RCN$RCN.se, digits=1)
	
		# Computing the mean initial cambial cell number (ICN) and its variation for each tree
		# ------------------------------------------------------------------------------------
		ICN2 <- aggregate(ICN1[, "nC"], by=list(Tree= ICN1$Tree), mean, na.rm=TRUE)
		names(ICN2) <- c("Tree", "ICN.mean")
		ICN2$ICN.mean <- round(ICN2$ICN.mean, digits=1)
	
		ICN3 <- aggregate(ICN1[, "nC"], by=list(Tree= ICN1$Tree), sd, na.rm=TRUE)
		names(ICN3) <- c("Tree", "ICN.sd")
		
		ICN4 <- aggregate(ICN1[, "c"], by=list(Tree= ICN1$Tree), sum, na.rm=TRUE)
		names(ICN4) <- c("Tree", "ICN.nb")
		
		ICNa <- merge(ICN2, ICN3, by="Tree")
		ICN <- merge(ICNa, ICN4, by="Tree")
		
		ICN$ICN.se <- ICN$ICN.sd /sqrt(ICN$ICN.nb)
		
		ICN$ICN.sd <- round(ICN$ICN.sd, digits=1)
		ICN$ICN.se <- round(ICN$ICN.se, digits=1)
		
	
		# Computing the mean final cambial cell number (FCN) and its variation for each tree
		# ----------------------------------------------------------------------------------
		FCN2 <- aggregate(FCN1[, "nC"], by=list(Tree=FCN1$Tree), mean, na.rm=TRUE)
		names(FCN2) <- c("Tree", "FCN.mean")
		FCN2$FCN.mean <- round(FCN2$FCN.mean, digits=1)
	
		FCN3 <- aggregate(FCN1[, "nC"], by=list(Tree= FCN1$Tree), sd, na.rm=TRUE)
		names(FCN3) <- c("Tree", "FCN.sd")
	
		FCN4 <- aggregate(FCN1[, "c"], by=list(Tree= FCN1$Tree), sum, na.rm=TRUE)
		names(FCN4) <- c("Tree", "FCN.nb")
		
		FCNa <- merge(FCN2, FCN3, by="Tree")
		FCN <- merge(FCNa, FCN4, by="Tree")
		
		FCN$FCN.se <- FCN$FCN.sd /sqrt(FCN$FCN.nb)
		
		FCN$FCN.sd <- round(FCN$FCN.sd, digits=1)
		FCN$FCN.se <- round(FCN$FCN.se, digits=1)
		
	} # end if (stat=="mean")
	
	else {
		
		# Computing the median tree-ring cell number (RCN) and its variation for each tree
		# --------------------------------------------------------------------------------
		RCN2 <- aggregate(RCN1[, "nT"], by=list(Tree=RCN1$Tree), median, na.rm=TRUE)
		names(RCN2) <- c("Tree", "RCN.median")
		
		RCN3 <- merge(RCN1, RCN2, by="Tree")
		RCN3$d <- abs(RCN3$nT - RCN3$RCN.median)	
		RCN4 <- aggregate(RCN3[, "d"], by=list(Tree=RCN1$Tree), median, na.rm=TRUE)
		names(RCN4) <- c("Tree", "RCN.mad")
		
		RCN5 <- aggregate(RCN1[, "c"], by=list(Tree= RCN1$Tree), sum, na.rm=TRUE)
		names(RCN5) <- c("Tree", "RCN.nb")

		RCNa <- merge(RCN2, RCN4, by="Tree")
		RCN <- merge(RCNa, RCN5, by="Tree")
		
		RCN$RCN.mae <- RCN$RCN.mad / sqrt(RCN$RCN.nb)
		RCN$RCN.mae.normcor <- 3/4 * RCN$RCN.mae
		
		RCN$RCN.mae <- round(RCN$RCN.mae, digits=1)
		RCN$RCN.mae.normcor  <- round(RCN$RCN.mae.normcor, digits=1)
		
	
		# Computing the initial cambial cell number (ICN) for each tree
		# -------------------------------------------------------------
		ICN2 <- aggregate(ICN1[, "nC"], by=list(Tree= ICN1$Tree), median, na.rm=TRUE)
		names(ICN2) <- c("Tree", "ICN.median")
	
		ICN3 <- merge(ICN1, ICN2, by="Tree")
		ICN3$d <- abs(ICN3$nC - ICN3$ICN.median)	
		ICN4 <- aggregate(ICN3[, "d"], by=list(Tree=ICN1$Tree), median, na.rm=TRUE)
		names(ICN4) <- c("Tree", "ICN.mad")
		
		ICN5 <- aggregate(ICN1[, "c"], by=list(Tree= ICN1$Tree), sum, na.rm=TRUE)
		names(ICN5) <- c("Tree", "ICN.nb")

		ICNa <- merge(ICN2, ICN4, by="Tree")
		ICN <- merge(ICNa, ICN5, by="Tree")
		
		ICN$ICN.mae <- ICN$ICN.mad / sqrt(ICN$ICN.nb)
		ICN$ICN.mae.normcor <- 3/4 * ICN$ICN.mae
		
		ICN$ICN.mae <- round(ICN$ICN.mae, digits=1)
		ICN$ICN.mae.normcor  <- round(ICN$ICN.mae.normcor, digits=1)
	
	
		# Computing the final cambial cell number (FCN) for each tree
		# -----------------------------------------------------------
		FCN2 <- aggregate(FCN1[, "nC"], by=list(Tree=FCN1$Tree), median, na.rm=TRUE)
		names(FCN2) <- c("Tree", "FCN.median")
	
		FCN3 <- merge(FCN1, FCN2, by="Tree")
		FCN3$d <- abs(FCN3$nC - FCN3$FCN.median)	
		FCN4 <- aggregate(FCN3[, "d"], by=list(Tree=FCN1$Tree), median, na.rm=TRUE)
		names(FCN4) <- c("Tree", "FCN.mad")
		
		FCN5 <- aggregate(FCN1[, "c"], by=list(Tree= FCN1$Tree), sum, na.rm=TRUE)
		names(FCN5) <- c("Tree", "FCN.nb")

		FCNa <- merge(FCN2, FCN4, by="Tree")
		FCN <- merge(FCNa, FCN5, by="Tree")
		
		FCN$FCN.mae <- FCN$FCN.mad / sqrt(FCN$FCN.nb)
		FCN$FCN.mae.normcor <- 3/4 * FCN$FCN.mae
		
		FCN$FCN.mae <- round(FCN$FCN.mae, digits=1)
		FCN$FCN.mae.normcor  <- round(FCN$FCN.mae.normcor, digits=1)
	
	} # end else
	
	
	# Merging data frame
	# ------------------
	RCN$Tree <- as.factor(RCN$Tree)
	ICN$Tree <- as.factor(ICN$Tree)
	FCN$Tree <- as.factor(FCN$Tree)
	
	ODF1 <- merge(RCN, ICN, by="Tree")
	ODF2 <- merge(ODF1, FCN, by="Tree")
	
	ODF <- ODF2[order(ODF2$Tree), ]
	
	
	# Returning the results
	# ---------------------
	return(ODF)

} # end function computeRingCellNb