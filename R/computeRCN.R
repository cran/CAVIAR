# **************************************************************************************
# FUNCTION computeRCN
# -------------------
#
# Author: Cyrille RATHGEBER - LERFoB UMR1092 - INRA Nancy
# Purpose: Computing characteristic ring cell numbers
# Version:
#		1. Compute RCN, ICN, FCNB
#
# Started: 25 March 2010
# Last modified: 29 March 2010
#
# **************************************************************************************

computeRCN <- function(data, stat="median") {
 # *************************************************************************************
 # computeRCN() function definition
 # --------------------------------
 #
 # Data processing function for computing final ring cell number (RCN)
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
	
	
	# Option "mean"
	# -------------
	if (stat=="mean") {
		
		# Computing the tree-ring cell number (RCN) for each tree
		# --------------------------------------------------------
		RCN1  <- DF[is.na(DF$nL) == FALSE & DF$nL < 0.1
			& is.na(DF$nM) == FALSE & DF$nM >= 1, c("Tree", "DY", "RF", "nM")]

		RCN2 <- aggregate(RCN1[, "nM"], by=list(Tree= RCN1$Tree), mean, na.rm=TRUE)
		names(RCN2) <- c("Tree", "RCN.mean")
		RCN2$RCN.mean <- round(RCN2$RCN.mean, digits=1)
	
		RCN3 <- aggregate(RCN1[, "nM"], by=list(Tree= RCN1$Tree), sd, na.rm=TRUE)
		names(RCN3) <- c("Tree", "RCN.sd")
		RCN3$RCN.sd <- round(RCN3$RCN.sd, digits=1)
	
		RCN <- merge(RCN2, RCN3, by="Tree")
	
	
		# Computing the initial cambial cell number (ICN) for each tree
		# -------------------------------------------------------------
		ICN1  <- DF[is.na(DF$nE) == FALSE & DF$nE < 0.1
			& is.na(DF$nM) == FALSE & DF$nM < 0.1, c("Tree", "DY", "RF", "nC")]

		ICN2 <- aggregate(ICN1[, "nC"], by=list(Tree= ICN1$Tree), mean, na.rm=TRUE)
		names(ICN2) <- c("Tree", "ICN.mean")
		ICN2$ICN.mean <- round(ICN2$ICN.mean, digits=1)
	
		ICN3 <- aggregate(ICN1[, "nC"], by=list(Tree= ICN1$Tree), sd, na.rm=TRUE)
		names(ICN3) <- c("Tree", "ICN.sd")
		ICN3$ICN.sd <- round(ICN3$ICN.sd, digits=1)
	
		ICN <- merge(ICN2, ICN3, by="Tree")
		
	
		# Computing the final cambial cell number (FCN) for each tree
		# -----------------------------------------------------------
		FCN1  <- DF[is.na(DF$nE) == FALSE & DF$nE < 0.1
			& is.na(DF$nM) == FALSE & DF$nM > 1, c("Tree", "DY", "RF", "nC")]

		FCN2 <- aggregate(FCN1[, "nC"], by=list(Tree=FCN1$Tree), mean, na.rm=TRUE)
		names(FCN2) <- c("Tree", "FCN.mean")
		FCN2$FCN.mean <- round(FCN2$FCN.mean, digits=1)
	
		FCN3 <- aggregate(FCN1[, "nC"], by=list(Tree= FCN1$Tree), sd, na.rm=TRUE)
		names(FCN3) <- c("Tree", "FCN.sd")
		FCN3$FCN.sd <- round(FCN3$FCN.sd, digits=1)
	
		FCN <- merge(FCN2, FCN3, by="Tree")
		
	} # end if (stat=="mean")
	
	else {
		
		# Computing the tree-ring cell number (RCN) for each tree
		# --------------------------------------------------------
		RCN1  <- DF[is.na(DF$nL) == FALSE & DF$nL < 0.1
			& is.na(DF$nM) == FALSE & DF$nM >= 1, c("Tree", "DY", "RF", "nM")]

		RCN <- aggregate(RCN1[, "nM"], by=list(Tree= RCN1$Tree), median, na.rm=TRUE)
		names(RCN) <- c("Tree", "RCN.median")
	
		RCN2 <- aggregate(RCN1[, "nM"], by=list(Tree= RCN1$Tree), quantile, probs=0.25, na.rm=TRUE)
		names(RCN2) <- c("Tree", "q25")
	
		RCN3 <- aggregate(RCN1[, "nM"], by=list(Tree= RCN1$Tree), quantile, probs=0.75, na.rm=TRUE)
		names(RCN3) <- c("Tree", "q75")
	
		RCN$RCN.mad <- round(((RCN3$q75 - RCN2$q25) / 2), digits=1)
	
	
		# Computing the initial cambial cell number (ICN) for each tree
		# -------------------------------------------------------------
		ICN1  <- DF[is.na(DF$nE) == FALSE & DF$nE < 0.1
			& is.na(DF$nM) == FALSE & DF$nM < 0.1, c("Tree", "DY", "RF", "nC")]

		ICN <- aggregate(ICN1[, "nC"], by=list(Tree= ICN1$Tree), median, na.rm=TRUE)
		names(ICN) <- c("Tree", "ICN.median")
		#ICN2$ICN.median <- round(ICN2$ICN.median, digits=1)
	
		ICN2 <- aggregate(ICN1[, "nC"], by=list(Tree= ICN1$Tree), quantile, probs=0.25, na.rm=TRUE)
		names(ICN2) <- c("Tree", "q25")
	
		ICN3 <- aggregate(ICN1[, "nC"], by=list(Tree= ICN1$Tree), quantile, probs=0.75, na.rm=TRUE)
		names(ICN3) <- c("Tree", "q75")
	
		ICN$ICN.mad <- round(((ICN3$q75 - ICN2$q25) / 2), digits=1)
	
	
		# Computing the final cambial cell number (FCN) for each tree
		# -----------------------------------------------------------
		FCN1  <- DF[is.na(DF$nE) == FALSE & DF$nE < 0.1
			& is.na(DF$nM) == FALSE & DF$nM > 1, c("Tree", "DY", "RF", "nC")]

		FCN <- aggregate(FCN1[, "nC"], by=list(Tree=FCN1$Tree), median, na.rm=TRUE)
		names(FCN) <- c("Tree", "FCN.median")
	
		FCN2 <- aggregate(FCN1[, "nC"], by=list(Tree= FCN1$Tree), quantile, probs=0.25, na.rm=TRUE)
		names(FCN2) <- c("Tree", "q25")
	
		FCN3 <- aggregate(FCN1[, "nC"], by=list(Tree= FCN1$Tree), quantile, probs=0.75, na.rm=TRUE)
		names(FCN3) <- c("Tree", "q75")
	
		FCN$FCN.mad <- round(((FCN3$q75 - FCN2$q25) / 2), digits=1)
	
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

} # end function computeRCN
