standardiseCellCount <- function(data, na.rm="FALSE") {

 # *************************************************************************************************
 # standardiseCellCount() function definition
 # ----------------------------------------
 #
 # Data processing function for creating a data frame where values are standardised by the
 #  measurements (number of cells or alternatively width) of the previous tree-ring
 #
 # Arguments:
 # 		data: data frame with imposed column names (see exemple)
 #		na.rm: logical indicating if lines where P is NA should be left as they are without applying any
 #			standardisation (FALSE by default) or removed from the record (TRUE)
 #
 # Output:
 #		a dataframe containing the standardised values
 #
 # Version: 1.0-0
 # Started: 25 January 2011
 # Last modifications: 25 January 2011
 # *************************************************************************************************

	IDF <- data
	
	# Computing the average value of the previous ring for each tree
	IDF$Tree <- as.factor(IDF$Tree)	
	Means <- aggregate(IDF[, c("P")], by=list(IDF$Tree), FUN=mean, na.rm=TRUE)
	names(Means) <- c("Tree", "P.mean")

	# Computing the correction coeficient for each sample
	DF <- merge(IDF, Means)
	
	if (na.rm=="TRUE") { DF$Coef <- DF$P / DF$P.mean }
	if (na.rm=="FALSE") { DF$Coef <- ifelse(is.finite(DF$P / DF$P.mean), DF$P / DF$P.mean, 1) }
	
	DF <- DF[is.finite(DF$Coef)=="TRUE", ]
	
	# Computing the standardised value for each sample
	ODF <- DF[, c(-10, -11, -12)]
	ODF$nC <- round(DF$nC / DF$Coef, digit=1)
	ODF$nE <- round(DF$nE / DF$Coef, digit=1)
	ODF$nL <- round(DF$nL / DF$Coef, digit=1)
	ODF$nM <- round(DF$nM / DF$Coef, digit=1)
	ODF$nLM <- round(DF$nLM / DF$Coef, digit=1)
	ODF$nELM <- round(DF$nELM / DF$Coef, digit=1)
	
	# Returning output dataframe
	# --------------------------
	ODF <- ODF[order(ODF$Tree), ]
	
	return(ODF)
	
} # End function standardiseCellCount