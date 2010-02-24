averageRF <- function(data) {
 # *************************************************************************************
 # averageRF() function definition
 # -------------------------------
 #
 # Data processing function for creating a data frame where values are averaged by radial
 #  files from the raw data table and where cumulative number of cells are computed
 #
 # Arguments:
 # 		data: data frame with imposed column names (see exemple)
 #
 # Output:
 #		dataframe of averaged values
 #
 # Version: 1.1-0
 # Started: 27 November 2008
 # Last modifications: 1 April 2010
 # *************************************************************************************

	IDF <- data
	
	# Converting Tree and DY as factors
	# ---------------------------------
	IDF$Tree <- as.factor(IDF$Tree)
	IDF$DY <- as.factor(IDF$DY)

	# Creating a new aggregated dataframe by radial files
	# ---------------------------------------------------
	ODF <- aggregate(IDF[, c("nC", "nE", "nL", "nM")],
		by=list(IDF$Tree, IDF$DY), FUN=mean, na.rm=TRUE)
	names(ODF) <- c("Tree", "DY", "nC", "nE", "nL", "nM")
	
	# Ordering the new aggregated data frame
	# --------------------------------------
	ODF <- ODF[order(ODF$Tree, ODF$DY), ]

	ODF$DY <- as.numeric(levels(ODF$DY))[ODF$DY]
	
	# Computing cumulative cell number
	# --------------------------------
	# nM: number of mature cells (not changed)
	# nLM: number of cells in lignification phase + mature cells
	# nELM: number of cells in expansion phase + cells in lignification phase + mature cells
	ODF$nLM <- ODF$nL + ODF$nM
	ODF$nELM <- ODF$nE + ODF$nL + ODF$nM

	# Rounding values
	# ---------------
	ODF[, c("nC", "nE", "nL", "nM", "nLM", "nELM")] <- round(ODF[, 
		c("nC", "nE", "nL", "nM", "nLM", "nELM")], digits=1)
	
	# Returning output dataframe
	# --------------------------
	return(ODF)
	
} # End function averageRF

