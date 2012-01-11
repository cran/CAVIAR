# **************************************************************************************
# FUNCTION readXLSCellCountTable
# ------------------------------
#
# Author: Cyrille RATHGEBER - LERFoB UMR1092 - INRA Nancy
# Purpose: Reading Excel spreadsheets and collecting them in a database-like data.frame
# Version:
#  	1.1. Using data
#	1.2. Taking into account previous ring measurement
#	1.3. Adding the year if asked (optional)
#   	1.4-1. Removing the year option and correcting a bug concerning P
#
# Started: 13 Octobre 2010
# Last modified: 11 January 2012
#
# **************************************************************************************


readXLSCellCountTable <- function(xls.file.name, sheets, sampling.dates, tree.names) {

 # *******************************************************************************************************
 # readXLSCellCountTable() function definition
 # -------------------------------------------
 #
 # Data manipulation function for converting Excel spreadsheet file to database-like data.frame
 # 
 # Arguments:
 # 		- xls.file.name: a character containing the path to the Excel file
 #		- sheets: a vector containing the numbers of the sheet to be read
 #		- sampling.dates: a data.frame containing the correspondence between sampling numbers
 #			and dates (in DOY)
 #		- tree.names: a vector containing the names of the trees
 #
 # Output:
 #		- a data.frame in database-like format
 #
 #
 # *******************************************************************************************************
	
	# Loading required library
	library(gdata)
	
	# Reading the first sheet of the input data file
	InDF  <- read.xls(xls.file.name, sheet=sheets[1])
	N <- names(InDF)
	
	# Initiating objects
	cDF <- data.frame()	

	# Case 1: cell count table without previous ring counting
	# -------------------------------------------------------

	if (is.na(match("P1", N)) == TRUE) {

		for (i in 1:length(sheets)) {

			# Reading the input data file
			InDF  <- read.xls(xls.file.name, sheet=sheets[i])
		
			# Merging cell count table and sampling dates
			IDF <- merge(sampling.dates, InDF)

			# Creating a "database-type" data.frame
			Tree <- as.factor(rep.int(tree.names[i], 3*nrow(IDF)))
			Sample <- as.factor(rep.int(IDF$Sample, 3))
			DY <- as.integer(rep.int(IDF$DY, 3))
			RF <- as.factor(c(rep.int(1, nrow(IDF)), rep.int(2, nrow(IDF)), rep.int(3, nrow(IDF))))
			nC <- as.integer(c(IDF$C1, IDF$C2, IDF$C3))
			nE <- as.integer(c(IDF$E1, IDF$E2, IDF$E3))
			nL <- as.integer(c(IDF$L1, IDF$L2, IDF$L3))
			nM <- as.integer(c(IDF$M1, IDF$M2, IDF$M3))

			aDF <- data.frame(Tree, Sample, DY, RF, nC, nE, nL, nM)

			# Ordering the final table
			bDF <- aDF[order(aDF$Sample), ]
	
			# Adding the new table to the database
			cDF <- rbind(cDF, bDF)

		} # End for

		# Removing blank lines
		# Constructing a function for testing blank lines
		removeTest <- function(X) {is.na(X[1]) & is.na(X[2]) & is.na(X[3]) &
			is.na(X[4])}		
		
		# Applying the test function to the data
		cDF$Remove <- apply(cDF[, 5:8], 1, removeTest)
		
		# Removing blank lines				
		dDF <- cDF[cDF$Remove == FALSE, ]
	
	} # End if
	

	# Case 2: cell count table with previous ring counting
	# ----------------------------------------------------

	if (is.na(match("P1", N)) == FALSE) {

		for (i in 1:length(sheets)) {

			# Reading the input data file
			InDF  <- read.xls(xls.file.name, sheet=sheets[i])
		
			# Merging cell count table and sampling dates
			IDF <- merge(sampling.dates, InDF)

			# Creating a "database-type" data.frame
			Tree <- as.factor(rep.int(tree.names[i], 3*nrow(IDF)))
			Sample <- as.factor(rep.int(IDF$Sample, 3))
			DY <- as.integer(rep.int(IDF$DY, 3))
			RF <- as.factor(c(rep.int(1, nrow(IDF)), rep.int(2, nrow(IDF)), rep.int(3, nrow(IDF))))
			nC <- as.integer(c(IDF$C1, IDF$C2, IDF$C3))
			nE <- as.integer(c(IDF$E1, IDF$E2, IDF$E3))
			nL <- as.integer(c(IDF$L1, IDF$L2, IDF$L3))
			nM <- as.integer(c(IDF$M1, IDF$M2, IDF$M3))
			P <- as.numeric(c(IDF$P1, IDF$P2, IDF$P3))

			aDF <- data.frame(Tree, Sample, DY, RF, nC, nE, nL, nM, P)
	
			# Ordering the final table
			bDF <- aDF[order(aDF$Sample), ]
	
			# Adding the new table to the database
			cDF <- rbind(cDF, bDF)

		} # End for

		# Removing blank lines
		# Constructing a function for testing blank lines
		removeTest <- function(X) {is.na(X[1]) & is.na(X[2]) & is.na(X[3]) &
			is.na(X[4])}		
		
		# Applying the test function to the data
		cDF$Remove <- apply(cDF[, 5:8], 1, removeTest)
		
		# Removing blank lines				
		dDF <- cDF[cDF$Remove == FALSE, ]
		
	} # End if
		
	# Returning the data-base like table 
	# -----------------------------------
	
	eDF <- dDF[, 1:(ncol(dDF)-1)]
	
	return(eDF)

} # End function readXLSCellCountTable
