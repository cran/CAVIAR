# *************************************************************************************
# Cyrille RATHGEBER - INRA Nancy - 14 July to 15 September 2009, 9 & 11 February 2010
# Purpose -> Testing if critical dates or durations are different between groups or not
# Versions:
# 		1. Two sample bootstrap test
#     	2. Test on median
#		3. Test on variance
#		4. Test with dispersion of critical dates
# Warnings:
#		-> Dispersion is not implemented for test based on centred statistics !
# ************************************************************************************

computeBootstrapTest <- function(y, z, y.sd=NA, z.sd=NA, stat="median", centring=FALSE, iter=1000, 	var.name="", out=FALSE, plot=FALSE, plot.label="") {
		
 # *************************************************************************************************
 # computeBootstrapTest function definition
 # ----------------------------------------
 #
 # Unilateral bootstrap permutation tests for testing equality of mean, median or variance between two
 # groups
 #
 # Arguments:
 #				- y: data vector for the first group
 #				- z: data vector for the second group
 #				- y.sd (optional): data vector for the first group with standard deviations
 #				- z.sd (optional): data vector for the second group with standard deviations
 #				- stat (optional): 	"mean" for mean comparison;
 #										"Student" for mean comparison using Student's t statistic;
 #										"median" for median comparison, default = median;
 #										"var" for variance comparison;
 #										"disp.mean" for dispersion around the mean comparison;
 #										"disp.median" for dispersion around the median comparison;
 #				- centring (optional): centring data (TRUE), or not (FALSE), default = FALSE
 #				- iter (optional): number of iteration for the bootstrap, default = 1,000
 #				- var.name (optional): name of the variable to be tested
 #				- out (optional): TRUE (or FALSE) if bootstrapped series are desired (or not),
 #					default = FALSE
 #				- plot (optional): TRUE (or FALSE) if plot is desired (or not), default = FALSE
 # Outputs:
 #			- list containing the results of the test
 # 				Observed statistic of group 1 is superior to the one of group 2.
 #				H0 : the statistic of reference population 1 is equal the one of reference population 2
 #				H1 : the statistic of reference population 1 is greater than the one of reference
 #					population 2
 #			- Histrogram of the bootrap distribution with the observed statistic
 #
 # *************************************************************************************************
		
	# Associated function
	# -------------------	
	
	C <- function(n, p) {
		# *********************************************************************************************
		# C() function definition
		# ---------------------------------------------------------------------------------------------
		# C(n, p) compute the number of possible combinations of p elements drawn simultaneously from 
		# an ensemble of n elements.
		# *********************************************************************************************
		
		factorial(n) / (factorial(p) * factorial(n - p))
	}
	
	
	# Removing missing values from x and y
	# ------------------------------------
	y <- y[is.na(y)==FALSE]	
	z <- z[is.na(z)==FALSE]
	
	# Constructing the common dataset
	# -------------------------------	
	x <- c(y, z)
	
	if (is.na(y.sd)==FALSE && is.na(z.sd)==FALSE) {
		Y <- numeric()
		Z <- numeric()
	
		for (i in 1: length(y)) {
			Y <- c(Y, rnorm(100, mean=y[i], sd=y.sd[i]))
			Z <- c(Z, rnorm(100, mean=z[i], sd=z.sd[i]))
		}
		x <- c(Y, Z)		
	}

	
	# Local variables initialization
	# ------------------------------
	
	V1greaterthanV2 <- TRUE
	
	ts.boot <- numeric(iter)	
	
	# Computing the number of possible combinations
	# ---------------------------------------------
	p <- length(y)
	q <- length(z)
	n <- p + q	
	a1 <- n + p - 1
	a2 <- n + q - 1
	b <- n - 1
	NbC <- C(a1, b) * C(a2, b) # Number of possible combinations considering the two groups
	
	
	# Simple mean comparison
	# ----------------------
	 if (stat == "mean") {
		
		test.stat <- "mean difference"
		xlab <- "Difference between group means"
		
		# observed studied statistic for group 1 and 2
		s1 <- mean(y)
		s2 <- mean(z)

		# Computing observed difference between means
		if (s1 > s2) {
			sa <- s1
			na <- length(y)
			sb <- s2
			nb <- length(z)
		}
		else {
			V1greaterthanV2 <- FALSE
			sa <- s2
			na <- length(z)
			sb <- s1
			nb <- length(y)
		}
		ts.obs <- sa - sb
		
		if (centring == TRUE) {
			# adjusting data of group 1 and 2
			y <- y - mean(y) + mean(x)
			z <- z - mean(z) + mean(x)
			x <- c(y, z)
		}
		
		# Computing bootstrapped differences between means
		for (i in 1:iter) {
			ts.boot[i] <- mean(sample(x, na, replace=TRUE)) - mean(sample(x, nb, replace=TRUE))
		}
	}
	
	
	# Studentized mean comparison
	# ---------------------------
	if (stat == "Student") {
		
		test.stat <- "Student t"
		xlab <- "Studentized difference between group means"
		
		# observed studied statistic for group 1 and 2
		s1 <- mean(y)
		s2 <- mean(z)

		# Computing observed difference between means
		if (s1 > s2) {
			sa <- s1
			va <- var(y)
			na <- length(y)
			sb <- s2
			vb <- var(z)
			nb <- length(z)
		}
		else {
			V1greaterthanV2 <- FALSE
			sa <- s2
			va <- var(z)
			na <- length(z)
			sb <- s1
			vb <- var(y)
			nb <- length(y)
		}
		ts.obs <- (sa - sb) / sqrt(va/na + vb/nb)
		
		if (centring == TRUE) {
			# adjusting data of group 1 and 2
			y <- y - mean(y) + mean(x)
			z <- z - mean(z) + mean(x)
			x <- c(y, z)
		}
		
		# Computing bootstrapped differences between means
		for (i in 1:iter) {
			a <- sample(x, na, replace=TRUE)
			b <- sample(x, nb, replace=TRUE)
			ts.boot[i] <- (mean(a) - mean(b)) / sqrt(var(a)/na + var(b)/nb)
		}
	}
	
	
	# Median comparison
	# -----------------
	if (stat == "median") {
		
		test.stat <- "median difference"
		xlab <- "Difference between group medians"

		# observed studied statistic for group 1 and 2
		s1 <- median(y)
		s2 <- median(z)
		
		# Computing observed difference between medians
		if (s1 > s2) {
			sa <- s1
			na <- length(y)
			sb <- s2
			nb <- length(z)			
		}
		else {
			V1greaterthanV2 <- FALSE
			sa <- s2
			na <- length(z)			
			sb <- s1
			nb <- length(y)			
		}
		ts.obs <- sa - sb
		
		# Computing bootstrapped differences between medians
		for (i in 1:iter) {
			ts.boot[i] <- median(sample(x, na, replace=TRUE)) - median(sample(x, nb, replace=TRUE))
		}
	}
	
		
	# Variance comparison
	# -------------------
	if (stat == "var") {
		
		test.stat <- "log of variance ratio (F)"
		xlab <- "Ratio between group variance"
		
		# observed studied statistic for group 1 and 2
		s1 <- var(y)
		s2 <- var(z)
		
		# Computing observed ratio between variance
		if (s1 > s2) {
			sa <- s1
			na <- length(y)
			sb <- s2
			nb <- length(z)		
		}
		else {
			V1greaterthanV2 <- FALSE
			sa <- s2
			na <- length(z)		
			sb <- s1
			nb <- length(y)	
		}
		ts.obs <- log(sa / sb)
		
		# Computing bootstrapped ratio between variances
		for (i in 1:iter) {
			ts.boot[i] <- log(var(sample(x, na, replace=TRUE)) / var(sample(x, nb, replace=TRUE)))
		}
	}
	
	
	# Dispersion comparison based on mean
	# -----------------------------------
	if (stat == "disp.mean") {
		
		test.stat <- "dispersion around the mean ratio"
		xlab <- "Ratio between group dispersion around the mean"
		
		# data transformation (Levene, 1960)
		y.trans <- abs(y - mean(y))
		z.trans <- abs(z - mean(z))
		
		# observed studied statistic for group 1 and 2
		s1 <- mean(y.trans)
		s2 <- mean(z.trans)
		
		# Computing observed ratio between dispersions
		if (s1 > s2) {
			sa <- s1
			na <- length(y.trans)
			sb <- s2
			nb <- length(z.trans)		
		}
		else {
			V1greaterthanV2 <- FALSE
			sa <- s2
			na <- length(z.trans)		
			sb <- s1
			nb <- length(y.trans)	
		}
		ts.obs <- sa / sb
		
		# Computing bootstrapped ratio between dispersions
		for (i in 1:iter) {
			ts.boot[i] <- mean(sample(x, na, replace=TRUE)) / mean(sample(x, nb, replace=TRUE))
		}
	}
	
	
	# Dispersion comparison based on median
	# -------------------------------------
	if (stat == "disp.median") {
		
		test.stat <- "dispersion around the median ratio"
		xlab <- "Ratio between group dispersion around the median"
		
		# data transformation (Levene, 1960)
		y.trans <- abs(y - median(y))
		z.trans <- abs(z - median(z))
		
		# observed studied statistic for group 1 and 2
		s1 <- mean(y.trans)
		s2 <- mean(z.trans)		
		
		# Computing observed ratio between dispersions
		if (s1 > s2) {
			sa <- s1
			na <- length(y.trans)
			sb <- s2
			nb <- length(z.trans)		
		}
		else {
			V1greaterthanV2 <- FALSE
			sa <- s2
			na <- length(z.trans)		
			sb <- s1
			nb <- length(y.trans)	
		}
		ts.obs <- sa / sb
		
		# Computing bootstrapped ratio between dispersions
		for (i in 1:iter) {
			ts.boot[i] <- mean(sample(x, na, replace=TRUE)) / mean(sample(x, nb, replace=TRUE))
		}
	} 
	
	
	# Computing the significance level
	# --------------------------------
	ts <- c(ts.boot, ts.obs)
	dr <- as.integer(round(log10(iter), digits=0))
	asl <- signif(length(ts[ts > ts.obs]) / length(ts), dr)
	
	
	# Writting the results
	# --------------------
	L1 <- paste("Bootstrapped test of significance  - variable: ", var.name)
	L2 <- paste("Studied statistic (s): ", stat, "; tested statistic: ", test.stat)
	L3 <- paste("s group 1  = ", sa, " ; s group 2 = ", sb, " ; ts obs. = ", ts.obs)
	
	L4 <- " "
	if (V1greaterthanV2 == TRUE) {
		L4 <- paste("H0: s group 1 = s group 2 ; H1: s group 1 > s group 2")
	}
	if (V1greaterthanV2 == FALSE) {
		L4 <- paste("H0: s group 1 = s group 2 ; H1: s group 1 < s group 2")
	}
	
	L5 <- paste("Nb. iterations = ", iter, "; ASL = ", asl)
	
	Summary <- c(L1, L2, L3, L4, L5)
	KeyValues <- c(ts.obs, iter, NbC, asl)
	names(KeyValues)  <- c("D.obs", "Nb Combi.", "Nb iter.", "ASL")

	if (out == FALSE) Output <- list(summary=Summary, ASL=asl)
	else Output  <- list(summary=Summary, ts.obs=ts.obs, iterations=iter, nb.combinations=NbC,
							ASL=asl, ts=ts)
	
	
	# Plotting the results
	# --------------------
	
	if (plot==TRUE) {
	
		# Drawing the simulated histogram
		H <- hist(ts, main=NULL, xlab=xlab)

		# Drawing the observation line
		abline(v=ts.obs, col="red")
		
		maxcount <- max(H$counts)
		text(ts.obs, maxcount, labels="Observed", pos=2, col="red", srt=90)
	
		# Writting plot title
		if (plot.label == "") {
			mtext(L1, side=3, line=3, adj=0.0, cex=1)
			mtext(L2, side=3, line=2, adj=0.0, cex=1)
			mtext(L4, side=3, line=1, adj=0.0, cex=1)
			mtext(L5, side=3, line=0, adj=0.0, cex=1)
		}
		else {
			mtext(plot.label, side=3, line=1, adj=0.0, cex=1)
		}

	}
	
	
	# Returning the results
	# ---------------------
	
	return(Output)
	
} # End function computeBootstrapTest
