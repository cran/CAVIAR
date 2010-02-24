fitGompertz <- function(x, y, plot=TRUE, ttl1, ttl2, 
	lblx="Date", lbly="Cumulated cell number (CCN)") {
 # *************************************************************************************
 # fitGompertz() function definition
 # ---------------------------------
 #
 # Data analysis function for fitting a Gompertz model to cumulative cell count data
 #
 # Arguments:
 # 		x: a numeric vector containing dates in day of year
 #		y: a numeric vector containing cell numbers
 #		plot: an optional logical indicating if a graph must be plotted (default) or not
 #		ttl1: an optional string of characters containing the plot title
 #		ttl2: an optional string of characters containing the plot subtitle
 #		lblx: an optional string of characters containing x-axis label
#		lbly: an optional string of characters containing y-axis label
 #
 # Output:
 #		dataframe of averaged values
 #
 # Version: 1.1-0
 # Started: 27 November 2008
 # Last modifications: 19 July 2010
 # *************************************************************************************
	
	# Gompertz fitting
	# ----------------	
	# Gompertz assymptote initiated by the maximal number of cells
	# b and k parameters initiated using value from Rossi et al. 2003
	Gompertz <- nls(y ~ a*exp(-exp(b - k*x)), start=list(a=max(y), b=7, k=0.04))
	a <- as.double(coef(Gompertz)["a"])
	b <- as.double(coef(Gompertz)["b"])
	k <- as.double(coef(Gompertz)["k"])
	
	# Computing the predicted value and the goodness-of-fit
	G <- predict(Gompertz, list(x=x))
	
	# Computing the R2 --> To be replaced by EF !!!
	R2 <- (cor(y, G))^2
	
	
	# Computing the biological parameters from the fitted equations
	# -------------------------------------------------------------
	# - t5: date at which 5% of the cells are produced
	# - tip: date of inflexion point
	# - t95: date at which 95% of the cells are produced
	# - Dt90: t95 - t5
	# - rx: maximal growth rate
	# - r90: mean growth rate computed between 5 and 95% of the produced cells
	t5 <- (b - log(-log(0.05)))/k
	tip <- b/k
	t95 <- (b - log(-log(0.95)))/k
	Dt90 <- log(log(0.05)/log(0.95))/k
	rmax <- (k*a)/exp(1)
	r90 <- (0.9*k*a)/(log(log(0.05)/log(0.95)))
	
	
	# Ploting the data and the fitted Gompertz curve if option plot=TRUE (default)
	# ----------------------------------------------------------------------------
	if(plot==TRUE) {
		
		# Setting the plot region
		ymax <- max(c(y, G), na.rm=TRUE)
		xmin <- min(x)
		xmax <- max(x)
		plot(x, y, type="n", xlim=c(xmin, xmax), ylim=c(0, ymax), ann=FALSE, axes=FALSE)

		# Drawing individual curves
		points(x, y, type="p", pch=19, col="blue")
		points(x, G, type="l", col="red")
	
		# Customising axes
		axis(1)
		mtext(lblx, side=1, line=2.5)
	
		axis(2)
		mtext(lbly, side=2, line=2.5)

		# Writting plot title
		mtext(ttl1, side=3, line=2.5, adj=0.0, cex=1.0)
		mtext(ttl2, side=3, line=1.5, adj=0.0, cex=1.0)
		mtext(paste("R-squared = ", round(R2, digit=2)), side=3, line=0.5, adj=0.0, cex=1.0)
	
		# Additional labels
		# -----------------
		adjx <- 2 * (xmax - xmin)/100

		# Assymptote
		abline(h=a, lty=2, col="green")
		text(xmin, a - 0.05*a, "A", col="green")

		# t5
		points(c(t5, t5), c(-10, 0.05*a), type="l", lty=3, col="green")
		mtext(expression(t[5]), side=1, line=-1.5, at=(t5+adjx) , col="green")
	
		# t95
		points(c(t95, t95), c(-10, 0.95*a), type="l", lty=3, col="green")
		mtext(expression(t[95]), side=1, line=-1.5, at=(t95+adjx), col="green")

		# Point of inflexion
		VLXi <- c(tip, tip)
		VLYi <- c(-10, a/exp(1))
		points(VLXi, VLYi, type="l", lty=2, col="green")
		mtext(expression(t[ip]), side=1, line=-1.5, at=(tip+adjx), col="green")

		# Mean speed
		arrows(t5, 0.05*a, t95, 0.95*a, length=0.15, angle=25, code=2, col="green")
		xt <- 2 * adjx + ((t5 + t95) / 2)
		text(xt, 0.5*a, expression(r[90]), col="green")

		# Maximal speed
		alpha <- atan(rmax)
		xs <- (rmax) * cos(alpha)
		ys <- sqrt((rmax)^2 - x^2)
		c <- 35
		arrows(tip, a/exp(1), (tip + c*xs), (a/exp(1) + c*ys), length=0.15, angle=25, code=2, 			col="green")
		text(tip - adjx, 0.5*a, expression(r[max]), col="green")
	
		# Legend
		xmin <- min(x)
		r <- ymax / 15
		points(xmin, ymax, pch=19, col="blue")
		text(xmin + 1, ymax, "Observations", pos=4)
		points(xmin, ymax - r, pch=15, col="red")
		text(xmin + 1, ymax - r, "Gompertz", pos=4)
	
		box()
		
	} # End option plot
	
	# Returning the results of the fitting
	# ------------------------------------
	list(a=a, b=b, k=k, R2=R2, Y.predicted=G, t5=t5, tip=tip, t95=t95, Dt90=Dt90, rmax=rmax, r90=r90)

} # End of function fit Gompertz

