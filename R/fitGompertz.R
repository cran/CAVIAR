fitGompertz <- function(data, asymptote=NULL,
	plot.fitting=FALSE, ttl="", lblx="Date", lbly="Total cumulated cell number",
	plot.OP=FALSE) {
 # *************************************************************************************
 # fitGompertz() function definition
 # ---------------------------------
 #
 # Data analysis function for fitting a Gompertz model to cumulative cell count data
 #
 # Arguments:
 # 		data: a dataframe containing tree index, dates in day of year and total number of cells
 #		asymptote: an optional dataframe containing tree index the value of the assymptote
 #		plot.Gompertz: an optional logical indicating if a graph of the fitting must be plotted (default) or not
 #		ttl1: an optional string of characters containing the plot title
 #		ttl2: an optional string of characters containing the plot subtitle
 #		lblx: an optional string of characters containing x-axis label
 #		lbly: an optional string of characters containing y-axis label
 #		plot.OP: an optional logical indicating if a verification graph must be plotted (default) or not
 #
 # Output:
 #		data.frames of values
 #
 # Version: 2.0-0
 # Started: 27 November 2008
 # Last modifications: 6 Decembre 2010
 # *************************************************************************************

	DF  <-  data
	DF$Tree <- as.factor(DF$Tree)
	
	# Initialising objects
	Results <- data.frame()
	
	for (t in 1:nlevels(DF$Tree)) {
		D <- DF$DY[DF$Tree==t]
		N <- DF$nELM[DF$Tree==t]
		
		# Gompertz fitting
		# ----------------	
		if (is.null(asymptote) == TRUE) {
			
			# Gompertz fitting
			# Gompertz asymptote initiated by the maximal number of cells
			# b and k parameters initiated using value from Rossi et al. 2003
			Gompertz <- nls(N ~ a*exp(-exp(b - k*D)), start=list(a=max(N), b=7, k=0.04))
		
			# Extracting the fitted parameters
			a <- as.double(coef(Gompertz)["a"])
			b <- as.double(coef(Gompertz)["b"])
			k <- as.double(coef(Gompertz)["k"])
		}
		else {
			if (nrow(asymptote) == nlevels(DF$Tree)) {
				
			# Gompertz fitting
			# Gompertz asymptote fixed by the user
			# b and k parameters initiated using value from Rossi et al. 2003
			A <- asymptote$A[asymptote$Tree==t]
			Gompertz <- nls(N ~ A*exp(-exp(b - k*D)), start=list(b=7, k=0.04))
		
			# Extracting the fitted parameters
			a <- A
			b <- as.double(coef(Gompertz)["b"])
			k <- as.double(coef(Gompertz)["k"])
			}			
			else {
				print("Error message: The number of asymptote values doesn't match the number of trees !")
			}
		}
	
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
		
		
		# Evaluating the goodness-of-fit
		# ------------------------------	
		# Computing the predicted value
		G <- predict(Gompertz, list(D=D))
	
		# Computing the R2
		R2 <- cor(N, G)^2
		
		# Computing modelling efficiency
		EF <- 1 - sum((N - G)^2) / sum((N - mean(N))^2)
		
		# Computing root mean squared deviation
		RMSD <- sqrt(1/(length(N) - 1) * sum((N - G)^2))	
		
		# Recording the results of the fitting
		# ------------------------------------
		R <- c(t, a, b, k, t5, tip, t95, Dt90, rmax, r90, R2, EF, RMSD)
		Results <- rbind(Results, R)
	
	
		# Ploting the data and the fitted Gompertz curve if option plot.fitting=TRUE (default)
		# -------------------------------------------------------------------------------------
		if (plot.fitting==TRUE) {
		
			# Setting the plot region
			ymax <- max(c(N, G), na.rm=TRUE)
			xmin <- min(D)
			xmax <- max(D)
			plot(D, N, type="n", xlim=c(xmin, xmax), ylim=c(0, ymax), ann=FALSE, axes=FALSE)

			# Drawing individual curves
			points(D, N, type="p", pch=19, col="blue")
			points(D, G, type="l", col="red")
	
			# Customising axes
			axis(1)
			mtext(lblx, side=1, line=2.5)
	
			axis(2)
			mtext(lbly, side=2, line=2.5)

			# Writting plot title
			mtext(ttl, side=3, line=2.5, adj=0.0, cex=1.0)
			mtext(paste("Tree ", t, sep=""), side=3, line=1.5, adj=0.0, cex=1.0)
			mtext(paste("R2 = ", round(R2, digit=2)), side=3, line=0.5, adj=0.0, cex=1.0)
	
			# Additional labels
			# -----------------
			adjx <- 2 * (xmax - xmin)/100

			# Asymptote
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
			ys <- sqrt(rmax^2 - xs^2)
			c <- 50
			arrows(tip, a/exp(1), (tip + c*xs), (a/exp(1) + c*ys), length=0.15,
				angle=25, code=2, col="green")
			text(tip - adjx, 0.5*a, expression(r[max]), col="green")
	
			# Legend
			xmin <- min(D)
			r <- ymax / 15
			points(xmin, ymax, pch=19, col="blue")
			text(xmin + 1, ymax, "Observations", pos=4)
			points(xmin, ymax - r, pch=15, col="red")
			text(xmin + 1, ymax - r, "Gompertz", pos=4)
	
			box()
		
		} # End option plot.fitting
		
		
		# Ploting the observed vs. predicted values if option plot.OP=TRUE (default)
		# --------------------------------------------------------------------------
		if (plot.OP==TRUE) {
					
			# Plotting the data
			plot(G, N, type="p", xlim=c(0, max(G, N)), ylim=c(0, max(G, N)), pch=19, ann=FALSE, axes=FALSE)
			
			# Plotting the xy and regression line
			abline(a=0, b=1, lty=3)
			
			reg1 <- lm(N ~ G)
			abline(reg1, lty=2)			
	
			# Customising axes
			axis(1)
			mtext("Predictions", side=1, line=2.5)
	
			axis(2)
			mtext("Observations", side=2, line=2.5)

			# Writting plot title
			mtext(paste("Verification plot for tree ", t, sep=""), side=3, line=1, adj=0.0,
				cex=1.0, font=2)
				
			# Verification tests
			mtext(paste("Obs = ", round(reg1$coef[2], digit=2), "Pre + ", round(reg1$coef[1], digit=0),
				sep=""), side=3, line=-2, at=0, adj=0, cex=1.0)
			
			r2 <- round(summary(reg1)$adj.r.squared, digit=2)
			mtext(paste("R2 = ", r2), side=3, line=-3, at=0, adj=0, cex=1.0)
			
			diff <- G - N
			reg2 <- lm(diff ~ G)
			s <- summary(reg2)
			t1 <- round(s$coefficients[1, 4], digit=2)
			t2 <- round(s$coefficients[2, 4], digit=2)
			
			mtext(paste("Significance of test a = 0: P = ", t1, sep=""), side=3, line=-4,
				at=0, adj=0, cex=1.0)
			mtext(paste("Significance of test b = 1: P = ", t1, sep=""), side=3, line=-5,
				at=0, adj=0, cex=1.0)
	
			box()
		
		} # End option plot.OP

	} # End for loop
	
	names(Results) <- c("Tree", "a", "b", "k", "t5", "tip", "t95", "Dt90", "rmax", "r90",
		"R2", "EF", "RMSD")
	
	return(Results)

} # End of function fitGompertz