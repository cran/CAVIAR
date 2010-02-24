\name{computeBootstrapTest}
\alias{computeBootstrapTest}

\title{computeBootstrapTest}

\description{
Computes bootstrap test for comparing 2 groups.
}

\usage{
computeBootstrapTest(y, z, y.sd = NA, z.sd = NA, stat = "median", centring = FALSE,
		     iter = 1000, var.name = "", out = FALSE, plot = FALSE, plot.label = "")
}


\arguments{
  \item{y}{a vector containing data for the first group}
  \item{z}{a vector containing data for the second group}
  \item{y.sd}{an optional vector containing standard deviations for the first group}
  \item{z.sd}{an optional vector containing standard deviations for the second group}
  \item{stat}{ the method to be used: \cr
	"mean" for mean comparison;\cr
	"Student" for mean comparison using Student's t statistic;\cr
	"median" for median comparison (default);\cr
	"var" for variance comparison;\cr
	"disp.mean" for dispersion around the mean comparison;\cr
	"disp.median" for dispersion around the median comparison.
  }
  \item{centring}{an optional logical indicating if data must be centred or not (default=FALSE)}
  \item{iter}{an optional integer indicating the number of iterations for the bootstrap (default=1000)}
  \item{var.name}{an optional character labeling the variable to be tested}
  \item{out}{an optional logical indicating if bootstrapped series must be outputted or not (default=FALSE)}
  \item{plot}{an optional logical indicating if a plot must be outputted or not (default=FALSE)}
  \item{plot.label}{an optional character containing the plot title}
}

\details{
Computes unilateral bootstrap permutation tests for testing equality of mean, median or variance between two groups. \cr
Several test statistics can be used: mean difference, Studentized mean difference, median difference, logarithm ratio of variance, dispersion around the mean ratio and dispersion around the median ratio. \cr
For test statistics based on the mean, data can be centred. \cr
Moreover this function can also be feed with two additional vectors containing standard deviations around the critical dates in order to take dispersion into account.
It can plot (option{plot}) an histogram of the bootstrap distribution of the tested statistic along with its observed value.
}

\value{
A list containing the results of the test:\cr
\item{summary}{Summary of the test}
\item{ASL}{Signification level of the test}
\item{ts.obs}{Observed value of the tested statistic (\option{out=TRUE})}
\item{iterations}{number of iterations (\option{out=TRUE})}
\item{nb.combinations}{number of possible combinations (\option{out=TRUE})}
\item{ts}{bootstrap series (\option{out=TRUE})}
}

\section{Warning }{
Dispersion is not implemented for tests based on centred statistics (\option{centring=TRUE})!
}

\note{
Version: 4.0 \cr
Last modifications: 19 July 2010
}

\author{
Cyrille Rathgeber - LERFoB, UMR1092 - INRA Nancy (France) \cr
Email: \email{cyrille.rathgeber@nancy.inra.fr}
}

\references{
Rathgeber C.B.K., Longuetaud F. , Mothe F., Cuny H., Le Moguedec G. 2010. Phenology of wood formation: data processing, analysis and visualisation using R. Accepted in \emph{Dendrochronologia}.
}

\seealso{
\code{\link{AMA2006}, \link{averageRF}, \link{computeCriticalDates}, \link{computeRCN}, \link{fitGompertz}, \link{plotWoodCalendar}}
}

\examples{
## Loading the example dataset:
data(AMA2006)

## Computing wood formation critical dates for firs and pines:
Fir.cdd <- computeCriticalDates(AMA2006[AMA2006$Sp=="ABAL", ])
Pine.cdd <- computeCriticalDates(AMA2006[AMA2006$Sp=="PISY", ])

## Testing if pines start xylogenesis before firs
computeBootstrapTest(Fir.cdd$bE, Pine.cdd$bE)

## Testing if pines start xylogenesis before firs with taking dispersion into account
computeBootstrapTest(Fir.cdd$bE, Pine.cdd$bE, Fir.cdd$bE.sd, Pine.cdd$bE.sd)
}

\keyword{htest}
\keyword{nonparametric}
\keyword{hplot}
