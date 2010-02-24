\name{fitGompertz}

\alias{fitGompertz}

\title{fitGompertz}

\description{
fits a Gompertz model to the observations.
}

\usage{
fitGompertz(x, y, plot=TRUE, ttl1, ttl2, lblx="Date", lbly="Cumulated cell number (CCN)")
}


\arguments{
   \item{x}{a numeric vector containing dates in day of year}
   \item{y}{a numeric vector containing cell numbers}
   \item{plot}{an optional logical indicating if a graph must be plotted (default) or not}
   \item{ttl1}{an optional string of characters containing the plot title}
   \item{ttl2}{an optional string of characters containing the plot subtitle}
   \item{lblx}{an optional string of characters containing x-axis label}
   \item{lbly}{an optional string of characters containing y-axis label}
}

\details{
This function fits a Gompertz model (!!! provide here the equation !!!) to the observations and records model parameters, predicted values and goodness-of-fit. \cr
Additionally it computes characteristic biological parameters.
}

\value{
	A list containing the results of the fitting:\cr
	\item{a}{parameter a of the model}
	\item{b}{parameter b of the model}
	\item{k}{parameter k of the model}
	\item{R2}{R-square}
	\item{Y.predicted}{values predicted by the Gompertz model}
	\item{t5}{date at which 5 percent of the cells are produced}
	\item{tip}{date at which the cell production is maximal (Gompertz curve inflection point)}
	\item{t95}{date at which 95 percent of the cells are produced}
	\item{Dt90}{time elapsed between t5 and t95}
	\item{rmax}{maximal rate of cell production}
	\item{r90}{mean rate of cell production computed between t5 and t95}
}

\note{
Version: 1.0 \cr
Development started: 27 November 2008
Last modifications: 19 July 2010
}

\author{
Cyrille Rathgeber - LERFoB, UMR1092 - INRA Nancy (France) \cr
Email: \email{cyrille.rathgeber@nancy.inra.fr}
}

\references{
Rathgeber C.B.K., Rossi S. , Bontemps J.-D., Longuetaud F. Influence of social status on intra-annual dynamics of xylem cell production in a silver fir plantation. In preparation.
}

\seealso{
\code{\link{AMA2006}, \link{averageRF}, \link{computeBootstrapTest}, \link{computeCriticalDates}, \link{computeRCN}, \link{plotWoodCalendar}}
}

\examples{
## Loading the example dataset:
data(AMA2006)

## Averaging raw data:
DF <- averageRF(AMA2006)

## Fitting Gompertz model to the first tree
D <- DF$DY[DF$Tree==1]
N <- DF$nELM[DF$Tree==1]
ttl1 <- paste("Tree ", 1)
ttl2 <- "Amance Forest, year 2006, silver firs"
lblx <- "Date"
lbly <- "nELM"
		
results <- fitGompertz(D, N, plot=TRUE, ttl1, ttl2, lblx, lbly)
}

\keyword{Statistical Models}
\keyword{Non-linear Regression}
\keyword{hplot}
