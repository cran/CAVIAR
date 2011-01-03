\name{fitGompertz}

\alias{fitGompertz}

\title{fitGompertz}

\description{
fits a Gompertz model to the observations.
}

\usage{
fitGompertz(data, asymptote=NULL, plot.fitting=FALSE, ttl="", lblx="Date", lbly="Total cumulated cell number", plot.OP=FALSE)
}


\arguments{
   \item{data}{a data.frame containing dates in day of year and total cell numbers for each tree to be treated}
   \item{asymptote}{a data.frame containing a value for the asymptote parameter for each tree}
   \item{plot.fitting}{an optional logical indicating if a graph of the fitting must be plotted (default) or not}
   \item{ttl}{an optional string of characters containing the plot title}
   \item{lblx}{an optional string of characters containing x-axis label}
   \item{lbly}{an optional string of characters containing y-axis label}
   \item{plot.OP}{an optional logical indicating if a verification graph must be plotted or not (default)}
}

\details{
This function fits a Gompertz model (!!! provide here the equation !!!) to the observations and records model parameters, predicted values and goodness-of-fit. \cr
Additionally it computes characteristic biological parameters.
}

\value{
	A data.frame containing the results of the fitting:\cr
	\item{a}{parameter a of the model}
	\item{b}{parameter b of the model}
	\item{k}{parameter k of the model}
	\item{t5}{date at which 5 percent of the cells are produced}
	\item{tip}{date at which the cell production is maximal (Gompertz curve inflection point)}
	\item{t95}{date at which 95 percent of the cells are produced}
	\item{Dt90}{time elapsed between t5 and t95}
	\item{rmax}{maximal rate of cell production}
	\item{r90}{mean rate of cell production computed between t5 and t95}
	\item{R2}{R-square}
	\item{EF}{Modelling Efficiency}
	\item{RMSD}{Root Mean Squared Deviation}
}

\note{
Version: 2.0-0 \cr
Development started: 27 November 2008
Last modifications: 6 Decembre 2010
}

\author{
Cyrille Rathgeber - LERFoB, UMR1092 - INRA Nancy (France) \cr
Email: \email{cyrille.rathgeber@nancy.inra.fr}
}

\references{
Rathgeber C.B.K., Rossi S. , Bontemps J.-D. Influence of social status on intra-annual dynamics of xylem cell production in a silver fir plantation. In preparation.
}

\seealso{
\code{\link{AMA2006}, \link{averageRadialFiles}, \link{computeBootstrapTest}, \link{computeCriticalDates}, \link{computeRingCellNb}, \link{readXLSCellCountTable}, \link{plotWoodFormationCalendar}, \link{plotWoodFormationDynamics}}
}

\examples{
## Loading the example dataset:
data(AMA2006)

## Averaging radial files
DF <- averageRadialFiles(AMA2006)

## Fitting Gompertz
R1 <- fitGompertz(data=DF, plot.fitting=TRUE, ttl="", lblx="Date", lbly="Total cumulated cell number", plot.OP=TRUE)

## Compute the ring cell number to constrain the asymptote of the Gompertz
RCN <- computeRingCellNb(AMA2006, stat="median")

asymp <- RCN[, c("Tree", "RCN.median")]
names(asymp) <- c("Tree", "A")

## Fitting Gompertz with fixed asymptotes
R2 <- fitGompertz(data=DF, asymptote=asymp, plot.fitting=TRUE, ttl="", lblx="Date", lbly="Total cumulated cell number", plot.OP=TRUE)

}

\keyword{Statistical Models}
\keyword{Non-linear Regression}
\keyword{hplot}
