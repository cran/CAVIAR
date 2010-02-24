\name{computeRCN}

\alias{computeRCN}

\title{computeRCN}

\description{
Computes characteristic ring cell numbers
}

\usage{
computeRCN(data, stat="median")
}


\arguments{
  \item{data}{a data frame containing cell count data. This data frame must contain at least 7 imposed columns with compulsory names: Tree, DY, RF, nC, nE, nL, nM (see \code{\link{AMA2006}})}
  \item{stat}{an optional string of characters indicating if the function must compute mean or median (default) as summary statistic}
}

\details{
This function computes characteristic tree ring cell numbers (median or mean depending on the option selected) as well as their associated measure of s (median absolute deviation or standard deviation) for each tree of the data frame. \cr
The characteristic tree ring cell numbers are: the final ring cell number (RCN), the initial number of dormant cambial cells (ICN) and the final number of dormant cambial cells (FCN)
}

\value{
A data frame of 6 columns containing the computed characteristic tree ring cell numbers mean or median (e.g. RCN.median or RCN.mean) with their associated measure of variations (e.g. ICN.mad or FCN.sd). \cr
}

\note{
Version: 1.0 \cr
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
\code{\link{AMA2006}, \link{averageRF}, \link{computeBootstrapTest}, \link{computeCriticalDates}, \link{fitGompertz}, \link{plotWoodCalendar}}
}

\examples{
## Loading the example dataset:
data(AMA2006)

## Computing characteristic ring cell numbers for firs and pines:
Fir.RCN <- computeRCN(AMA2006[AMA2006$Sp=="ABAL", ])
Pine.RCN <- computeRCN(AMA2006[AMA2006$Sp=="PISY", ])
}

\keyword{datagen}
