\name{computeRingCellNb}

\alias{computeRingCellNb}

\title{computeRingCellNb}

\description{
Computes characteristic ring cell numbers
}

\usage{
computeRingCellNb(data, stat="median")
}


\arguments{
  \item{data}{a data frame containing cell count data. This data frame must contain at least 6 imposed columns with compulsory names: Tree, DY, nC, nE, nL, nM (see \code{\link{AMA2006}})}
  \item{stat}{an optional string of characters indicating if the function must compute mean or median (default) as summary statistic}
}

\details{
This function computes characteristic tree ring cell numbers (median or mean depending on the option selected) as well as their associated measure of s (median absolute deviation or standard deviation) for each tree of the data frame. \cr
The characteristic tree ring cell numbers are: the final ring cell number (RCN), the initial number of dormant cambial cells (ICN) and the final number of dormant cambial cells (FCN)
}

\value{
A data frame containing the computed characteristic tree ring cell numbers mean or median (e.g. RCN.median or RCN.mean) with their associated measure of variations (e.g. ICN.mad or FCN.sd) and signal-to-noise ratio (SNR). \cr
}

\note{
Version: 1.4-1 \cr
Development started: 25 March 2010 \cr
Last modifications: 27 January 2011
}

\author{
Cyrille Rathgeber - LERFoB, UMR1092 - INRA Nancy (France) \cr
Email: \email{cyrille.rathgeber@nancy.inra.fr}
}

\references{
Rathgeber C.B.K., Longuetaud F., Mothe F., Cuny H. & Le Moguedec G. 2011. Phenology of wood formation: Data processing, analysis and visualisation using R (Package CAVIAR).  \emph{Dendrochronologia} \bold{29}: 139-149. \cr

Rathgeber C.B.K., Rossi S.& Bontemps J.-D. 2011. Cambial activity related to tree size in a mature silver-fir plantation. \emph{Annals of Botany} \bold{108}: 429-438.
}

\seealso{
\code{\link{AMA2006}, \link{averageRadialFiles}, \link{computeBootstrapTest}, \link{computeCriticalDates}, \link{readXLSCellCountTable}, \link{fitGompertz}, \link{plotWoodFormationCalendar}, \link{plotWoodFormationDynamics}, \link{standardiseCellCount}}
}

\examples{
## Loading the example dataset:
data(AMA2006)

## Computing characteristic ring cell numbers for dominant firs and pines:
Fir.RCN <- computeRingCellNb(AMA2006[AMA2006$Sp=="ABAL" & AMA2006$Crown=="Dom", ])
Pine.RCN <- computeRingCellNb(AMA2006[AMA2006$Sp=="PISY", ])

}

\keyword{datagen}
