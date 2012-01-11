\name{averageRadialFiles}

\alias{averageRadialFiles}

\title{averageRadialFiles}

\description{
Average raw radial file data to compute mean cell numbers by trees
}

\usage{
averageRadialFiles(data)
}


\arguments{
  \item{data}{a data frame containing cell count data. This data frame must contain at least 7 imposed columns with compulsory names: Tree, DY, RF, nC, nE, nL, nM   (see \code{\link{AMA2006}})}
}

\details{
This function takes the raw cell count data (nC, NE, nL, nM and optionally P) recorded by radial file and averages them by trees. \cr
Additionally, it computes the sums of mature and maturing cells (nLM = nL + nM) and mature, maturing and enlarging cells (nELM = nE + nL + nM). \cr
}

\value{
A data frame of at least 8 columns containing the computed values: Tree, DY, nC, nE, nL, nM, nLM, nELM. \cr
}

\note{
Version: 2.1-1 \cr
Development started: 27 November 2008 \cr
Last modifications: 25 January 2011
}

\author{
Cyrille Rathgeber - LERFoB, UMR1092 - INRA Nancy (France) \cr
Email: \email{cyrille.rathgeber@nancy.inra.fr}
}

\references{
Rathgeber C.B.K., Longuetaud F., Mothe F., Cuny H. & Le Moguedec G. 2011. Phenology of wood formation: Data processing, analysis and visualisation using R (Package CAVIAR).  \emph{Dendrochronologia} \bold{29}: 139-149. \cr
}

\seealso{
\code{\link{AMA2006}, \link{computeBootstrapTest}, \link{computeCriticalDates}, \link{computeRingCellNb}, \link{readXLSCellCountTable}, \link{fitGompertz}, \link{plotWoodFormationCalendar}, \link{plotWoodFormationDynamics}, \link{standardiseCellCount}}
}

\examples{
## Loading the example dataset:
data(AMA2006)

## Averaging raw data:
AMA2006.mean <- averageRadialFiles(AMA2006)
}

\keyword{datagen}
