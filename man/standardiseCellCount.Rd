\name{standardiseCellCount}

\alias{standardiseCellCount}

\title{standardiseCellCount}

\description{
Standardise average cell count values (nC, nE, nL, nM) using previous ring cell count or measurements (P)
}

\usage{
standardiseCellCount(data, na.rm="FALSE")
}


\arguments{
  \item{data}{a data frame of cell counting with 6 compulsory columns (Tree, DY, nC, nE, nL, nM, P) located by their names, typically output from \code{\link{averageRadialFiles}} function}
  \item{na.rm}{logical indicating if lines where P is missing should be left as they are without applying any standardisation (FALSE by default) or should be removed from the record (TRUE)}
}

\details{
This function standardizes average cell count values (nC, nE, nL, nM) using previous ring cell count or measurements (P). ! Insert formula here ! \cr
The function requires that the input data frame contains 6 compulsory columns: Tree, DY, nC, nE, nL, nM and P, typically output from \code{\link{averageRadialFiles}} function. \cr
}

\value{
A data frame similar to the input data frame but containing standardized values for the columns: nC, nE, nL, nM \cr
}

\note{
Version: 1.0-0 \cr
Development started: 25 January 2011
Last modifications: 26 January 2011
}

\author{
Cyrille Rathgeber - LERFoB, UMR1092 - INRA Nancy (France) \cr
Email: \email{cyrille.rathgeber@nancy.inra.fr}
}

\references{
Rathgeber C.B.K., Rossi S.& Bontemps J.-D. 2011. Cambial activity related to tree size in a mature silver-fir plantation. \emph{Annals of Botany} \bold{108}: 429-438.
}

\seealso{
\code{\link{AMA2006}, \link{averageRadialFiles}, \link{computeBootstrapTest}, \link{computeCriticalDates}, \link{computeRingCellNb}, \link{fitGompertz}, \link{plotWoodFormationCalendar}, \link{plotWoodFormationDynamics}, \link{readXLSCellCountTable}}
}

\examples{
## Loading the example dataset:
data(AMA2006)

## Averaging raw data for dominant trees:
AMA2006.mean <- averageRadialFiles(AMA2006[AMA2006$Crown=="Dom", ])

## Standardising data:
AMA2006.std <- standardiseCellCount(AMA2006.mean)
}

\keyword{datagen}
