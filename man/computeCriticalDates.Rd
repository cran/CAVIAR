\name{computeCriticalDates}
\alias{computeCriticalDates}

\title{computeCriticalDates}

\description{
Computes wood formation critical dates and durations and their standard deviations
}

\usage{
computeCriticalDates(data, plot=TRUE)
}


\arguments{
  \item{data}{a data frame containing cell count data. This data frame must contain at least 7 imposed columns with compulsory names: Tree, DY, RF, nC, nE, nL, nM (see \code{\link{AMA2006}})}
  \item{plot}{an optional logical indicating if a PDF plot must be outputted or not (default=TRUE)}
}

\details{
This function computes wood formation critical dates and durations from a raw data table containing cell count (for each radial file) for the 4 development phases (cambial, enlarging, maturing & mature) per sampling dates (see AMA2006). \cr
Critical dates as well as their associated standard deviations are computed using logistic regressions \code{\link{glm}}, the main date corresponds to the 50\% probability of the phase being active while the confidence interval limits corresponds to the 2.5 and 97.5\% probabilities. The standard deviation is obtained by dividing the 95\% confidence interval by 2. \cr 
From these critical dates, the function computes 3 critical durations with their associated standard deviations (assuming independence between dates). \cr
By default (\option{plot=TRUE}), the function outputs a PDF plot entitled \file{Critical dates computation report... .pdf}, illustrating the computation of each critical date for a quick visual checking.
}

\value{
A data frame containing the computed critical dates and durations with their associated standard deviations. \cr
This data frame is composed of 5 columns for the critical dates, which are the beginning of enlarging phase (bE), the beginning of maturing phase (bL), the beginning of mature phase (bM), the cessation of the enlarging phase (cE), the cessation of the maturing phase (cL); 3 columns for the critical durations, which are: the duration of the enlarging phase (dE = cE - bE), the duration of the maturing phase (dL = cL - bL), and the total duration of the xylogenesis (dX = cL - bE); as well as 8 columns for their associated standard deviations noted xX.sd (e.g. bE.sd, cE.sd, dE.sd).
}

\note{
Version: 4.1-1 \cr
Development started: 21 October 2009 \cr
Last modifications: 21 December 2012
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
\code{\link{AMA2006}, \link{averageRadialFiles}, \link{computeBootstrapTest}, \link{computeRingCellNb}, \link{readXLSCellCountTable}, \link{fitGompertz}, \link{plotWoodFormationCalendar}, \link{plotWoodFormationDynamics}, \link{standardiseCellCount}}
}

\examples{
## Loading the example dataset:
data(AMA2006)

## Computing wood formation critical dates and durations for dominant firs and pines:
Fir.cdd <- computeCriticalDates(AMA2006[AMA2006$Sp=="ABAL" & AMA2006$Crown=="Dom", ])
Pine.cdd <- computeCriticalDates(AMA2006[AMA2006$Sp=="PISY", ])
}

\keyword{datagen}
\keyword{hplot}
