\name{plotWoodFormationDynamics}
\alias{plotWoodFormationDynamics}

\title{plotWoodFormationDynamics}

\description{
Plots intra-annual dynamics of wood formation
}

\usage{
plotWoodFormationDynamics(data, main = " ", sub=" ", x.axis.labels = "DOY", wrap = FALSE)
}

\arguments{
  \item{data}{a data frame of cell counting with 4 compulsory columns (nC, nE, nL, nM) and 8 optional columns (nC.inf, nC.sup, nE.inf, nE.sup, nL.inf, nL.sup, nM.inf, nM.sup) located by their names, typically output from \code{\link{averageRadialFiles}} function}
  \item{main}{an optional character containing the plot title}
  \item{sub}{an optional character containing the plot subtitle}
  \item{x.axis.labels}{an optional character for selecting between day of year (DOY, by default) or sample number (Sample)}
  \item{wrap}{an optional logical indicating if an envelope must be drawn around the main curves (if the 8 optional columns are provided!)}
}

\details{
Cell counting can be plotted for each individual tree in order to allow a visual checking of data consistencies using \option{wrap = FALSE}. \cr
Cell counting can then be plotted for a population of tree taking advantage of the \option{wrap = TRUE}. \cr
}

\value{
none
}

\note{
Version: 1.1-1 \cr
Development started: 27 November 2008 \cr
Last modifications: 4 November 2010
}

\author{
Cyrille Rathgeber - LERFoB, UMR1092 - INRA Nancy (France) \cr
Email: \email{cyrille.rathgeber@nancy.inra.fr}
}

\references{
Rathgeber C.B.K., Longuetaud F., Mothe F., Cuny H. & Le Moguedec G. 2011. Phenology of wood formation: Data processing, analysis and visualisation using R (Package CAVIAR).  \emph{Dendrochronologia} \bold{29}: 139-149. \cr
}

\seealso{
\code{\link{AMA2006}, \link{averageRadialFiles}, \link{computeBootstrapTest}, \link{computeCriticalDates}, \link{computeRingCellNb}, \link{readXLSCellCountTable}, \link{fitGompertz}, \link{plotWoodFormationCalendar}, \link{standardiseCellCount}}
}

\examples{
## Loading the AMA2006 dataset:
data(AMA2006)

## Averaging raw data for dominant trees:
AMA2006.mean <- averageRadialFiles(AMA2006[AMA2006$Crown=="Dom", ])

## Visualising intra-annual dynamics of wood formation:
## Plotting data for one tree
plotWoodFormationDynamics(AMA2006.mean[AMA2006$Tree==1, ])

## Plotting data for a group of trees
# Computing the median
FirMed <- aggregate(AMA2006.mean[, c("nC", "nE", "nL", "nM")], by=list(AMA2006.mean$DY), FUN=median, na.rm=TRUE)
names(FirMed) <- c("DY", "nC", "nE", "nL", "nM")

# Computing the first quartile
FirInf <- aggregate(AMA2006.mean[, c("nC", "nE", "nL", "nM")], by=list(AMA2006.mean$DY), FUN=quantile, probs=0.25, na.rm=TRUE)
names(FirInf) <- c("DY", "nC.inf", "nE.inf", "nL.inf", "nM.inf")

# Computing the third quartile
FirSup <- aggregate(AMA2006.mean[, c("nC", "nE", "nL", "nM")], by=list(AMA2006.mean$DY), FUN=quantile, probs=0.75, na.rm=TRUE)
names(FirSup) <- c("DY", "nC.sup", "nE.sup", "nL.sup", "nM.sup")

# Creating the final data frame
Temp <- merge(FirMed, FirInf)
Fir <- merge(Temp, FirSup)

plotWoodFormationDynamics(Fir, main = "Firs, Amance Forest, year 2006", wrap = TRUE)
}

\keyword{hplot}
