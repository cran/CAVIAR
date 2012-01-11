\name{AMA2006}
\docType{data}
\alias{AMA2006}

\title{AMA2006 cell count table}

\description{
This dataset was obtained from a cambial activity and wood formation monitoring study. Data were collected weekly in 2006 for 15 silver firs and 5 Scots pines grown in the same stand, located in Amance forest near Nancy (northeast France). \cr
Column "Tree" indicates the tree reference number, "Sp" indicates the tree species code, "Crown" Indicates tree crown class, "Sample" indicates the sampling campaign reference number, "DY" indicates the date of the sampling in day of year, RF indicates the radial file reference number. "nC" indicates the number of cells in the cambial zone, "nE" indicates the number of cells in the enlarging zone, "nL" indicates the number of cells in the maturing zone and "nM" indicates the number of mature cells. "P" is mesure (in number of cells, expression(mu) or mm) of the previous tree-ring width.
}

\usage{AMA2006}

\format{A data frame containing 11 columns and 1775 rows}

\source{
Cyrille Rathgeber - LERFoB, UMR1092 - INRA Nancy (France) \cr
Email: \email{cyrille.rathgeber@nancy.inra.fr}
}

\references{
Rathgeber C.B.K., Longuetaud F., Mothe F., Cuny H. & Le Moguedec G. 2011. Phenology of wood formation: Data processing, analysis and visualisation using R (Package CAVIAR).  \emph{Dendrochronologia} \bold{29}: 139-149. \cr

Rathgeber C.B.K., Rossi S.& Bontemps J.-D. 2011. Cambial activity related to tree size in a mature silver-fir plantation. \emph{Annals of Botany} \bold{108}: 429-438.
}

\seealso{
\code{\link{averageRadialFiles}, \link{computeBootstrapTest}, \link{computeCriticalDates}, \link{computeRingCellNb}, \link{readXLSCellCountTable}, \link{fitGompertz}, \link{plotWoodFormationCalendar}, \link{plotWoodFormationDynamics}, \link{standardiseCellCount}}
}

\note{
Version: 4.1-1 \cr
Development started: 4 March 2007
Last modifications: 13 Janvier 2012
}

\keyword{datasets}
