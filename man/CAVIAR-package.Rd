\name{CAVIAR-package}
\alias{CAVIAR-package}
\alias{CAVIAR}
\docType{package}

\title{Cambial activity and wood formation: data processing, visualisation and analysis using R}

\description{
Process, analyses and plots data from cambial activity and wood formation monitoring studies
}

\details{
\tabular{ll}{
Package: \tab CAVIAR \cr
Type: \tab Package \cr
Version: \tab 1.3-1 \cr
Date: \tab 2010-10-22 \cr
License: \tab GPL (>= 2) \cr
LazyLoad: \tab yes \cr
URL: \tab \url{http://www.nancy.inra.fr/foret_bois_lerfob/personnels_de_l_umr/scientifiques/rathgeber_cyrille}
}

CAVIAR is a package that allow processing, visualisation and analysis of data coming from cambial activity and wood formation monitoring studies. \cr
The developed functions require input data files as semicolon separated value text files (i.e. CSV format). \cr
The columns names are imposed but not the order. That means that the input table may contain a free number of columns, as long as the compulsory columns, identified by their headers, are provided to the function. \cr
Missing values are coded "NA" in the data table.
}

\author{
Cyrille Rathgeber - LERFoB, UMR1092 - INRA Nancy (France) \cr
Email: \email{cyrille.rathgeber@nancy.inra.fr}
}

\references{
Rathgeber C.B.K., Longuetaud F., Mothe F., Cuny H., Le Moguedec G. 2010. Phenology of wood formation: data processing, analysis and visualisation using R. Accepted in \emph{Dendrochronologia}.
}

\seealso{
\code{\link{AMA2006}, \link{averageRadialFiles}, \link{computeBootstrapTest}, \link{computeCriticalDates}, \link{computeRingCellNb}, \link{readXLSCellCountTable}, \link{fitGompertz}, \link{plotWoodFormationCalendar}, \link{plotWoodFormationDynamics}}
}

\note{
Development started: 27 Novembre 2008
Last modifications: 06 Decembre 2010
}

\examples{}

\keyword{package}

