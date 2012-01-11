\name{readXLSCellCountTable}

\alias{readXLSCellCountTable}

\title{readXLSCellCountTable}

\description{
Convert Excel spreadsheet file to database-like text file
}

\usage{
readXLSCellCountTable(xls.file.name, sheets, sampling.dates, tree.names)
}


\arguments{
  \item{xls.file.name}{a character containing the path to the Excel file}
  \item{sheets}{a vector containing the numbers of the sheet to be read}
  \item{sampling.dates}{a data.frame containing the correspondence between sampling numbers and dates (in DOY)}
  \item{tree.names}{a vector containing a key to identify records like the names of the trees}
}

\details{
This function takes the selected sheets of an Excel file and convert them into one single database-like data.frame. \cr
The function requires that the sheets contain 13 compulsory columns: Sample, C1, C2, C3, E1, E2, E3, L1, L2, L3, M1, M2, M3. \cr
Additionally the sheets can contain the 3 supplementary columns: P1, P2, P3. \cr
Other columns are discarded.  \cr
}

\value{
A data frame of at least 8 columns containing the values: Tree, Sample, DY, RF, nC, nE, nL, nM, and optionally P. \cr
}

\note{
Version: 1.3 \cr
Development started: 13 October 2010
Last modifications: 12 January 2012
}

\author{
Cyrille Rathgeber - LERFoB, UMR1092 - INRA Nancy (France) \cr
Email: \email{cyrille.rathgeber@nancy.inra.fr}
}

\references{
Rathgeber C.B.K., Longuetaud F., Mothe F., Cuny H. & Le Moguedec G. 2011. Phenology of wood formation: Data processing, analysis and visualisation using R (Package CAVIAR).  \emph{Dendrochronologia} \bold{29}: 139-149. \cr
}

\seealso{
\code{\link{AMA2006}, \link{averageRadialFiles}, \link{computeBootstrapTest}, \link{computeCriticalDates}, \link{computeRingCellNb}, \link{fitGompertz}, \link{plotWoodFormationCalendar}, \link{plotWoodFormationDynamics}, \link{standardiseCellCount}}
}

\examples{
## Setting the parameters of the function
library(gdata)

xls.file.name <- file.path(.path.package('CAVIAR'), 'extdata', 'AMA2006.xls')

sheets <- 2:21
DY <- c(seq(82, 320, 7), 340)
Sample <- 1:36
sampling.dates <- data.frame(Sample, DY)

tree.names <- 1:20

## Converting the Excel file
converted.data.frame <- readXLSCellCountTable(xls.file.name=xls.file.name, sheets= sheets, sampling.dates=sampling.dates, tree.names=tree.names)
}

\keyword{datagen}
