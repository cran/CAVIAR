\name{readXLSCellCountTable}

\alias{readXLSCellCountTable}

\title{readXLSCellCountTable}

\description{
Convert Excel spreadsheet file to database-like text file
}

\usage{
readXLSCellCountTable(xls.file.name, sheets, sampling.dates, index)
}


\arguments{
  \item{xls.file.name}{a character containing the path to the Excel file}
  \item{sheets}{a vector containing the numbers of the sheet to be read}
  \item{sampling.dates}{a data.frame containing the correspondence between sampling numbers and dates (in DOY)}
  \item{index}{a vector containing a key to identify records like the names of the trees}
}

\details{
This function takes the selected sheets of an Excel file and convert them into one single database-like data.frame. \cr
The function requires that the sheets contain 13 compulsory columns: Sample, C1, C2, C3, E1, E2, E3, L1, L2, L3, M1, M2, M3 (see \code{\link{AMA2006}}). \cr
Additionally the sheets can contain the 3 supplementary columns: P1, P2, P3 (see AMA2007). \cr
Other columns are discarded.  \cr
}

\value{
A data frame of at least 8 columns containing the values: Index, Sample, DY, RF, nC, nE, nL, nM \cr
}

\note{
Version: 1.1-2 \cr
Development started: 13 October 2010
Last modifications: 4 January 2010
}

\author{
Cyrille Rathgeber - LERFoB, UMR1092 - INRA Nancy (France) \cr
Email: \email{cyrille.rathgeber@nancy.inra.fr}
}

\references{

}

\seealso{
\code{\link{AMA2006}, \link{averageRadialFiles}, \link{computeBootstrapTest}, \link{computeCriticalDates}, \link{computeRingCellNb}, \link{fitGompertz}, \link{plotWoodFormationCalendar}, \link{plotWoodFormationDynamics}}
}

\examples{
## Setting the parameters of the function
library(gdata)
xls.file.name <- file.path(.path.package('CAVIAR'), 'data', 'AMA2006.xls')
sheets <- 1:10

DY <- seq(93, 303, 7)
Sample <- 1:31
sampling.dates <- data.frame(Sample, DY)

tree.names <- 1:10

## Converting the Excel file
converted.data.frame <- readXLSCellCountTable(xls.file.name=xls.file.name, sheets= sheets, sampling.dates=sampling.dates, index=tree.names)
}

\keyword{datagen}
