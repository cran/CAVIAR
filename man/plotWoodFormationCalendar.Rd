\name{plotWoodFormationCalendar}
\alias{plotWoodFormationCalendar}

\title{plotWoodFormationCalendar}

\description{
Plots wood formation calendar
}

\usage{
plotWoodFormationCalendar(data, main = " ", sub = " ", xlim=NULL, plotype = 4)
}

\arguments{
  \item{data}{a data frame of critical dates and durations with 8 compulsory columns located by their names (bE, bL, bM, cE, cL, dE, dL, dX), typically output from \code{\link{computeCriticalDates}} function}
  \item{main}{an optional character containing the plot title}
  \item{sub}{an optional character containing the plot subtitle}
  \item{xlim}{an optional vector indicating the x-axis minimum and maximum values}
  \item{plotype}{an optional integer indicating the type of plot
	1: individual critical dates plot,\cr
	2: group critical dates plot,\cr
	3: group critical durations plot,\cr
	4: group critical dates and durations plots (default).
  }
}

\details{
Critical dates can be first plotted for each individual tree (5 trees per page) in order to allow a visual checking of date consistencies using \option{plotype = 1}. \cr
Consistent dates can then be used to draw a wood formation calendar using \option{plotype = 2} for dates and \option{plotype = 3} for durations. \cr
Xylem development phase (enlarging: E, maturing: L, mature: M) onsets and cessations are represented by diamond-crossed-by-a-line marks. The left end of the line represents the minimum, the left end of the diamond the first quartile, the middle of the diamond the median, the right end of the diamond the third quartile and the right end of the line the maximum. \option{plotype = 4} plots dates and durations at the same time. \cr
NB. These representations are specially designed for five trees datasets. They also work when more trees are available, but then the graph represents a summary of the dataset and not all the data individually.
}

\value{
A data frame of the critical dates and durations median and median absolute deviation.
}

\note{
Version: 4.1 \cr
Last modifications: 03 January 2011
}

\author{
Cyrille Rathgeber - LERFoB, UMR1092 - INRA Nancy (France) \cr
Email: \email{cyrille.rathgeber@nancy.inra.fr}
}

\references{
Rathgeber C.B.K., Longuetaud F., Mothe F., Cuny H. & Le Moguedec G. 2011. Phenology of wood formation: Data processing, analysis and visualisation using R (Package CAVIAR).  \emph{Dendrochronologia} \bold{29}: 139-149. \cr
}

\seealso{
\code{\link{AMA2006}, \link{averageRadialFiles}, \link{computeBootstrapTest}, \link{computeCriticalDates}, \link{computeRingCellNb}, \link{readXLSCellCountTable}, \link{fitGompertz}, \link{plotWoodFormationDynamics}, \link{standardiseCellCount}}
}

\examples{
## Loading the AMA2006 dataset:
data(AMA2006)

## Computing wood formation critical dates for dominant firs and pines:
Fir.cdd <- computeCriticalDates(AMA2006[AMA2006$Sp=="ABAL" & AMA2006$Crown=="Dom", ])
Pine.cdd <- computeCriticalDates(AMA2006[AMA2006$Sp=="PISY", ])

# Computing medians and medians absolute deviation for firs critical dates and durations
Fir.medians <- plotWoodFormationCalendar(Fir.cdd)

## Visualising wood formation critical dates and durations on the screen for firs:
# Dates by trees individually
plotWoodFormationCalendar(Fir.cdd, plotype=1)
# Dates by group of trees
plotWoodFormationCalendar(Fir.cdd, plotype=2)
# Duration by group of trees
plotWoodFormationCalendar(Fir.cdd, plotype=3)

## Plotting wood formation calendar for firs and pines in a PDF file
# General settings
pdf(file="Wood formation calendar.pdf")
layout(matrix(c(1, 2, 3, 4), nrow=2, ncol=2, byrow=TRUE))
par(mar=c(4, 4, 4, 1))

# Drawing the first plot - Firs critical dates
plotWoodFormationCalendar(Fir.cdd, xlim= c(100,305) , sub="a. Fir critical dates", plotype=2)

# Drawing the second plot - Firs critical durations
plotWoodFormationCalendar(Fir.cdd, xlim= c(100, 210), sub="b. Fir critical durations", plotype=3)

# Drawing the third plot - Pines critical dates
plotWoodFormationCalendar(Pine.cdd, xlim= c(100,305), sub="c. Pine critical dates", plotype=2)

# Drawing the fourth plot - Pines critical durations
plotWoodFormationCalendar(Pine.cdd, xlim= c(100, 210), sub="d. Pine critical durations", plotype=3)

dev.off()
}

\keyword{univar}
\keyword{hplot}
