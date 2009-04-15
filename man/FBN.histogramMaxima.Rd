\name{FBN.histogramMaxima}
\alias{FBN.histogramMaxima}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Local maxima of the histogram of SNP microarray data }
\description{
  Finds at most 6 of the local maxima of the histogram of the \code{inputData} representing the SNP microarray data. 
  Firstly, the function finds all local maxima onto the histogram, 
  and finally removes those that are closer than \code{minSpan}. 
  The histogram is estimated with equi-spaced breaks (also the default), 
  defined by \code{breaksData} (see the documentation of \code{\link{hist}} for more details). 
}
\usage{
FBN.histogramMaxima(inputData, minSpan, breaksData)
FBN.histogramMaxima(inputData, minSpan = 0.1, breaksData = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{inputData}{ A vector of SNP microarray values for which the local maxima of the histogram are desired }
  \item{minSpan}{ The minimum distance separating consecutive local maxima. If \code{minSpan} is negative, then a smoothing filter is applied on the histogram to reduce its noise, and therefore, estimates the most significative maxima. }
  \item{breaksData}{ One of: \itemize{ 
	\item{}{a vector giving the breakpoints between histogram cells,}
	\item{}{a single number giving the number of cells for the histogram,}
	\item{}{a character string naming an algorithm to compute the number of cells (see the Details section of \code{\link{hist}}),}
	\item{}{a function to compute the number of cells.}
	}}
}
\details{
 This function has been designed based on SNP microarray data and FISH resolution. As in FISH analyses high numbers of signals do not allow a correct identification of discrete CNs, an empirical number of no more than 6 local maxima is therefore imposed. These maxima are used to initialize the k-means algorithm to determine the CN clusters. 
}
\value{
	Returns a vector containing at most 6 values of the \code{inputData} in which the histogram shows local maxima.
}
\references{  }
\author{ Adrian Andronache \email{adi.andronache@gmail.com} \cr
Luca Agnelli \email{luca.agnelli@gmail.com}
}
\note{  

}
\seealso{\code{\link{hist}}, \code{\link{FBNormalization}} }
\examples{
require(stats)
require(graphics)

x = c(rnorm(1000, 1, .2), rnorm(1000, 2, .2))
y = FBN.histogramMaxima(x, minSpan = .1)
h = hist(x)
par(new = TRUE)
plot(y,vector(mode=mode(y), length = length(y)), xlim = c(min(h$breaks), 
 	max(h$breaks)), ylim = c(0,max(h$counts)), xlab = NA, ylab = NA, 
 	col = 'red' )


x = c(1, 2, 2, 3, 4, 5, 5, 5, 6, 7, 8, 9, 10, 10, 10, 11)
y = FBN.histogramMaxima(x, minSpan = 3, breaksData = 100)
h = hist(x, 100)
par(new = TRUE)
plot(y,vector(mode=mode(y), length = length(y)), xlim = c(min(h$breaks), 
 	max(h$breaks)), ylim = c(0,max(h$counts)), xlab = NA, ylab = NA, 
 	col = 'red' )
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
