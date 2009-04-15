\name{FBN.valueCenter}
\alias{FBN.valueCenter}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ The SNP normalization function }
\description{
  Normalization of the raw SNP microarray values, by multiplication (on linear scale) or addition (in log scale) 
  of all the raw SNP values with the normalization factor. 
  The normalization factor is estimated such that it brings the \code{normalizingValue} of the raw 
  SNP values onto the \code{nominalValueCN}. 
}
\usage{
FBN.valueCenter(inputData, normalizingValue, nominalValueCN, 
	logScale)
FBN.valueCenter(inputData = NULL, normalizingValue = NULL, 
	nominalValueCN = 2, logScale = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{inputData}{ The vector of raw SNP values, as they come out from, e.g. Circular Binary Segmentation in \code{DNAcopy} package from Bioconductor  }
  \item{normalizingValue}{ The value representing the center of the cluster identified as having a certain CN }
  \item{nominalValueCN}{ The nominal value representing a certain CN on which the \code{normalizingValue} has to be brought. }
  \item{logScale}{ A logical value, specifying wether the data is on linear (\code{FALSE}) or logarithmic scale (\code{TRUE}).  }
}
\details{
  The \code{nominalValueCN} is a real value representing the CN, e.g. CN\eqn{=2} has a \code{nominalValueCN} of 2, 
  but all other CN\eqn{=n} (\eqn{n} != 2) will have a \code{nominalValueCN} different from \eqn{n}. 
  Such \code{nominalValueCN} is identified by the \code{\link{FBN.kmeans}} function. 
}
\value{
	Returns a vector containing the normalized values of the \code{inputData}
}
\references{ }
\author{ Adrian Andronache \email{adi.andronache@gmail.com} \cr
Luca Agnelli \email{luca.agnelli@gmail.com}
}
\note{ 
}
\seealso{ \code{\link{FBN.kmeans}}, \code{\link{FBNormalization}}  }
\examples{
require(stats)
require(graphics)
x = c(rnorm(1000, 1, .1), rnorm(1000, 1.5, .1))
y = FBN.valueCenter(x, normalizingValue = 1, nominalValueCN = 2, 
	logScale = FALSE)
par(mfrow = c(2, 1), new = FALSE)
h = hist(x)
par(new = TRUE)
plot(1, 0, col = 'red', xlim = c(min(h$breaks), max(h$breaks)), 
	ylim = c(0,max(h$counts)), xlab = NA, ylab = NA)
par(new = FALSE)
h = hist(y)
par(new = TRUE)
plot(2, 0, col = 'red', xlim = c(min(h$breaks), max(h$breaks)), 
	ylim = c(0,max(h$counts)), xlab = NA, ylab = NA)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ FBN }
%\keyword{ FBN.valueCenter }% __ONLY ONE__ keyword per line
