\name{medianFilter}
\alias{medianFilter}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ 1D Median Filter }
\description{
  1-dimensinal median filter with a specified \code{windowSize} of the \code{inputData}
}
\usage{
medianFilter(inputData, windowSize)
medianFilter(inputData = NULL, windowSize = 3)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{inputData}{ The vector of values to be filtered }
  \item{windowSize}{ The half-size of the filtering window (default \code{windowSize} = 3) }
}
\details{
  Classical implementation of a median filter, using a sliding window. By default, the half-size of the sliding window is set to 3 unless otherwise specified. 
}
\value{
  The output data has the same size of the input data. If the window half-size is smaller or equal to 1, then the input data is passed directly to the output data.
}
\references{ }
\author{ Adrian Andronache \email{adi.andronache@gmail.com} \cr
Luca Agnelli \email{luca.agnelli@gmail.com}
}
\note{
}
\seealso{  }
\examples{
x <- medianFilter(c(0, 0, 0, 1, 1, 1, 0, 0, 1, 0))
x <- medianFilter(c(0, 0, 0, 1, 0, 0, 0, 0, 1, 0), windowSize = 5)
}
