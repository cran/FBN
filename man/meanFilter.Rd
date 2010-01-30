\name{meanFilter}
\alias{meanFilter}
\title{ 1D Mean Filter }
\description{
1-dimensional mean filter with a specified \code{windowSize} of the \code{inputData}}
\usage{
meanFilter(inputData, windowSize)
meanFilter(inputData = NULL, windowSize = 3)
}
\arguments{
  \item{inputData}{ The vector of values to be filtered }
  \item{windowSize}{ The half-size of the filtering window (default \code{windowSize} = 3) }
}
\details{
  Classical implementation of a mean filter, using a sliding window. By default, the half-size of the sliding window is set to 3 unless otherwise specified. 
}
\value{
  The output data has the same size of the input data. If the window half-size is smaller or equal to 1, then the input data is passed directly to the output data.
}
\author{ Adrian Andronache \email{adi.andronache@gmail.com} \cr
Luca Agnelli \email{luca.agnelli@gmail.com}
}
\examples{
x <- meanFilter(c(0, 0, 0, 1, 1, 1, 0, 0, 1, 0))
x <- meanFilter(c(0, 0, 0, 1, 0, 0, 0, 0, 1, 0), windowSize = 5)
}
