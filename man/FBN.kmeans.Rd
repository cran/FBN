\name{FBN.kmeans}
\alias{FBN.kmeans}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ K-Means clustering of SNP microarray data }
\description{
  Performs a k-means clustering of SNP microarray data. Returns clusters of values as being putatively characterized by different CN. 
}
\usage{
FBN.kmeans(inputData, minSpan, breaksData)
FBN.kmeans(inputData = NULL, minSpan = 0.2, breaksData = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{inputData}{ A vector of values containig the SNP microarray data }
  \item{minSpan}{ The minimum distance separating consecutive local maxima that are to be detected on 
  	the histogram of the \code{inputData}. These maxima are used to initialize the k-means clustering process.
  	For more details concerning the local maxima detection, check the documentation of \code{\link{FBN.histogramMaxima}}}
  \item{breaksData}{ One of: \itemize{
	\item{}{a vector giving the breakpoints between histogram cells,}
	\item{}{a single number giving the number of cells for the histogram,}
	\item{}{a character string naming an algorithm to compute the number of cells (see Details section of \code{\link{hist}}),}
	\item{}{a function to compute the number of cells.}
	}}
}
\details{
  This fuction takes as input the vector of raw SNP microarray values, and perform a k-means clustering 
  trying to identify the groups of raw values characterized by different CNs. The clustering process is initialized 
  with the local maxima detected on the histogram of the input data (see the documentation of \code{\link{FBN.histogramMaxima}}).
  To increase the robustness of the clustering process and to remove possible small or noisy clusters,
  a double filtering is done: firstly, removing those clusters populated by less than 1\% values from the entire \code{inputData}, and then, due to putative noisy histograms, merging those clusters whose centers are closer than 0.2 in nominal values.   
}
\value{
  An object of class \code{\link{kmeans}}
}
\author{ Adrian Andronache \email{adi.andronache@gmail.com} \cr
	Luca Agnelli \email{luca.agnelli@gmail.com}
}
\examples{
	
require(stats)
require(graphics)
x = c(rnorm(1000, 1, .2), rnorm(1000, 2, .2))
y = FBN.kmeans(x, minSpan = .001)
h = hist(x)
par(new = TRUE)
plot(y$centers,vector(mode=mode(y$centers), length = length(y$centers)), 
	xlim = c(min(h$breaks), max(h$breaks)), ylim = c(0,max(h$counts)), 
	xlab = NA, ylab = NA, col = 'red' )
	
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
