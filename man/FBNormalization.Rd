\name{FBNormalization}
\alias{FBNormalization}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ FISH Based Normalization and Copy Number inference of 
	SNP microarray data }
\description{
  Normalizes the data from a file containing the raw values of the SNP probes of microarrray data
  by using the FISH probes and their corresponding CNs.
}
\usage{
FBNormalization(rawDataFileName, fishProbesFileName, 
	normDataFileName, debugFlag, plotFlag, plotAndSaveFlag)
FBNormalization(rawDataFileName = NULL, fishProbesFileName = NULL, 
	normDataFileName = NULL, debugFlag = FALSE, plotFlag = FALSE, plotAndSaveFlag = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rawDataFileName}{ The file containig the raw values of the SNP probes. It should be .txt format, 
  	tab delimited, containing on each column the SNP probes data of the different samples. 
  	The structure of the input file should be as follows: \itemize{
   		\item{}{ a column containing the identification name of the SNP probes, entitled 'SNP',}
  		\item{}{ a column containing the chromosome number of the SNP probes, entitled 'chrom',}
  		\item{}{ a column containing the physical positon of the SNP probes, entitled 'physical.position',}
  		\item{}{ a column containing the cytoband of the SNP probe, entitled 'cytoband',}
  		\item{}{ 'n' columns containing the values of the SNP probes, entitled with the custom identification name of the samples.}
  	}}
  \item{fishProbesFileName}{ The file containing the FISH probes information. It should be .txt format, 
  	tab delimited, containing a separate row for each available FISH probe and on each column the various 
  	CN data revealed by FISH of the different samples. The structure of the input file should be as follows: \itemize{
		\item{}{ each row should represent one FISH probe}
		\item{}{ one column containing the name of the BAC clone, entitled 'BACclone';}
		\item{}{ one column containing the chromosome number where the FISH probe is located, entitled 'chromosome';}
		\item{}{ one column containing the cytoband location of the Fish probe, entitled 'cytoband';}
		\item{}{ one column containing the start position of the physical localization of the FISH probe, entitled 'start.loc';}
		\item{}{ one column containing the end position of the physical localization of the FISH probe, entitled 'end.loc';}
		\item{}{ 'n' columns containing the values of the CN coresponding to the FISH probe, entitled with the custom identification name of the samples.}
	}}
  \item{normDataFileName}{ The name of the output file where the normalized data will be written. Also, the function writes
  	a file starting with the same name as \code{normDataFileName} and ending with \emph{'\_thresholds.txt'}, containing the
  	thresholds values that can be used to infer CNs to the normalized data. }
  \item{debugFlag}{ Logical value, specifying wether the function should run in debug mode. If \code{TRUE}, it plots the various histograms and the detected centers of the clusters in the different samples. }
  \item{plotFlag}{ Logical value, specifying wether the function should plot the histograms of the individual data. If \code{TRUE}, it plots the various histograms and the detected centers of the clusters in the different samples. }
  \item{plotAndSaveFlag}{ Logical value, specifying wether the function should plot and save the various plots of the histograms. }
}
\details{
	For further detalis, see Supporting Information in \cite{Agnelli L et al. (2009), "A SNP Microarray and FISH-Based Procedure to Detect Allelic Imbalances in Multiple Myeloma: an Integrated Genomics Approach Reveals a Wide Gene Dosage Effect", \emph{Genes Chrom Cancer}}
}
\value{
	No value is returned. The function writes on the disk the output files as previously described. 
}
\author{ Adrian Andronache \email{adi.andronache@gmail.com} \cr
Luca Agnelli \email{luca.agnelli@gmail.com}
}
\examples{

## set path to FBN package data directory
rawDataFileName = './../data/hmcls.txt'
fishProbesFileName = './../data/FISHprobes.txt'
normDataFileName= 'hmcls_NORM.txt'
FBNormalization(rawDataFileName, fishProbesFileName, normDataFileName, 
	debugFlag = FALSE)
}
