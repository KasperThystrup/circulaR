% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/S4_generics.R
\name{bam.file<-}
\alias{bam.file<-}
\title{Replacement method for the BAM file path}
\usage{
bam.file(object) <- value
}
\arguments{
\item{object}{circSample object}
}
\value{
New class object with replaced BAM file path
}
\description{
circSample class replacement methods for replacing the chimeric file path
The replacment value must be valid to enable replacement
}
\examples{
# Replacing BAM file path of a sample object
#bam.file(sample.object) <- new.path

}
