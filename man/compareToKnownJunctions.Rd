% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/S4_generics.R
\name{compareToKnownJunctions}
\alias{compareToKnownJunctions}
\title{Compare backsplice junction identified by STAR to junction sites involved in linear splicing}
\usage{
compareToKnownJunctions(object, known.junctions, cores = 1L, ...)
}
\arguments{
\item{object}{circSample or circExperiment object}

\item{known.junctions}{A vector of (valid) chromosome names that is used to restrict the chimeric read data}

\item{cores}{Only supported on unix! Number of cores used for parallel processing for circExperiment objets.}
}
\value{
List of linear splice junction count summaries (SJ.tab.out) for all samples stored in the circExperiment-object.
}
\description{
The function takes a circSample or circExperiment object and compares backsplice read data for circSample(s) to a database of known junction sites (generated by \code{\link{constructSJDB,ANY-method}}).
}
\examples{
#kj <- constructSJDB(annotationDB = db)
#readBSJdata(object, known.junctions = kj, cores = 4)

}
