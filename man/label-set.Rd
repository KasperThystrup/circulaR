% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/S4_generics.R
\name{label<-}
\alias{label<-}
\title{Replacement method for sample label}
\usage{
label(object) <- value
}
\arguments{
\item{object}{circSample or circExperiment object}
}
\value{
New class object with replaced sample label(s)
}
\description{
circSample and circExperiment class replacement methods for replacing the sample label(s)
The replacment value must be valid to enable replacement
}
\examples{
# Replacing label of a sample object
#label(sample.object) <- .new.label

# Replacing label in each sample of an experiment
#label(experiment.object) <- new.label.vector

}
