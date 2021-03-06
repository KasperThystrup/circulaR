% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/S4_generics.R
\name{gb}
\alias{gb}
\title{Accessor for genome build}
\usage{
gb(object)
}
\arguments{
\item{object}{circSample or circExperiment object}
}
\value{
Character vector with the genome build(s)
}
\description{
circSample and circExperiment class accessors for showing the genome build(s) of circSample object(s)
}
\examples{
# The genome build of a sample object
#gb(sample.object)

# The genome build in each sample of an experiment
#gb(experiment.object)

}
