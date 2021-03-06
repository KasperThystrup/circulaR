% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/S4_generics.R
\name{bsj.reads}
\alias{bsj.reads}
\title{Accessor for the backsplice junction read data}
\usage{
bsj.reads(object, ...)
}
\arguments{
\item{object}{circSample or circExperiment object}

\item{returnAs}{(only valid for circSample) Whether to return data as a list of samples or compiled into a single long table with a sample.id column}
}
\value{
Datatable(s) of filtered backsplice junction read data
}
\description{
circSample and circExperiment class accessors for showing the filtered backsplice junction read data of a circSample object(s)
}
\examples{
# A datatable of the filtered backsplice junction read data in a sample object
#bsj.reads(sample.object)

# A datatable of the filtered backsplice junction read data in each sample of an experiment
#bsj.reads(experiment.object)

}
