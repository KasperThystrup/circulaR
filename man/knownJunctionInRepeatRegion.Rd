% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helper_functions.r
\name{knownJunctionInRepeatRegion}
\alias{knownJunctionInRepeatRegion}
\title{Function to identify chimeric reads eligible for shifting}
\usage{
knownJunctionInRepeatRegion(cbs, junc = "donor")
}
\arguments{
\item{cbs}{Tibble of candidate backspliced reads}
}
\value{
Boolean vector
}
\description{
Function to identify chimeric reads eligible for shifting
}
