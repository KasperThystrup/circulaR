% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/circRNA_annotation.r
\name{addKnownJunctions}
\alias{addKnownJunctions}
\title{Function to add information on annotated junctions and transcript termini to candidate backsplice sites}
\usage{
addKnownJunctions(cbs, kj)
}
\arguments{
\item{cbs}{Table of candidate backspliced reads (produced by getCandidateBackspliceSites).}

\item{kj}{GRanges object of known junctions (produced by constructSJDB).}
}
\value{
Table of candidate backspliced reads with information on closest donor and acceptor.
}
\description{
This function takes a table of candidate backspliced reads and finds the closest known junction.
}
\examples{
\dontrun{
# Build database of known junctions
ah <- AnnotationHub()
ahdb <- query(ah, pattern = c("Homo Sapiens", "EnsDb", 92))[[1]]
sjdb <- constructSJDB(txdb)

df <- readChimFile(file = "~/projects/circulaR/misc/sub.Chimeric.out.junction") ## Replace with 'data(ExampleJunctions)'
df <- addKnownJunctions(cbs = df, kj = sjdb)
}
}
