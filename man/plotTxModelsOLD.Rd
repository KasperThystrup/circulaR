% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_functions.r
\name{plotTxModelsOLD}
\alias{plotTxModelsOLD}
\title{Plot transcripts stored in a TxDB.}
\usage{
plotTxModelsOLD(
  db = NULL,
  g = NULL,
  tid = NULL,
  addChevrons = T,
  rangex = NULL
)
}
\arguments{
\item{db}{A TxDb containing transcript models from Ensembl.}

\item{g}{If the TxDb contains more than one gene, then use g to specify what to plot.}

\item{addChevrons}{Boolean indicating wheter to add chevrons to intronic regions to indicate direction of the transcript.}

\item{rangex}{The range in genomic coordinates for the region of interest. This is important for correct alignment of the tx model plot with splice junction plot.}

\item{ex}{GrangesList with exonsBy tx.}

\item{int}{GrangesList with intronsByTranscript information.}
}
\value{
Don't know?
}
\description{
Plot the transcript models of a gene.
}
\examples{
ensembl <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
annot <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name", "ensembl_transcript_id", "external_transcript_name", "transcript_tsl", "chromosome_name", "start_position", "end_position", "strand"),
  filters = "external_gene_name",
  values = "FIRRE",
  mart = ensembl
)
txdb <- makeTxDbFromBiomart(transcript_ids = annot$ensembl_transcript_id)#, host = "uswest.ensembl.org")
plotTranscripts(tdb = txdb)
Plot transcript models.

Using data stored in a TxDb to visualize the transcript models.

}
