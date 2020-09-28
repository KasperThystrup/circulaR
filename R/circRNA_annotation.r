#' Function to add information on annotated junctions and transcript termini to candidate backsplice sites
#'
#' This function takes a table of candidate backspliced reads and finds the closest known junction.
#' @param cbs Table of candidate backspliced reads (produced by getCandidateBackspliceSites).
#' @param kj GRanges object of known junctions (produced by constructSJDB).
#'
#' @return Table of candidate backspliced reads with information on closest donor and acceptor.
#' @examples
#' \dontrun{
#' # Build database of known junctions
#' ah <- AnnotationHub()
#' ahdb <- query(ah, pattern = c("Homo Sapiens", "EnsDb", 92))[[1]]
#' sjdb <- constructSJDB(txdb)
#'
#' df <- readChimFile(file = "~/projects/circulaR/misc/sub.Chimeric.out.junction") ## Replace with 'data(ExampleJunctions)'
#' df <- addKnownJunctions(cbs = df, kj = sjdb)
#' }
#' @importFrom dplyr filter bind_cols bind_rows inner_join left_join right_join
#' @importFrom GenomicRanges GRanges nearest
#' @importFrom magrittr %>%
#' @importFrom tibble as_tibble
#' @export
addKnownJunctions <- function(cbs, kj){
  # # Checks
  # if(!is.null(err <- checkVariables(obj = cbs, expect_class = c("tbl_df", "tbl", "data.frame") , vari = "cbs"))) stop(err)
  # if(!is.null(err <- checkVariables(obj = kj, expect_class = "GRanges", vari = "kj"))) stop(err)
  # if(!all(colnames(mcols(kj)) %in% c("type", "g_id", "tx_id", "jID"))){stop("Known junctions appears to have the wrong format.")}
  cbs$queryHits <- 1:nrow(cbs)

  # Make GRanges objects of identified donors and acceptors
  donors <- junctionsToGRanges(seqnames = cbs$X1, start = cbs$X2, width = 1, strand = cbs$X3)
  acceptors <- junctionsToGRanges(seqnames = cbs$X4, start = cbs$X5, width = 1, strand = cbs$X6)

  nearestDonors <- findNearest(juncCand = donors, kj, type = "donor", cbs)
  nearestAcceptors <- findNearest(juncCand = acceptors, kj, type = "acceptor", cbs)

  # Determining ambiguous hits
  ambiguousHits <- determineAmbiguousHits(acceptors = nearestAcceptors, donors = nearestDonors)

  totalReads <- length(unique(nearestAcceptors$queryHits))
  ambi <- length(unique(ambiguousHits))
  message(" Ambiguous: ", ambi, " (", round(ambi*100/totalReads, digits = 2),"%)")

  if (ambi == 0) {
    output <- cbs %>% left_join(nearestDonors, by = "queryHits") %>% left_join(nearestAcceptors, by = "queryHits")
  } else {
    # Trying to resolve ambiguities
    resolvedHits <- rescueAmbiguous(acceptors = nearestAcceptors, donors = nearestDonors, ambiguousHits)

    res <- nrow(resolvedHits)
    message(" Rescued: ", res, " (", round(res*100/ambi, digits = 2),"%)")

    output <- rbind(
      subset(cbs, !queryHits %in% ambiguousHits) %>%
        left_join(nearestDonors, by = "queryHits") %>%
        left_join(nearestAcceptors, by = "queryHits"),
      right_join(cbs, resolvedHits, by = "queryHits")
    )
  }

  return(output[,!grepl("queryHits|subjectHits", colnames(output))])
}


#' Annotate circRNAs (bsID) by identifying overlapping genes
#'
#' @param bsid Vector of bsIDs to annotate.
#' @param TXDB Annotation database (from AnnotationHub!!!)
#' @param t Type of overlap
#' @export
#' @importFrom GenomicFeatures genes
#' @importFrom GenomicRanges findOverlaps
#' @importFrom AnnotationDbi select
#' @importFrom dplyr left_join tibble
annotateByOverlap <- function(bsids = NULL, db = ahdb, t = "within"){
  db.type <- class(db)
  if(!db.type %in% c("TxDb", "EnsDb")){stop("Only TxDb or AhDb databases are supported.")}
  colN <- switch(
    db.type,
    "EnsDb" = c(
      gene_id = "GENEID", genename = "GENENAME", "GENEBIOTYPE",
      seqname = "SEQNAME", start = "GENESEQSTART", end = "GENESEQEND", strand = "SEQSTRAND"
    ),
    "TxDb" = c(
      #seqname = "TXCHROM", strand = "TXSTRAND",
      gene_id = "GENEID"
    )
  )
  getGenes <- switch( ### It seems not to matter, both pkgs handles the Genomic Ranges object the same way AFAIK!
    db.type,
    "EnsDb" = ensembldb::genes,
    "TxDb" = GenomicFeatures::genes
  )

  output <- tibble(
    bsID = bsids
  )

  message("Converting bsID input to GRanges ... ", appendLF = F)
  bsid.gr <- bsid2gr(bsids)
  bsid.gr <- bsid.gr - 1 # reduce interval to match with exon seq
  message("Done")

  message("Finding overlap to TXDB - ", appendLF = F)
  g.gr <- getGenes(db)
  ol <- GenomicRanges::findOverlaps(bsid.gr, g.gr, type = t, select = "all")
  if(length(ol) == 0){stop("No overlap found. If t = 'within', you could try to rerun using t = 'any'.")}
  hits <- tibble(
    bsID = bsids[queryHits(ol)],
    GENEID = g.gr$gene_id[subjectHits(ol)]
  )
  message("Done")

  message("Getting annotation on genes ... ", appendLF = F)
  # 1st level Annot
  # Gene id and gene symbol
  annot.1st <- AnnotationDbi::select(db, keys = g.gr$gene_id[subjectHits(ol)], columns = colN, keytype = "GENEID") ### AnnotationDBI::select or dplyr::select?
  hits <- dplyr::left_join(hits, annot.1st, by = "GENEID")
  message("Done")

  output <- dplyr::left_join(output, hits, by = "bsID")

  # message("Constructing output")
  # output <- data.frame(
  #   bsID = bsid.gr$bsID[queryHits(ol)],
  #   ensembl_gene_id = g.gr$gene_id[subjectHits(ol)],
  #   stringsAsFactors = F
  # )
  # output <- merge(output, annot)
  # output <- output %>% group_by(bsID) %>% summarise(ensembl_gene_id = paste(ensembl_gene_id, collapse = "|"), external_gene_name = paste(external_gene_name, collapse = "|"))
  # unknown <- bsids[!bsids %in% output$bsID]
  # if(length(unknown) != 0){
  #   output <- rbind(
  #     output,
  #     data.frame(bsID = unknown, ensembl_gene_id = NA, external_gene_name = NA)
  #   )
  # }
  # message("Done")

  return(output)
}


