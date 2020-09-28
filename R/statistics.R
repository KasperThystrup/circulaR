# bsjExonStats ####

#' @title Calculate exon statistics
#' @name bsjExonStats
#'
#' @description \code{bsjExonStats} determines which circRNAs covers exons of the transcript models provided through the \code{annot} argument,
#'  and calculates overall statistics for the circRNAs in the input \code{circSample} or \code{circExperiment} object.
#'  The circRNA size are calculated from splicing exons, and covered exons are counted. overlapping exons in a given transcript model is merged,
#'  and in case of a circRNA covering multiple transcritp model, the largest size is calculated and selected.
#'  Be patient: It currently take minutes to complete when ncores = 5.
#'
#' @param object A circSample or circExperiment object
#' @param annot The annotation object, either as an \code{EnsDb} object or a pre defined exons by transcript \code{CompressedGRangesList} object, generated from running \code{exonsBy(x = annot, by = "tx")}
#' @param ncores The number of threads (must be a whole number)
#'
#' @author Kasper Thystrup Karstensen
#'
#' @return
#' @export
#'
#' @examples
#' # Execute with circSample object using default threading (4 cores)
#' ahdb <- AnnotationHub()[["AH64923"]]  #  Ensembl Homo sapiens release-94
#' bsjExonStats(object = circ_smpl, annot = ahdb)
#'
#' # Exons by transcripts can be used as annotation object
#' ahdb <- AnnotationHub()[["AH64923"]]
#' ex_tx <- exonsBy(x = ahdb, by = "tx")
#' bsjExonStats(object = circ_smpl, annot = ex_tx)
#'
#' # Execute with circExperiment object using 10 cores
#' bsjExonStats(object = circ_exp, annot = ahdb, ncores = 10L)
#'
#' @import IRanges
#' @import GenomicRanges
#' @importFrom ensembldb exonsBy
#' @importFrom pbmcapply pbmclapply
#' @importFrom tidyr tibble
#'
setGeneric(name = "bsjExonStats",
           def = function(object, ...)
             standardGeneric("bsjExonStats"))

#' @rdname bsjExonStats

setMethod(
  f = "bsjExonStats",
  signature = "circSample",
  definition = function(object,
                        annot,
                        ncores = 4,
                        ...) {
    # Check arguments
    if (!is.null(err <- checkVariables(obj = ncores,
                                       expect_class = c("integer",
                                                        "numeric"),
                                       vari = "ncores"))) stop(err)

    if (is.numeric(ncores))
      ncores <- as.integer(ncores)

    # Check annotation object
    if (class(annot) == "EnsDb"){
      annot <- ensembldb::exonsBy(x = annot, by = "tx")
    } else if (class(annot) != "CompressedGRangesList") {
      stop("The class of the annotation object: ", class(annot), " is not supported with this function!")
    }

    # Fetch all circRNA candidates
    smpl <- bsj.counts(object, returnAs = "gr")

    # Determine which circRNA that covers whole exons
    circs.idx <- IRanges::findOverlaps(query = annot,
                                             subject = smpl,
                                             type = "within") %>%
      subjectHits

    exon_stats <- pbmcapply::pbmclapply(X = seq_along(circs.idx), FUN = function(j) {
      idx <- circs.idx[j]
      circ <- smpl[idx]

      covered.transcripts <- IRanges::subsetByOverlaps(x = annot,
                                                       ranges = circ,
                                                       type = "within")  %>%
        IRanges::reduce()

      # Handle circRNAs covering multiple transcript models
      if (length(covered.transcripts) > 1) {
        largest.txs <- list()

        # Determine the total covered exonic sizes and index by transcript models
        for (k in seq_along(covered.transcripts)) {
          covered.exons.tmp <- covered.transcripts[[k]]
          exon.size.tmp <- IRanges::subsetByOverlaps(x = covered.exons.tmp,
                                                     ranges = circ,
                                                     type = "within")

          # Record total exonic size and transcript model index
          largest.txs[k] <- IRanges::width(exon.size.tmp) %>% sum
        }

        # Select the transcript model which haves the larges total exonic size
        covered.exons <- covered.transcripts[[which.max(largest.txs)]]

      } else if (length(covered.transcripts) == 1) {
        covered.exons <- unlist(covered.transcripts)
      } else {
        stop(paste("A selected circRNA did not overlap the given exons, despite being selected for overlapping exons previously.",
                   "Perhabs the annotation object or a circSample / circExperiment object is corrupt:",
                   annot, object, sep = "\n"))
      }

      exon.count <- length(covered.exons)
      exon.size <- IRanges::width(covered.exons) %>%
        sum

      output <- tidyr::tibble(bsID = circ$bsID,
                              exonic_size = exon.size,
                              exon_count = exon.count)

      return(output)
    }, mc.cores = ncores) %>% do.call(rbind, .)

    return(exon_stats)
  })

#' @rdname bsjExonStats

setMethod(
  f = "bsjExonStats",
  signature = "circExperiment",
  definition = function(object,
                        annot,
                        ...) {

    if (class(annot) == "EnsDb") {
      # Determine exon by transcript models
      annot <- exonsBy(x = annot, by = "tx")
    }

    exon_stats <- lapply(X = samples(object),
                         FUN = function(sample.list) {
                           smpl <- unlist(sample.list)
                           exon_stats <- bsjExonStats(object = smpl,
                                                      annot = annot,
                                                      ...)
                           return(exon_stats)
                         }) %>% do.call(rbind, .)


    return(exon_stats)
  })



#' @description \code{plotExonStats} Plots crude histograms of the overall statistics calculated for exons from transcript models, covered by detected circRNAs.
#'
#' @param stats \code{plotExonStats}: Output from \code{bsjExonStats}
#' @param size_limit \code{plotExonStats}: Collapse data points where the exonic size is higher thatn the given limit. When NULL (default) data is not collapsed.
#' @param count_limit \code{plotExonStats}: Collapse data points where the exonic count is higher thatn the given limit. When NULL (default) data is not collapsed.
#' @param ...
#'
#' @return Generates histograms over the exon statistics
#' @export
#'
#' @examples
#' stats <- bsjExonStats(object = circ_exp, annot = ahdb, ncores = 10L)
#' # Collapse circRNAs with exon sizes > 2000 and exon counts > 10, into a single column.
#' plotExonStats(stats, size_limit = 2000, count_limit = 10)
#'
#' # Return a list of ggplot object
#' p <- plotExonStats(stats)
#'
#' @rdname bsjExonStats
setGeneric(
  name = "plotExonStats",
  def = function(stats, size_limit = NULL, count_limit = NULL, ...)
    standardGeneric("plotExonStats"))

setMethod(
  f = "plotExonStats",
  signature = "tbl_df",
  definition = function(stats, ...) {

    if (!is.null(size_limit)) {
      limit <- as.integer(size_limit)

      stats <- dplyr::mutate(stats,
                             exonic_size = ifelse(exonic_size > limit,
                                                  limit,
                                                  exonic_size))
    }

    if (!is.null(count_limit)) {
      limit <- as.integer(count_limit)

      stats <- dplyr::mutate(stats,
                             exon_count = ifelse(exon_count > limit,
                                                 limit,
                                                 exon_count))
    }


    exon_sizes.plot <- ggplot2::ggplot(data = stats, mapping = aes(x = exonic_size)) +
      ggplot2::stat_bin(binwidth = 100, boundary = 0.5, ...)

    exon_count.plot <- ggplot2::ggplot(data = stats, mapping = aes(x = exon_count)) +
      ggplot2::stat_bin(binwidth = 1, boundary = 0.5, ...)

    plots <- list("exon_sizes" = exon_sizes.plot, "exon_counts" = exon_count.plot)
    return(plots)
  })


# NEXT ####

