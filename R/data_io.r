#' Read candidate backsplice sites from Chimeric.out.junction file
#'
#' \code{getCandidateBackspliceSites} reads entire Chimeric.out.junction file and returns reads consistent with backsplicing.
#' Candidate backsplice junctions are identified by applying the following filters
#' \enumerate{
#'   \item Donor and acceptor are on the same chromosome
#'   \item Donor and acceptor are on the same strand
#'   \item The acceptor is upstream of donor (On plus strand: acceptor position < donor position ; On minus strand: donor position < acceptor position)
#'   \item Distance between donor and acceptor is no more than 100.000 bases (by default, can be changed)
#' }
#'
#' @param fn Path to a Chimeric.out.junction file.
#' @param backSpliceDist The maximum allowed distance between donor and acceptor sites.
#' @param seqAcrossSj Boolean indicating whether to return all backspliced reads or just reads where one of the mates span the candidate backsplice site.
#' @param libType Character string indicating which stranded library protocol was used to generate the data, must be 'fr-firststrand' for e.g. dUTP method or 'fr-secondstrand'.
#' @param verbose Inlcude basic statstics on the chimeric readsm ust be logical (default FALSE).
#' @return tibble of reads covering backsplice sites
#' @export
getCandidateBackspliceSites <- function(fn, backSpliceDist = 1e5, seqAcrossSj = TRUE, libType = "fr-firststrand", verbose = FALSE){ ### name "seqAcrossSj" is not describtive, should it be something like "all_candidates"??
  message("File: ", fn)
  # Ensure that the function arguments are correct
  if (!is.null(err <- checkVariables(obj = backSpliceDist, expect_class = c("numeric", "integer"), vari = "backSpliceDist"))) stop(err)
  if (!is.null(err <- checkVariables(obj = seqAcrossSj, expect_class = "logical", vari = "seqAcrossSj"))) stop(err)
  if (!libType %in% c("fr-firststrand", "fr-secondstrand")) stop("Currently, only support for Illumina stranded seq data.") # To be develloped futher in the near future!
  if (!is.null(err <- checkVariables(obj = verbose, expect_class = "logical", vari = "verbose"))) stop(err)

  # Ensure that files are correct format
  if (!is.null(err <- checkInputFiles(files = fn, single_file = TRUE))) stop(err)
  if (!is.null(err <- checkChimFile(file = fn))) stop(err)

  # Import chimeric reads
  message(" - loading file... ", appendLF = F)
  chim <- readChimFile(fn)
  message("Done")

  # Create filters
  message(" - filtering data ... ", appendLF = F)
  filters <- generateCandidateFilters(data = chim, backSpliceDist = backSpliceDist, seqAcrossSj)
  message("Done")

  candidates <- chim[filters, ]

  # For some library preparation types we need to move the identifed junction to the opposite strand
  if(libType == "fr-firststrand"){
    message(" - swapping strands ... ", appendLF = F)
    candidates <- swapStrands(data = candidates)
    message("Done")
  }


  ### Do we still want to print messages on stats? My suggestion is to have a Verbose mode!
  if (verbose){
    tot <- nrow(chim)
    bsj <- nrow(candidates)
    msg <- paste0(" - basic statistics:\n",
                  "   - total chimeric reads: ", tot,
                  "\n   - chimeric reads consistent with backsplicing: ", bsj, " (", round(bsj*100/tot, digits = 2), "%)\n"#,
                  #paste0("- candidate backsplice sites on chromosomes: ", paste(BiocGenerics::sort(unique(candidates$X1)), collapse = "; "))
                  )
    message(msg)
  }
  return(candidates)
}

# getCandidateBackspliceSites <- function(file, backSpliceDist = 1e5, seqAcrossSj = T, libType = "fr-firststrand"){ ### name "seqAcrossSj" is not describtive, should it be something like "all_candidates"??
#   message("Processing ", basename(file))
#   # Check if files exists
#   if(!file.exists(file)){stop("File doest not exist: ", file)}
#   if(length(file) != 1){stop("Can only process one file. Use lapply to process multiple files (lapply(fns, getCandidateBackspliceSites).")}
#   if(!grepl("Chimeric.out.junction", basename(file))){warning("Non-standard filename!")}
#
#   if(!checkChimFile(file)){stop("The input file appears to have the wrong format.")} ### Same error msg as in wrapper function? -> maybe change check function to be interruptive?
#   if(libType != "fr-firststrand"){stop("Currently, only support for Illumina stranded seq data.")}
#
#   # Read the data
#   message("Loading ... ", appendLF = F) ### Add filename to msg??
#   db <- readChimFile(file)
#   message("Done")
#
#   # Create filters
#   message("Filtering data ... ", appendLF = F)
#   same.chr <- db$X1 == db$X4
#   same.str <- db$X3 == db$X6
#   acceptorUpstreamOfDonor <- (db$X3 == "-" & (db$X5 > db$X2)) | (db$X3 == "+" & (db$X2 > db$X5))
#   dist.below.cutoff <- abs(db$X2-db$X5) <= backSpliceDist ### Cutoff = "at most" (instead of "less than")
#   return.these <- same.chr & same.str & acceptorUpstreamOfDonor & dist.below.cutoff ### Curious: How does this work?
#
#   if(seqAcrossSj){
#     return.these <- return.these & db$X7 >= 0 ### What are type(0) junction types??
#   }
#   message("Done")
#
#   output <- db[return.these, ]
#
#   if(libType == "fr-firststrand"){ # Reverse strand + ac/do
#     message("Swapping strands ... ", appendLF = F)
#     output$X3 <- output$X6 <- ifelse(output$X3 == "+", "-", "+") # Swap strands
#     tmp <- output$X5
#     output$X5 <- output$X2
#     output$X2 <- tmp
#     message("Done")
#   }
#
#   #tot <- nrow(db)
#   #bsj <- nrow(output)
#   #message("Statistics")
#   #message("- total chimeric reads: ", tot)
#   #message("- chimeric reads consistent with backsplicing: ", bsj, " (", round(bsj*100/tot, digits = 2), "%)")
#   #message("- candidate backsplice sites on chromosomes: ", paste(BiocGenerics::sort(unique(output$X1)), collapse = "; "))
#   return(output)
# }
#


# generateCandidateFilters <- function(data, backSpliceDist = 1e5, seqAcrossSj = TRUE){
#   same.chr <- data$X1 == data$X4
#   same.str <- data$X3 == data$X6
#   acceptorUpstreamOfDonor <- (data$X3 == "-" & (data$X5 > data$X2)) | (data$X3 == "+" & (data$X2 > data$X5))
#   dist.below.cutoff <- abs(data$X2-data$X5) <= backSpliceDist
#
#   filters <- same.chr & same.str & acceptorUpstreamOfDonor & dist.below.cutoff
#   if(seqAcrossSj){
#     filters <- filters & data$X7 > -1
#   }
#
#   return(filters)
# }
