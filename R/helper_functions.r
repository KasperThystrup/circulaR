#' Extract junction motif(s) from genome
#'
#' This function takes a vector of bsIDs and a BSgenome-object and constructs junction-motifs of the backsplice junctions.
#'
#' @param bsids An character vector with bsIDs.
#' @param g A BS genome object. Default Hsapiens.
#' @return Character vector with the junction motifs.
#' @importFrom BSgenome providerVersion Views
#' @importFrom tibble add_column
#' @export
junctionMotif <- function(bsids, g = Hsapiens){
  # Check input
  if(!class(g) == "BSgenome"){stop("Please provide valid BS genome object.")}
  if(class(bsids) != "character" & !all(regexpr(":", bsids, fixed = T) == 5)){stop("Does not look like valid bsIDs.")}
  org <- bsids %>% strsplit(":") %>% sapply(function(x)x[1]) %>% unique
  if(length(org) > 1){stop("Backsplice IDs from multiple organisms detected. Currently, only a single organism is supported.")}
  if(BSgenome::providerVersion(g) != org){stop("Make sure that genome build identifier in bsID and the supplied BS genome object are identical.")}

  d.gr <- bsid2junc(bsids, junc.type = "donor") %>% GenomicRanges::resize(., 2, fix = "start")
  a.gr <- bsid2junc(bsids, junc.type = "acceptor") %>% GenomicRanges::resize(., 2, fix = "end")

  output <- paste(
    as.character(BSgenome::Views(g, d.gr)),
    as.character(BSgenome::Views(g, a.gr)),
    sep = "/"
  )

  return(output)
}


#' Wrapper for easy reading of Chimeric.out.junction file
#'
#' This function uses read_tsv from the tidyverse to read the 14-coulmn file.
#'
#' @param file Path to file.
#' @inheritDotParams readr::read_tsv n_max
#' @examples
#' tmp <- readChimFile(fn = "Chimeric.out.junction")
#' @return Tibble data frame with data on chimeric reads
#' @importFrom readr read_tsv cols col_character col_integer
#' @export
readChimFile <- function(file, n_max = Inf, ...){

  # Ensure that the file is readable
  data <- readr::read_tsv(
    file, col_names = F, n_max = n_max,
    col_types = readr::cols(
      X1 = readr::col_character(),
      X2 = readr::col_integer(),
      X3 = readr::col_character(),
      X4 = readr::col_character(),
      X5 = readr::col_integer(),
      X6 = readr::col_character(),
      X7 = readr::col_integer(),
      X8 = readr::col_integer(),
      X9 = readr::col_integer(),
      X10 = readr::col_character(),
      X11 = readr::col_integer(),
      X12 = readr::col_character(),
      X13 = readr::col_integer(),
      X14 = readr::col_character()
    ),  progress = F,
    ...
  )
  return(data)
}

#' Summarize the read count for each backsplice junction
#'
#' @param df Table (tibble) of chimeric reads
#' @param asGRanges Boolean indication if results should be returned as a GRanges
#' @export
countSumChimFile <- function(df, asGRanges = F){
  if (!is.null(err <- checkVariables(obj = df, expect_class = c("tbl_df", "tbl", "data.frame"), vari = "df"))) stop(err)
  if(base::ncol(df) == 14 & all(colnames(df) %in% paste0("X", 1:14))){warning("Looks like the data have not been through the circulaR-pipeline.")}

  if(!"bsID" %in% colnames(df)){
    message("bsID column not found, adding it now.")
    df$bsID <- constructBsId(df) ###! Previosuly paste("hg38", df$X1, df$X2, df$X5, df$X3, sep = ":")
  }

  # Count amount of reads that covers unique backsplice sites
  df <- df %>% dplyr::group_by(bsID) %>% dplyr::summarise(count = dplyr::n())

  # Convert to GRAnges
  if(asGRanges){
    df$chromosome_name <- strsplit(df$bsID, ":") %>% sapply(function(x)x[2])
    df$donor <- strsplit(df$bsID, ":") %>% sapply(function(x)x[3]) %>% as.numeric
    df$acceptor <- strsplit(df$bsID, ":") %>% sapply(function(x)x[4]) %>% as.numeric
    df$strand <- strsplit(df$bsID, ":") %>% sapply(function(x)x[5])

    df$start <- apply(df[,c("donor", "acceptor")],1,min)
    df$end <- apply(df[,c("donor", "acceptor")],1,max)
    df <- GenomicRanges::GRanges(df)
  }
  return(df)
}

#' Convert bsID to BED like format.
#'
#' @export
bsid2bed <- function(id, UCSCstyle = F){
  id.split <- strsplit(id, ":")
  sapply(id.split, function(x){
    paste0(
      ifelse(UCSCstyle, paste0("chr", x[2]), x[2]), "|",
      min(as.integer(x[3:4])), "-", max(as.integer(x[3:4]))-1, "|",
      x[5]
    )
  })
}

#' Convert bsID to GRanges object
#'
#' @export
bsid2gr <- function(x){
  tmp <- strsplit(x, ":")

  gr <- GenomicRanges::GRanges(
    seqnames = sapply(tmp, function(x)x[2]),
    ranges = IRanges(
      start = sapply(tmp, function(x)min(as.numeric(x[3:4]))),
      end = sapply(tmp, function(x)max(as.numeric(x[3:4])))
    ),
    strand = sapply(tmp, function(x)x[5]),
    bsID = x
  )
  return(gr)
}

#' Convert bsID to GRanges object of acceptor or donor positions
#'
#' @export
bsid2junc <- function(bsids=NULL, junc.type = "acceptor"){
  if(!junc.type %in% c("acceptor", "donor")){stop("junc.type must be either acceptor or donor.")}
  if(junc.type == "acceptor") i<-4
  if(junc.type == "donor") i<-3
  tmp <- strsplit(bsids, ":")

  gr <- GenomicRanges::GRanges(
    seqnames = sapply(tmp, function(x)x[2]),
    ranges = IRanges(
      start = as.numeric(sapply(tmp, function(x)x[i])),
      width = 1
    ),
    strand = sapply(tmp, function(x)x[5]),
    bsID = bsids
  )
  return(gr)
}

#' #' Construct bsID from genomic coordinates
#' #'
#' #' @param gb Genome build
#' #' @export
#' constructBsId <- function(df, gb){
#'   bsID <- paste(gb, df$X1, df$X2, df$X5, df$X3, sep = ":")
#'   return(bsID)
#' }




#' Calculate distance covered by the alignment from a CIGAR string
#'
#' ### Do we wish to export?? (only if used outside of "wrapper functions")
#'
#' @param cig The CIGAR string.
#' @return value indicating number of genomic nucleotides covered by the alignment.
#' @examples
#' bla bla
#' parseCIGAR(cig = "77M107N24M-135p4M107N77M20S")
#' @export
parseCIGAR <- function(cig, returnLength = T) {
  if(length(cig) > 1){warning("Supplied multiple CIGAR strings, only using first!")}

  L <- as.numeric(strsplit(cig, "[A-Zp]", perl=T)[[1]])
  C <- strsplit(cig, "-?\\d+", perl=T)[[1]][-1] # Remove initial empty slot
  if(length(L) != length(C)){stop("Parsing of CIGAR string ", cig, " went wrong.")}
  if(!returnLength){
    names(L) <- C
    return(L)
  } else {
    align.length <- 0
    for(i in 1:length(L)){
      if(C[i] %in% c("S", "I")) next() # Soft mapped or insertions should be disregarded
      align.length <- align.length + L[i]
    }
    return(align.length)
  }
}


#' Function to return covered chromosomal positions from a start position and CIGAR string
#'
#' @param first.pos Genomic position where the CIGAR string starts (from 11th or 13th column in the Chimeric.out.junction file)
#' @param cigar CIGAR string (from 12th or 14th column in the Chimeric.out.junction file)
#' @export
cigar2covr <- function(first.pos = NULL, cigar = NULL){
  cigar.groups <- parseCIGAR(cigar, returnLength = F)
  if(names(cigar.groups)[1] == "S") cigar.groups <- cigar.groups[-1]
  if(names(cigar.groups)[length(cigar.groups)] == "S") cigar.groups <- cigar.groups[-length(cigar.groups)]

  index <- lapply(seq_along(cigar.groups), function(i){
    if(names(cigar.groups)[i] == "M" & i == 1){
      return(first.pos:(first.pos+cigar.groups[i]-1))
    }
    if(names(cigar.groups)[i] == "M" & i > 1) {
      st <- first.pos + sum(cigar.groups[1:(i-1)])
      ed <- st + cigar.groups[i] - 1
      return(st:ed)
    }
  }) %>% do.call(c, .)

  return(index)
}

#' Function to identify chimeric reads eligible for shifting
#'
#' @param cbs Tibble of candidate backspliced reads
#' @return Boolean vector
#' @export
knownJunctionInRepeatRegion <- function(cbs, junc = "donor"){
  if(! junc %in% c("donor", "acceptor", "both")){stop("junc must be either 'donor', 'acceptor' or 'both'")}

  donor.index <- cbs$X3 == "+" & sign(cbs$shiftDonorToNearestJ) == -1 & abs(cbs$shiftDonorToNearestJ) <= cbs$X8 | # Donor, on plus strand, kj and repeat upstream
    cbs$X3 == "+" & sign(cbs$shiftDonorToNearestJ) == 1 & abs(cbs$shiftDonorToNearestJ) <= cbs$X9 | # Donor, on plus strand, kj and repeat downstream
    cbs$X3 == "-" & sign(cbs$shiftDonorToNearestJ) == 1 & abs(cbs$shiftDonorToNearestJ) <= cbs$X8 | # Donor, on minus strand, kj and repeat upstream
    cbs$X3 == "-" & sign(cbs$shiftDonorToNearestJ) == -1 & abs(cbs$shiftDonorToNearestJ) <= cbs$X9 # Donor, on minus strand, kj and repeat downstream

  acceptor.index <- cbs$X3 == "+" & sign(cbs$shiftAcceptorToNearestJ) == -1 & abs(cbs$shiftAcceptorToNearestJ) <= cbs$X8 | # Acceptor, on plus strand, kj and repeat upstream
    cbs$X3 == "+" & sign(cbs$shiftAcceptorToNearestJ) == 1 & abs(cbs$shiftAcceptorToNearestJ) <= cbs$X9 | # Acceptor, on plus strand, kj and repeat downstream
    cbs$X3 == "-" & sign(cbs$shiftAcceptorToNearestJ) == 1 & abs(cbs$shiftAcceptorToNearestJ) <= cbs$X8 | # Acceptor, on minus strand, kj and repeat upstream
    cbs$X3 == "-" & sign(cbs$shiftAcceptorToNearestJ) == -1 & abs(cbs$shiftAcceptorToNearestJ) <= cbs$X9 # Acceptor, on minus strand, kj and repeat downstream

  if(junc == "donor"){
    return(donor.index)
  }
  if(junc == "acceptor"){
    return(acceptor.index)
  }
  if(junc == "both"){
    return(donor.index & acceptor.index)
  }
}


#' Function to calculate coverage of backspliced reads
#'
#' After subsetting table of candidate backsplice sites to a single backsplice sites
#' we can calculate the coverage.
#' @param df Tibble of candidate backspliced reads
#' @param asGranges Boolean indicating whether result should be returned as GRanges obejct, otherwise a table is returned
#' @export
calcCoverage <- function(df, asGRanges = T){
  if (nrow(df) != 0) {
    cv <- lapply(1:nrow(df), function(i){
      seg1 <- cigar2covr(first.pos = df$X11[i], cigar = df$X12[i])
      seg2 <- cigar2covr(first.pos = df$X13[i], cigar = df$X14[i])
      return(c(seg1,seg2))
    }) %>% unlist %>% table

    gr <- GenomicRanges::GRanges(
      seqnames = unique(df$X1),
      range = IRanges(
        start = as.numeric(names(cv)),
        width = 1
      ),
      strand = unique(df$X3)
    )
    gr$cov <- as.numeric(cv)

    if (asGRanges) {
      return(gr)
    } else {
      return(cv)
    }
  } else {
    return(NULL)
  }
}

#' Prune chromosomes
#'
#' Define which chromosomes to use from genomic annotation object, default `NULL` does not exclude any chromosomes
#'
#' @param db An object containing a Seqinfo class, for chromosomal values.
#' @param chromosomes A vector of chromosomes to keep, must be `integer`, `character` or `NULL` (default). Alterantive provide the value "standard" No chromosomes are removed if set to `NULL`.
#'
#' @return Same as input, but with certin chromosomes removed.
#' @examples
#' @export
pruneChromosomes <- function(db, chromosomes = NULL){
  if (!is.null(err <- checkVariables(obj = chromosomes, expect_class = c("NULL", "character", "integer"), vari = "chromosomes"))) stop(err)
  message("Removing chromosomes.")

  # Check if chromosomes are provided
  if(is.empty(chromosomes)){
    message(" - 'chromosomes' has not been defined, no chromosomes will be removed!")
  } else {
    dbType <- class(db)[1]

    message("- Subsetting chromosomes ...", appendLF = F)

    if (dbType == "tbl_df"){
      db <- subset(db, X1 %in% chromosomes)


    } else if (dbType == "GenomicFeatures"){
      if (!all(chromosomes %in% seqlevels(db))) stop("The following provided chromosomes, was not found in the database and will be removed:\n",
                                                        paste0(chromosomes[!chromosomes %in% seqlevels(db)]), collapse = "\n")
      db <- restoreSeqlevels(db)
      seqlevels(db) <- BiocGenerics::sort(unique(chromosomes))
    }
  }
  message(" Done")
  return(db)
}


### Latest and Undocumented functions

#' Simple wrapper for generating a GRanges object, it supports the same input as 'GenomicRanges::Ranges', i.e. addition of metadata columns.
#' As a minimum a chromosome name and aligned strand must be provided, in addition two of either must be provided of following; a 'start', 'end', or 'width'
#'
junctionsToGRanges <- function(seqnames, start = NULL, end = NULL, width = NULL, strand, ...){
  # Checks
  if (is.null(start) & is.null(end) | is.null(start) & is.null(width) | is.null(end) & is.null(width)) stop("Too few variables are provided for generating a GRanges object.\n",
                                                                                                            "Define only two of these variables: start, end & width")

  if (is.null(start)) {
    range <- IRanges(end = end, width = width)
  } else if (is.null(end)) {
    range <- IRanges(start = start, width = width)
  } else if (is.null(width)) {
    range <- IRanges(start = start, end = end)
  } else if (end - start != width + 1 | start + width != end - 1 | end - width != start + 1) {
    stop("Something went wrong with defining IRanges. Please ensure that 'start', 'end', or 'width' is propperly defined!")
  }


  grObject <- GenomicRanges::GRanges(
    seqnames = seqnames,
    ranges = range,
    strand = strand,
    ... = ...

  )

  return(grObject)
}


#' Candidate backsplice data are split into donors and acceptors as an GRanges object,
#' the annotation data are screened, to determine which individual donor and acceptor site maps to a nearby known junction.
#' The candidate backsplice data that overlaps or are near a known junction site are annotated and returned.
#' Following annotation are used for donors and acceptors:
#' * queryHit: The dataset postioin of a given backsplice candidate read.
#' * jID: The closest known junction unqiue ID (provided during 'constructSJDB')
#' * type: The closest known junction type (acceptor or donor)
#' * strand: The closest known junction strand position
#' * ToNearestJ: Minimum distance of the read donor or acceptor position to the nearest known donor or acceptor site.
#' @importFrom BiocGenerics start
#' @importFrom GenomicRanges strand
findNearest <- function(juncCand, kj, type, cbs){
  # Determine which known junctions that are closest to junction candidates
  nearestHits <- GenomicRanges::nearest(juncCand, kj, select = "all", ignore.strand = FALSE)

  # Esnure that all junction candidates are included
  if (length(unique(queryHits(nearestHits))) != nrow(cbs)) stop(paste0("One or more ", type, " reads are missing a nearest junction!"))

  # Generate junction annotation data
  near <- tibble::as_data_frame(as.matrix(nearestHits))
  closest.jID <- kj[near$subjectHits]$jID
  closest.type <- kj[near$subjectHits]$type
  closest.strand <- as.character(GenomicRanges::strand(kj[near$subjectHits]))
  ToNearestJ <- BiocGenerics::start(kj)[near$subjectHits] - BiocGenerics::start(juncCand)[near$queryHits]

  if (type == "donor") {
    near$donor.closest.jID <- closest.jID
    near$donor.closest.type <- closest.type
    #near$donor.closest.strand <- closest.strand
    near$shiftDonorToNearestJ <- ToNearestJ
  } else if (type == "acceptor") {
    near$acceptor.closest.jID <- closest.jID
    near$acceptor.closest.type <- closest.type
    #near$acceptor.closest.strand <- closest.strand
    near$shiftAcceptorToNearestJ <- ToNearestJ
  } else {
    stop("Wrong format used on type. Must be 'donor' or 'acceptor' (lowercase sensitive)!")
  }

  return(near)
}

#' The function used to identify nearest known junction ('IRange::nearest') reports the gap between two genomic positions and not the distance.
#' This means that a position (x) can be reported as nearest to three positions (x-1, x, x+1). If such ambiguities arise by using the nearest function,
#' we need to determine if the results can be disambiguated. This function will return reads where:
#' closest.donor.type == donor & closest.acceptor.type == acceptor
#' AND
#' shiftAcceptorToNearestJ == 0 | shiftDonorToNearestJ == 0
#' AND
#' donor/acceptor pairs that are closest to known junctions
rescueAmbiguous <- function(acceptors, donors, ambiguousHits){
  # Ensure that there are ambiguous reads in the dataset
  if(length(ambiguousHits) != 0){
    ambiguousJuncs <- dplyr::full_join(
      subset(acceptors, queryHits %in% ambiguousHits),
      subset(donors, queryHits %in% ambiguousHits),
      by = "queryHits"
    )
    rescuedJuncs <- ambiguousJuncs %>% dplyr::group_by(queryHits) %>%
      dplyr::filter(
        donor.closest.type == "donor" & acceptor.closest.type == "acceptor" &
          abs(shiftAcceptorToNearestJ) == min(abs(shiftAcceptorToNearestJ)) &
          abs(shiftDonorToNearestJ) == min(abs(shiftDonorToNearestJ)) &
          (shiftAcceptorToNearestJ == 0 | shiftDonorToNearestJ == 0)
      )

  return(rescuedJuncs)
  }
}


#' Analyses the provided donors and acceptors to determine if there are any junctions that cannot be unequivocally determined by the nearest function
determineAmbiguousHits <- function(acceptors, donors){
  message("Screening for ambiguities in closest junction.")
  duplicatedDonors <- unique(donors$queryHits[BiocGenerics::duplicated(donors$queryHits)])  ###!CHECK: BiocGenerics:: added, if not working, remove dependency!
  duplicatedAcceptors <- unique(acceptors$queryHits[BiocGenerics::duplicated(acceptors$queryHits)])  ###!CHECK: BiocGenerics:: added, if not working, remove dependency!
  ambiguousHits <- BiocGenerics::sort(unique(c(duplicatedDonors, duplicatedAcceptors)))

  #message("- ambiguities for closest junction found in ", length(ambiguousHits), " reads (", signif(length(ambiguousHits)*100/nrow(cbs), digits = 2), "% of input).")

  return(ambiguousHits)
}


#' FURTHER FUNCTIONALITY should be implemented if we wish to determine lariats and intronic reads -> further functionality are currently commented out!!
#' Determines which known backsplice donor and acceptor junction candidates can be rescued and indexes these.
indexJunctions <- function(donorSubset, acceptorSubset, type){

  # Generate indexes for Donor ambiguities
  if (tolower(type) == "donor"){
    index.exonic <- donorSubset$donor.closest.type == "donor" &
      donorSubset$shiftDonorToNearestJ == 0

    # index.intronic <- donorSubset$donor.closest.type == "acceptor" &
    #   donorSubset$shiftDonorToNearestJ == ifelse(unique(donorSubset$donor.closest.strand) == "-", -1, 1) &
    #   acceptorSubset$acceptor.closest.type == "donor" &
    #   acceptorSubset$shiftAcceptorToNearestJ == ifelse(unique(acceptorSubset$acceptor.closest.strand) == "-", 1, -1)
    #
    # index.lariat <- donorSubset$donor.closest.type == "acceptor" &
    #   donorSubset$shiftDonorToNearestJ == ifelse(unique(donorSubset$donor.closest.strand) == "-", 1, -1) &
    #   acceptorSubset$acceptor.closest.type == "donor" &
    #   abs(acceptorSubset$shiftAcceptorToNearestJ) > 1

    # Generate indexes for Acceptor ambiguities
  } else if (tolower(type) == "acceptor"){
    index.exonic <- acceptorSubset$acceptor.closest.type == "acceptor" &
      acceptorSubset$shiftAcceptorToNearestJ == 0

    # index.intronic <- donorSubset$donor.closest.type == "acceptor" &
    #   donorSubset$shiftDonorToNearestJ == ifelse(unique(donorSubset$donor.closest.strand) == "-", -1, 1) &
    #   acceptorSubset$acceptor.closest.type == "donor" &
    #   acceptorSubset$shiftAcceptorToNearestJ == ifelse(unique(acceptorSubset$acceptor.closest.strand) == "-", 1, -1)
    #
    # index.lariat <- donorSubset$donor.closest.type == "acceptor" &
    #   abs(donorSubset$shiftDonorToNearestJ) > 1 &
    #   acceptorSubset$acceptor.closest.type == "donor" &
    #   abs(acceptorSubset$shiftAcceptorToNearestJ) == ifelse(unique(donorSubset$donor.closest.strand) == "-", 1, -1)

    # Generate indexes for both Donor and Acceptor ambiguities
  } else if (tolower(type) == "both"){
    indexDonor <- donorSubset$donor.closest.type == "donor" & donorSubset$shiftDonorToNearestJ == 0
    indexAcceptor <- acceptorSubset$acceptor.closest.type == "acceptor" & acceptorSubset$shiftAcceptorToNearestJ == 0

    # Yield results early
    index <- list(indexDonor, indexAcceptor)
    return(index)

  } else {
    stop("The variable 'type' in the function 'indexJuncions' had a wrong value. \n",
         "'type' must be set to either 'donor', 'acceptor' or 'both'!")
  }

  # Ensure that ambiguities doens't have multiple categories
  #if(any(index.exonic) & any(index.intronic) & any(index.lariat)) stop(paste0("Read appears to be both intronic, exonic, and/or lariat!")) # More specificity!!
  index <- index.exonic#|index.intronic|index.lariat

  return(index)
}


#' Compares the full backsplice candidate set with reads that are ambiguously mapped to the known junctions,
#' and deselcts all ambiguous reads. The ambiguous reads that can be rescued are then defined and included.
#' Finally the selected backsplice sites are fully annotated and returned.
resolveReads <- function(cbs, ambiguousHits, toRescue, donors, acceptors){
  resolved <- cbs
  resolved$queryHits <- 1:nrow(cbs)
  resolved <- subset(resolved, !queryHits %in% ambiguousHits)

  if(!is.null(toRescue)){
    toRescue <- unique(c(toRescue, resolved$queryHits))
    acceptors <- subset(acceptors, queryHits %in% toRescue)
    donors <- subset(donors, queryHits %in% toRescue)
  }

  resolved <- dplyr::left_join(resolved, acceptors[,c("queryHits", "acceptor.closest.jID", "acceptor.closest.type", "shiftAcceptorToNearestJ")], by = "queryHits")
  resolved <- dplyr::left_join(resolved, donors[,c("queryHits", "donor.closest.jID", "donor.closest.type", "shiftDonorToNearestJ")], by = "queryHits")

  return(resolved)
}


#' Add information on usage of donor and acceptor sites for linear splicing
#'
#' @importFrom GenomicRanges findOverlaps GRanges
#' @export
addLinCounts <- function(cbj = cbj.sum, sj = sj.data){
  # Split into plus and minus strand
  cbj.plus <- cbj[grepl("\\+$", cbj$bsID),]
  cbj.minus <- cbj[grepl("\\-$", cbj$bsID),]


  ## On plus
  # acceptor
  ac.plus <- GenomicRanges::findOverlaps(
    bsid2junc(bsids = cbj.plus$bsID, junc.type = "acceptor"),
    sj,
    type = "end"
  ) %>% data.frame
  # donor
  do.plus <- GenomicRanges::findOverlaps(
    bsid2junc(bsids = cbj.plus$bsID, junc.type = "donor"),
    sj,
    type = "start"
  ) %>% data.frame

  ac.plus$bsID <- cbj.plus$bsID[ac.plus$queryHits]
  ac.plus <- ac.plus %>% dplyr::group_by(bsID) %>% dplyr::summarise(lin.acceptor.count = sum(sj$read.count.unique[subjectHits]))

  do.plus$bsID <- cbj.plus$bsID[do.plus$queryHits]
  do.plus <- do.plus %>% dplyr::group_by(bsID) %>% dplyr::summarise(lin.donor.count = sum(sj$read.count.unique[subjectHits]))

  ## On plus
  # acceptor
  ac.minus <- GenomicRanges::findOverlaps(
    bsid2junc(bsids = cbj.minus$bsID, junc.type = "acceptor"),
    sj,
    type = "start"
  ) %>% data.frame
  # donor
  do.minus <- GenomicRanges::findOverlaps(
    bsid2junc(bsids = cbj.minus$bsID, junc.type = "donor"),
    sj,
    type = "end"
  ) %>% data.frame

  ac.minus$bsID <- cbj.minus$bsID[ac.minus$queryHits]
  ac.minus <- ac.minus %>% dplyr::group_by(bsID) %>% dplyr::summarise(lin.acceptor.count = sum(sj$read.count.unique[subjectHits]))

  do.minus$bsID <- cbj.minus$bsID[do.minus$queryHits]
  do.minus <- do.minus %>% dplyr::group_by(bsID) %>% dplyr::summarise(lin.donor.count = sum(sj$read.count.unique[subjectHits]))

  ac <- dplyr::bind_rows(ac.plus, ac.minus)
  do <- dplyr::bind_rows(do.plus, do.minus)

  output <- dplyr::left_join(data.frame(cbj), ac, by = "bsID") %>% dplyr::left_join(., do, by = "bsID")
  return(GenomicRanges::GRanges(output))
}

#' Infer the exons structure of a circular RNA
#'
#' The genomic region covered by a backsplice junction is used to identify overlapping transcript models
#' and subsequently, whether the nucleotides adjacent to the donor and acceptor positions overlap with known exon boundaries.
#' If both positions overlap with known exon termini, the resulting exon structure is returned. Otherwise, NULL is returned.
#' @param bsid A valid backsplice ID
#' @param eBt A GRangesList of exons-by-transcript, generated using the exonsBy-function.
#' @return GRanges object or list of GRanges obejcts giving the inferred exon structure of a circRNA identified by the supplied backsplice ID.
#' @importFrom GenomicRanges strand
#' @importFrom dplyr inner_join
#' @importFrom magrittr %>%
#' @importFrom IRanges subsetByOverlaps
#' @importFrom S4Vectors mcols
#' @export
inferExons <- function(bsid, eBt){
  if(!is.null(err <- checkVariables(eBt, "CompressedGRangesList", "eBt"))){stop(err)}
  if(!class(bsid) == "character" & length(strsplit(bsid, ":")[[1]]) != 5){stop("bsID appears to have wrong format.")}
  if(length(bsid) != 1){stop("Please input 1 bsID.")}

  gr <- bsid2gr(bsid)-1 # Need to convert donor/acceptor positions to first and last exon positions
  if(as.logical(GenomicRanges::strand(gr) == "*")){
    stop("Unstranded data not supported in this version.")
  }

  txs <- IRanges::subsetByOverlaps(eBt, gr)
  # Identify the transcripts that have exon structure that overlaps with the backsplice junction
  valid.models <- lapply(txs, function(x){
    first <- IRanges::subsetByOverlaps(
      GenomicRanges::resize(x, fix = "start", width = 1),
      GenomicRanges::resize(gr, fix = "start", width = 1)
    )
    last <- IRanges::subsetByOverlaps(
      GenomicRanges::resize(x, fix = "end", width = 1),
      GenomicRanges::resize(gr, fix = "end", width = 1)
    )
    if(length(first) == 1 & length(last) == 1){
      output <- x[first$exon_rank:last$exon_rank]
      output$exon_rank <- output$exon_rank - first$exon_rank + 1
      return(output)
    }
  }) %>% .[!sapply(., is.null)]

  # Identify unique models
  output <- lapply(valid.models, function(x){
    S4Vectors::mcols(x) <- NULL
    return(x)
  }) %>% unique

  if(length(output) == 1){
    names(output) <- paste0(bsid, paste0("_variant",1:length(output)))
  } else if(length(output) > 1){
    names(output) <- paste0(bsid, paste0("_variant",1:length(output)))
  } else {
    output <- NULL
  }


  return(output)
}

#' Read and parse Log.final.out file
#'
#' @param file Path to file.
#' @inheritDotParams readr::read_tsv n_max
#' @examples
#' tmp <- readChimFile(fn = "Chimeric.out.junction")
#' @return Tibble data frame with data on chimeric reads
#' @importFrom magrittr %>%
#' @importFrom readr read_lines
#' @export
STARlogParser <- function(file, return = "All"){
  # Read the log file
  output <- readr::read_lines(file) %>%
    gsub("^\\s+","", .) %>%
    strsplit("\\s\\|\t") %>%
    sapply(function(x){
      if(length(x) == 2) {
        return(x)
      } else {
        return(NULL)
      }
    }) %>%
    do.call(rbind, .)
  index <- ifelse(
    output[,1] %in% c("Started job on", "Started mapping on", "Finished on"),
    "time",
    "stats"
  )
  output <- split(data.frame(output, stringsAsFactors = F), index)
  output$time$X2 <- strptime(as.character(output$time$X2), format = "%b %d %H:%M:%S")

  output$stats$X2 <- sub("%", "", output$stats$X2) %>% as.numeric()
  return(output$stats)
}


####### Experimental code to automatically detect a lot of required information
#' Guess whether data are paired ended or not by parsing the BAM file header
#'
#' @param bf path to bam file
#'
#' @return Logical vector of length 1, indicating if sample is paired ended (TRUE), single ended (FALSE), or not able to determine (NULL)
#' @export
#'
#' @examples
#' guessPE(path.to.bam.file)
guessPE <- function(bf){
  hd <- Rsamtools::scanBamHeader(bf)

  # Command line
  index <- grepl("@CO", names(hd[[1]]$text))
  lookForFastq <- gregexpr("\\.[fastq|fq](.gz)?", hd[[1]]$text[index])[[1]]

  if(length(lookForFastq) == 1){
    #message("Looks like single ended.")
    return(FALSE)}
  if(length(lookForFastq) == 2){
    #message("Looks like paired ended.")
    return(TRUE)
  }
  if(!length(lookForFastq) %in% c(1,2)){message("Could not determine if sample is SE or PE."); return(NULL)}
}


#' Parse a directory to identify samples for processing in the circulaR pipeline
#'
#' @param dir The path to directory where STAR output files are stored
#' @examples
#' ...
#' @return Unknown
#' @export
#' @importFrom tibble tibble
parseDataDir <- function(dir, PE = NULL, frfs = NULL, org = NULL, geno = NULL){
  # Checks
  if(!dir.exists(dir)){stop("Directory does not exist.")}

  all.files <- list.files(dir, recursive = T, full.names = T)

  # Determine number of samples based on how many Chimeric.out.junction files are present
  chim.index <- grepl("Chimeric.out.junction$", all.files)
  num.samples <- sum(chim.index)
  sample.prefix <- all.files[chim.index] %>% sub("Chimeric.out.junction", "", ., fixed = T)

  output <- tibble::tibble(
    id = basename(sample.prefix),
    chim.file = all.files[grepl("Chimeric.out.junction$", all.files)],
    bam.file = all.files[grepl("Aligned.*\\.out\\.bam$", all.files)],
    log.file = all.files[grepl("Log.final.out$", all.files)],
    sj.file = all.files[grepl("SJ.out.tab$", all.files)],
    PE = NA,
    frfs = NA,
    org = NA,
    geno = NA
  )

  # Add info on whether data are PE or SE
  if(is.null(PE)){
    message("Determining if samples are SE or PE.")
    output$PE <- sapply(output$bam.file, guessPE)
  }
  if(length(PE) == 1 & is.logical(PE)){
    message("Setting all samples to PE=", PE)
    output$PE <- PE
  }
  if(length(PE) == nrow(output) & is.logical(PE)){
    message("Setting PE as specified.")
    output$PE <- PE
  }
  if(!is.null(PE) & (!length(PE) %in% c(1, nrow(output)) | !is.logical(PE))){
    stop("Please provide correct PE input")
  }

  #Add info on whether samples are firstread.firststrand
  if(is.null(frfs)){
    stop("Please specify if the sequencing data are 'first read first strand'.")
  }
  if(length(frfs) == 1 & is.logical(frfs)){
    message("Setting all samples to frfs=", frfs)
    output$frfs <- frfs
  }
  if(length(frfs) == nrow(output) & is.logical(frfs)){
    message("Setting frfs as specified.")
    output$frfs <- frfs
  }
  if(!is.null(frfs) & (!length(frfs) %in% c(1, nrow(output)) | !is.logical(frfs))){
    stop("Please provide correct frfs")
  }

  #Add organism
  if(is.null(org)){
    stop("Please specify organism.")
    # Eventually, determine by parsing BAM header
  }
  if(length(org) == 1){
    message("Setting all samples to org=", org)
    output$org <- org
  }
  if(length(org) == nrow(output)){
    message("Setting org as specified.")
    output$org <- org
  }
  if(!is.null(org) & (!length(org) %in% c(1, nrow(output)))){
    stop("Please provide correct org")
  }

  #Add genome build
  if(is.null(geno)){
    stop("Please specify genome build")
    # Eventually, determine by parsing BAM header
  }
  if(length(geno) == 1){
    message("Setting all samples to geno=", geno)
    output$geno <- geno
  }
  if(length(geno) == nrow(output)){
    message("Setting geno as specified.")
    output$geno <- geno
  }
  if(!is.null(geno) & (!length(geno) %in% c(1, nrow(output)))){
    stop("Please provide correct geno")
  }

  # Check if files are present
  output$chim.file <- ifelse(file.exists(output$chim.file), output$chim.file, NA)
  output$bam.file <- ifelse(file.exists(output$chim.file), output$bam.file, NA)
  output$log.file <- ifelse(file.exists(output$chim.file), output$log.file, NA)
  output$sj.file <- ifelse(file.exists(output$chim.file), output$sj.file, NA)

  return(output)
}


#' Determine whether reads in a bsj.read-table uses known splice junctions
#'
#' @param df BSJ data table
#' @examples
#' ...
#' @return A vector that is as long as there are rows in the input.
#' @export
knownJunctionSites <- function(df){
  knownDonor <- df$shiftDonorToNearestJ == 0
  knownAcceptor <- df$shiftAcceptorToNearestJ == 0

  output <- rep(NA, nrow(df))
  output[knownDonor & knownAcceptor] <- "[D,A]" #"[known,known]" #"donor,acceptor"
  output[knownDonor & !knownAcceptor] <- "[D,_]" #"[known,?]" #"donor,NA"
  output[!knownDonor & knownAcceptor] <- "[_,A]" #"[?,known]" #"NA,acceptor"
  output[!knownDonor & !knownAcceptor] <- "[_,_]" #"[?,?]" #"NA,NA"
  return(output)
}
