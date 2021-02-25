#' Check if both mate pairs are consistent with backsplice junction.
#'
#' Deprecated function; not particularly intuitive
#'
#' @param df Path to a SJ.out.tab file or a data.frame.
#' @param endTol How many nucleotides on the "wrong" side of the junction can be tolerated.
#' @return data.frame with added column of boolean values indicating whether read pair is OK.
#' @examples
#' bsj.data <- checkPEreads(df = chi.data)
checkPEReads_old <- function(df = NULL, endTol = 5){
  # Check that input looks like candidate backsplice sites filtered from Chimeric.out.junction
  if(ncol(df) < 14 | !any(class(df) %in% c("data.frame","tbl")) | is.null(df)){stop("Too few columns, not consistent with Chimeric.out.junction data.")}
  if(any(df[, "X1"] != df[, "X4"]) ){stop("Acceptor and donor found on different chromosomes. Have the data been filtered?")}
  if(any(df[, "X3"] != df[, "X6"]) ){stop("Acceptor and donor found on different strands. Have the data been filtered?")}

  output <- ifelse(
    df[, "X3"] == "+",
    df[, "X11"] + endTol > df[, "X5"] & df[, "X13"] + sapply(df[, "X14"], parseCIGAR, returnLength = T) - endTol <= df[, "X2"], # YES
    df[, "X13"] + endTol > df[, "X2"] & df[, "X11"] + sapply(df[, "X12"], parseCIGAR, returnLength = T) - endTol <= df[, "X5"] # NO
  )

  return(output[,1])
}


#' Function to categorize backsplice sites
#'
#' Read are classified into various classes based on distance/overlap with known junctions. Need to make
#' function strand insensitive.
#'
#' @param cbs Table (tibble) of candidate baclsplice sites
#' @export
# categorizeBSJ <- function(cbs, stranded = T){
#   totalBSJ <- cbs$X10
#
#   # Reads that need no further processing
#   both.OK <- subset(
#     cbs,
#     shiftDonorToNearestJ == 0 & grepl("donor", donor.closest.type) &
#       shiftAcceptorToNearestJ == 0 & grepl("acceptor", acceptor.closest.type)
#   )$X10 # Both overlap. Exon containing circRNAs
#
#   one.OK <- subset(cbs, (shiftDonorToNearestJ == 0 & shiftAcceptorToNearestJ != 0 & grepl("donor", donor.closest.type)) | (shiftDonorToNearestJ != 0 & shiftAcceptorToNearestJ == 0 & grepl("acceptor", acceptor.closest.type)))$X10
#
#   intronic.circles <- subset(cbs, shiftDonorToNearestJ == -1 & shiftAcceptorToNearestJ == 1)$X10
#
#   # intron lariats. "The RT often mutates the branch-point A to any other nucleotide"
#   intron.lariats <- c(
#     subset(cbs, shiftAcceptorToNearestJ == -1 & grepl("donor", acceptor.closest.type))$X10,
#     subset(cbs, shiftDonorToNearestJ == 1 & grepl("acceptor", donor.closest.type))$X10
#   )
#
#   # Reads that can be shifted
#   shiftable.reads <- c(
#     ## DonorShift
#     subset(cbs, X3 == "-" & sign(shiftDonorToNearestJ) == -1 & abs(shiftDonorToNearestJ) <= X9)$X10, # sign(0) = 0, hence this subset will not get any reads where donor overlaps with known junc!
#     subset(cbs, X3 == "-" & sign(shiftDonorToNearestJ) == 1 & abs(shiftDonorToNearestJ) <= X8)$X10,
#     subset(cbs, X3 == "+" & sign(shiftDonorToNearestJ) == -1 & abs(shiftDonorToNearestJ) <= X8)$X10,
#     subset(cbs, X3 == "+" & sign(shiftDonorToNearestJ) == 1 & abs(shiftDonorToNearestJ) <= X9)$X10,
#
#     ## Acceptor shift
#     subset(cbs, X3 == "-" & sign(shiftAcceptorToNearestJ) == -1 & abs(shiftAcceptorToNearestJ) <= X9)$X10,
#     subset(cbs, X3 == "-" & sign(shiftAcceptorToNearestJ) == 1 & abs(shiftAcceptorToNearestJ) <= X8)$X10,
#     subset(cbs, X3 == "+" & sign(shiftAcceptorToNearestJ) == -1 & abs(shiftAcceptorToNearestJ) <= X8)$X10,
#     subset(cbs, X3 == "+" & sign(shiftAcceptorToNearestJ) == 1 & abs(shiftAcceptorToNearestJ) <= X9)$X10
#   )
#   shiftable.reads <- dplyr::setdiff(shiftable.reads, c(both.OK, one.OK, intronic.circles, intron.lariats)) # Remove reads that have already been accounted for.
#
#   category <- rep(NA, length(totalBSJ))
#   category[totalBSJ %in% both.OK] <- "Both OK"
#   category[totalBSJ %in% one.OK] <-  "One OK"
#   category[totalBSJ %in% intronic.circles] <- "Intronic circle"
#   category[totalBSJ %in% intron.lariats] <- "Intron lariat"
#   category[totalBSJ %in% shiftable.reads] <- "shiftable"
#   category[is.na(category)] <- "non-shiftable"
#   return(category)
# }
categorizeBSJ <- function(cbs, stranded = T){
  totalBSJ <- cbs$X10

  # Reads that need no further processing
  index.donorOK <- cbs$shiftDonorToNearestJ == 0
  index.acceptorOK <- cbs$shiftAcceptorToNearestJ == 0
  if(stranded){
    index.donorOK <- index.donorOK & grepl("donor", cbs$donor.closest.type)
    index.acceptorOK <- index.acceptorOK & grepl("acceptor", cbs$acceptor.closest.type)
  }

  index.bothOK <- index.donorOK & index.acceptorOK
  index.oneOK <- index.donorOK & !index.acceptorOK | !index.donorOK & index.acceptorOK

  both.OK <- totalBSJ[index.bothOK] # Both overlap. Exon containing circRNAs
  one.OK <- totalBSJ[index.oneOK]

  if(stranded){
    index.intronicCirc <- cbs$shiftDonorToNearestJ == -1 & cbs$shiftAcceptorToNearestJ == 1
  } else {
    index.intronicCirc <- NULL
  }
  intronic.circles <- totalBSJ[index.intronicCirc]
  if(stranded){
    index.lariats <- cbs$shiftAcceptorToNearestJ == -1 & grepl("donor", cbs$acceptor.closest.type) |
      cbs$shiftDonorToNearestJ == 1 & grepl("acceptor", cbs$donor.closest.type)
  } else {
    index.lariats <- NULL
  }
  intron.lariats <- totalBSJ[index.lariats]

  # if(stranded){
  index.shiftableDonor <- cbs$X3 == "-" & sign(cbs$shiftDonorToNearestJ) == -1 & abs(cbs$shiftDonorToNearestJ) <= cbs$X9 |
    cbs$X3 == "-" & sign(cbs$shiftDonorToNearestJ) == 1 & abs(cbs$shiftDonorToNearestJ) <= cbs$X8 |
    cbs$X3 == "+" & sign(cbs$shiftDonorToNearestJ) == -1 & abs(cbs$shiftDonorToNearestJ) <= cbs$X8 |
    cbs$X3 == "+" & sign(cbs$shiftDonorToNearestJ) == 1 & abs(cbs$shiftDonorToNearestJ) <= cbs$X9
  index.shiftableAcceptor <- cbs$X3 == "-" & sign(cbs$shiftAcceptorToNearestJ) == -1 & abs(cbs$shiftAcceptorToNearestJ) <= cbs$X9 |
    cbs$X3 == "-" & sign(cbs$shiftAcceptorToNearestJ) == 1 & abs(cbs$shiftAcceptorToNearestJ) <= cbs$X8 |
    cbs$X3 == "+" & sign(cbs$shiftAcceptorToNearestJ) == -1 & abs(cbs$shiftAcceptorToNearestJ) <= cbs$X8 |
    cbs$X3 == "+" & sign(cbs$shiftAcceptorToNearestJ) == 1 & abs(cbs$shiftAcceptorToNearestJ) <= cbs$X9
  index.shiftable <- index.shiftableDonor | index.shiftableAcceptor
  # } else {
  #   index.shiftable <- NA
  # }
  shiftable.reads <- totalBSJ[index.shiftable]

  category <- rep(NA, length(totalBSJ))
  category[totalBSJ %in% both.OK] <- "Both OK"
  category[totalBSJ %in% one.OK] <-  "One OK"
  category[totalBSJ %in% intronic.circles] <- "Intronic circle"
  category[totalBSJ %in% intron.lariats] <- "Intron lariat"
  category[totalBSJ %in% shiftable.reads] <- "shiftable"
  category[is.na(category)] <- "non-shiftable"
  return(category)
}

#' Function to shift alignment to overlap with known junctions
#'
#' @param cbs Table (tibble) of candidate backspliced reads
#' @export
shiftAlignment <- function(cbs){
  cbs$SHIFTED <- F

  ## Identify the different kinds shifts ##
  dont.touch.these <- cbs$shiftAcceptorToNearestJ == 0 | cbs$shiftDonorToNearestJ == 0

  # Group 1: Shifting that will bring both donor and acceptor in alignment
  rescue.both <- knownJunctionInRepeatRegion(cbs, junc = "both") &
    cbs$shiftAcceptorToNearestJ == cbs$shiftDonorToNearestJ &
    cbs$donor.closest.type == "donor" &
    cbs$acceptor.closest.type == "acceptor" &
    !dont.touch.these
  # Group 2: Shifting that will bring donor into alignment
  rescue.donor.only <- knownJunctionInRepeatRegion(cbs, junc = "donor") &
    !knownJunctionInRepeatRegion(cbs, junc = "acceptor") &
    cbs$donor.closest.type == "donor" &
    !dont.touch.these
  # group 3: Shifting that will bring acceptor into alignment
  rescue.acceptor.only <- knownJunctionInRepeatRegion(cbs, junc = "acceptor") &
    !knownJunctionInRepeatRegion(cbs, junc = "donor") &
    cbs$donor.closest.type == "acceptor" &
    !dont.touch.these
  message("Reads that can be shifted: ")
  message("- Both OK: ", length(which(rescue.both)))
  message("- Donor OK: ", length(which(rescue.donor.only)))
  message("- Acceptor OK: ", length(which(rescue.acceptor.only)))
  message("Only shifting 'Both OK'.")

  rescue.index <- rescue.both #| rescue.donor.only | rescue.acceptor.only
  shiftSize <- cbs$shiftDonorToNearestJ
  #shiftSize[rescue.acceptor.only] <- cbs$shiftAcceptorToNearestJ[rescue.acceptor.only]

  ## All three groups in one go ##
  cbs$X2[rescue.index] <- cbs$X2[rescue.index] + shiftSize[rescue.index]
  cbs$X5[rescue.index] <- cbs$X5[rescue.index] + shiftSize[rescue.index]
  cbs$X8[rescue.index] <- cbs$X8[rescue.index] + ifelse(cbs$X3[rescue.index] == "+", 1, -1) * shiftSize[rescue.index]
  cbs$X9[rescue.index] <- cbs$X9[rescue.index] + ifelse(cbs$X3[rescue.index] == "+", -1, 1) * shiftSize[rescue.index]

  cbs$shiftAcceptorToNearestJ[rescue.index] <- cbs$shiftAcceptorToNearestJ[rescue.index] - shiftSize[rescue.index]
  cbs$shiftDonorToNearestJ[rescue.index] <- cbs$shiftDonorToNearestJ[rescue.index] - shiftSize[rescue.index]

  # Shift/modify the CIGARs
  # It is always the numerically minimum value of X11/X13 that needs to be shifted (by shiftSize)!
  x11.min <- cbs$X11 < cbs$X13
  cbs$X11[rescue.index & x11.min] <- cbs$X11[rescue.index & x11.min] + shiftSize[rescue.index & x11.min]
  cbs$X13[rescue.index & !x11.min] <- cbs$X13[rescue.index & !x11.min] + shiftSize[rescue.index & !x11.min]

  # Modify the CIGARs
  ## plus
  index <- rescue.index & cbs$X3 == "+"
  fb <- lapply(cbs$X12[index], parseCIGAR, returnLength = F)
  sb <- lapply(cbs$X14[index], parseCIGAR, returnLength = F)
  ss <- shiftSize[index]
  fb.shifted <- lapply(seq_along(fb), function(i){
    output <- fb[[i]]
    slots <- 1:2
    if(!all(names(output)[slots] == c("S", "M"))){stop("Something is wrong: ", i)}
    output[slots] <- output[slots] + c(-1,1)*ss[i]
    return(output)
  }) %>% sapply(function(x){paste0(x, names(x), collapse = "")})
  sb.shifted <- lapply(seq_along(sb), function(i){
    output <- sb[[i]]
    slots <- (length(output)-1):length(output)
    if(!all(names(output)[slots] == c("M", "S"))){stop("Something is wrong: ", i)}
    output[slots] <- output[slots] + c(1,-1)*ss[i]
    return(output)
  }) %>% sapply(function(x){paste0(x, names(x), collapse = "")})

  cbs$X12[index] <- fb.shifted
  cbs$X14[index] <- sb.shifted
  cbs$SHIFTED[index] <- T

  ## minus
  index <- rescue.index & cbs$X3 == "-"
  fb <- lapply(cbs$X12[index], parseCIGAR, returnLength = F)
  sb <- lapply(cbs$X14[index], parseCIGAR, returnLength = F)
  ss <- shiftSize[index]

  fb.shifted <- lapply(seq_along(fb), function(i){
    output <- fb[[i]]
    slots <- (length(output)-1):length(output)
    if(!all(names(output)[slots] == c("M", "S"))){stop("Something is wrong: ", output)}
    output[slots] <- output[slots] + c(1,-1)*ss[i]
    return(output)
  }) %>% sapply(function(x){paste0(x, names(x), collapse = "")})
  sb.shifted <- lapply(seq_along(sb), function(i){
    output <- sb[[i]]
    slots <- 1:2
    if(!all(names(output)[slots] == c("S", "M"))){stop("Something is wrong: ", output)}
    output[slots] <- output[slots] + c(-1,1)*ss[i]
    return(output)
  }) %>% sapply(function(x){paste0(x, names(x), collapse = "")})

  cbs$X12[index] <- fb.shifted
  cbs$X14[index] <- sb.shifted
  cbs$SHIFTED[index] <- T

  ## Correct the bsIDs
  # cbs$bsID <- constructBsId(object)
  cbs$bsID <- paste("GRCh38", cbs$X1, cbs$X2, cbs$X5, cbs$X3, sep = ":")

  ## Correct the junction types
  jt <- junctionMotif(bsids = cbs$bsID, g = Hsapiens)
  cbs$X7 <- 0
  cbs$X7[jt == "GT/AG"] <- 1
  cbs$X7[jt == "CT/AC"] <- 2

  return(cbs)
}

#' Process Chimeric.out.junction file(s) to identify and quantify backsplice sites
#'
#' Process one Chimeric.out.junction file to produce files with information on all reads covering candidate backsplice sites and quantitative summary.
#'
#' @param path A character vector (length 1 or more) with paths to Chimeric.out.junction files.
#' @importFrom dplyr filter
#' @export
wrapperFunction <- function(path, stranded = T, ...){ ### log = TRUE
  # Check that files exist and appear to be valid Chimeric.out.junction files
  if(!all(file.exists(path))){stop("One or more files do not exist.")}
  if(!all(sapply(path, checkChimFile))){stop("One or more files do appear to be valid Chimeric.out.junction files.")} ### invalid instead of valid??

  for(file in path){
    # Get number of lines in chim.file
    zz <- file(file, open="rb")
    chimTotal <- 0L
    while (length(tmp <- readBin(zz, "raw", n = 1e4)) > 0) {
      chimTotal <- chimTotal + sum(tmp == as.raw(10L))
    }
    close(zz)

    candidate.bsj <- getCandidateBackspliceSites(file, backSpliceDist = 1e5, seqAcrossSj = T)
    bsjTotal <- nrow(candidate.bsj)

    candidate.bsj$PEreadsOK <- checkPEReads(df = candidate.bsj)
    candidate.bsj <- dplyr::filter(candidate.bsj, PEreadsOK)
    bsjGood <- nrow(candidate.bsj)

    if(stranded){
      candidate.bsj <- addKnownJunctions(cbs = candidate.bsj, kj = knownJunctions, ...)
    }
    if(!stranded){
      candidate.bsj <- addKnownJunctions_unstranded(cbs = candidate.bsj, kj = knownJunctions, ...)
    }
    rm.index <- is.na(candidate.bsj$shiftAcceptorToNearestJ) | is.na(candidate.bsj$shiftDonorToNearestJ)
    nearestJunctionAmbiguity <- length(which(rm.index))
    candidate.bsj <- candidate.bsj[!rm.index,]

    candidate.bsj$category <- categorizeBSJ(cbs = candidate.bsj)
    candidate.bsj <- shiftAlignment(cbs = candidate.bsj)

    cats <- rep(0, 6)
    names(cats) <- c("Both OK", "One OK", "Intronic circle", "Intron lariat", "shiftable", "non-shiftable")
    cats[names(table(candidate.bsj$category))] <- table(candidate.bsj$category)

    # Summarize counts
    bsj.tab <- countSumChimFile(df = candidate.bsj)
    # construct logfile and write to disk
    logfile <- data.frame(
      part = c("Total chimeric", "Total candidate bsj", "Good candidate bsj", "Uncertain nearest junction", names(cats)),
      values = c(chimTotal, bsjTotal, bsjGood, nearestJunctionAmbiguity, cats),
      type = rep(c("total", "cbs"), c(4,6))
    )

    fn1 <- sub("Chimeric.out.junction", "Backsplice.out.junction", file)
    fn2 <- sub("Chimeric.out.junction", "BSJ.out.tab", file)
    fn3 <- sub("Chimeric.out.junction", "bsj.log", file)

    write.table(candidate.bsj, file = fn1, quote = F, sep = "\t", row.names = F)
    write.table(bsj.tab, file = fn2, quote = F, sep = "\t", row.names = F)
    write.table(logfile, file = fn3, quote = F, sep = "\t", row.names = F) ### Made optional??
  }
}

#' Reverse strand information and donor/acceptor order, in the chimeric data
#'
#' @return The input data where the strand information, and the donor/acceptor order has been reversed.
#' @export
swapStrands <- function(data){
  # First swap strands
  data$X3 <- ifelse(data$X3 == "+", "-", "+")
  data$X6 <- ifelse(data$X6 == "+", "-", "+")

  # Swap donor and acceptor positions
  donors <- data[,c("X1", "X2", "X3")]
  data[,c("X1", "X2", "X3")] <- data[,c("X4", "X5", "X6")]
  data[,c("X4", "X5", "X6")] <- donors

  # Change the junction type
  data$X7 <- data$X7 %>% ifelse(. == 1, 9, .) %>% ifelse(. == 2, 1, .) %>% ifelse(. == 9,2,.)

  # Swap repeats to left / right
  tmp <- data$X8
  data$X8 <- data$X9
  data$X9 <- tmp

  return(data)
}
