#' A error message generator, for checking classes of variables used in functions.
#'
#' @param obj Variable object
#' @param expect_class A string value of the expected class type
#' @param vari A string value of the function variable name
#' @return Null if no error occured, otherwise returns a custom error message.
#' @examples
#' # Will produce an error
#' a <- function(x){
#'   if (!is.null(err <- checkVariables(obj = x, req_class = "numerical", vari = "x"))) stop(err)
#'   return(x + 10)
#' }
#' # Will produce a costum error
#' a("one")
#' # Will execute the entire function
#' a(1)
#' @export
checkVariables <- function(obj, expect_class, vari){
  if (all(class(obj) %in% expect_class)){
    err <- NULL
  } else {
    err <- paste0("The variable: '", vari, "' was of the class: '", class(obj), "'.\n",
                  "This variable must be of the following class(es): '", paste0(expect_class, collapse = " OR "), "'")
  }
  return(err)
}

#' Check if amount of input files are correct, and notify if any files are missing.
#'
#' @param files A character vector of one or more file names
#' @param single_file A logical variable, defining whether the function supports one or multiple files
#' @param dev Developper mode, not yet implemented, may be removed!
#' @export
checkInputFiles <- function(files, single_file){

  # Ensure correct format of function variables
  if (!is.null(err <- checkVariables(obj = files, expect_class = "character", vari = "files"))) stop(err)
  if (!is.null(err <- checkVariables(obj = single_file, expect_class = "logical", vari = "single_file"))) stop(err)
  #if (!is.null(err <- checkVariables(obj = dev, expect_class = "logical", vari = "dev"))) stop(err)

  err <- NULL
  # Ensure that the amount of files are corect
  if (single_file & (length(files) != 1)){
    err <- paste0("This function can only process one file at a time.\n",
         "To process multiple files you can use lapply():\n",
         "(lapply(files, getCandidateBackspliceSites)")

  } else {
    # Check if files exist
    missing_files <- FALSE
    for (fn in files){
      if (!file.exists(fn)){
        missing_files <- TRUE
        message("The file '", fn, "' does not exist!")
      }
    }

    # Provide error message
    if (missing_files){
      err <- paste0("One or more of the files does not exist.\n",
                    "Check file names and make sure to use the full file path!")
    }
  }

  # Ensure that no error has been issued
  return(err)
}


#' Function to test if Chimeric.out.junction file can be read without problems
#'
#' Attempts to read the first 1000 lines of the file. If no errors/warnings are
#' encountered then it is assumed that the entire file is OK.
#'
#' @param fn Path to file
#' @return Boolean value, any errors/warnings when reading fiurst 1000 lines.
#' @importFrom readr problems
#' @export
checkChimFile <- function(file, ...){

  # Check if naming convention are maintained for the chimeric junctions file
  if (!grepl("Chimeric.out.junction", basename(file))){
    warning("The file name '", file, "' derivates for the standard name convention 'Chimeric.out.junction'!")
  }

  # Determine if any problems occurs when readin the file
  err <- NULL
  test <- readr::problems(readChimFile(file = file, n_max = 1000, ...))
  if(nrow(test)>0){
    err <- paste0("Unable to read the chimeric junctions file: ", file, "\n",
                  "Please check if file is corrupted!")
  }

  return(err)
}


#' Check if both mates of a read pair are consistent with the identified backsplice junction.
#'
#' This function tests if the two fragments of a chimeric read fall within the
#' range defined defined by the backsplice junctions. Using the genome coordinates
#' of the acceptor and donor (columns 2 and 5 in Chimeric.out.junction), the
#' function checks if the region covered by the two fragments of the chimeric reads
#' (i.e. column 11 plus the width of the CIGAR string in column 12 and column 13
#' plus the width of the CIGAR string in column 14) are within the region spanned
#' by the acceptor and donor sites. By setting the end tollerance (endTol) you can allow for some
#' coverage outside the circularized region.
#'
#' @param df Chimeric.out.junction table returned by getCandidateBackspliceSites
#' @param endTol How many nucleotides on the "wrong" side of the junction can be tolerated.
#' @export
#' @importFrom IRanges IRanges %within%
#' @return Logical vector
checkPEReads <- function(df, endTol = 5){
  if(!is.null(err <- checkVariables(obj = endTol, expect_class = c("numeric", "integer"), vari = "endTol"))) stop(err)
  if(!all(df$X3 == df$X6)){stop("Some junctions have donor and acceptor on opposite strands.\n",
                                "Please read Chimeric.out.junction file using function: getCandidateBackspliceSites")}


  # # Construct GRanges of region between backsplice sites
  # circ.region <- GRanges(
  #   seqnames = df$X1,
  #   ranges = IRanges(
  #     start = ifelse(df$X2 < df$X5, df$X2, df$X5),
  #     end = ifelse(df$X2 < df$X5, df$X5, df$X2)
  #   ),
  #   strand = df$X3
  # )
  # circ.region <- circ.region + endTol # Expand the region according to the end tolerance
  #
  # # Construct GRanges objects of the regions covered by the two fragments of the chimeric read
  # cigarRegion1 <- GRanges(
  #   seqnames = df$X1,
  #   ranges = IRanges(
  #     start = df$X11,
  #     end = df$X11 + sapply(df$X12, parseCIGAR, returnLength = T) - 1
  #   ),
  #   strand = df$X3
  # )
  #
  # cigarRegion2 <- GRanges(
  #   seqnames = df$X1,
  #   ranges = IRanges(
  #     start = df$X13,
  #     end = df$X13 + sapply(df$X14, parseCIGAR, returnLength = T) - 1
  #   ),
  #   strand = df$X3
  # )
  #
  # # Test which chimeric reads have both fragments within the circularized region
  # goodIndex <- cigarRegion1 %within% circ.region & cigarRegion2 %within% circ.region

  goodIndex <- ifelse(
    df$X3 == "-",
    ((df$X11 + endTol) >= df$X2) & (df$X13 + sapply(df$X14, parseCIGAR) - endTol <= df$X5), # YES
    ((df$X13 + endTol) >= df$X5) & (df$X11 + sapply(df$X12, parseCIGAR) - endTol <= df$X2) # NO
  )
  return(goodIndex)
}
