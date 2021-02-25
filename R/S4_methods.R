################################ Accessor methods #################################


setMethod(f = "sample.id",
          signature = "circSample",
          definition = function(object){
            return(object@sample.id)
          })


setMethod(f = "sample.id",
          signature = "circExperiment",
          definition = function(object){
            return(sapply(samples(object), sample.id))
          })



setMethod(f = "organism",
          signature = "circSample",
          definition = function(object){
            return(object@organism)
          })


setMethod(f = "organism",
          signature = "circExperiment",
          definition = function(object){
            return(sapply(samples(object), organism))
          })


setMethod(f = "gb",
          signature = "circSample",
          definition = function(object){
            return(object@gb)
          })


setMethod(f = "gb",
          signature = "circExperiment",
          definition = function(object){
            return(sapply(samples(object), gb))
          })


setMethod(f = "chim.file",
          signature = "circSample",
          definition = function(object){
            return(object@chim.file)
          })


setMethod(f = "chim.file",
          signature = "circExperiment",
          definition = function(object){
            return(sapply(samples(object), chim.file))
          })


setMethod(f = "bam.file",
          signature = "circSample",
          definition = function(object){
            return(object@bam.file)
          })


setMethod(f = "bam.file",
          signature = "circExperiment",
          definition = function(object){
            return(sapply(samples(object), bam.file))
          })



setMethod(f = "log.file",
          signature = "circSample",
          definition = function(object){
            return(object@log.file)
          })


setMethod(f = "log.file",
          signature = "circExperiment",
          definition = function(object){
            return(sapply(samples(object), log.file))
          })


setMethod(f = "sj.file",
          signature = "circSample",
          definition = function(object){
            return(object@sj.file)
          })


setMethod(f = "sj.file",
          signature = "circExperiment",
          definition = function(object){
            return(sapply(samples(object), sj.file))
          })


setMethod(f = "count.file",
          signature = "circSample",
          definition = function(object){
            return(object@count.file)
          })


setMethod(f = "count.file",
          signature = "circExperiment",
          definition = function(object){
            return(sapply(samples(object), count.file))
          })


setMethod(f = "firstread.firststrand",
          signature = "circSample",
          definition = function(object){
            return(object@firstread.firststrand)
          })


setMethod(f = "firstread.firststrand",
          signature = "circExperiment",
          definition = function(object){
            return(sapply(samples(object), firstread.firststrand))
          })



setMethod(f = "paired.end",
          signature = "circSample",
          definition = function(object){
            return(object@paired.end)
          })


setMethod(f = "paired.end",
          signature = "circExperiment",
          definition = function(object){
            return(sapply(samples(object), paired.end))
          })


setMethod(f = "bsj.reads",
          signature = "circSample",
          definition = function(object){
            return(object@bsj.reads)
          })


setMethod(f = "bsj.reads",
          signature = "circExperiment",
          definition = function(object, returnAs = "list"){
            if(! returnAs %in% c("list", "table")){stop("Can only return circSamples as 'list' or 'table'.")}
            if(returnAs == "table"){
              output <- lapply(
                seq_along(sample.id(object)),
                function(i){
                  o <- bsj.reads(samples(object)[[i]])
                  o$sample.id <- sample.id(object)[i]
                  return(o)
                }
              ) %>% bind_rows()
            } else {
              output <- lapply(samples(object), bsj.reads)
            }
            return(output)
          })


setMethod(f = "bsj.counts",
          signature = "circSample",
          definition = function(object, returnAs = "list") {
            if (returnAs == "list"){
              return(object@bsj.counts)
            } else if (returnAs == "gr") {
              return(GRanges(object@bsj.counts))
            }
          })


setMethod(f = "bsj.counts",
          signature = "circExperiment",
          definition = function(object, returnAs = "list", use.names = T){
            output <- lapply(samples(object), bsj.counts, returnAs)
            if (use.names){
              names(output) <- sample.id(object)
            }
            if (returnAs == "list") {
              return(output)
            } else if (returnAs == "gr") {
              return(GRangesList(output))
            }
          })


setMethod(f = "lsj.counts",
          signature = "circSample",
          definition = function(object){
            return(object@lsj.counts)
          })


setMethod(f = "lsj.counts",
          signature = "circExperiment",
          definition = function(object){
            return(lapply(samples(object), lsj.counts))
          })


setMethod(f = "label",
          signature = "circSample",
          definition = function(object){
            return(object@label)
          })


setMethod(f = "label",
          signature = "circExperiment",
          definition = function(object){
            return(sapply(samples(object), label))
          })

##NOTE Documentation must be placed here!##
#' Accessor for experiment samples
#'
#' circExperiment class accessors for showing the circSample objects
#'
#' @param object circExperiment object
#'
#' @return List of samples stored in the circExperiment object
#'
#' @examples
#' # The samples of an experiment
#' samples(experiment.object)
#'
#' @importFrom Biobase samples
#' @export
setMethod(f = "samples",
          signature = "circExperiment",
          definition = function(object){
            return(object@samples)
          })

setMethod(f = "exp.mat",
          signature = "circExperiment",
          definition = function(object, format, normFactors){
            if(!format %in% c("long", "wide")){stop("Format should be 'long' or 'wide'.")}
            if(!is.null(normFactors) & !length(normFactors) == length(samples(object))){stop("Number of samples and number of normalization factors do not match!")}
            if(!is.null(normFactors) & !all(names(normFactors) %in% sample.id(object))){stop("normFactor should be a named vector where names correspond to sample IDs in the circExperiment-object.")}

            output <- bsj.counts(object) %>%
              lapply(GenomicRanges::mcols) %>%
              lapply(as.data.frame) %>%
              dplyr::bind_rows(.id = "sample")

            if(!is.null(normFactors)){
              output <- output %>% group_by(sample) %>% mutate(count = count/normFactors[sample])
            }

            if(format == "wide"){
              output <- output %>% tidyr::spread(sample, count)
            }
            return(tibble::as_tibble(output))
          })



setMethod(f = "path",
          signature = "circExperiment",
          definition = function(object){
            return(object@path)
          })


setMethod(f = "name",
          signature = "circExperiment",
          definition = function(object){
            return(object@name)
          })

############################### Replace methods ###############################

setMethod(f = "sample.id<-",
          signature = "circSample",
          definition = function(object, value){
            object@sample.id <- value

            if (validObject(object)){
              return(object)
            }
          }
)


setMethod(f = "sample.id<-",
          signature = "circExperiment",
          definition = function(object, value){

            if (length(value) != length(samples(object))){
              stop("Number of samples do not match number of sample IDs!")
            } else if (any(duplicated(value))){
              stop("Sample IDs must be unique!")
            }

            s <- samples(object)
            samples(object) <- lapply(seq_along(s), function(i){
              .s <- s[[i]]
              sample.id(.s) <- value[[i]]
              return(.s)
            })

            return(object)
          }
)

#' @importFrom BiocGenerics organism
setMethod(f = "organism<-",
          signature = "circSample",
          definition = function(object, value){
            object@organism <- value

            if (validObject(object)){
              return(object)
            }
          }
)


setMethod(f = "organism<-",
          signature = "circExperiment",
          definition = function(object, value){

            if (length(value) == 1){
              value = rep(x = value, length(samples(object)))
              message("Recycling '",value,"' for all samples.")
            } else if (length(value) != length(samples(object))){
              stop("Number of samples do not match number of organism names.")
              #stop("There are ", length(value), "organisms out of ", length(samples(object)), " samples!")
            }

            s <- samples(object)
            samples(object) <- lapply(seq_along(s), function(i){
              .s <- s[[i]]
              organism(.s) <- value[[i]]
              return(.s)
            })

            return(object)
          }
)


setMethod(f = "gb<-",
          signature = "circSample",
          definition = function(object, value){
            object@gb <- value

            if (validObject(object)){
              return(object)
            }
          })


setMethod(f = "gb<-",
          signature = "circExperiment",
          definition = function(object, value){

            if (length(value) == 1){
              value = rep(x = value, length(samples(object)))
              message("Recycling '",value,"' for all samples.")
            } else if (length(value) != length(samples(object))){
              stop("Number of samples do not match number of genome build names.")
              #stop("There are ", length(value), "genome builds out of ", length(samples(object)), " samples!")
            }

            s <- samples(object)
            samples(object) <- lapply(seq_along(s), function(i){
              .s <- s[[i]]
              gb(.s) <- value[[i]]
              return(.s)
            })

            return(object)
          }
)


setMethod(f = "chim.file<-",
          signature = "circSample",
          definition = function(object, value){
            object@chim.file <- value

            if (validObject(object)){
              return(object)
            }
          })


setMethod(f = "bam.file<-",
          signature = "circSample",
          definition = function(object, value){
            object@bam.file <- value

            if (validObject(object)){
              return(object)
            }
          })


setMethod(f = "log.file<-",
          signature = "circSample",
          definition = function(object, value){
            object@log.file <- value

            if (validObject(object)){
              return(object)
            }
          })


setMethod(f = "sj.file<-",
          signature = "circSample",
          definition = function(object, value){
            object@sj.file <- value

            if (validObject(object)){
              return(object)
            }
          })


setMethod(f = "count.file<-",
          signature = "circSample",
          definition = function(object, value){
            object@count.file <- value

            if (validObject(object)){
              return(object)
            }
          })


setMethod(f = "firstread.firststrand<-",
          signature = "circSample",
          definition = function(object, value){
            object@firstread.firststrand <- value

            if (validObject(object)){
              return(object)
            }
          })


setMethod(f = "firstread.firststrand<-",
          signature = "circExperiment",
          definition = function(object, value){

            if (length(value) == 1){
              value = rep(x = value, length(samples(object)))
              message("Recycling '",value,"' for all samples.")
            } else if (length(value) != length(samples(object))){
              stop("Number of samples do not match number of firstread.firststrand values.")
              #stop("There are ", length(value), "firstread.firststrand values out of ", length(samples(object)), " samples!")
            }

            s <- samples(object)
            samples(object) <- lapply(seq_along(s), function(i){
              .s <- s[[i]]
              firstread.firststrand(.s) <- value[[i]]
              return(.s)
            })

            return(object)
          }
)


setMethod(f = "paired.end<-",
          signature = "circSample",
          definition = function(object, value){
            object@paired.end <- value

            if (validObject(object)){
              return(object)
            }
          })


setMethod(f = "paired.end<-",
          signature = "circExperiment",
          definition = function(object, value){

            if (length(value) == 1){
              value = rep(x = value, length(samples(object)))
              message("Recycling '",value,"' for all samples.")
            } else if (length(value) != length(samples(object))){
              stop("Number of samples do not match number of paired.end values.")
              #stop("There are ", length(value), "paired.end values out of ", length(samples(object)), " samples!")
            }

            s <- samples(object)
            samples(object) <- lapply(seq_along(s), function(i){
              .s <- s[[i]]
              paired.end(.s) <- value[[i]]
              return(.s)
            })

            return(object)
          }
)


setMethod(f = "bsj.reads<-",
          signature = "circSample",
          definition = function(object, value){

            # Validate BSJ.reads format
            if (!is.empty(value) & !all(paste0("X", 1:14) %in% colnames(value))){
              stop("Input BSJ read data seems to have the wrong format - check that input is a Chimeric.out.junction file and that format has not changed (works with STAR2.6.1d)!")
            }

            object@bsj.reads <- value

            if (validObject(object)){
              return(object)
            }
          })


setMethod(f = "bsj.reads<-",
          signature = "circExperiment",
          definition = function(object, value){

            if (length(value) != length(samples(object))){
              stop("Number of samples do not match number of BSJ read datatables.")
              #stop("There are ", length(value), "BSJ read datatables out of ", length(samples(object)), " samples!")
            }

            s <- samples(object)
            samples(object) <- lapply(seq_along(s), function(i){
              .s <- s[[i]]
              bsj.reads(.s) <- value[[i]]
              return(.s)
            })

            return(object)
          }
)

setMethod(f = "bsj.counts<-",
          signature = "circSample",
          definition = function(object, value){
            object@bsj.counts <- value

            if (validObject(object)){
              return(object)
            }
          })


setMethod(f = "bsj.counts<-",
          signature = "circExperiment",
          definition = function(object, value){

            if (length(value) != length(samples(object))){
              stop("Number of samples do not match number of BSJ count datatables.")
              #stop("There are ", length(value), "BSJ count datatables out of ", length(samples(object)), " samples!")
            }

            s <- samples(object)
            samples(object) <- lapply(seq_along(s), function(i){
              .s <- s[[i]]
              bsj.counts(.s) <- value[[i]]
              return(.s)
            })

            return(object)
          }
)


setMethod(f = "lsj.counts<-",
          signature = "circSample",
          definition = function(object, value){
            object@lsj.counts <- value

            if (validObject(object)){
              return(object)
            }
          })


setMethod(f = "lsj.counts<-",
          signature = "circExperiment",
          definition = function(object, value){

            if (length(value) != length(samples(object))){
              stop("Number of samples do not match number of LSJ count datatables.")
              #stop("There are ", length(value), "LSJ count datatables out of ", length(samples(object)), " samples!")
            }

            s <- samples(object)
            samples(object) <- lapply(seq_along(s), function(i){
              .s <- s[[i]]
              lsj.counts(.s) <- value[[i]]
              return(.s)
            })

            return(object)
          }
)


setMethod(f = "label<-",
          signature = "circSample",
          definition = function(object, value){
            object@label <- value

            if (validObject(object)){
              return(object)
            }
          })


setMethod(f = "label<-",
          signature = "circExperiment",
          definition = function(object, value){
            if (length(value) != length(samples(object))){
              stop("Number of samples do not match number of labels.")
            } else if (any(duplicated(value))){
              stop("Labels must be unique!")
            }
            samples(object) <- lapply(seq_along(samples(object)), function(i){
              .s <- samples(object)[[i]]
              label(.s) <- value[[i]]
              return(.s)
            })
            return(object)
          }
)


setMethod(f = "samples<-",
          signature = "circExperiment",
          definition = function(object, value){
            object@samples <- value

            if (validObject(object)){
              return(object)
            }
          })


#' @importFrom BiocGenerics path
setMethod(f = "path<-",
          signature = "circExperiment",
          definition = function(object, value){
            object@path <- value

            if (validObject(object)){
              return(object)
            }
          })


setMethod(f = "name<-",
          signature = "circExperiment",
          definition = function(object, value){
            object@name <- value

            if (validObject(object)){
              return(object)
            }
          })


############################### Populate objects ###############################

setMethod(
  f = "locateSamples", signature = "circExperiment",
  definition = function(object, organism, gb, firstread.firststrand, paired.end){
    message("Parsing data directory to identify samples")

    all.files <- list.files(path(object), recursive = T, full.names = T)

    # Determine number of samples based on how many Chimeric.out.junction files are present
    chim.index <- grepl("Chimeric.out.junction$", all.files)
    num.samples <- sum(chim.index)
    sample.prefix <- all.files[chim.index] %>% gsub(".Chimeric.out.junction", "", ., fixed = T)

    # split files
    all.files <- lapply(seq_along(sample.prefix), function(i){
      all.files[grepl(sample.prefix[i], all.files)]
    })

    # Generate unique sample names
    sn.split <- strsplit(sample.prefix, "/+")
    dupli <- T
    i<-0
    while(dupli){
      sn <- lapply(sn.split, function(x)x[(length(x)-i):length(x)]) %>% sapply(function(y)paste0(y, collapse="_"))
      i<-i+1
      dupli <- any(duplicated(sn))
    }

    # Construct table of samples
    sampleTable <- tibble::tibble(
      id = sn,
      chim.file = sapply(all.files, function(x)x[grepl("Chimeric.out.junction$", x)]),
      # bam.file = sapply(all.files, function(x)x[grepl("Aligned.*\\.out\\.bam$", x)]),
      log.file = sapply(all.files, function(x)x[grepl("Log.final.out$", x)]),
      sj.file = sapply(all.files, function(x)x[grepl("SJ.out.tab$", x)]),
      count.file = sapply(all.files, function(x)x[grepl("ReadsPerGene.out.tab$", x)]),
      all.required.files = NA,
      PE = NA,
      frfs = NA,
      org = NA,
      geno = NA
    )

    # Check if all required files exist
    sampleTable$chim.file <- ifelse(file.exists(sampleTable$chim.file), sampleTable$chim.file, NA)
    # sampleTable$bam.file <- ifelse(file.exists(sampleTable$bam.file), sampleTable$bam.file, NA)
    sampleTable$log.file <- ifelse(file.exists(sampleTable$log.file), sampleTable$log.file, NA)
    sampleTable$sj.file <- ifelse(file.exists(sampleTable$sj.file), sampleTable$sj.file, NA)
    sampleTable$count.file <- ifelse(file.exists(sampleTable$count.file), sampleTable$count.file, NA)
    sampleTable$all.required.files <- !apply(is.na(sampleTable[,c("chim.file", "log.file", "sj.file", "count.file")]), 1, any)

    ## Characterize samples ##
    # PE
    # if (is.null(paired.end)){
    #   message(" Detecting if reads are single end or paired end ...", appendLF = F)
    #   guess <- sapply(sampleTable$bam.file, guessPE)
    #   if (is.null(guess)){
    #     stop("Could not parse bam header.")
    #   } else {
    #     sampleTable$PE <- guess
    #     message(" Success")
    #   }
    # } else if ((length(paired.end) == 1 | length(paired.end) == nrow(sampleTable)) & is.logical(paired.end)){
    #
    if ((length(paired.end) == 1 | length(paired.end) == nrow(sampleTable)) & is.logical(paired.end)){
      sampleTable$PE <- paired.end
    } else {
      stop("'paired.end' should be a logical vector of length 1 or identical to number of samples.")
    }

    # First read first strand
    if (is.null(firstread.firststrand)){
      stop("Currently, automatic detection of lib prep method is not supported. Please provide valid input.")
    } else if ((length(firstread.firststrand) == 1 | length(firstread.firststrand) == nrow(sampleTable)) & is.logical(firstread.firststrand)){
      sampleTable$frfs <- firstread.firststrand
    } else {
      stop("'firstread.firststrand' should be a logical vector of length 1 or identical to number of samples.")
    }

    if ((length(organism) == 1 & length(gb) == 1) | (length(organism) == nrow(sampleTable) & length(gb) == nrow(sampleTable))){
      sampleTable$org <- organism
      sampleTable$geno <- gb
    } else if (is.null(organism) | is.null(gb)){
      stop("'organism' & 'gb' should be character vectors of length 1 or identical to number of samples.")
    } else {
      stop("please provide input!")
    }

    all.samples <- lapply(seq_along(sampleTable$id), function(i){
      circSample(
        sample.id = sampleTable$id[i],
        organism = sampleTable$org[i],
        gb = sampleTable$geno[i],
        chim.file = sampleTable$chim.file[i],
        # bam.file = sampleTable$bam.file[i],
        log.file = sampleTable$log.file[i],
        sj.file = sampleTable$sj.file[i],
        count.file = sampleTable$count.file[i],
        firstread.firststrand = sampleTable$frfs[i],
        paired.end = sampleTable$PE[i]
      )
    })

    samples(object) <- all.samples
    return(object)
  }
)

setMethod(
  f = "readBSJdata", signature = "circSample",
  definition = function(object, ...){
    message("Importing data from: ", sample.id(object))
    #message(chromosomes)
    if(!(is.numeric(maxGenomicDist) | is.integer(maxGenomicDist)) | !maxGenomicDist > 0){stop("Provide valid maxGenomicDist.")}
    if(!is.logical(onlySpanning)){stop("onlySpanning should be TRUE or FALSE")}

    # message("Chromosomes: ", chromosomes)
    # message("maxGenomicDist: ", maxGenomicDist)
    # message("onlySpanning: ", onlySpanning)

    # Read input chimdata
    chim <- readChimFile(chim.file(object))[,paste0("X", 1:14)]

    # Filter for backsplices junctions only
    same.chr <- chim$X1 == chim$X4
    same.str <- chim$X3 == chim$X6
    acceptorUpstreamOfDonor <- (chim$X3 == "-" & (chim$X5 > chim$X2)) | (chim$X3 == "+" & (chim$X2 > chim$X5))
    dist.below.cutoff <- abs(chim$X2-chim$X5) <= maxGenomicDist
    filters <- same.chr & same.str & acceptorUpstreamOfDonor & dist.below.cutoff

    if(onlySpanning){
      filters <- filters & chim$X7 > -1
    }
    chim <- chim[filters, ]

    if (!is.null(chromosomes)){
      chim <- subset(chim, X1 %in% chromosomes)
    }

    # Determine if data are first read first strand
    if (firstread.firststrand(object)){
      # Swap strand specific information (library preparation protocol dependant)
      chim <- swapStrands(chim)
    }

    if (paired.end(object) & removeBadPairs){
      chim <- chim[checkPEReads(chim),]
      chim$PEok <- TRUE
    } else if (paired.end(object)){
      chim$PEok <- checkPEReads(chim)
    } else {
      chim$PEok <- NA
    }

    # Add logical culumn intended for read filtering
    chim$include.read <- T # Per default, all reads are considered.

    bsj.reads(object) <- chim
    return(object)
  }
)


setMethod(
  f = "readBSJdata", signature = "circExperiment",
  definition = function(object, ...){
    message("Importing data from expriment: ", name(object))

    s <- samples(object)

    if (Sys.info()["sysname"] == "Windows" & cores > 1){
      warning("Sorry, Microsoft Windows doesn't support multicore computing with mclapply!")
      cores = 1L
    }
    s <- mclapply(s, function(o){readBSJdata(object = o, chromosomes = chromosomes, maxGenomicDist = maxGenomicDist, onlySpanning = onlySpanning, removeBadPairs = removeBadPairs)}, mc.cores = cores)

    samples(object) <- s
    return(object)
  }
)


setMethod(
  f = "readLSJdata", signature = "circSample",
  definition = function(object, chromosomes = NULL, ...){
    message("Importing SJ data from: ", sample.id(object))
    # Read SJ data
    data <- readr::read_tsv(
      sj.file(object), col_names = F,
      col_types = readr::cols(
        X1 = readr::col_character(),
        X2 = readr::col_integer(),
        X3 = readr::col_integer(),
        X4 = readr::col_integer(),
        X5 = readr::col_integer(),
        X6 = readr::col_integer(),
        X7 = readr::col_integer(),
        X8 = readr::col_integer(),
        X9 = readr::col_integer()
      ),  progress = F
    )
    colnames(data) <- c("seqnames", "start", "end", "strand", "int.motif", "annotated", "read.count.unique", "read.count.multi", "max.overhang")
    data$strand <- stringi::stri_trans_char(data$strand, "120", "+-*")

    if (!is.null(chromosomes)){
      data <- subset(data, seqnames %in% chromosomes)
    }
    lsj.counts(object) <- GenomicRanges::GRanges(data)
    return(object)
  }
)


setMethod(
  f = "readLSJdata", signature = "circExperiment",
  definition = function(object, chromosomes, cores){
    message("Importing LSJ data from all samples in expriment: ", name(object))

    s <- samples(object)

    if (Sys.info()["sysname"] == "Windows" & cores > 1){
      warning("Sorry, Microsoft Windows doesn't support multicore computing with mclapply!")
      cores = 1L
    }

    s <- mclapply(s, function(o){readLSJdata(o, chromosomes = chromosomes)}, mc.cores = cores)

    samples(object) <- s
    return(object)
  }
)

setMethod(f = "addFilter",
          signature = "circSample",
          definition = function(object, filter, mode){
            if(class(mode) != "character" | !(mode %in% c("strict", "last"))){stop("'mode' should be either 'strict' or 'last'")}
            if(class(filter) != "logical"){stop("Please provide logical vector")}
            if(length(filter) != nrow(bsj.reads(object))){stop("The supplied filter does not match the number of bsj.reads in sample ", sample.id(object))}
            bsj <- bsj.reads(object)
            if(mode == "last"){
              bsj$include.read <- filter
            } else if (mode == "strict") {
              bsj$include.read <- bsj$include.read & filter
            }
            bsj.reads(object) <- bsj

            return(object)
          })

setMethod(f = "addFilter",
          signature = "circExperiment",
          definition = function(object, filter, mode){
            if(class(mode) != "character" | !(mode %in% c("strict", "last"))){stop("'mode' should be either 'strict' or 'last'")}
            if(class(filter) != "list" & !all(sapply(filter, class) == "logical")){stop("Please provide a list of logical vectors.")}
            if(length(filter) != length(sample.id(object))){stop("Error adding filter to experiment: ", name(object), "\n\t- Number of samples: ", length(samples(object)), "\n\t- Number of filters: ", length(filter))}

            s <- samples(object)
            s <- lapply(seq_along(s), function(i){
              addFilter(s[[i]], filter[[i]], mode)
            })
            samples(object) <- s

            return(object)
          })

############################### Annotation ###############################

setMethod(
  f = "compareToKnownJunctions", signature = "circSample",
  definition = function(object, known.junctions){

    message("Comparing to known junctions ...")
    data <- bsj.reads(object)

    # Ensure that junction data previosuly was generated
    if (is.empty(data)){
      stop("Junction data has not been loaded, please run `readBSJdata()` first!")
    }

    # Determine junction reads that overlaps known junctions
    candidates <- addKnownJunctions(cbs = data, kj = known.junctions)

    # # Shift identified backsplice candidates
    # candidates <- shiftAlignment(cbs = candidates)

    bsj.reads(object) <- candidates
    object <- constructBsId(object = object)

    message("... Done")

    return(object)
  }
)


setMethod(
  f = "compareToKnownJunctions", signature = "circExperiment",
  definition = function(object, known.junctions, cores){

    # message("Comparing to known junctions ...", appendLF = F)

    if (Sys.info()["sysname"] == "Windows" & cores > 1){
      warning("Sorry, Microsoft Windows doesn't support multicore computing with mclapply!")
      cores = 1L
    }

    s <- mclapply(samples(object), compareToKnownJunctions, known.junctions, mc.cores = cores)
    #s <- lapply(samples(object), compareToKnownJunctions, known.junctions)
    samples(object) <- s

    # message(" Done")

    return(object)
  }
)


setMethod(f = "generateJunctionMotifs",
          signature = "circSample",
          definition = function(object, genome_seq) {
            bsj <- bsj.reads(object)

            motifs <- dplyr::pull(bsj, bsID)
            motifs <- junctionMotif(bsids = motifs, g = genome_seq)

            bsj.reads(object) <- add_column(bsj, motif = motifs)

            return(object)
          }
)


setMethod(f = "generateJunctionMotifs",
          signature = "circExperiment",
          definition = function(object, genome_seq, cores) {
            s <- samples(object)

            if (Sys.info()["sysname"] == "Windows" & cores > 1){
              warning("Sorry, Microsoft Windows doesn't support multicore computing with mclapply!")
              cores = 1L
            }

            s <- mclapply(s, function(o){generateJunctionMotifs(object = o, genome_seq)}, mc.cores = cores)
            samples(object) <- s
            return(object)
          }
)


setMethod(f = "adjustAlignment",
          signature = "circSample",
          definition = function(object){
            data <- bsj.reads(object)
            # Ensure that junction data previosuly was generated
            if (is.empty(data)){
              stop("Junction data has not been loaded, please run `readBSJdata()` first!")
            }

            if (!all(c("shiftDonorToNearestJ", "shiftAcceptorToNearestJ") %in% colnames(bsj.reads(s1)))){
              stop("Run the first part of the pipeline before running this function")
            }

            # Shift identified backsplice candidates
            candidates <- shiftAlignment(cbs = data)

            bsj.reads(object) <- candidates

            return(object)
          }
)


setMethod(f = "adjustAlignment",
          signature = "circExperiment",
          definition = function(object, cores){
            message("Adjusting alignments for samples in expriment: ", name(object))
            s <- samples(object)

            if (Sys.info()["sysname"] == "Windows" & cores > 1){
              warning("Sorry, Microsoft Windows doesn't support multicore computing with mclapply!")
              cores = 1L
            }

            s <- mclapply(s, function(o){adjustAlignment(o)}, mc.cores = cores)
            samples(object) <- s
            return(object)
          }
)


setMethod(
  f = "summarizeBSJreads", signature = "circSample",
  definition = function(object, applyFilter, ...){
    if (is.empty(bsj.reads(object))){stop("Please load chimeric data first.")}

    if (!"bsID" %in% colnames(bsj.reads(object))){
      message("bsID column not found, adding it now.")
      ## ...
    }

    # Count amount of reads that covers unique backsplice sites
    df <- bsj.reads(object)
    if(applyFilter){
      df <- df[df$include.read,]
    }

    df <- df %>% dplyr::group_by(bsID) %>% dplyr::summarise(count = dplyr::n())

    # Convert to GRanges
    output <- bsid2gr(df$bsID)
    output$count <- df$count

    bsj.counts(object) <- output
    return(object)
  }
)


setMethod(
  f = "summarizeBSJreads", signature = "circExperiment",
  definition = function(object, cores, applyFilter, colUsedForNorm = "X4", ...){
    s <- samples(object)

    if (Sys.info()["sysname"] == "Windows" & cores > 1){
      warning("Sorry, Microsoft Windows doesn't support multicore computing with mclapply!")
      cores = 1L
    }

    s <- mclapply(s, function(o){summarizeBSJreads(o, applyFilter = applyFilter)}, mc.cores = cores)

    # # Read the linear count data and calculate the normalization factor
    # normFactors <- lapply(count.file(o),
    #   read_tsv,
    #   col_names = F,
    #   skip = 4,
    #   col_types = cols(
    #     X1 = col_character(),
    #     X2 = col_double(),
    #     X3 = col_double(),
    #     X4 = col_double()
    #   )
    # ) %>%
    #   lapply(function(x)x[,c("X1", colUsedForNorm)]) %>%
    #   plyr::join_all(by = "X1", type = "left")
    # normFactors <- DESeq2::estimateSizeFactorsForMatrix(counts = normFactors[,-1])
    #
    # s <- mclapply(seq_along(normFactors), function(i){
    #   tmp <- bsj.counts(s[[i]])
    #   tmp$count.norm <- tmp$count / normFactors[i]
    #   return(tmp)
    # })

    samples(object) <- s
    return(object)
  }
)


setMethod(
  f = "circulaR", signature = "circSample",
  definition = function(object, annotationDB){
    # Make annotation database
    kj <- constructSJDB(annotationDB, chromosomes = chromosomes(object))

    # Read sample file
    object <- readBSJdata(object)

    # Annotate known junctons
    object <- compareToKnownJunctions(object, known.junctions)

    # Return all identified junctions, known and new candidates
    return(object)
  }
)


setMethod(
  f = "circulaR", signature = "circExperiment",
  definition = function(object, annotationDB, cores){

    message("Locating and importing sample data for experiment: ", name(object), " ... ", appendLF = FALSE)
    object <- locateSamples(object, organism = "Homo sapiens", gb = "hg38", firstread.firststrand = T)
    object <- readBSJdata(object)
    message(" Done")

    message("Constructing database of known splice junction... ", appendLF = FALSE)
    kj <- constructSJDB(annotationDB, chromosomes = chromosomes(object))
    message(" Done")


    message("Comparing identified splice junction candidates with known splice junctions... ", appendLF = FALSE)
    object <- compareToKnownJunctions(object, known.junctions = kj)
    message(" Done")


    message("Running circulaR pipeline... ", appendLF = FALSE)

    if (Sys.info()["sysname"] == "Windows" & cores > 1){
      warning("Sorry, Microsoft Windows doesn't support multicore computing with mclapply!")
      cores = 1L
    }

    s <- mclapply(samples(object), circulaR, mc.cores = cores)
    samples(object) <- s
    message(" Done")

    return(object)
  }
)


setMethod(f = "constructBsId",
          signature = "circSample",
          definition = function(object){
            if (is.empty(bsj.reads(object))){
              stop("Backsplice read data is missing, please run 'readBSJdata()' first!")
            }
            bsid <- paste(gb(object), bsj.reads(object)$X1, bsj.reads(object)$X2, bsj.reads(object)$X5, bsj.reads(object)$X3, sep = ":")

            bsj.reads(object)$bsID <- bsid

            return(object)
          })


setMethod(f = "constructBsId",
          signature = "circExperiment",
          definition = function(object, cores){
            if (Sys.info()["sysname"] == "Windows" & cores > 1){
              warning("Sorry, Microsoft Windows doesn't support multicore computing with mclapply!")
              cores = 1L
            }

            bsj.reads(object) <- parallel::mclapply(X = samples(object), FUN = constructBsId, mc.cores = cores)

            return(object)
          })

################################ Stats #################################

setMethod(f = "alignmentStats",
          signature = "circSample",
          definition = function(object){
            output <- STARlogParser(log.file(object))
            names(output) <- c("stat", sample.id(object))
            return(tibble::as_tibble(output))
          }
)


setMethod(f = "alignmentStats",
          signature = "circExperiment",
          definition = function(object, out_type){
            output <- lapply(log.file(object), STARlogParser)
            names(output) <- sample.id(object)
            output <- lapply(seq_along(output), function(x){
              setNames(output[[x]], c("stat", names(output)[x]))
            })
            output <- plyr::join_all(output, by = "stat", type = "left")

            if (out_type == "long"){
              output <- pivot_longer(
                data = output, cols = -stat,
                names_to = "sample", values_to = "value"
              )
            } else if (out_type != "wide") {
              stop('The variable "out_type" could not be interpretted.\n',
                   'Please set it to either "wide" or "long"!')
            }
            return(tibble::as_tibble(output))
          }
)


setMethod(
  f = "bsjStats", signature = "circSample",
  definition = function(object, ...){
    # message("Compiling BSJ statistics: ", sample.id(object))

    # Determine the total number of chimeric reads
    # zz <- file(chim.file(object), open="rb")
    # chimTotal <- 0L
    # while (length(tmp <- readBin(zz, "raw", n = 1e4)) > 0) {
    #   chimTotal <- chimTotal + sum(tmp == as.raw(10L))
    # }
    # close(zz)
    chimTotal <- alignmentStats(object) %>% subset(stat == "Number of chimeric reads")
    chimTotal <- as.numeric(chimTotal[1,2])

    # Determine the number of backspliced reads
    bsjTotal <- nrow(bsj.reads(object))

    # Determine junction types
    bsjJuncTypes <- table(bsj.reads(object)$X7)

    # construct output tibble
    output <- tibble::tibble(
      stat = c(
        "Number of chimeric reads", "Number of BSJ reads",
        "BSJ reads %", "Number of splices: GT/AG",
        "Number of splices: CT/AC", "Number of splices: other"
      ),
      value = c(chimTotal, bsjTotal, bsjTotal*100/chimTotal, bsjJuncTypes["1"], bsjJuncTypes["2"], bsjJuncTypes["0"])
    )
    names(output)[2] <- sample.id(object)
    return(output)
  }
)


setMethod(f = "bsjStats",
          signature = "circExperiment",
          definition = function(object, out_type) {
            output <- lapply(samples(object), bsjStats)
            output <- plyr::join_all(output, by = "stat", type = "left")
            if (out_type == "long") {
              output <- tidyr::pivot_longer(
                data = output, cols = -stat,
                names_to = "sample", values_to = "value"
              )
              } else if (out_type != "wide") {
                stop('The variable "out_type" could not be interpretted.\n',
                     'Please set it to either "wide" or "long"!')
              }

            return(tibble::as_tibble(output))
          }
)


################################ Vizualization #################################

setMethod(
  f = "vizJunctions", signature = "circSample",
  definition = function(object, ...) {

    if (circulaR::is.empty(bsj.counts(object)))
      stop("Please load chimeric data.")
    if (circulaR::is.empty(lsj.counts(object)))
      stop("Please load linear splice data.")
    if (is.null(symbol) & is.null(range))
      stop("Please tell me what to plot.")
    if (class(db) != "EnsDb")
      stop("Please provide an EnsDb-object.")
    if (!is.null(xlim) & !(class(xlim) %in% c("integer", "numeric") & length(xlim) == 2))
      stop("Please provide correctly formatted xlim.")
    if (!device %in% c("bmp", "jpeg", "png", "tiff", "pdf", "svg"))
      stop("Please provide a supported device!")

    if ((length(symbol) >= 1 & length(range) >= 1) & length(symbol) != length(range)){
      stop("The dimensions of defined ranges and symbols does not match!")
    } else {
      dimensions <- max(c(length(symbol), length(range)))
    }

    for (i in 1:dimensions) {
      s <- symbol[i]
      r <- range[i]
      # Check if symbol is a valid genename in the suppled db
      if (!is.null(s)) {

        # Ensure that the correct gene symbol have been determined
        if (!s %in% ensembldb::genes(db)$gene_name) {
          stop("Cannot find ", s, " in the supplied database.")
        }

        # Construct genomic range corresponding to gene
        symbol_r <- ensembldb::genes(db, filter = ~ "gene_name" == s)

        if (is.null(r)) {
          r <- symbol_r
          if (isTRUE(chromosomes)) {
            message("Vizulising backsplice and linear splice junctions using standard chromosomes")
            r <- GenomeInfoDb::keepStandardChromosomes(x = r, pruning.mode = "coarse")
          } else if (class(chromosomes) %in% c("character", "numeric")) {
            r <- GenomeInfoDb::keepSeqlevels(x = r, pruning.mode = "coarse", value = chromosomes)
          }
        } else {
          overlaps <- suppressWarnings(IRanges::subsetByOverlaps(x = symbol_r, ranges = r))
          if (length(overlaps) == 0)
            stop("The specified range does not overlap with the specified symbol!")

        }
      } else {
        if (length(reduce(r)) != 1)
          stop("Looks like multiple genomic regions.")

        # Construct genomic range corresponding to gene
        s <- IRanges::subsetByOverlaps(ensembldb::genes(db), r)$gene_name

        if (length(s) != 1)
          stop("None or multiple gene symbols where identified. Please specify the gene symbol!")
      }

      gr <- suppressWarnings(ensembldb::getGeneRegionTrackForGviz(
        db,
        filter = ~ "gene_name" == s))
      if (isTRUE(chromosomes)) {
        gr <- GenomeInfoDb::keepStandardChromosomes(x = gr, pruning.mode = "coarse")
      } else if (class(chromosomes) %in% c("character", "numeric")) {
        gr <- GenomeInfoDb::keepSeqlevels(x = gr, value = chromosomes, pruning.mode = "coarse")
      } else {
        # Ensure that seqlevels matches
        gr <- GenomeInfoDb::keepSeqlevels(x = gr, value = seqlevels(r), pruning.mode = "coarse")
        r <- GenomeInfoDb::keepSeqlevels(x = r, value = seqlevels(gr), pruning.mode = "coarse")
      }

      gr.tr <- Gviz::GeneRegionTrack(gr, name = paste0(s, " models"))

      # Get the splicing data
      lsj <- suppressWarnings(subsetByOverlaps(lsj.counts(object), r))
      lsj <- lsj[strand(lsj) == strand(r)]

      bsj <- suppressWarnings(subsetByOverlaps(bsj.counts(object), r))

      if (!is.null(xlim)) {
        start(r) <- min(xlim)
        end(r) <- max(xlim)
      }

      if (onlyJunctionsWithinRange) {
        lsj <- suppressWarnings(subsetByOverlaps(lsj, r, type = "within"))
        bsj <- suppressWarnings(subsetByOverlaps(bsj, r, type = "within"))
      }

      if (!is.null(validExonModels)) {
        ex.gr <- lapply(seq_along(validExonModels), function(i) {
          output <- validExonModels[[i]]
          output$groupid <- strsplit(names(validExonModels)[i],"_")[[1]][2]
          return(output)
        })
        ex.tr <- Gviz::AnnotationTrack(do.call(c, ex.gr), name = "Exon Models", shape = "box", group = do.call(c, ex.gr)$groupid)
      }

      # Construct the tracks
      sj.tr <- SJTrack(sj = lsj)
      bsj.tr <- BSJTrack(bsj = bsj)

      tr.list <- c(bsj.tr, gr.tr, sj.tr)

      # Setup tracks for plotting
      if (as.character(strand(r)) == "+") {
        if (!is.null(validExonModels)) {
          # Insert before lastposition
          tr.list <- c(
            tr.list[1:(length(tr.list) - 1)],
            ex.tr,
            tr.list[length(tr.list)]
          )
        }
        if (!is.null(deducedCoverage)) {
          tr.list <- c(
            tr.list[1:(length(tr.list) - 1)],
            deducedCoverage,
            tr.list[length(tr.list)]
          )
        }
      }
      if (as.character(strand(r)) == "-") {
        if (!is.null(validExonModels)) {
          # Insert at 2nd position
          tr.list <- c(
            tr.list[1],
            ex.tr,
            tr.list[2:length(tr.list)]
          )
        }
        if (!is.null(deducedCoverage)) {
          # Insert at 2nd position
          tr.list <- c(
            tr.list[1],
            deducedCoverage,
            tr.list[2:length(tr.list)]
          )
        }
      }

      header <- paste0(paste(s, collapse = "; "), ", ", seqnames(r), ":", start(r), "-", end(r), ":", strand(r))
      si <- rep(1, length(tr.list))
      si[sapply(tr.list, names) == "Exon Models"] <- 0.5

      if (is.null(path)) {
        Gviz::plotTracks(tr.list, from = start(r), to = end(r), sizes = si, main = header, cex.main = 1, fontfamily = "Arial")
      } else {
        suffix <- paste0(sample.id(object), "_", s, "_", r, ".", device)
        if (is.null(prefix)) {
          fn <- file.path(path, suffix)
        } else {
          fn <- file.path(path, paste0(prefix, suffix))
        }
        do.call(device, args = list(filename = fn, ...))
        Gviz::plotTracks(tr.list, from = start(r), to = end(r), sizes = si, main = header, cex.main = 1, fontfamily = "Arial")
        dev.off()
      }
    }
  }
)

setMethod(
  f = "vizJunctions", signature = "circExperiment",
  definition = function(object, cores, ...) {

    if (is.empty(bsj.counts(object)))
      stop("Please load chimeric data.")
    if (is.empty(lsj.counts(object)))
      stop("Please load linear splice data.")
    if (is.null(symbol) & is.null(range))
      stop("Please tell me what to plot.")
    if (class(db) != "EnsDb")
      stop("Please provide an EnsDb-object.")
    if (!is.null(xlim) & !(class(xlim) %in% c("integer", "numeric") & length(xlim) == 2))
      stop("Please provide correctly formatted xlim.")
    if (!device %in% c("bmp", "jpeg", "png", "tiff", "pdf", "svg"))
      stop("Please provide a supported device!")
    if (is.null(path))
      stop("Please provide an output path!")

    message("Vizualizing linear and circular splice coverage of expriment: ", name(object))

    s <- samples(object)

    if (Sys.info()["sysname"] == "Windows" & cores > 1) {
      warning("Sorry, Microsoft Windows doesn't support multicore computing with mclapply!")
      cores = 1L
    }

    for (o in samples(object)){
      vizJunctions(object = o, symbol = symbol, range = range, db = db, xlim = xlim, device = device, path = path, cores = cores)
    }
  }
)



############################### General functions ###############################

setMethod(f = "is.empty",
          definition = function(object){
            empty <- FALSE
            if (is.null(object)){
              empty <- TRUE
            } else if (identical(object, character())){
              empty <- TRUE
            } else if (identical(object, numeric())){
              empty <- TRUE
            } else if (identical(object, integer())){
              empty <- TRUE
            } else if (identical(object, logical())){
              empty <- TRUE
            } else if (length(object) == 0){
              empty <- TRUE
            }
            # } else if (is.na(object)){  ### causes length > 1 WARININGs also i flist of multiple values has one NA, then all list is NA!
            #   empty <- TRUE
            # }

            return(empty)
          })


setMethod(f = "constructSJDB",
          definition = function(annotationDB, force){
            message("Retrieving known splice junctions. ", appendLF = F)

            # Get information on supplied DB
            metadata <- RSQLite::dbReadTable(annotationDB@ensdb, "metadata")

            # construct path to local version
            p2db <- file.path(
              "~/.circulaR/",
              gsub(" ", "_", subset(metadata, name == "Organism")$value),
              subset(metadata, name == "Db type")$value,
              subset(metadata, name == "ensembl_version")$value
            )

            if (file.exists(p2db) & !force){
              message(" Reading from cache ...", appendLF = FALSE)
              known.junctions <- readRDS(file = file.path(p2db, "knownJuncs.RData.gz"))
              message(" Done")
            } else {
              message(" Building splice junction database ...", appendLF = FALSE)

              # Determine all known junction
              junctions <- suppressMessages(getKnownJunctions(annotationDB))

              # Combine splice junctions and transcript ends
              message("  Compiling information on unique junction sites (This will take some time!)... ", appendLF = F)
              known.junctions <- junctions %>% as.data.frame %>%
                dplyr::group_by(seqnames, start, end, strand) %>%
                dplyr::summarise(
                  type = paste(unique(type), collapse = "|"),
                  g_id = paste(unique(ensembl_gene_id), collapse = "|"),
                  tx_id = paste(unique(ensembl_transcript_id), collapse = "|")
                )
              message(" Done")

              known.junctions$jID <- paste0("J",1:nrow(known.junctions))
              known.junctions <- GenomicRanges::GRanges(known.junctions)

              # Writing to local file
              message(" Saving to local cache ...", appendLF = FALSE)
              if (dir.exists(dirname(p2db)) & force){
                unlink(p2db)
                dir.create(path = p2db, recursive = TRUE)
              }
              if (!dir.exists(p2db)){
                dir.create(path = p2db, recursive = TRUE)
              }
              saveRDS(object = known.junctions, file = file.path(p2db, "knownJuncs.RData.gz"), compress = "gzip")
              message(" Done.")
            }

            return(known.junctions)
          })


setMethod(f = "getKnownJunctions",
          signature = "TxDb",
          definition = function(annotationDB){

            colN <- c(seqname = "TXCHROM",
                      strand = "TXSTRAND",
                      exon_id = "EXONID",
                      start = "EXONSTART",
                      end = "EXONEND",
                      gene_id = "GENEID")

            message("Retrieving all transcript IDs from DB ... ", appendLF = F)
            all.tx <- GenomicFeatures::transcripts(annotationDB)$tx_name
            message(" Done")

            message("Extracting known splice junctions & transcript termini ... ", appendLF = F)
            annot <- suppressMessages(
              AnnotationDbi::select(annotationDB,
                                    keys = all.tx,
                                    columns = colN,
                                    keytype = "TXNAME")
            )

            known.junctions <- c(
              GRanges(seqnames = annot[,colN["seqname"]],
                      ranges = IRanges(start = annot[,colN["start"]]-1, width = 1),
                      strand = annot[,colN["strand"]],
                      type = ifelse(annot[,colN["strand"]] == "+" | annot[,colN["strand"]] == 1, "acceptor", "donor"),
                      ensembl_gene_id = annot[,colN["gene_id"]],
                      ensembl_transcript_id = annot$TXNAME),
              GRanges(seqnames = annot[,colN["seqname"]],
                      ranges = IRanges(start = annot[,colN["end"]]+1, width = 1),
                      strand = annot[,colN["strand"]],
                      type = ifelse(annot[,colN["strand"]] == "+" | annot[,colN["strand"]] == 1, "donor", "acceptor"),
                      ensembl_gene_id = annot[,colN["gene_id"]],
                      ensembl_transcript_id = annot$TXNAME)
            )
            message(" Done")

            return(known.junctions)
          }
)


setMethod(f = "getKnownJunctions",
          signature = "EnsDb",
          definition = function(annotationDB){
            colN <- c(seqname = "SEQNAME",
                      strand = "SEQSTRAND",
                      exon_id = "EXONID",
                      start = "EXONSEQSTART",
                      end = "EXONSEQEND",
                      gene_id = "GENEID")

            message("Retrieving all transcript IDs from DB ... ", appendLF = F)
            all.tx <- ensembldb::transcripts(annotationDB)$tx_name
            message(" Done")

            message("Extracting known splice junctions & transcript termini ... ", appendLF = F)
            annot <- suppressMessages(
              AnnotationDbi::select(annotationDB,
                                    keys = all.tx,
                                    columns = colN,
                                    keytype = "TXNAME")
            )

            known.junctions <- c(
              GRanges(seqnames = annot[,colN["seqname"]],
                      ranges = IRanges(start = annot[,colN["start"]]-1, width = 1),
                      strand = annot[,colN["strand"]],
                      type = ifelse(annot[,colN["strand"]] == "+" | annot[,colN["strand"]] == 1, "acceptor", "donor"),
                      ensembl_gene_id = annot[,colN["gene_id"]],
                      ensembl_transcript_id = annot$TXNAME),
              GRanges(seqnames = annot[,colN["seqname"]],
                      ranges = IRanges(start = annot[,colN["end"]]+1, width = 1),
                      strand = annot[,colN["strand"]],
                      type = ifelse(annot[,colN["strand"]] == "+" | annot[,colN["strand"]] == 1, "donor", "acceptor"),
                      ensembl_gene_id = annot[,colN["gene_id"]],
                      ensembl_transcript_id = annot$TXNAME)
            )
            message(" Done")

            return(known.junctions)
          }
)

setMethod(f = "bsReadCoverage",
          signature = "circSample",
          definition = function(object, bsid){
            df <- subset(bsj.reads(object), bsID == bsid)

            output <- calcCoverage(df, asGRanges = T)
            return(output)
          }
)

setMethod(
  f = "c2l.ratio", signature = "circSample",
  definition = function(object, ...){

    # Find linear spliced reads, that use same donor/acceptor
    bsid.gr <- bsid2gr(bsid)

    first <- IRanges::subsetByOverlaps(IRanges::resize(lsj.counts(object), width = 1, fix = "end"), IRanges::resize(bsid.gr, width = 1, fix = "start"))
    second <- IRanges::subsetByOverlaps(IRanges::resize(lsj.counts(object), width = 1, fix = "start"), IRanges::resize(bsid.gr, width = 1, fix = "end"))

    if(onlyAnnotated){
      first <- subset(first, annotated != 0)
      second <- subset(second, annotated != 0)
    }

    lin.counts <- mean(c(sum(first$read.count.unique), sum(second$read.count.unique)), na.rm = T)
    circ <- subset(bsj.counts(object), bsID == bsid)

    #if(is.na(lin.counts)){message("No linear splicing!")}
    #if(length(circ) == 0){message("No circular counts!")}

    r <- circ$count / lin.counts

    output <- ifelse(length(r) == 0, NA, log2(r))

    return(output)
  }
)
