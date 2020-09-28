############################### circulaR Sample class definition ###############################

#' @title S4 class for storing circRNA relevant information.
#' @name circSample
#'
#' @slot sample.id The id of the sample, specified by the input sample file naming (e.g. `aligned/[sample]/Chimeric.out.junction`, where `[sample]` are the applied id.), `sample.id``. Must be unique!
#' @slot organism The organism of which the smaple originates from.
#' @slot gb The specific genome build used for alining the sample.
#' @slot chim.file File location for the sample `Chimeric.out.junction` alignment output.
#' @slot bam.file File locaion for the sample `Aligned.X.out.bam` alignment output, where X specifies a sorting option during alignment by the STAR algorithm. #[PROPPER REFS?]
#' @slot log.file File location for the sample `Log.final.out` alignment run log.
#' @slot sj.file File location for the sample `SJ.out.tab` alignment output.
#' @slot firstread.firststrand Library preparation protocol, relative to library strandedness. If a first_read-first_strand  protocol is used (e.g. dUTP), then `firstread.firststrand` is `TRUE`, if first_read-second_strand or unstranded then it is `FALSE`.
#' @slot paired.end Flag used for describing if whether the sample sequnecing data are paired-end or single-end. Must be logical!
#' @slot bsj.reads Annotated and shifted candidate junctions
#' @slot bsj.counts The number of reads covering a backsplice junction for each detected basckplice junction CANDIDATES?.
#' @slot lsj.counts The input linear splice junction data.
#' @slot label An optional name for a given sample, if no `label` is provided, `sample.id`` will be used. Need to check if labels are unique when constructing an experiment!!!!
#'
#' @return circSample-object
#' @examples
#' circSample(
#'     sample.id = "A549",
#'     organism = "Homo sapiens",
#'     gb = "hg38",
#'     chim.file = "Chimeric.out.junction",
#'     bam.file = "Aligned.sortedByCoord.out.bam",
#'     log.file = "Log.final.out",
#'     sj.file = "SJ.out.tab",
#'     firstread.firststrand = TRUE,
#'     paired.end = TRUE
#' )
#' @export circSample
#' @exportClass circSample
circSample <- setClass(
  Class = "circSample",
  slots = c(
    sample.id = "character",
    organism = "character",
    gb = "character",
    chim.file = "character",
    bam.file = "character",
    log.file = "character",
    sj.file = "character",
    count.file = "character",
    firstread.firststrand = "logical",
    paired.end = "logical",
    bsj.reads = "data.frame",
    bsj.counts = "GRanges",
    lsj.counts = "GRanges",
    label = "character"
  )
)


setValidity(Class = "circSample",
            method = function(object){
              msg <- NULL
              valid <- TRUE

              # Define basic slots
              core.slots <- c("sample.id",
                              "organism",
                              "gb",
                              "chim.file",
                              "log.file",
                              "sj.file",
                              "count.file",
                              "firstread.firststrand",
                              "paired.end")

              # Note basic slots which have wrong dimensions
              invalids <- c()
              empty.slots <- lapply(core.slots, function(s){
                if (length(slot(object, s)) < 1){
                  invalids <- c(invalids, s)

                  return(invalids)
                }

              }) %>% do.call(rbind, .) %>% as.vector

              invalids <- c()
              long.slots <- lapply(core.slots, function(s){
                if (length(slot(object, s)) > 1){
                  invalids <- c(invalids, s)

                  return(invalids)
                }

              }) %>% do.call(rbind, .) %>% as.vector

              # Generate error message for empty basic slots
              if (!is.empty(empty.slots)){
                msg <- c(msg, paste(empty.slots, "must be specified!"))
                valid <- FALSE

              # Generate error message for basic slots with multiple values
              } else if (!is.empty(long.slots)){
                msg <- c(msg, paste(long.slots, "have more than one value!"))
                valid <- FALSE
              }

              # Same goes for the label slot
              if (length(label(object)) > 1){
                msg <- c(msg, "label has more than one value!")
                valid <- FALSE
              }

              if (valid) TRUE else msg
            }
)


setMethod(f = "show",
          signature = "circSample",
          definition = function(object){
            # ### Show label instead of sample.id if it is not empty...
            # msg <- ifelse(test = is.empty(label(object)),
            #               yes = paste0(class(object),
            #                            ": ", sample.id(object), "\n"),
            #                no = paste0(class(object),
            #                            ": ", label(object), "\n"))
            # cat(msg)
            cat(class(object), ": ", sample.id(object), "\n", sep = "")


            cat(" Organism: ", organism(object), "\n", sep = "")
            cat(" Genome: ", gb(object), "\n", sep = "")
            cat(ifelse (paired.end(object)," Paired end reads", " Single end read"), "\n", sep = "")

            cat("Candidate backsplice data:\n")
            if (!is.empty(bsj.reads(object))){
              n <- nrow(bsj.reads(object))
              f <- length(which(bsj.reads(object)$include.read))
              cat(" Chimeric reads detected:", n, paste0("(included by current filter: ", round(f*100/n, digits=1),"%)"), "\n")
            } else {
              cat(" Chimeric reads detected: <data not loaded>\n")
            }

            if (!is.empty(bsj.counts(object))){
              cat(" Unique backsplice junctions:", length(bsj.counts(object)), "\n")
            } else {
              cat(" Unique backsplice junctions: <data not summarized>\n")
            }

            cat("Linear splice data:\n")
            if (!is.empty(lsj.counts(object))){
              cat(" Linear splice junctions detected:", length(lsj.counts(object)), "\n")
            } else {
              cat(" Linear splice junctions detected: <data not loaded>\n")
            }

            cat("\n")
          }
)


############################### circulaR Experiment class definition ###############################

#' S4 class for representing an experiment consisting of multiple samples (i.e. circSample objects).
#'
#' @slot path character Path to the directory containing output from the STAR aligner.
#' @slot samples list A list of circSamples.
#' @slot name character A free text field holding the name of the experiment.
#'
#' @return circExperiment-object
#'
#' @importFrom parallel mclapply
#' @examples
#' circExperiment(
#'     path = "path_to_data",
#'     name = "TEST experiment"
#' )
#'
#' @export circExperiment
#' @exportClass circExperiment

circExperiment <- setClass(
  Class = "circExperiment",
  slots = c(
    path = "character",
    samples = "list",
    name = "character"
  )
)


setValidity(Class = "circExperiment",
            method = function(object){
              msg <- NULL
              valid <- TRUE

              # Check required slots.
              if (is.empty(path(object))){
                msg <- c(msg, "Data directory must be specified!")
                valid <- FALSE
              } else if (!dir.exists(path(object))){
                msg <- c(msg, "Data directory does not exist!")
                valid <- FALSE
              }

              if (valid) TRUE else msg

            }
)


setMethod(f = "show",
          signature = "circExperiment",
          definition = function(object){
            cat(class(object), ": ", name(object), "\n", sep = "")
            cat(" Path: ", path(object), "\n", sep = "")

            if(is.empty(samples(object))){
              cat("\nNumber of samples: <undetermined>\n")
            } else {
              sn <- sapply(samples(object), sample.id)

              cat("\nNumber of samples: ", length(sn), "\n")
              if (length(sn) > 0){
                cat(" Sample names:\t", ifelse(length(sn)>5, paste0(c(sn[1:5], "..."), collapse = ", "), paste0(sn, collapse = ", ")), "\n", sep = "")
              }
            }

            cat("\n")
          }
)

