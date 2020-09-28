################################ Accessor methods #################################

#' Accessor for sample.id
#'
#' circSample and circExperiment class accessors for showing the sample.id(s) of cicSample object(s)
#'
#' @param object circSample or circExperiment object
#'
#' @return Character vector with sample.id(s)
#'
#' @examples
#' # The sample.id of a sample object
#' sample.id(sample.object)
#'
#' # The sample.id in each sample of an experiment
#' sample.id(experiment.object)
#'
#' @export
setGeneric(name = "sample.id",
           def = function(object)
             standardGeneric("sample.id"))


#' Accessor for organism
#'
#' circSample and circExperiment class accessors for showing the organism(s) of cicSample object(s)
#'
#' @param object circSample or circExperiment object
#'
#' @return Character vector with organism(s)
#'
#' @examples
#' # The organism of a sample object
#' organism(sample.object)
#'
#' # The organism in each sample of an experiment
#' organism(experiment.object)
#' @importFrom GenomeInfoDb organism
#' @export
setGeneric(name = "organism",
           def = function(object)
             standardGeneric("organism"))


#' Accessor for genome build
#'
#' circSample and circExperiment class accessors for showing the genome build(s) of circSample object(s)
#'
#' @param object circSample or circExperiment object
#'
#' @return Character vector with the genome build(s)
#'
#' @examples
#' # The genome build of a sample object
#' gb(sample.object)
#'
#' # The genome build in each sample of an experiment
#' gb(experiment.object)
#'
#' @export
setGeneric(name = "gb",
           def = function(object)
             standardGeneric("gb"))


#' Accessor for the chimeric file path
#'
#' circSample and circExperiment class accessors for showing the chimeric file path(s) of circSample object(s)
#'
#' @param object circSample or circExperiment object
#'
#' @return Character vector with the chimeric file path(s)
#'
#' @examples
#' # The chimeric file path of a sample object
#' chim.file(sample.object)
#'
#' # The chimeric file path in each sample of an experiment
#' chim.file(experiment.object)
#'
#' @export
setGeneric(name = "chim.file",
           def = function(object)
             standardGeneric("chim.file"))


#' Accessor for the BAM file path
#'
#' circSample and circExperiment class accessors for showing the chimeric file path(s) of circSample object(s)
#'
#' @param object circSample or circExperiment object
#'
#' @return Character vector with the chimeric file path(s)
#'
#' @examples
#' # The BAM file path of a sample object
#' bam.file(sample.object)
#'
#' # The BAM file path of each sample in each sample of an experiment
#' bam.file(experiment.object)
#'
#' @export
setGeneric(name = "bam.file",
           def = function(object)
             standardGeneric("bam.file"))


#' Accessor for the STAR alignment Log file path
#'
#' circSample and circExperiment class accessors for showing the STAR alignment Log file path(s) of circSample object(s)
#'
#' @param object circSample or circExperiment object
#'
#' @return Character vector with the STAR alignment Log file path(s)
#'
#' @examples
#' # The STAR alignment Log file path of a sample object
#' log.file(sample.object)
#'
#' # The STAR alignment Log file path in each sample of an experiment
#' log.file(experiment.object)
#'
#' @export
setGeneric(name = "log.file",
           def = function(object)
             standardGeneric("log.file"))


#' Accessor for the Linear splice junctions file path
#'
#' circSample and circExperiment class accessors for showing the Linear splice junctions file path(s) of circSample object(s)
#'
#' @param object circSample or circExperiment object
#'
#' @return Character vector with the Linear splice junctions file path(s)
#'
#' @examples
#' # The Linear splice junctions file path of a sample object
#' sj.file(sample.object)
#'
#' # The Linear splice junctions file path in each sample of an experiment
#' sj.file(experiment.object)
#'
#' @export
setGeneric(name = "sj.file",
           def = function(object)
             standardGeneric("sj.file"))


#' Accessor for the linear read count data file path
#'
#' circSample and circExperiment class accessors for showing the linear read count data file path(s) of circSample object(s)
#'
#' @param object circSample or circExperiment object
#'
#' @return Character vector with the linear read count data file path(s)
#'
#' @examples
#' # The linear read count data file path of a sample object
#' count.file(sample.object)
#'
#' # The Linear splice junctions file path in each sample of an experiment
#' count.file(experiment.object)
#'
#' @export
setGeneric(name = "count.file",
           def = function(object)
             standardGeneric("count.file"))


#' Accessor for the strandedness protocol of sample library preparation
#'
#' circSample and circExperiment class accessors for showing whether the sample are prepared using e.g. dUTP- or unstranded protocol
#'
#' @param object circSample or circExperiment object
#'
#' @return Logical vector with a boolean value stating if the library preparationion protocol yielded firstread-firststrand reads
#'
#' @examples
#' # The logical value whether the protocol yielded firstread-firststrand reads in a sample object
#' firstread.firststrand(sample.object)
#'
#' # The logical value whether the protocol yielded firstread-firststrand reads in each sample of an experiment
#' firstread.firststrand(experiment.object)
#'
#' @export
setGeneric(name = "firstread.firststrand",
           def = function(object)
             standardGeneric("firstread.firststrand"))


#' Accessor for the library preparation protocol realtive to end-pairing
#'
#' circSample and circExperiment class accessors for showing whether the sample are prepared as paired-end or single-end
#'
#' @param object circSample or circExperiment object
#'
#' @return Logcal vector with a boolean value stating if the library preparationion protocol yielded paired-end reads
#'
#' @examples
#' # The logical value whether the protocol yielded paired-end reads in a sample object
#' paried.end(sample.object)
#'
#' # The logical value whether the protocol yielded paired-end reads in each sample of an experiment
#' paried.end(experiment.object)
#'
#' @export
setGeneric(name = "paired.end",
           def = function(object)
             standardGeneric("paired.end"))


#' Accessor for the backsplice junction read data
#'
#' circSample and circExperiment class accessors for showing the filtered backsplice junction read data of a circSample object(s)
#'
#' @param object circSample or circExperiment object
#' @param returnAs (only valid for circSample) Whether to return data as a list of samples or compiled into a single long table with a sample.id column
#'
#' @return Datatable(s) of filtered backsplice junction read data
#'
#' @examples
#' # A datatable of the filtered backsplice junction read data in a sample object
#' bsj.reads(sample.object)
#'
#' # A datatable of the filtered backsplice junction read data in each sample of an experiment
#' bsj.reads(experiment.object)
#'
#' @export
setGeneric(name = "bsj.reads",
           def = function(object, ...)
             standardGeneric("bsj.reads"))


#' Accessor for the backsplice junction count data
#'
#' circSample and circExperiment class accessors for showing the backsplice junction count data
#'
#' @param object circSample or circExperiment object
#' @param returnAs Determines if output should be list (default) or GRanges (GRangesList for circExperiment objects)
#' @param use.names Boolean. If TRUE, the list being returned will be named using sample IDs
#'
#' @return List of datatable(s) with backsplice junction count data
#'
#' @examples
#' # A datatable of the backsplice junction count data in a sample object
#' bsj.counts(sample.object)
#'
#' # A datatable of the backsplice junction count data in each sample of an experiment
#' bsj.counts(experiment.object)
#'
#' @export
setGeneric(name = "bsj.counts",
           def = function(object, returnAs = "list", ...)
             standardGeneric("bsj.counts"))


#' Accessor for summarized data on linear splice sites
#'
#' circSample and circExperiment class accessors for showing the linear splice junction count data
#'
#' @param object circSample or circExperiment object
#'
#' @return Datatable(s) of linear splice junction count data
#'
#' @examples
#' # A datatable of the linear splice junction count data in a sample object
#' lsj.counts(sample.object)
#'
#' # A datatable of the linear splice junction count data in each sample of an experiment
#' lsj.counts(experiment.object)
#'
#' @export
#' @import dplyr
setGeneric(name = "lsj.counts",
           def = function(object)
             standardGeneric("lsj.counts"))


#' Accessor for sample label
#'
#' circSample and circExperiment class accessors for showing the sample label(s)
#'
#' @param object circSample or circExperiment object
#'
#' @return Character vector with sample label(s)
#'
#' @examples
#' # The label of a sample object
#' label(sample.object)
#'
#' # The label in each sample of an experiment
#' label(experiment.object)
#'
#' @export
setGeneric(name = "label",
           def = function(object)
             standardGeneric("label"))


#' Accessor for main directory path
#' circExperiment class accessors for showing the experiment main path
#'
#' @param object circExperiment object
#'
#' @return Character string with the main directory path
#'
#' @examples
#' # The main path to the experiment
#' path(experiment.object)
#'
#' @import BiocGenerics
#' @export
setGeneric(name = "path",
           def = function(object)
             standardGeneric("path"))


#' Accessor for experiment name
#' circExperiment class accessors for showing the experiment name
#'
#' @param object circExperiment object
#'
#' @return Character string with the experiment name
#'
#' @examples
#' # The experiment name
#' name(experiment.object)
#'
#' @export
setGeneric(name = "name",
           def = function(object)
             standardGeneric("name"))

############################### Replace methods ###############################

#' Replacement method for sample.id
#'
#' circSample and circExperiment class replacement methods for replacing the sample.id(s) of cicSample object(s)
#' The replacment value must be valid to enable replacement
#'
#' @param object circSample or circExperiment object
#'
#' @return New class object with replaced sample.id(s)
#'
#' @examples
#' # Replacing sample.id of a sample object
#' sample.id(sample.object) <- sample.id
#'
#' # Replacing sample.id in each sample of an experiment
#' sample.id(experiment.object) <- sample.id.vector
#'
#' @export
setGeneric(name = "sample.id<-",
           def = function(object, value)
             standardGeneric("sample.id<-"))


#' Replacement method for genome build
#'
#' circSample and circExperiment class replacement methods for replacing the genome build(s) of circSample object(s)
#' The replacment value must be valid to enable replacement
#'
#' @param object circSample or circExperiment object
#'
#' @return New class object with replaced genome build(s)
#'
#' @examples
#' # Replacing genome build of a sample object
#' gb(sample.object) <- new.gb
#'
#' # Replacing genome build in each sample of an experiment
#' gb(experiment.object) <- new.gb.vector
#'
#' @export
setGeneric(name = "gb<-",
           def = function(object, value)
             standardGeneric("gb<-"))


#' Replacement method for the chimeric file path
#'
#' circSample class replacement methods for replacing the chimeric file path
#' The replacment value must be valid to enable replacement
#'
#' @param object circSample object
#'
#' @return New class object with replaced chimeric file path
#'
#' @examples
#' # Replacing chimeric file path of a sample object
#' chim.file(sample.object) <- new.path
#'
#' @export
setGeneric(name = "chim.file<-",
           def = function(object, value)
             standardGeneric("chim.file<-"))


#' Replacement method for the BAM file path
#'
#' circSample class replacement methods for replacing the chimeric file path
#' The replacment value must be valid to enable replacement
#'
#' @param object circSample object
#'
#' @return New class object with replaced BAM file path
#'
#' @examples
#' # Replacing BAM file path of a sample object
#' bam.file(sample.object) <- new.path
#'
#' @export
setGeneric(name = "bam.file<-",
           def = function(object, value)
             standardGeneric("bam.file<-"))


#' Replacement method for the STAR alignment Log file path
#'
#' circSample class replacement methods for replacing the STAR alignment Log file path
#' The replacment value must be valid to enable replacement
#'
#' @param object circSample bject
#'
#' @return New class object with replaced STAR alignment Log file path
#'
#' @examples
#' # Replacing STAR alignment Log file path of a sample object
#' log.file(sample.object) <- new.path
#'
#' @export
setGeneric(name = "log.file<-",
           def = function(object, value)
             standardGeneric("log.file<-"))


#' Replacement method for the Linear splice junctions file path
#'
#' circSample class replacement methods for replacing the Linear splice junctions file path
#' The replacment value must be valid to enable replacement
#'
#' @param object circSample object
#'
#' @return New class object with replaced Linear splice junctions file path
#'
#' @examples
#' # Replacing Linear splice junctions file path of a sample object
#' sj.file(sample.object) <- new.path
#'
#' @export
setGeneric(name = "sj.file<-",
           def = function(object, value)
             standardGeneric("sj.file<-"))

#' Replacement method for the linear read count data file path
#'
#' circSample class replacement methods for replacing the linear read count data file path
#' The replacment value must be valid to enable replacement
#'
#' @param object circSample object
#'
#' @return New class object with replaced linear read count data file path
#'
#' @examples
#' # Replacing linear read count data file path of a sample object
#' count.file(sample.object) <- new.path
#'
#' @export
setGeneric(name = "count.file<-",
           def = function(object, value)
             standardGeneric("count.file<-"))


#' Replacement method for the strandedness protocol of sample library preparation
#'
#' circSample and circExperiment class replacement methods for replacing the boolean value stating if the library preparationion protocol yielded firstread-firststrand reads
#' The replacment value must be valid to enable replacement
#'
#' @param object circSample or circExperiment object
#'
#' @return New class object with replaced logical firstread-firststrand vector
#'
#' @examples
#' # Replacing boolean value of whether reads are firstread-firststrand in a sample object
#' firstread.firststrand(sample.object) <- new.logical
#'
#' # Replacing boolean value whether reads are firstread-firststrand in each sample of an experiment
#' firstread.firststrand(experiment.object) <- new.logical.vector
#'
#' @export
setGeneric(name = "firstread.firststrand<-",
           def = function(object, value)
             standardGeneric("firstread.firststrand<-"))


#' Replacement method for the library preparation protocol realtive to end-pairing
#'
#' circSample and circExperiment class replacement methods for replacing the boolean value stating if the library preparationion protocol yielded paired-end reads
#' The replacment value must be valid to enable replacement
#'
#' @param object circSample or circExperiment object
#'
#' @return New class object with replaced logical paired-end vector
#'
#' @examples
#' # Replacing boolean value whether reads are paired-end reads in a sample object
#' paried.end(sample.object) <- new.logical
#'
#' # Replacing logical value whether reads are paired-end reads in each sample of an experiment
#' paried.end(experiment.object) <- new.logical.vector
#'
#' @export
setGeneric(name = "paired.end<-",
           def = function(object, value)
             standardGeneric("paired.end<-"))


#' Replacement method for the backsplice junction read data
#'
#' circSample and circExperiment class replacement methods for replacing the filtered backsplice junction read data of a circSample object(s)
#' The replacment value must be valid to enable replacement
#'
#' @param object circSample or circExperiment object
#'
#' @return New class object with replaced filtered backsplice junction read data
#'
#' @examples
#' # A datatable of the filtered backsplice junction read data in a sample object
#' bsj.reads(sample.object) <- new.datatable
#'
#' # A datatable of the filtered backsplice junction read data in each sample of an experiment
#' bsj.reads(experiment.object) <- new.datatable.list
#'
#' @export
setGeneric(name = "bsj.reads<-",
           def = function(object, value)
             standardGeneric("bsj.reads<-"))


#' Replacement method for the backsplice junction count data
#'
#' circSample and circExperiment class replacement methods for replacing the backsplice junction count data
#' The replacment value must be valid to enable replacement
#'
#' @param object circSample or circExperiment object
#'
#' @return New class object with replaced backsplice junction count data
#'
#' @examples
#' # A datatable of the backsplice junction count data in a sample object
#' bsj.counts(sample.object) <- new.datatable
#'
#' # A datatable of the backsplice junction count data in each sample of an experiment
#' bsj.counts(experiment.object) <- new.datatable.list
#'
#' @export
setGeneric(name = "bsj.counts<-",
           def = function(object, value)
             standardGeneric("bsj.counts<-"))


#' Replacement method for summarized data on linear splice sites
#'
#' circSample and circExperiment class replacement methods for replacing the linear splice junction count data
#' The replacment value must be valid to enable replacement
#'
#' @param object circSample or circExperiment object
#'
#' @return New class object with replaced linear splice junction count data
#'
#' @examples
#' # A datatable of the linear splice junction count data in a sample object
#' lsj.counts(sample.object) <- new.datatable
#'
#' # A datatable of the linear splice junction count data in each sample of an experiment
#' lsj.counts(experiment.object) <- new.datatable.list
#'
#' @export
setGeneric(name = "lsj.counts<-",
           def = function(object, value)
             standardGeneric("lsj.counts<-"))


#' Replacement method for sample label
#'
#' circSample and circExperiment class replacement methods for replacing the sample label(s)
#' The replacment value must be valid to enable replacement
#'
#' @param object circSample or circExperiment object
#'
#' @return New class object with replaced sample label(s)
#'
#' @examples
#' # Replacing label of a sample object
#' label(sample.object) <- .new.label
#'
#' # Replacing label in each sample of an experiment
#' label(experiment.object) <- new.label.vector
#'
#' @export
setGeneric(name = "label<-",
           def = function(object, value)
             standardGeneric("label<-"))



#' Replacement method for samples
#'
#' circExperiment class replacement methods for replacing circSample's
#' The replacment value must be valid to enable replacement
#'
#' @param object circExperiment object
#'
#' @return New class object with replaced circSample's
#'
#' @examples
#' # Replacing samples in an experiment
#' samples(experiment.object) <- new.circSample.list
#'
#' @export
setGeneric(name = "samples<-",
           def = function(object, value)
             standardGeneric("samples<-"))


#' Replacement method for experiment name
#'
#' circExperiment class replacement methods for replacing the experiment name
#' The replacment value must be valid to enable replacement
#'
#' @param object circExperiment object
#'
#' @return New class object with replaced experiment name
#'
#' @examples
#' # Replacing the experiment name
#' name(experiment.object) <- new.name
#'
#' @export
setGeneric(name = "name<-",
           def = function(object, value)
             standardGeneric("name<-"))


############################### Populate objects ###############################

#' Method for automatically detecting and setting up all samples in an experiment.
#'
#' @return circExperiment object with populated sample-slot.
#'
#' @examples
#' locateSamples(experiment.object)
#'
#' @importFrom magrittr %>%
#' @importFrom tibble tibble
#' @export
setGeneric(name = "locateSamples",
           def = function(object, organism = "Homo sapiens", gb = "hg38", firstread.firststrand = T, paired.end = T)
             standardGeneric("locateSamples"))


#' Read the chimeric data of a circSample or circExperiment.
#'
#' The function takes a circSample or circExperiment object and reads the chimeric data for samples described in these objects.
#'
#' @param object circSample or circExperiment object
#' @param chromosomes A vector of (valid) chromosome names that is used to restrict the chimeric read data
#' @param cores Only supported on unix! Number of cores used for parallel processing for circExperiment objets.
#' @param maxGenomicDist Maximal allowed size between donor and acceptor (defaults to 100.000 nt).
#' @param onlySpanning Return only BSJ spanning or all (spanning and encompassing) reads.
#'
#' @return Same as input object, but populated slot for bsj data.
#' @export
#'
#' @examples
#' readBSJdata(testExperiment)
setGeneric(name = "readBSJdata",
           def = function(object, chromosomes = c(1:22, "X", "Y"), cores = 1L, maxGenomicDist = 1e5, onlySpanning = TRUE, removeBadPairs = T, ...)
             standardGeneric("readBSJdata"))


#' Read the data on reads spanning linear splice junctions
#'
#' The function takes a circSample or circExperiment object and reads the chimeric data for samples described in these objects (SJ.out.tab file).
#'
#' @param object circSample or circExperiment object
#' @param chromosomes A vector of (valid) chromosome names that is used to restrict the chimeric read data
#' @param cores Only supported on unix! Number of cores used for parallel processing for circExperiment objets.
#'
#' @return Same as input object, but populated slot for lsj counts.
#'
#' @examples
#' readLSJdata(object)
#'
#' @importFrom stringi stri_trans_char
#' @importFrom readr read_tsv
#' @importFrom GenomicRanges GRanges
#' @export
setGeneric(name = "readLSJdata",
           def = function(object, chromosomes = c(1:22, "X", "Y"), cores = 1L, ...)
             standardGeneric("readLSJdata"))


############################### Annotation ###############################

#' Compare backsplice junction identified by STAR to junction sites involved in linear splicing
#'
#' The function takes a circSample or circExperiment object and compares backsplice read data for circSample(s) to a database of known junction sites (generated by \code{\link{constructSJDB,ANY-method}}).
#'
#' @param object circSample or circExperiment object
#' @param known.junctions A vector of (valid) chromosome names that is used to restrict the chimeric read data
#' @param cores Only supported on unix! Number of cores used for parallel processing for circExperiment objets.
#'
#' @examples
#' kj <- constructSJDB(annotationDB = db)
#' readBSJdata(object, known.junctions = kj, cores = 4)
#'
#' @return List of linear splice junction count summaries (SJ.tab.out) for all samples stored in the circExperiment-object.
#' @importFrom dplyr full_join
#' @export
setGeneric(name = "compareToKnownJunctions",
           def = function(object, known.junctions, cores = 1L, ...)
             standardGeneric("compareToKnownJunctions"))


#' Generate junction motifs
#'
#' Generates splice motif sequences and adds them t the backsplice junction datatable
#'
#' @param object A circSample or circExperiment object
#' @param genome_seq A BSGenome object used for determining the sequences for splice motifs
#'
#' @return an updated circSample or circExperiment object
#' @export
#'
#' @importFrom tibble add_column
#'
#' @examples
setGeneric(name = "generateJunctionMotifs",
           def = function(object, genome_seq, cores = 1L)
             standardGeneric("generateJunctionMotifs"))


#' Read alignment to nearest genes
#'
#' Reads circSample bsj.reads, shifts shiftable reads, and replaces bsj.read data with the shifted data
#' @param object circSample or circExperiment object
#' @param cores Only supported on unix! Number of cores used for parallel processing for circExperiment objets.
#' @param ...
#'
#' @return A new object with shifted bsj.read data
#'
#' @examples
#' adjustAlignment(object)
#'
#'
setGeneric(name = "adjustAlignment",
           def = function(object, cores = 1L, ...)
             standardGeneric("adjustAlignment"))


#' Summarize number of reads per backsplice junction.
#'
#' @param object circSample or circExperiment object
#' @param cores Only supported on unix! Number of cores used for parallel processing for circExperiment objets.
#' @param applyFilter Boolean indicating whether to applyFilter (default is TRUE).
#'
#' @return Object with same class as input, but with bsj.count slot populated with data
#'
#' @examples
#' summarizeBSJreads(object)
#'
#' @importFrom dplyr group_by summarise
#' @importFrom GenomicRanges GRanges
#' @importFrom DESeq2 estimateSizeFactorsForMatrix
#' @importFrom plyr join_all
#' @export
setGeneric(name = "summarizeBSJreads",
           def = function(object, cores = 1L, applyFilter = TRUE, ...)
             standardGeneric("summarizeBSJreads"))

#' Add filter to bsj.read-data.
#'
#' @param object circSample or circExperiment object
#' @param filter A logical vector (circSample) or list of logical vectors (circExperiment) indicating which bsj.reads should be included in downstream analyses.
#' @param mode Character, "strict" sets new filter value to FALSE if either the previous \cide{include.read} or new filter value is \code{FALSE}. "last" overwrites values with the new regardles of previous value.
#'
#' @return Object with same class as input.
#'
#' @examples
#' addFilter(object, filter)
#'
#' @export
setGeneric(name = "addFilter",
           def = function(object, filter, mode = "strict", ...)
             standardGeneric("addFilter"))

#' circulaR pipeline wrapper
#'
#' Performs all nescesary steps for filtering, identifying and generate information on backsplice junction candidates
#'
#' @param object circSample or circExperiment object
#' @param annotationDB A TxDb or EnsDb object
#' @param chromosomes Character vector on which chromosomes to keep, while discardin all other chromosomes. The default value (NULL)
#' @param cores Only supported on unix! Number of cores used for parallel processing for circExperiment objets.
#' @param ...
#'
#' @return circSample or circExperiment with all data slots populated with processed data
#'
#' @examples
#' circulaR(circExperiment)
#'
#' @export
setGeneric(name = "circulaR",
           def = function(object, annotationDB, cores = 1L, ...)
             standardGeneric("circulaR"))


#' Construct unique backsplice junction IDs
#'
#' Generate unique readID from information of the candidate backsplice junction read data.
#'
#' @param object circSample or circExperiment object
#' @param cores Only supported on unix! Number of cores used for parallel processing for circExperiment objets.
#'
#' @return A new circSample or circExperiment object, with a bsID column in the bsj.read datatable.
#'
#' @examples
#' constructBsId(object)
#'
#' @export
setGeneric(name = "constructBsId",
           def = function(object, cores = 1L){
             standardGeneric("constructBsId")
           })

################################ Stats #################################

#' Extract alignment stats
#'
#' Parse STAR alignment log files, and extract alignment stats.
#'
#' @param object circSample or circExperiment object
#'
#' @return Datatable with STAR log file information
#'
#' @examples
#' alignmentStats(object)
#'
#' @export
setGeneric(name = "alignmentStats",
           def = function(object, out_type = "wide")
             standardGeneric("alignmentStats"))


#' Make summary statistics of the chimeric/backspliced reads.
#'
#' @param object circSample or circExperiment object
#'
#' @return tibble with summary statistics.
#'
#' @examples
#' bsjStats(object)
#'
#' @importFrom plyr join_all
#' @import tidyr
#' @export
setGeneric(name = "bsjStats",
           def = function(object, out_type = "wide", ...)
             standardGeneric("bsjStats"))


################################ Vizualization #################################

#' Visualize the linear and backsplice junction of a specific gene or genomic region in a specific sample.
#'
#' @param object circSample
#' @param sampleIndex Relevant if we should support circExperiment. Numeric value corresponding to the sample to plot.
#' @param sampleName Relevant if we should support circExperiment. Name of the sample to plot.
#' @param symbol Character vector of length 1 with the symbol (genename) of the gene to plot.
#' @param r GenomicRange corresponding to the region to plot.
#'
#' @return A plot
#'
#' @examples
#' vizJunctions(testSample, symbol = "FIRRE")
#' @importFrom dplyr group_by summarise
#' @importFrom GenomicRanges GRanges
#' @importFrom ensembldb genes getGeneRegionTrackForGviz
#' @importFrom IRanges subsetByOverlaps
#' @export
setGeneric(name = "vizJunctions",
           def = function(object, symbol = NULL, range = NULL, db, xlim = NULL,
                          onlyJunctionsWithinRange = FALSE,
                          validExonModels = NULL, deducedCoverage = NULL,
                          chromosomes = TRUE, path = NULL,
                          prefix = NULL, device = "png", cores = 1L, ...)
             standardGeneric("vizJunctions"))



############################### General functions ###############################
#' Check for empty values
#'
#' Check if the object has an empty value, does not include checking empty strings.
#' @param object The object that are to be checked.
#'
#' @return A logical value that states whether the object is considered empty.
#' @export
#'
#' @examples
#' # Test that turns out TRUE
#' a <- NULL
#' b <- character()
#' is.empty(a)
#' is.empty(b)
#' is.empty(list())  # is TRUE as length = 0
#' is.empty(c())
#'
#' # Test that turns out FALSE
#' is.empty("")
#' is.empy(NA)
setGeneric(name = "is.empty",
           def = function(object)
             standardGeneric("is.empty"))


#' Construct a database of known splice junctions from a database of annotated genes.
#'
#' This function takes an annotation database, either downloaded from AnnotationHub
#' or generated using makeTxDbFromBiomart (GenomicFeatures) and extracts all junction sites and transcript termini.
#'
#' @param annotationDB A TxDb or EnsDb object
#' @param force Boolean (default = FALSE) indicating whether or not to overwrite any cached version splice junction database.
#'
#' @export
#' @importFrom GenomicRanges GRanges
#' @importFrom RSQLite dbReadTable
#' @importFrom magrittr %>%
#' @importFrom dplyr group_by summarise
setGeneric(name = "constructSJDB",
           def = function(annotationDB, force = FALSE)
             standardGeneric("constructSJDB"))


#' Generate database of known splice junction sites from a TxDb- or EnsDb-database
#'
#' @param annotationDB TxDb or EnsDb object
#'
#' @importFrom GenomicFeatures transcripts
#' @importFrom GenomicRanges GRanges
#' @importFrom ensembldb transcripts
#' @importFrom AnnotationDbi select
#' @importFrom GenomicRanges GRanges
#' @export
setGeneric(name = "getKnownJunctions",
           def = function(annotationDB)
             standardGeneric("getKnownJunctions"))

#' Determine coverage of backspliced reads.
#'
#' @param object circSample
#' @param bsid Backsplice ID
#'
#' @importFrom GenomicFeatures transcripts
#' @importFrom GenomicRanges GRanges
#' @importFrom ensembldb transcripts
#' @importFrom AnnotationDbi select
#' @importFrom GenomicRanges GRanges
#' @export
setGeneric(name = "bsReadCoverage",
           def = function(object, bsid)
             standardGeneric("bsReadCoverage"))

#' Accessor for expression data from circExperiment
#'
#' circExperiment class accessors for expression data from a circExperiment
#'
#' @param object circExperiment object
#' @param format Should result be returned as a "long" (fast) or "wide" (slow) table
#'
#' @return A tibble with expression values
#'
#' @examples
#' # The samples of an experiment
#' exprs(circExperiment)
#'
#' @importFrom magrittr %>%
#' @importFrom GenomicRanges mcols
#' @importFrom dplyr bind_rows
#' @importFrom tidyr spread
#'
#' @export
setGeneric(name = "exp.mat",
           def = function(object, format = "long", normFactors = NULL, ...)
             standardGeneric("exp.mat"))

#' Function to calculate the circular to linear ratio for a given bsID
#'
#' This function identifies linear splice sites that use donor or acceptor from a given bsID
#' and calculate the circ:lin ratio
#'
#' @param object circSample
#' @param bsid Backsplice ID
#' @param onlyAnnotated Boolean value indicating whether to use only known linear splice junctions
#'
#' @importFrom IRanges subsetByOverlaps resize

#' @export
setGeneric(
  name = "c2l.ratio",
  def = function(object, bsid, onlyAnnotated = TRUE, ...)
    standardGeneric("c2l.ratio")
)

