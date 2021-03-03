#' Plot transcripts stored in a TxDB.
#'
#' Plot the transcript models of a gene.
#'
#' @param ex GrangesList with exonsBy tx.
#' @param int GrangesList with intronsByTranscript information.
#' @param addChevrons Boolean indicating wheter to add chevrons to intronic regions to indicate direction of the transcript.
#' @param rangex The positions of the leftmost and rightmost nucleotides to plot.
#' @return Don't know?
#' @examples
#' ensembl <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
#' annot <- getBM(
#'   attributes = c("ensembl_gene_id", "external_gene_name", "ensembl_transcript_id", "external_transcript_name", "transcript_tsl", "chromosome_name", "start_position", "end_position", "strand"),
#'   filters = "external_gene_name",
#'   values = "FIRRE",
#'   mart = ensembl
#' )
#' txdb <- makeTxDbFromBiomart(transcript_ids = annot$ensembl_transcript_id)#, host = "uswest.ensembl.org")
#' plotTranscripts(tdb = txdb)
# plotTranscripts <- function(ex = NULL, int = NULL, addChevrons = T, rangex = NULL) {
#   if(!grepl("GRangesList", class(ex)) & !grepl("GRangesList", class(int))){stop("Expecting GRangesLists object describing ex and int.")}
#   if(!length(ex) == length(int)){stop("Number ex and int don't seem to match.")}
#
#   if(is.null(rangex)){stop("Need to explicitly define range of genomic coordinates.")}
#
#   ntranscripts <- length(ex)
#   nex <- lapply(ex, length)
#
#   # Make top viewport. A grid with one column and as many rows as transcripts to plot.
#   top.vp <- grid::viewport(layout = grid::grid.layout(ntranscripts, 1))
#   grid::pushViewport(top.vp)
#
#   # Make viewports for each intron in each transcript
#   for(t in 1:ntranscripts){
#     # Make viewport for transcript
#     vp <- grid::viewport(layout.pos.col = 1, layout.pos.row = t, name = names(ex)[t], xscale = rangex)
#     grid::pushViewport(vp)
#
#     int.widths <- BiocGenerics::width(int[[t]])
#     for(i in 1:length(int[[t]])){
#       if(length(int.widths) == 0){next()}
#       if(int.widths[i] < 5)next()
#       # Make viewport for intron
#       x.boundaries <- c(BiocGenerics::start(int[[t]][i])-0.5, BiocGenerics::end(int[[t]][i])+0.5)
#
#       vp <- viewport(
#         x = grid::unit(x.boundaries[1], "native"),
#         y = 0.4,
#         width = grid::unit(diff(x.boundaries), "native"),
#         height = 0.2,
#         just = c("left", "bottom"),
#         name = paste0(names(ex)[t], "_intron_",i)
#       )
#
#       grid::pushViewport(vp)
#       #grid::grid.rect(gp=grid::gpar(col=colors()[sample(x = 2:length(colors()), size = 1)]))
#       grid::grid.lines(
#         x = grid::unit(c(0,1), "npc"),
#         y = grid::unit(0.5, "npc")
#       )
#
#       if(addChevrons){# Plot chevrons
#         # Determine number of chevrons to plot
#         size.of.intron.viewport <- grid::convertunit(grid::unit(x.boundaries[2], "native") - grid::unit(x.boundaries[1], "native"), unitTo = "cm")
#         size.of.plot.area <- grid::convertunit(grid::unit(rangex[2], "native") - grid::unit(rangex[1], "native"), unitTo = "cm")
#         rel.size <- as.numeric(size.of.intron.viewport)/as.numeric(size.of.plot.area)
#         nc <- ceiling(rel.size / 0.04)
#
#         if(rel.size > 0.05){ # Exclude small int
#           for(che in 1:(nc-1)){
#             rel.x <- che/nc
#             grid::grid.lines(
#               x = grid::unit(c(rel.x, rel.x+0.001), "npc"),
#               y = grid::unit(c(0.5), "npc"),
#               arrow = grid::arrow(angle = 30, length = grid::unit(0.1, "inches"), ends = ifelse(all(GenomicRanges::strand(int[[t]]) == "+"), "last", "first"), type = "open")
#             )
#           }
#         }
#       }
#       grid::upViewport()
#     }
#
#     # Plot the ex
#     grid::grid.rect(y = rep(0.5, nex[[t]]),
#               x = grid::unit(BiocGenerics::start(ex[[t]])-0.5, "native"),
#               width = grid::unit(apply(cbind(BiocGenerics::start(ex[[t]]), BiocGenerics::end(ex[[t]])),1,diff)+1, "native"),
#               height = 0.6, just = "left", gp = grid::gpar(fill = "steelblue")
#     )
#     grid::upViewport()
#   }
#   grid::popViewport()
# }

#' Plot transcript models.
#'
#' Using data stored in a TxDb to visualize the transcript models.
#'
#' @param db A TxDb containing transcript models from Ensembl.
#' @param g If the TxDb contains more than one gene, then use g to specify what to plot.
#' @param addChevrons Boolean indicating wheter to add chevrons to intronic regions to indicate direction of the transcript.
#' @param rangex The range in genomic coordinates for the region of interest. This is important for correct alignment of the tx model plot with splice junction plot.
#' @export
#' @importFrom AnnotationDbi select
#' @importFrom magrittr %>%
#' @importFrom grid arrow convertUnit gpar grid.layout grid.lines grid.rect popViewport pushViewport unit upViewport viewport
plotTxModelsOLD <- function(db = NULL, g = NULL, tid = NULL, addChevrons = T, rangex = NULL) { # Add option to define transcripts to plot, tx = NULL
  if(is.null(db)){stop("Need a database of transcripts.")}
  if(!class(db) %in% c("TxDb")){stop("Database of transcripts should be made using GenomicFeatures.")}
  if(is.null(g)){stop("Please indicate which gene to plot.")}
  if(length(g) > 1){stop("Can only plot a single gene.")}
  g.db <- suppressMessages(AnnotationDbi::select(db, g, c("TXNAME", "TXSTART", "TXEND"), keytype = "GENEID"))
  if(nrow(g.db) == 0){stop("Cannot find gene in database.")}
  if(is.null(rangex)){stop("Need to explicitly define range of genomic coordinates.")}

  # Get transcript names from db
  tx <- suppressMessages(AnnotationDbi::select(db, keys(db), "TXNAME", keytype = "GENEID")) %>% subset(GENEID %in% g) %>% .$TXNAME
  if(!is.null(tid) & !all(tid %in% tx)){
    stop("Some of the supplied Ensembl transcript IDs are not found in the local database.")
  }

  # exons and intron coordinates
  if(!is.null(tid)){
    tx <- tid
  }
  ex <- exonsBy(db, by = "tx", use.names = T)[tx]
  int <- intronsByTranscript(db, use.names = T)
  int <- int[names(int) %in% tx]

  # Check if rangex is much smaller than gene region. If so, turn off chevrons
  if(diff(rangex) / diff(c(min(g.db$TXSTART), max(g.db$TXEND))) < 0.1){
    message("Tight zoom, changing addChevrons to FALSE")
    addChevrons <- FALSE
  }

  ntranscripts <- length(ex)
  nex <- lapply(ex, length)

  # Make top viewport. A grid with one column and as many rows as transcripts to plot.
  top.vp <- grid::viewport(layout = grid::grid.layout(ntranscripts, 1))
  grid::pushViewport(top.vp)

  # Make viewports for each intron in each transcript
  for(t in 1:ntranscripts){
    # Make viewport for transcript
    vp <- grid::viewport(layout.pos.col = 1, layout.pos.row = t, name = names(ex)[t], xscale = rangex)
    grid::pushViewport(vp)

    if(names(ex)[t] %in% names(int)){ # If there is an intron in the transcript, then plot it!
      int.widths <- BiocGenerics::width(int[[t]])
      for(i in 1:length(int[[t]])){
        if(length(int.widths) == 0){next()}
        if(int.widths[i] < 5)next()
        # Make viewport for intron
        x.boundaries <- c(BiocGenerics::start(int[[t]][i])-0.5, BiocGenerics::end(int[[t]][i])+0.5)

        vp <- viewport(
          x = grid::unit(x.boundaries[1], "native"),
          y = 0.4,
          width = grid::unit(diff(x.boundaries), "native"),
          height = 0.2,
          just = c("left", "bottom"),
          name = paste0(names(ex)[t], "_intron_",i)
        )

        grid::pushViewport(vp)
        #grid::grid.rect(gp=grid::gpar(col=colors()[sample(x = 2:length(colors()), size = 1)]))
        grid::grid.lines(
          x = grid::unit(c(0,1), "npc"),
          y = grid::unit(0.5, "npc")
        )

        if(addChevrons){# Plot chevrons
          # Determine number of chevrons to plot
          size.of.intron.viewport <- grid::convertunit(grid::unit(x.boundaries[2], "native") - grid::unit(x.boundaries[1], "native"), unitTo = "cm")
          size.of.plot.area <- grid::convertunit(grid::unit(rangex[2], "native") - grid::unit(rangex[1], "native"), unitTo = "cm")
          rel.size <- as.numeric(size.of.intron.viewport)/as.numeric(size.of.plot.area)
          nc <- ceiling(rel.size / 0.04)

          if(rel.size > 0.05){ # Exclude small int
            for(che in 1:(nc-1)){
              rel.x <- che/nc
              grid::grid.lines(
                x = grid::unit(c(rel.x, rel.x+0.001), "npc"),
                y = grid::unit(c(0.5), "npc"),
                arrow = grid::arrow(angle = 30, length = grid::unit(0.1, "inches"), ends = ifelse(all(GenomicRanges::strand(int[[t]]) == "+"), "last", "first"), type = "open")
              )
            }
          }
        }
        grid::upViewport()
      }
    }

    # Plot the ex
    grid::grid.rect(y = rep(0.5, nex[[t]]),
              x = grid::unit(BiocGenerics::start(ex[[t]])-0.5, "native"),
              width = grid::unit(apply(cbind(BiocGenerics::start(ex[[t]]), BiocGenerics::end(ex[[t]])),1,diff)+1, "native"),
              height = 0.6, just = "left", gp = grid::gpar(fill = "steelblue")
    )
    grid::upViewport()
  }
  grid::popViewport()
}

#' Plot transcript models.
#'
#' Using data stored in an EnsDB from AnnotationHub.
#'
#' @param db A db containing transcript models from AnnotationHub/Ensembl.
#' @param g Specify which gene to plot (Ensembl Gene ID).
#' @param showDirection Boolean indicating wheter to add chevrons to intronic regions to indicate direction of the transcript.
#' @param rangex The range in genomic coordinates for the region of interest. This is important for correct alignment of the tx model plot with splice junction plot.
#' @export
#' @importFrom AnnotationDbi select
#' @importFrom magrittr %>%
#' @importFrom grid arrow convertUnit gpar grid.layout grid.lines grid.rect popViewport pushViewport unit upViewport viewport
plotTxModels <- function(db = NULL, g = NULL, tid = NULL, showDirection = T, rangex = NULL) { # Add option to define transcripts to plot, tx = NULL
  if(is.null(db)){stop("Need a database of transcripts.")}
  if(!class(db) %in% c("EnsDb")){stop("Database of transcripts should be EnsDb from AnnotationHub.")}
  if(is.null(g)){stop("Please indicate which gene to plot.")}
  if(length(g) > 1){stop("Can only plot a single gene.")}
  g.db <- suppressMessages(AnnotationDbi::select(db, g, c("TXNAME", "TXSEQSTART", "TXSEQEND"), keytype = "GENEID"))
  if(nrow(g.db) == 0){stop("Cannot find gene in database.")}
  if(is.null(rangex)){stop("Need to explicitly define range of genomic coordinates.")}

  # Get transcript names from db
  tx <- suppressMessages(AnnotationDbi::select(db, keys(db), "TXNAME", keytype = "GENEID")) %>% subset(GENEID %in% g) %>% .$TXNAME
  if(!is.null(tid) & !all(tid %in% tx)){
    stop("Some of the supplied Ensembl transcript IDs are not found in the local database.")
  }

  # exons and intron coordinates
  if(!is.null(tid)){
    tx <- tid
  }
  ex <- AnnotationDbi::select(db, g, c("TXNAME", "TXSUPPORTLEVEL", "EXONID", "EXONIDX", "EXONSEQSTART", "EXONSEQEND", "SEQSTRAND"), keytype = "GENEID")
  ex <- subset(ex, TXNAME %in% tx)

  # Check if rangex is much smaller than gene region. If so, turn off chevrons
  if(diff(rangex) / diff(c(min(ex$EXONSEQSTART), max(ex$EXONSEQEND))) < 0.1){
    message("Tight zoom, changing addChevrons to FALSE")
    showDirection <- FALSE
  }

  ntranscripts <- ex$TXNAME %>% unique %>% length
  nex <- ex %>% dplyr::group_by(TXNAME) %>% dplyr::summarise(no.ex = length(EXONID))

  # Make top viewport. A grid with one column and as many rows as transcripts to plot.
  top.vp <- grid::viewport(layout = grid::grid.layout(ntranscripts, 1))
  grid::pushViewport(top.vp)

  # Make viewports for each intron in each transcript
  for(t in 1:ntranscripts){
    tx <- unique(ex$TXNAME)[t]
    # Make viewport for transcript
    vp <- grid::viewport(layout.pos.col = 1, layout.pos.row = t, name = tx, xscale = rangex)
    grid::pushViewport(vp)
    tmp.tx <- subset(ex, TXNAME == tx)

    if(nrow(tmp.tx) >= 2){ # If there are 2 or more exons, plot the intron!
      # Derive introns from the tx information
      tmp.int <- data.frame(
        INTRONSEQSTART = tmp.tx$EXONSEQEND[-nrow(tmp.tx)] + 1,
        INTRONSEQEND = tmp.tx$EXONSEQSTART[-1] - 1
      )
      # Go through all introns and plot
      for(i in 1:nrow(tmp.int)){
        # If intron is too narrow, then skip
        if(abs(diff(as.numeric(tmp.int[i,]))) < 5) next()

        # Make viewport for intron
        x.boundaries <- as.double(tmp.int[i,]) + c(-0.5, 0.5)

        vp <- viewport(
          x = grid::unit(x.boundaries[1], "native"),
          y = 0.4,
          width = grid::unit(diff(x.boundaries), "native"),
          height = 0.2,
          just = c("left", "bottom"),
          name = paste0(tx, "_intron_",i)
        )

        grid::pushViewport(vp)
        #grid::grid.rect(gp=grid::gpar(col=colors()[sample(x = 2:length(colors()), size = 1)]))
        grid::grid.lines(
          x = grid::unit(c(0,1), "npc"),
          y = grid::unit(0.5, "npc")
        )

        if(showDirection){# Plot chevrons
          # Determine number of chevrons to plot
          size.of.intron.viewport <- grid::convertunit(grid::unit(x.boundaries[2], "native") - grid::unit(x.boundaries[1], "native"), unitTo = "cm")
          size.of.plot.area <- grid::convertunit(grid::unit(rangex[2], "native") - grid::unit(rangex[1], "native"), unitTo = "cm")
          rel.size <- as.numeric(size.of.intron.viewport)/as.numeric(size.of.plot.area)
          nc <- ceiling(rel.size / 0.04)

          if(rel.size > 0.05){ # Exclude small int
            for(che in 1:(nc-1)){
              rel.x <- che/nc
              grid::grid.lines(
                x = grid::unit(c(rel.x, rel.x+0.001), "npc"),
                y = grid::unit(c(0.5), "npc"),
                arrow = grid::arrow(angle = 30, length = grid::unit(0.1, "inches"), ends = ifelse(all(ex$SEQSTRAND == 1), "last", "first"), type = "open")
              )
            }
          }
        }
        grid::upViewport()
      }
    }

    # Plot the ex
    grid::grid.rect(y = rep(0.5, nrow(tmp.tx)),
              x = grid::unit(tmp.tx$EXONSEQSTART-0.5, "native"),
              width = grid::unit(apply(tmp.tx[,c("EXONSEQSTART", "EXONSEQEND")],1,diff)+1, "native"),
              height = 0.6, just = "left", gp = grid::gpar(fill = "steelblue")
    )
    grid::upViewport()
  }
  grid::popViewport()
}



#' Plot linear splice sites stored in a GRanges object.
#'
#' @param df yada yada yada
#' @param trimJ Boolean, should junctions with donor/acceptor outside the range defined in 'rangex' be excluded.
#' @param line.scaling Factor for overall scaling of line widths.
#' @param rangex see plotTranscripts
#'
#' @importFrom BiocGenerics start end
#' @importFrom grid gpar grid.bezier grid.rect grid.text popViewport pushViewport unit viewport
#' @export
plotSJ <- function(df = NULL, rangex, line.scaling = 0.5, trimJ = T) {
  if(is.null(df)){stop("Please input some data.")}
  if(class(df) != "GRanges"){stop("Please provide a GRanges object!")}
  if(trimJ){
    df <- df[!(BiocGenerics::start(df) < min(rangex) | BiocGenerics::end(df) > max(rangex)),]
  }
  # Remove row with no counts in unique mapping
  df <- df[df$read.count.unique != 0,]

  if(length(df) == 0){
    grid::pushViewport(viewport())
    grid::grid.rect(gp = grid::gpar(col = "grey"))
    grid::grid.text("No linear splice junctions in this region")
  } else {
    # Scale height to length of junction. Height should be in the range [0.1 , 0.95].
    # The min value should be determined by the number of junctions to plot
    scale.range <- c(ifelse(length(df)<10, 1-length(df)*0.1, .1), 0.95)
    df$height <- BiocGenerics::width(df)
    df$height <- (df$height-min(df$height)) * diff(scale.range) / (max(df$height)-min(df$height)) + min(scale.range)

    tmp <- GenomicRanges::strand(df)
    tmp <- tmp[tmp!="*"]

    if(all(tmp == "-")){
      df$height <- 1-df$height
      y.start = 0.95
    } else {
      y.start = 0.05
    }

    # Scale line width to number of reads. Initially, I'll hardcode absolute values
    # 1-10:   1
    # 11-25:  3
    # 26-75:  5
    # 76-150: 7
    # 151-400:9
    # 401+:   11
    lw <- cut(df$read.count.unique, breaks = c(1,11,26,71,151,401, Inf), include.lowest = T)
    levels(lw) <- seq(from = 1, to = 11, by = 2)
    df$w <- as.integer(as.character(lw))

    grid::pushViewport(grid::viewport(xscale = rangex))

    for (i in 1:length(df)) {
      x <- grid::unit(rep(c(BiocGenerics::start(df[i]), BiocGenerics::end(df[i])), each=2), "native")
      y <- grid::unit(c(y.start, df$height[i], df$height[i], y.start), "npc")

      grid::grid.bezier(x, y, gp=grid::gpar(lwd = df$w[i], lex = line.scaling, col = "blue"))
    }
    grid::popViewport()
  }
}

#' Plot backsplice junctions stored in a GRanges object.
#'
#' @param df yada yada yada
#' @param line.scaling Factor for overall scaling of line widths.
#' @param rangex see plotTranscripts
#' @importFrom grid gpar grid.rect grid.text grid.xspline popViewport pushViewport unit viewport
#' @export
plotBSJ <- function(df = NULL, rangex, line.scaling = 0.5) {
  if(is.null(df)){stop("Please input some data.")}
  if(class(df) != "GRanges"){stop("Please provide a GRanges object!")}

  if(length(df) == 0){
    grid::pushViewport(viewport())
    grid::grid.rect(gp = grid::gpar(col = "grey"))
    grid::grid.text("No backsplice junctions in this region")
  } else {
    # Scale height to length of junction. Height should be in the range [0.15 , 0.85].
    if(length(df) > 1){
      scale.range <- c(ifelse(length(df) < 8, 0.85 - length(df)*0.1, 0.15), 0.85)
      df$height <- BiocGenerics::width(df)
      df$height <- (df$height-min(df$height)) * diff(scale.range) / (max(df$height)-min(df$height)) + min(scale.range)
    } else {
      df$height = 0.85
    }

    if(all(GenomicRanges::strand(df) == "+")){
      df$height <- 1-df$height
      y.start = 0.95
    } else {
      y.start = 0.05
    }

    # Scale line width to number of reads. Initially, I'll hardcode absolute values
    # 1-10:   1
    # 11-25:  3
    # 26-75:  5
    # 76-150: 7
    # 151-400:9
    # 401+:   11
    lw <- cut(df$count, breaks = c(1,11,26,71,151,401, Inf), include.lowest = T)
    levels(lw) <- seq(from = 1, to = 11, by = 2)
    df$w <- as.integer(as.character(lw))

    grid::pushViewport(grid::viewport(xscale = rangex))

    for (i in 1:length(df)) {
      x <- grid::unit(rep(c(BiocGenerics::start(df[i]), BiocGenerics::end(df[i])), each=2), "native")
      y <- grid::unit(c(y.start, df$height[i], df$height[i], y.start), "npc")

      grid::grid.xspline(x, y, shape=c(0, 1, 1, 0), open=TRUE, gp=grid::gpar(lwd = df$w[i], lex = line.scaling, col = "red"))
      #grid::grid.bezier(x, y, gp=grid::gpar(lwd = df$w[i], lex = line.scaling))
    }
    grid::popViewport()
  }
}

#' Zoomed in view of the donor and acceptor sites
#'
#' @param df Data.frame of candidate backspliced reads
#' @param rid Name of read to focus on
#' @param size Number of nucleotides adjacent to the junctions
#' @param l.fontsize Sice of the font used for legend
#' @importFrom grid gpar grid.clip grid.layout grid.newpage grid.points grid.rect grid.text pushViewport unit upViewport viewport
#' @importFrom BiocGenerics sort
#' @importFrom GenomicRanges GRanges
#' @export
plotBSJrepeats <- function(df = candidate.bsj, rid = NULL, size = 10, l.fontsize = 9, unstranded = F) {
  r <- subset(df, X10 == rid)
  xlim <- BiocGenerics::sort(as.integer(r[,c("X2", "X5")]))
  ex <- subsetByOverlaps(EXONS, GenomicRanges::GRanges(seqnames = r$X1, ranges = IRanges(start = xlim[1], end = xlim[2]), strand = ifelse(unstranded, "*", r$X3)))
  int <- INTRONS[names(INTRONS) %in% names(ex)]

  r.rep <- r$X8
  l.rep <- r$X9

  if(r$X3 == "+"){
    if(r.rep == 0 & l.rep == 0){
      l.extra <- "no adjacent repeats"
      l.pos <- NULL
      r.pos <- NULL
    }
    if(r.rep != 0 & l.rep == 0){
      l.extra <- paste0(r.rep, " nt repeat to the right")
      l.pos <- NULL
      r.pos <- (size+1):(size+r.rep) # rel to acceptor, add 1 for donor coordinates
    }
    if(r.rep == 0 & l.rep != 0){
      l.extra <- paste0(l.rep, " nt repeat to the left")
      l.pos <- (size-l.rep+1):(size) # rel to acceptor, add 1 for donor coordinates
      r.pos <- NULL
    }
    if(r.rep != 0 & l.rep != 0){
      l.extra <- paste0(r.rep, " nt repeat to the right & ", l.rep, " nt repeat to the left")
      l.pos <- (size-l.rep+1):(size) # rel to acceptor, add 1 for donor coordinates
      r.pos <- (size+1):(size+r.rep) # rel to acceptor, add 1 for donor coordinates
    }
  }

  if(r$X3 == "-"){
    if(r.rep == 0 & l.rep == 0){
      l.extra <- "no adjacent repeats"
      l.pos <- NULL
      r.pos <- NULL
    }
    if(r.rep != 0 & l.rep == 0){
      l.extra <- paste0(r.rep, " nt repeat to the right")
      l.pos <- NULL
      r.pos <- (size-r.rep+2):(size+1) # rel to acceptor, subtract one for donor coordinates
    }
    if(r.rep == 0 & l.rep != 0){
      l.extra <- paste0(l.rep, " nt repeat to the left")
      l.pos <- (size+2):(size+l.rep+1) # rel to acceptor, subtract one for donor coordinates
      r.pos <- NULL
    }
    if(r.rep != 0 & l.rep != 0){
      l.extra <- paste0(r.rep, " nt repeat to the right & ", l.rep, " nt repeat to the left")
      l.pos <- (size+2):(size+l.rep+1) # rel to acceptor, subtract one for donor coordinates
      r.pos <- (size-r.rep+2):(size+1) # rel to acceptor, subtract one for donor coordinates
    }
  }



  # Focus on acceptor
  ac.region <- (r$X2-size):(r$X2+size)
  do.region <- (r$X5-size):(r$X5+size)

  # Setup plot area
  vp.ac <- grid::viewport(x = 0, y = 1, w = 0.49, h = 1, just = c("left", "top"), name = "acceptor")
  vp.do <- grid::viewport(x = 0.51, y = 1, w = 0.49, h = 1, just = c("left", "top"), name = "donor")

  grid::grid.newpage()
  ### Acceptor
  grid::pushViewport(vp.ac)
  grid::grid.rect()

  vp.junc <- grid::viewport(x = 0, y = 1, w = 1, h = 0.4, just = c("left", "top"), name = "juncs", xscale = range(ac.region))
  vp.tx   <- grid::viewport(x = 0, y = 0.6, w = 1, h = 0.4, just = c("left", "top"), name = "tx", xscale = range(ac.region))
  vp.seq <- grid::viewport(x = 0, y = 0.2, w = 1, h = 0.2, just = c("left", "top"), name = "seq", xscale = range(ac.region))

  # Plot DNA strands
  s <- Hsapiens[[r$X1]][ac.region]
  sc <- complement(s)

  grid::pushViewport(vp.seq)
  vp.s <- grid::viewport(layout = grid::grid.layout(nrow = 2, ncol = length(s)))
  grid::pushViewport(vp.s)
  for(i in 1:length(s)){
    grid::pushViewport(grid::viewport(layout.pos.row = 1, layout.pos.col = i))
    if(i %in% l.pos & r$X3 == "+"){
      grid::grid.rect(height = 0.8, gp = grid::gpar(col = "red", fill = "red"))
    }
    if(i %in% r.pos & r$X3 == "+"){
      grid::grid.rect(height = 0.8, gp = grid::gpar(col = "green", fill = "green"))
    }
    grid::grid.text(label = as.character(s[i]), gp = grid::gpar(col = ifelse(r$X3 == "+", "black", "grey85")))
    grid::upViewport()
    grid::pushViewport(grid::viewport(layout.pos.row = 2, layout.pos.col = i))
    if(i %in% l.pos & r$X3 == "-"){
      grid::grid.rect(height = 0.8, gp = grid::gpar(col = "red", fill = "red"))
    }
    if(i %in% r.pos & r$X3 == "-"){
      grid::grid.rect(height = 0.8, gp = grid::gpar(col = "green", fill = "green"))
    }
    grid::grid.text(label = as.character(sc[i]), gp = grid::gpar(col = ifelse(r$X3 == "+", "grey85", "black")))
    grid::upViewport()
  }
  grid::upViewport()
  grid::upViewport()

  grid::pushViewport(vp.tx)
  grid::grid.rect(gp = grid::gpar(col = "grey"))
  grid::grid.clip()
  plotTranscripts(ex = ex, int = int, rangex = range(ac.region), addChevrons = F)
  # for(i in 1:length(int)){ # Add acceptor and donor sites
  #   grid::pushViewport(grid::viewport(layout.pos.col=1, layout.pos.row=i, xscale = range(ac.region)))
  #   grid::grid.points(
  #     x = tmp<-grid::unit(c(BiocGenerics::start(int[[i]]), BiocGenerics::end(int[[i]])), units = "native"),
  #     y = grid::unit(rep(0.5, length(tmp)), "npc"),
  #     pch="x")
  #   grid::upViewport()
  # }
  grid::upViewport()

  grid::pushViewport(vp.junc)
  grid::grid.clip()
  grid::grid.points(x = grid::unit(r$X2, "native"), y = grid::unit(0.3, "npc"), pch=16)
  #grid::upViewport()

  # Plot fragments
  et <- r$X12
  to <- r$X14

  grid::grid.rect(
    x = grid::unit(r$X11-0.5, "native"),
    y = grid::unit(0.1, "npc"),
    width = grid::unit(parseCIGAR(r$X12, returnLength = T), "native"),
    height = grid::unit(0.1, "npc"), just = "left", gp = grid::gpar(fill = "chartreuse")
  )

  # Plot second fragment
  grid::grid.rect(
    x = grid::unit(r$X13-0.5, "native"),
    y = grid::unit(0.3, "npc"),
    width = grid::unit(parseCIGAR(r$X14, returnLength = T), "native"),
    height = grid::unit(0.1, "npc"), just = "left", gp = grid::gpar(fill = "steelblue")
  )


  #grid::pushViewport(vp.junc)
  l <- paste0("Backsplice acceptor\n", l.extra, "\nreadID: ",rid)
  grid::grid.text(x=0.5, y = 0.8, label = l, gp = grid::gpar(fontsize = l.fontsize))
  grid::upViewport()

  grid::upViewport()

  ### Donor
  grid::pushViewport(vp.do)
  grid::grid.rect()

  vp.junc <- grid::viewport(x = 0, y = 1, w = 1, h = 0.4, just = c("left", "top"), name = "juncs", xscale = range(do.region))
  vp.tx   <- grid::viewport(x = 0, y = 0.6, w = 1, h = 0.4, just = c("left", "top"), name = "tx", xscale = range(do.region))
  vp.seq <- grid::viewport(x = 0, y = 0.2, w = 1, h = 0.2, just = c("left", "top"), name = "seq", xscale = range(do.region))

  # Plot DNA strands
  s <- Hsapiens[[r$X1]][do.region]
  sc <- complement(s)

  grid::pushViewport(vp.seq)
  vp.s <- grid::viewport(layout = grid::grid.layout(nrow = 2, ncol = length(s)))
  grid::pushViewport(vp.s)
  for(i in 1:length(s)){
    grid::pushViewport(grid::viewport(layout.pos.row = 1, layout.pos.col = i))
    if(i %in% (l.pos+ifelse(r$X3 == "+", 1, -1)) & r$X3 == "+"){
      grid::grid.rect(height = 0.8, gp = grid::gpar(col = "red", fill = "red"))
    }
    if(i %in% (r.pos+ifelse(r$X3 == "+", 1, -1)) & r$X3 == "+"){
      grid::grid.rect(height = 0.8, gp = grid::gpar(col = "green", fill = "green"))
    }
    grid::grid.text(label = as.character(s[i]), gp = grid::gpar(col = ifelse(r$X3 == "+", "black", "grey85")))
    grid::upViewport()
    grid::pushViewport(grid::viewport(layout.pos.row = 2, layout.pos.col = i))
    if(i %in% (l.pos+ifelse(r$X3 == "+", 1, -1)) & r$X3 == "-"){
      grid::grid.rect(height = 0.8, gp = grid::gpar(col = "red", fill = "red"))
    }
    if(i %in% (r.pos+ifelse(r$X3 == "+", 1, -1)) & r$X3 == "-"){
      grid::grid.rect(height = 0.8, gp = grid::gpar(col = "green", fill = "green"))
    }
    grid::grid.text(label = as.character(sc[i]), gp = grid::gpar(col = ifelse(r$X3 == "+", "grey85", "black")))
    grid::upViewport()
  }
  grid::upViewport()
  grid::upViewport()

  grid::pushViewport(vp.tx)
  grid::grid.rect(gp = grid::gpar(col = "grey"))
  grid::grid.clip()
  plotTranscripts(ex = ex, int = int, rangex = range(do.region), addChevrons = F)
  grid::upViewport()

  grid::pushViewport(vp.junc)
  grid::grid.clip()
  grid::grid.points(x = grid::unit(r$X5, "native"), y = grid::unit(0.1, "npc"), pch=1)
  #grid::upViewport()

  # Plot fragments
  et <- r$X12
  to <- r$X14

  grid::grid.rect(
    x = grid::unit(r$X11-0.5, "native"),
    y = grid::unit(0.1, "npc"),
    width = grid::unit(parseCIGAR(r$X12, returnLength = T), "native"),
    height = grid::unit(0.1, "npc"), just = "left", gp = grid::gpar(fill = "chartreuse")
  )

  # Plot second fragment
  grid::grid.rect(
    x = grid::unit(r$X13-0.5, "native"),
    y = grid::unit(0.3, "npc"),
    width = grid::unit(parseCIGAR(r$X14, returnLength = T), "native"),
    height = grid::unit(0.1, "npc"), just = "left", gp = grid::gpar(fill = "chartreuse")
  )

  #grid::pushViewport(vp.junc)
  l <- paste0("Backsplice donor\n", l.extra, "\nreadID: ",rid)
  grid::grid.text(x=0.5, y = 0.8, label = l, gp = grid::gpar(fontsize = l.fontsize))
  grid::upViewport()

  grid::upViewport()
}

#' Make homemade sashimi plot
#'
#' Plot transcript models and show linear and backsplice jjunctions.
#'
#' @param DB TxDb of transcripts
#' @param sj Linear splice junction data (output from readSJ)
#' @param bsj Backsplice junction data (output from readBSJ)
#' @param region GRanges object describing the region of interest
#' @param gid Ensembl gene id to plot grid.newpage
#' @importFrom grid downViewport pushViewport upViewport viewport
#' @export
plotAll <- function(DB = txdb, sj = sj.data, bsj = bsj.data, region = region.of.interest, gid = g, ...){
  st <- as.character(unique(GenomicRanges::strand(bsj)))

  vp.upper <- grid::viewport(x = 0, y = 1, w = 1, h = 1/3, just = c("left", "top"), name = "upper")
  vp.tx <- grid::viewport(x = 0, y = 2/3, w = 1, h = 1/3, just = c("left", "top"), name = "tx")
  vp.lower <- grid::viewport(x = 0, y = 1/3, w = 1, h = 1/3, just = c("left", "top"), name = "lower")

  # Initialize viewports
  grid::grid.newpage()
  grid::pushViewport(vp.upper)
  grid::upViewport()
  grid::pushViewport(vp.tx)
  grid::upViewport()
  grid::pushViewport(vp.lower)
  grid::upViewport()


  xlim <- c(BiocGenerics::start(region), BiocGenerics::end(region))
  #g <- IRanges::subsetByOverlaps(genes(db), region)

  # Plot transcripts
  grid::downViewport("tx")
  #grid::grid.rect(gp = grid::gpar(col = "grey"))
  plotTxModels(db = DB, g = gid, rangex = xlim)
  grid::upViewport()

  # If transcript is on plus strand, then plot linear junction in upper panel and backsplice junctions in lower panel, else vice versa
  o <- switch(as.character(GenomicRanges::strand(region)), "+" = c("upper", "lower"), "-" = c("lower", "upper"))

  grid::downViewport(o[1])
  plotSJ(df = IRanges::subsetByOverlaps(sj, region), rangex = xlim)
  grid::upViewport()

  grid::downViewport(o[2])
  plotBSJ(df = IRanges::subsetByOverlaps(bsj, region), rangex = xlim)
  grid::upViewport()
}


#' Function to plot coverage as barplot.
#'
#' @param df GRanges object returned by calcCoverage
#' @param rangex The range in genomic coordinates for the region of interest. This is important for correct alignment of the tx model plot with splice junction plot.
#'
#' @importFrom grid grid.lines pushViewport unit upViewport viewport
#' @export
plotCoverage <- function(df = NULL, rangex){
  rangey = c(0,max(df$cov)*1.1)

  # Make top viewport
  grid::pushViewport(grid::viewport(xscale = rangex, yscale = rangey))
  for(i in 1:length(df)){
    grid::grid.lines(x = grid::unit(rep(BiocGenerics::start(df)[i],2), units = "native"), y = grid::unit(c(0, df$cov[i]), unit = "native"))
  }
  grid::upViewport()
}


#' Create a track of linear splice junction
#'
#' Function to create a CustomTrack of linear splice junctions for plotting using Gviz.
#' @param sj Data on linear splice sites. Should be a GRanges object returned by readSJout. Ideally, the data should be restricted to the genomic reange of a single gene (e.g. GenomicFeatures::subsetByOverlaps(sj, region.of.intereset, type = "within")).
#'
#' @importFrom Gviz CustomTrack
#' @importFrom grid gpar grid.bezier grid.text unit
#' @export
SJTrack <- function(sj = NULL){
  ctrack <- Gviz::CustomTrack(
    plottingFunction = function(GdObject, prepare, df = sj) {
      if(length(df) == 0){
        if(!prepare){
          grid::grid.text("No linear splice junctions in this region")
        }
        return(invisible(GdObject))
      } else {
        if(!prepare){
          # Scale height to length of junction. Height should be in the range [0.1 , 0.95].
          # The min value should be determined by the number of junctions to plot
          scale.range <- c(ifelse(length(df)<10, 1-length(df)*0.1, .1), 0.95)
          df$height <- BiocGenerics::width(df)
          if(length(df) == 1){
            df$height <- 0.95
          } else {
            df$height <- (df$height-min(df$height)) * diff(scale.range) / (max(df$height)-min(df$height)) + min(scale.range)
          }

          df$height <- 1-df$height
          y.start = 0.95

          ###! Rotate LSJ track relative to strandedness

          # tmp <- GenomicRanges::strand(df)
          # tmp <- tmp[tmp!="*"]
          #
          # if(all(tmp == "-")){
            # df$height <- 1-df$height
            # y.start = 0.95
          # } else {
          #   y.start = 0.05
          # }

          # Scale line width to number of reads.
          # ## Harcode values
          # # 1-10:   1
          # # 11-25:  3
          # # 26-75:  5
          # # 76-150: 7
          # # 151-400:9
          # # 401+:   11
          # lw <- cut(df$read.count.unique, breaks = c(1,11,26,71,151,401, Inf), include.lowest = T)
          # levels(lw) <- seq(from = 1, to = 11, by = 2)
          # df$w <- as.integer(as.character(lw))
          ## log2 transform
          df$w <- log(df$read.count.unique+1.2, 2)

          for (i in 1:length(df)) {
            x <- grid::unit(rep(c(BiocGenerics::start(df[i]), BiocGenerics::end(df[i])), each=2), "native")
            y <- grid::unit(c(y.start, df$height[i], df$height[i], y.start), "npc")

            grid::grid.bezier(x, y, gp=grid::gpar(lwd = df$w[i], lex = 0.5, col = "blue"))
          }
        }
        return(invisible(GdObject))
      }
    },
    name = "LSJ"
  )
  return(ctrack)
}

#' Create track of backsplice sites
#'
#' Function to create a CustomTrack of backsplice junctions for plotting using Gviz
#' @param sj Data on backsplice sites. Should be a GRanges object returned by readBSJout or countSumChimFile.
#'
#' @importFrom Gviz CustomTrack
#' @importFrom grid gpar grid.text grid.xspline unit
#' @export
BSJTrack <- function(bsj = NULL){
  ctrack <- Gviz::CustomTrack(
    plottingFunction = function(GdObject, prepare, df = bsj) {
      if(length(df) == 0){
        if(!prepare){
          grid::grid.text("No backsplice junctions in this region")
        }
        return(invisible(GdObject))
      } else {
        if(!prepare){
          # Scale height to length of junction. Height should be in the range [0.15 , 0.85].
          if(length(df) > 1){
            scale.range <- c(ifelse(length(df) < 8, 0.85 - length(df)*0.1, 0.15), 0.85)
            df$height <- BiocGenerics::width(df)
            df$height <- (df$height-min(df$height)) * diff(scale.range) / (max(df$height)-min(df$height)) + min(scale.range)
          } else {
            df$height = 0.85
          }

          y.start = 0.05

          ###! Rotate BSJ track relative to strandedness

          # if(all(GenomicRanges::strand(df) == "+")){
          #   df$height <- 1-df$height
          #   y.start = 0.95
          # } else {
          #   y.start = 0.05
          # }

          # Scale line width to number of reads.
          # ## Harcode values
          # # 1-10:   1
          # # 11-25:  3
          # # 26-75:  5
          # # 76-150: 7
          # # 151-400:9
          # # 401+:   11
          # lw <- cut(df$read.count.unique, breaks = c(1,11,26,71,151,401, Inf), include.lowest = T)
          # levels(lw) <- seq(from = 1, to = 11, by = 2)
          # df$w <- as.integer(as.character(lw))
          ## log2 transform

          df$w <- log(df$count+1.2, 2)

          for (i in 1:length(df)) {
            x <- grid::unit(rep(c(BiocGenerics::start(df[i]), BiocGenerics::end(df[i])), each=2), "native")
            y <- grid::unit(c(y.start, df$height[i], df$height[i], y.start), "npc")

            grid::grid.xspline(x, y, shape=c(0, 1, 1, 0), open=TRUE, gp=grid::gpar(lwd = df$w[i], lex = 0.5, col = "red"))
          }
        }
        return(invisible(GdObject))
      }
    },
    name = "BSJ"
  )
  return(ctrack)
}


