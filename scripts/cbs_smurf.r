#!/usr/bin/env Rscript

# cbs_smurf.r: xxxxxx
#
# Copyright (C) 2019 Rish Prabakar and Andrew D Smith
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

library("DNAcopy")

LowessGc <- function(jtkx, jtky) {
  jtklow <- lowess(jtkx, log(jtky), f=0.05)
  jtkz <- approx(jtklow$x, jtklow$y, jtkx)
  return(exp(log(jtky) - jtkz$y))
}

FormatProgressMessage <- function(segs, idx) {
  fields <- c('chrom', 'loc.start', 'loc.end', 'num.mark', 
              'seg.mean', 'seg.start', 'seg.end')
  fld.format <- c('%s', '%d', '%d', '%d', '%f', '%d', '%d')
  mesg <- ''
  for (i in 1:length(fields)) {
    mesg <- paste(mesg, sprintf(fld.format[i], segs[idx, fields[i]]))
  }
  return(mesg)
}

### The remove segment function removes segments that are too short,
### and merges the removed segment with either the left or right
### neighbor, depending on which one seems more similar in terms of
### the mean of that segment.
RemoveSegment <- function(rs.short, rs.seg.num, ratio.data, sd.undo) {

  append.left <- TRUE
  check.sd.undo <- FALSE

  if (rs.seg.num == 1) { # current is first segment
    append.left <- FALSE
  } else if (rs.seg.num == nrow(rs.short)) {
    append.left <- TRUE
  } else {
    right.idx <- rs.seg.num + 1
    left.idx <- rs.seg.num - 1

    if (rs.short[right.idx, "chrom"] != rs.short[rs.seg.num, "chrom"]) {
      append.left <- TRUE
    } else {
      if (rs.short[left.idx, "chrom"] != rs.short[rs.seg.num, "chrom"]) {
        append.left <- FALSE
      } else {
        if (abs(rs.short[left.idx, "seg.mean"] - rs.short[rs.seg.num, "seg.mean"]) <
            abs(rs.short[right.idx, "seg.mean"] - rs.short[rs.seg.num, "seg.mean"])) {
          append.left <- TRUE
          check.sd.undo <- TRUE
        } else {
          append.left <- FALSE
          check.sd.undo <- TRUE
        }
      }
    }
  }
  apnd.idx <- 0
  if (append.left) {
    apnd.idx <- rs.seg.num - 1
  } else {
    apnd.idx <- rs.seg.num + 1
  }

  segs <- rs.short
  if (append.left) {
    segs[apnd.idx, "loc.end"] <- segs[rs.seg.num, "loc.end"]
    segs[apnd.idx, "seg.end"] <- segs[rs.seg.num, "seg.end"]
  } else {
    segs[apnd.idx, "loc.start"] <- segs[rs.seg.num, "loc.start"]
    segs[apnd.idx, "seg.start"] <- segs[rs.seg.num, "seg.start"]
  }

  segs[apnd.idx, "num.mark"] <- segs[apnd.idx, "num.mark"] + segs[rs.seg.num, "num.mark"]
  segs[apnd.idx, "seg.mean"] <- mean(log2(ratio.data$lowratio[segs[apnd.idx, "seg.start"]:segs[apnd.idx, "seg.end"]]))

  cat('append', FormatProgressMessage(segs, apnd.idx), '\n');

  segs <- segs[-rs.seg.num, ]
  segs$segnum <- seq(1:nrow(segs))

  if (check.sd.undo) {
    cur.sd <- -1
    if (append.left) {
      left.idx <- apnd.idx
      right.idx <- apnd.idx + 1
    } else {
      left.idx <- apnd.idx - 2
      right.idx <- apnd.idx - 1
    }
    ##cur.sd <- sd(ratio.data[segs$seg.start[left.idx]:segs$seg.start[right.idx], "lowratio"])
    cur.sd <- mad(diff(ratio.data[, "lowratio"])) / sqrt(2)

    if (abs(segs$seg.mean[left.idx] -
            segs$seg.mean[right.idx]) < (sd.undo * cur.sd)) {

      cat('left', FormatProgressMessage(segs, left.idx), '\n');
      cat('right', FormatProgressMessage(segs, right.idx), '\n');

      ##  remove breakpoint
      segs[left.idx, "loc.end"] <- segs[right.idx, "loc.end"]
      segs[left.idx, "seg.end"] <- segs[right.idx, "seg.end"]
      segs[left.idx, "num.mark"] <- segs[left.idx, "num.mark"] + segs[right.idx, "num.mark"]
      segs[left.idx, "seg.mean"] <- mean(log2(ratio.data$lowratio[segs[left.idx, "seg.start"]:segs[right.idx, "seg.end"]]))
      segs <- segs[-right.idx, ]
      segs$segnum <- seq(1:nrow(segs))
    }
  }

  return(segs)
}


SdUndoAll <- function (sd.short, ratio.data, sd.undo) {

  segs <- sd.short
  cur.sd <- mad(diff(ratio.data[, "lowratio"])) / sqrt(2)

  while (TRUE) {

    chrom <- segs$chrom
    chrom.shift <- c(chrom[-1], chrom[1])

    ### RISH: This doesn't seem to be working
    ## breakpoints <- which(chrom == chrom.shift)
    ### RISH: Fixed by adding as.numeric
    breakpoints <- which(as.numeric(chrom) == chrom.shift)
    cat("SdUndoAll intrachrom breakpoints", length(breakpoints), "\n")

    if (length(breakpoints) < 1) {
      break
    }

    breakpoints.shift <- breakpoints + 1

    undo.breakpoints <- breakpoints[which(abs(segs$seg.mean[breakpoints] -
                                              segs$seg.mean[breakpoints.shift]) < cur.sd * sd.undo)]

    cat("SdUndoAll undo breakpoints", length(undo.breakpoints), "\n")

    if (length(undo.breakpoints) < 1) {
      break
    }

    undo.breakpoints.shift <- undo.breakpoints + 1

    undo.df <- segs[undo.breakpoints, ]
    undo.df$seg.mean.diff <- abs(segs$seg.mean[undo.breakpoints] -
                                 segs$seg.mean[undo.breakpoints.shift])

    min.index <- which.min(undo.df$seg.mean.diff)

    left.idx <- undo.df$segnum[min.index]
    right.idx <- left.idx + 1

    cat("SdUndoAll left", FormatProgressMessage(segs, left.idx), "\n");
    cat("SdUndoAll right", FormatProgressMessage(segs, right.idx), "\n");

    segs[left.idx, "loc.end"] <- segs[right.idx, "loc.end"]
    segs[left.idx, "seg.end"] <- segs[right.idx, "seg.end"]
    segs[left.idx, "num.mark"] <- segs[left.idx, "num.mark"] + segs[right.idx, "num.mark"]
    segs[left.idx, "seg.mean"] <- mean(log2(ratio.data$lowratio[segs[left.idx, "seg.start"]:segs[right.idx, "seg.end"]]))
    segs <- segs[-right.idx, ]
    segs$segnum <- seq(1:nrow(segs))
  }

  return(segs)
}


PlotSegment <- function(cur.ratio, cur.ratio.bad, sample.name) {

  chr <- cur.ratio.bad$chrom
  chr.shift <- c(chr[-1], chr[length(chr)])

  vlines <- c(1, cur.ratio$abspos[which(chr != chr.shift) + 1], cur.ratio$abspos[nrow(cur.ratio)])
  hlines <- c(0.5, 1.0, 1.5, 2.0)
  chr.text <- c(1:22, "X", "Y")
  vlines.shift <- c(vlines[-1], 4 * 10^9)
  chr.at <- vlines + (vlines.shift - vlines) / 2
  x.at <- c(0, 0.5, 1, 1.5, 2, 2.5, 3) * 10^9
  x.labels <- c("0", "0.5", "1.0", "1.5", "2.0", "2.5", "3.0")
  y.at <- c(0.005, 0.020, 0.100, 0.500, 1.000, 2.000, 10, 100)
  y.labels <- c("0.005", "0.020", "0.100", "0.5", "1", "2", "10", "100")

  pdf(paste(sample.name, ".5k.wg.nobad.pdf", sep=""), height=3.5, width=6, useDingbats=FALSE)
  par(pin=c(5.0, 1.75))
  plot(x=cur.ratio.bad$abspos,
       y=cur.ratio.bad$lowratio,
       log="y", main=sample.name,
       xaxt="n", xlab="Genome Position Gb",
       yaxt="n", ylab="Ratio", col="#517FFF", cex=0.01)

  ## axis(1, at=x.at, labels=x.labels)
  axis(2, at=y.at, labels=y.labels)
  ## lines(x=cur.ratio.bad$abspos, y=cur.ratio.bad$lowratio, col="#CCCCCC")
  ## points(x=cur.ratio.bad$abspos, y=cur.ratio.bad$seg.mean.LOWESS, col="#0000AA")
  lines(x=cur.ratio.bad$abspos, y=cur.ratio.bad$seg.mean.LOWESS, col="red", cex=1.0)
  ## abline(h=hlines)
  ## abline(v=vlines)
  abline(v=vlines, lwd=0.1, col="grey")
  ## mtext(chr.text, at = chr.at)
  dev.off()
}


CbsSegment01 <- function(varbin.gc, bad.bins.file,
                         varbin.data, sample.name,
                         alpha, nperm, 
                         undo.sd, min.width) {

  ## Load the gc content information
  gc <- read.table(varbin.gc, header=T)

  ## Convert the chromosome names into numeric ids sorted by
  ## chromosome number with sex chroms at the end
  chrom.numeric <- substring(gc$bin.chrom, 4)
  chrom.numeric[which(gc$bin.chrom == "chrX")] <- "23"
  chrom.numeric[which(gc$bin.chrom == "chrY")] <- "24"
  chrom.numeric <- as.numeric(chrom.numeric)

  ## Load the data in the form of bin counts for genomic
  ## positions. The "ratio" column present in the input file will be
  ## over-written shortly.
  cur.ratio <- read.table(varbin.data, header=F)
  names(cur.ratio) <- c("chrom", "chrompos", "abspos", "bincount", "ratio")
  cur.ratio$chrom <- chrom.numeric

  ## Take the fractional counts after using Laplaces correction
  a <- cur.ratio$bincount + 1
  cur.ratio$ratio <- a / mean(a)
  cur.ratio$gc.content <- gc$gc.content
  cur.ratio$lowratio <- LowessGc(cur.ratio$gc.content, cur.ratio$ratio)

  ## Load the "bad" bins, which are pre-determined as having problems
  ## due to technical issues with the genome or the sequencing, etc.
  bad.bins <- read.table(bad.bins.file, header=F, as.is=T, stringsAsFactors=F)
  cur.ratio.bad <- cur.ratio[-bad.bins[, 1], ]

  set.seed(25)

  ## Do the work of the copy number analysis using CNA, and
  ## immediately apply smoothing to the result.
  cna.result <- smooth.CNA(CNA(log2(cur.ratio.bad$lowratio),
                               gc$chrom.arm[-bad.bins[, 1]],
                               cur.ratio.bad$chrompos,
                               data.type="logratio",
                               sampleid=sample.name))

  ## Obtain segments from the smoothed CNA result, using only the
  ## short form
  segs <- segment(cna.result, alpha=alpha, nperm=nperm,
                  undo.splits="sdundo", undo.SD=undo.sd, min.width=2)[[2]]

  ## RISH: This can probably be removed by ordering the CNV object
  sortcol <- segs$chrom
  sortcol <- gsub("chr", "", sortcol)
  sortcol <- gsub("p", "", sortcol)
  sortcol <- gsub("q", "", sortcol)
  segs <- segs[order(as.numeric(sortcol)), ]

  #####  NEW STUFF  also check min.width=2 above
  work.segs <- segs
  work.segs$segnum <- c()
  work.segs$seg.start <- c()
  work.segs$seg.end <- c()
  prev.end <- 0
  for (i in 1:nrow(segs)) {
    work.segs$seg.start[i] <- prev.end + 1
    curr.end <- prev.end + segs$num.mark[i]
    work.segs$seg.end[i] <- curr.end
    work.segs$segnum[i] <- i
    prev.end <- curr.end
  }

  discard.segs <- TRUE
  while (discard.segs) {
    work.segs.ord <- work.segs[order(work.segs$num.mark, abs(work.segs$seg.mean)), ]
    if (work.segs.ord[1, "num.mark"] < min.width) {
      work.segs <- RemoveSegment(work.segs, work.segs.ord[1, "segnum"], cur.ratio.bad, undo.sd)
    } else {
      discard.segs <- FALSE
    }
  }
  work.segs <- SdUndoAll(work.segs, cur.ratio.bad, undo.sd)
  segs <- work.segs
  #####  END NEW STUFF

  m <- matrix(data=0, nrow=nrow(cur.ratio.bad), ncol=1)
  prev.end <- 0
  for (i in 1:nrow(segs)) {
    cur.start <- prev.end + 1
    cur.end <- prev.end + segs$num.mark[i]
    m[cur.start:cur.end, 1] <- 2^segs$seg.mean[i]
    prev.end <- cur.end
  }
  cur.ratio.bad$seg.mean.LOWESS <- m[, 1]

  return(list(ratio=cur.ratio, ratio.bad=cur.ratio.bad, segs=segs))
}

main <- function() {

  # CONSTANTS
  kAlphaValue  <- 0.02
  kNPermutations  <- 1000
  kStandardDev <- 0.5
  kMinWidth <- 4

  args = commandArgs(trailingOnly=TRUE)
  if (length(args) != 4) {
    stop("cbs_smurf.r <varbin-file> <sample-name> ",
         "<gc-content-file> <bad-bin-file>", call.=FALSE)
  }

  varbin.file <- args[1]
  sample.name <- args[2]
  gc.file <- args[3]
  bad.bins.file <- args[4]

  cbs.seg = CbsSegment01(varbin.gc=gc.file, bad.bins.file=bad.bins.file,
                         varbin.data=varbin.file, sample.name=sample.name,
                         alpha=kAlphaValue,
                         nperm=kNPermutations,
                         undo.sd=kStandardDev,
                         min.width=kMinWidth)
  cur.ratio = cbs.seg$ratio
  cur.ratio.bad = cbs.seg$ratio.bad
  segs = cbs.seg$segs
  # plot segment
  PlotSegment(cur.ratio, cur.ratio.bad, sample.name)
  # save results
  write.table(cur.ratio.bad, sep="\t",
              file=paste(sample.name, ".hg19.5k.nobad.varbin.data.txt", sep=""),
              quote=F, row.names=F)
  write.table(segs, sep="\t",
              file=paste(sample.name, ".hg19.5k.nobad.varbin.short.txt", sep=""),
              quote=F, row.names=F)
}

main()
