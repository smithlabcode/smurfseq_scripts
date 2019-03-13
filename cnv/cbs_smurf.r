#!/usr/bin/env Rscript

library("DNAcopy")

lowess.gc <- function(jtkx, jtky) {
  jtklow <- lowess(jtkx, log(jtky), f=0.05)
  jtkz <- approx(jtklow$x, jtklow$y, jtkx)
  return(exp(log(jtky) - jtkz$y))
}

format.progress.message <- function(segs, idx) {
  fields <- c('chrom', 'loc.start', 'loc.end',
              'num.mark', 'seg.mean', 'seg.start', 'seg.end')
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
remove.segment <- function(rsShort, rsSegnum, ratioData, sd.undo) {

  appendLeft <- TRUE
  checkSdundo <- FALSE

  if (rsSegnum == 1) { # current is first segment
    appendLeft <- FALSE
  }
  else if (rsSegnum == nrow(rsShort)) {
    appendLeft <- TRUE
  }
  else {
    right.idx <- rsSegnum + 1
    left.idx <- rsSegnum - 1

    if (rsShort[right.idx, "chrom"] != rsShort[rsSegnum, "chrom"]) {
      appendLeft <- TRUE
    }
    else {
      if (rsShort[left.idx, "chrom"] != rsShort[rsSegnum, "chrom"]) {
        appendLeft <- FALSE
      }
      else {
        if (abs(rsShort[left.idx, "seg.mean"] - rsShort[rsSegnum, "seg.mean"]) <
            abs(rsShort[right.idx, "seg.mean"] - rsShort[rsSegnum, "seg.mean"])) {
          appendLeft <- TRUE
          checkSdundo <- TRUE
        }
        else {
          appendLeft <- FALSE
          checkSdundo <- TRUE
        }
      }
    }
  }
  apnd.idx <- 0
  if (appendLeft) {
    apnd.idx <- rsSegnum - 1
  }
  else {
    apnd.idx <- rsSegnum + 1
  }

  segs <- rsShort
  if (appendLeft) {
    segs[apnd.idx, "loc.end"] <- segs[rsSegnum, "loc.end"]
    segs[apnd.idx, "seg.end"] <- segs[rsSegnum, "seg.end"]
  }
  else {
    segs[apnd.idx, "loc.start"] <- segs[rsSegnum, "loc.start"]
    segs[apnd.idx, "seg.start"] <- segs[rsSegnum, "seg.start"]
  }

  segs[apnd.idx, "num.mark"] <- segs[apnd.idx, "num.mark"] + segs[rsSegnum, "num.mark"]
  segs[apnd.idx, "seg.mean"] <- mean(log2(ratioData$lowratio[segs[apnd.idx, "seg.start"]:segs[apnd.idx, "seg.end"]]))

  cat('append', format.progress.message(segs, apnd.idx), '\n');

  segs <- segs[-rsSegnum, ]
  segs$segnum <- seq(1:nrow(segs))

  if (checkSdundo) {
    thisSd <- -1
    if (appendLeft) {
      left.idx <- apnd.idx
      right.idx <- apnd.idx + 1
    }
    else {
      left.idx <- apnd.idx - 2
      right.idx <- apnd.idx - 1
    }
    ##thisSd <- sd(ratioData[segs$seg.start[left.idx]:segs$seg.start[right.idx], "lowratio"])
    thisSd <- mad(diff(ratioData[, "lowratio"])) / sqrt(2)

    if (abs(segs$seg.mean[left.idx] -
            segs$seg.mean[right.idx]) < (sd.undo * thisSd)) {

      cat('left', format.progress.message(segs, left.idx), '\n');
      cat('right', format.progress.message(segs, right.idx), '\n');

      ##  remove breakpoint
      segs[left.idx, "loc.end"] <- segs[right.idx, "loc.end"]
      segs[left.idx, "seg.end"] <- segs[right.idx, "seg.end"]
      segs[left.idx, "num.mark"] <- segs[left.idx, "num.mark"] + segs[right.idx, "num.mark"]
      segs[left.idx, "seg.mean"] <- mean(log2(ratioData$lowratio[segs[left.idx, "seg.start"]:segs[right.idx, "seg.end"]]))
      segs <- segs[-right.idx, ]
      segs$segnum <- seq(1:nrow(segs))
    }
  }

  return(segs)
}


sdundo.all <- function (sdShort, ratioData, sd.undo) {

  segs <- sdShort
  thisSd <- mad(diff(ratioData[, "lowratio"])) / sqrt(2)

  while (TRUE) {

    chrom <- segs$chrom
    chrom.shift <- c(segs$chrom[-1], segs$chrom[1])

### RISH: This doesn't seem to be working
    ## breakpoints <- which(chrom == chrom.shift)
### RISH: Fixed by adding as.numeric
    breakpoints <- which(as.numeric(chrom) == chrom.shift)
    cat("sdundo.all intrachrom breakpoints", length(breakpoints), "\n")

    if (length(breakpoints) < 1) {
      break
    }

    breakpoints.shift <- breakpoints + 1

    undo.breakpoints <- breakpoints[which(abs(segs$seg.mean[breakpoints] -
                                              segs$seg.mean[breakpoints.shift]) < thisSd*sd.undo)]

    cat("sdundo.all undo breakpoints", length(undo.breakpoints), "\n")

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

    cat("sdundo.all left", format.progress.message(segs, left.idx), "\n");
    cat("sdundo.all right", format.progress.message(segs, right.idx), "\n");

    segs[left.idx, "loc.end"] <- segs[right.idx, "loc.end"]
    segs[left.idx, "seg.end"] <- segs[right.idx, "seg.end"]
    segs[left.idx, "num.mark"] <- segs[left.idx, "num.mark"] + segs[right.idx, "num.mark"]
    segs[left.idx, "seg.mean"] <- mean(log2(ratioData$lowratio[segs[left.idx, "seg.start"]:segs[right.idx, "seg.end"]]))
    segs <- segs[-right.idx, ]
    segs$segnum <- seq(1:nrow(segs))
  }

  return(segs)
}

cbs.segment01 <- function(indir, outdir,
                          varbin.gc, bad.bins.file,
                          varbin.data, sample.name,
                          alt.sample.name, alpha,
                          nperm, undo.SD, min.width) {

  gc <- read.table(varbin.gc, header=T)

  chrom.numeric <- substring(gc$bin.chrom, 4)
  chrom.numeric[which(gc$bin.chrom == "chrX")] <- "23"
  chrom.numeric[which(gc$bin.chrom == "chrY")] <- "24"
  chrom.numeric <- as.numeric(chrom.numeric)

  thisRatio <- read.table(paste(indir, varbin.data, sep="/"), header=F)
  names(thisRatio) <- c("chrom", "chrompos", "abspos", "bincount", "ratio")
  thisRatio$chrom <- chrom.numeric
  a <- thisRatio$bincount + 1
  thisRatio$ratio <- a / mean(a)
  thisRatio$gc.content <- gc$gc.content
  thisRatio$lowratio <- lowess.gc(thisRatio$gc.content, thisRatio$ratio)

  bad <- read.table(bad.bins.file, header=F, as.is=T, stringsAsFactors=F)
  thisRatioNobig <- thisRatio[-bad[, 1], ]

  set.seed(25)
  CNA.object <- CNA(log(thisRatioNobig$lowratio, base=2),
                    gc$chrom.arm[-bad[, 1]],
                    thisRatioNobig$chrompos,
                    data.type="logratio",
                    sampleid=sample.name)
  smoothed.CNA.object <- smooth.CNA(CNA.object)
  segment.smoothed.CNA.object <- segment(smoothed.CNA.object,
                                         alpha=alpha, nperm=nperm,
                                         undo.splits="sdundo",
                                         undo.SD=undo.SD, min.width=2)
  segs <- segment.smoothed.CNA.object[[2]]

  ## RISH: This can probably be removed by ordering the CNV object
  sortcol <- segs$chrom
  sortcol <- gsub("chr", "", sortcol)
  sortcol <- gsub("p", "", sortcol)
  sortcol <- gsub("q", "", sortcol)
  segs <- segs[order(as.numeric(sortcol)), ]

#####  NEW STUFF  also check min.width=2 above

  work.segs <- segs
  work.segs$segnum <- 0
  work.segs$seg.start <- 0
  work.segs$seg.end <- 0
  prevEnd <- 0
  for (i in 1:nrow(segs)) {
    thisStart <- prevEnd + 1
    thisEnd <- prevEnd + segs$num.mark[i]
    work.segs$seg.start[i] <- thisStart
    work.segs$seg.end[i] <- thisEnd
    work.segs$segnum[i] <- i
    prevEnd = thisEnd
  }

  discard.segs <- TRUE
  while (discard.segs) {
    work.segs.ord <- work.segs[order(work.segs$num.mark, abs(work.segs$seg.mean)), ]
    if (work.segs.ord[1, "num.mark"] < min.width) {
      work.segs <- remove.segment(work.segs, work.segs.ord[1, "segnum"], thisRatioNobig, undo.SD)
    }
    else {
      discard.segs <- FALSE
    }
  }
  work.segs <- sdundo.all(work.segs, thisRatioNobig, undo.SD)
  segs <- work.segs

#####  END NEW STUFF

  m <- matrix(data=0, nrow=nrow(thisRatioNobig), ncol=1)
  prevEnd <- 0
  for (i in 1:nrow(segs)) {
    thisStart <- prevEnd + 1
    thisEnd <- prevEnd + segs$num.mark[i]
    m[thisStart:thisEnd, 1] <- 2^segs$seg.mean[i]
    prevEnd = thisEnd
  }
  thisRatioNobig$seg.mean.LOWESS <- m[, 1]

  chr <- thisRatioNobig$chrom
  chr.shift <- c(chr[-1], chr[length(chr)])

  vlines <- c(1, thisRatio$abspos[which(chr != chr.shift) + 1], thisRatio$abspos[nrow(thisRatio)])
  hlines <- c(0.5, 1.0, 1.5, 2.0)
  chr.text <- c(1:22, "X", "Y")
  vlines.shift <- c(vlines[-1], 4*10^9)
  chr.at <- vlines + (vlines.shift - vlines) / 2
  x.at <- c(0, 0.5, 1, 1.5, 2, 2.5, 3) * 10^9
  x.labels <- c("0", "0.5", "1.0", "1.5", "2.0", "2.5", "3.0")
  y.at <- c(0.005, 0.020, 0.100, 0.500, 1.000, 2.000, 10, 100)
  y.labels <- c("0.005", "0.020", "0.100", "0.5", "1", "2", "10", "100")

  pdf(paste(outdir, "/", sample.name, ".5k.wg.nobad.pdf", sep=""), height=3.5, width=6, useDingbats=FALSE)
  par(pin=c(5.0,1.75))
  plot(x=thisRatioNobig$abspos,
       y=thisRatioNobig$lowratio,
       log="y", main=paste(sample.name, alt.sample.name),
       xaxt="n", xlab="Genome Position Gb",
       yaxt="n", ylab="Ratio", col="#517FFF", cex=0.01)

  ## axis(1, at=x.at, labels=x.labels)
  axis(2, at=y.at, labels=y.labels)
  ## lines(x=thisRatioNobig$abspos, y=thisRatioNobig$lowratio, col="#CCCCCC")
  ## points(x=thisRatioNobig$abspos, y=thisRatioNobig$seg.mean.LOWESS, col="#0000AA")
  lines(x=thisRatioNobig$abspos, y=thisRatioNobig$seg.mean.LOWESS, col="red", cex=1.0)
  ## abline(h=hlines)
  ## abline(v=vlines)
  abline(v=vlines, lwd=0.1, col="grey")
  ## mtext(chr.text, at = chr.at)
  dev.off()

  write.table(thisRatioNobig, sep="\t",
              file=paste(outdir, "/", sample.name,
                         ".hg19.5k.nobad.varbin.data.txt", sep=""),
              quote=F, row.names=F)
  write.table(segs, sep="\t",
              file=paste(outdir, "/", sample.name,
                         ".hg19.5k.nobad.varbin.short.txt", sep=""),
              quote=F, row.names=F)
}

main <- function() {

  alpha.value  <- 0.02
  n.permutations  <- 1000
  standard.dev <- 0.5
  min.width <- 4

  args = commandArgs(trailingOnly=TRUE)
  if (length(args) != 4) {
    stop("cbs_smurf.r <varbin-file> <sample-name> ",
         "<gc-content-file> <bad-bin-file>", call.=FALSE)
  }

  varbin_file <- args[1]
  sample_name <- args[2]
  gc_file <- args[3]
  badbins_file <- args[4]

  cbs.segment01(indir=".", outdir=".",
                varbin.gc=gc_file, bad.bins.file=badbins_file,
                varbin.data=varbin_file, sample.name=sample_name,
                alt.sample.name="",
                alpha=alpha.value,
                nperm=n.permutations,
                undo.SD=standard.dev,
                min.width=min.width)
}

main()
