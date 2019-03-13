#!/usr/bin/env Rscript

library("DNAcopy")

lowess.gc <- function(jtkx, jtky) {
  jtklow <- lowess(jtkx, log(jtky), f=0.05)
  jtkz <- approx(jtklow$x, jtklow$y, jtkx)
  return(exp(log(jtky) - jtkz$y))
}

remove.segment <- function( rsShort, rsSegnum, ratioData, sd.undo ) {

  appendLeft <- TRUE
  checkSdundo <- FALSE

  if (rsSegnum == 1) {
    appendLeft <- FALSE
  }
  else {
    if (rsSegnum == nrow(rsShort)) {
      appendLeft <- TRUE
    }
    else {
      rightIndex <- rsSegnum + 1
      leftIndex <- rsSegnum - 1

      if (rsShort[rightIndex, "chrom"] != rsShort[rsSegnum, "chrom"]) {
        appendLeft <- TRUE
      }
      else {
        if (rsShort[leftIndex, "chrom"] != rsShort[rsSegnum, "chrom"]) {
          appendLeft <- FALSE
        }
        else {
          if (abs(rsShort[leftIndex, "seg.mean"] - rsShort[rsSegnum, "seg.mean"]) <
              abs(rsShort[rightIndex, "seg.mean"] - rsShort[rsSegnum, "seg.mean"])) {
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
  }
  apndIdx <- 0
  if (appendLeft) {
    apndIdx <- rsSegnum - 1
  }
  else {
    apndIdx <- rsSegnum + 1
  }

  tempShort <- rsShort
  if (appendLeft) {
    tempShort[apndIdx, "loc.end"] <- tempShort[rsSegnum, "loc.end"]
    tempShort[apndIdx, "seg.end"] <- tempShort[rsSegnum, "seg.end"]
  }
  else {
    tempShort[apndIdx, "loc.start"] <- tempShort[rsSegnum, "loc.start"]
    tempShort[apndIdx, "seg.start"] <- tempShort[rsSegnum, "seg.start"]
  }

  tempShort[apndIdx, "num.mark"] <- tempShort[apndIdx, "num.mark"] + tempShort[rsSegnum, "num.mark"]
  tempShort[apndIdx, "seg.mean"] <- mean(log(ratioData$lowratio[tempShort[apndIdx, "seg.start"]:tempShort[apndIdx, "seg.end"]], base=2))

  cat("append",
      tempShort[apndIdx, "chrom"],
      tempShort[apndIdx, "loc.start"],
      tempShort[apndIdx, "loc.end"],
      tempShort[apndIdx, "num.mark"],
      tempShort[apndIdx, "seg.mean"],
      tempShort[apndIdx, "seg.start"],
      tempShort[apndIdx, "seg.end"], "\n")

  tempShort <- tempShort[-rsSegnum, ]
  tempShort$segnum <- seq(1:nrow(tempShort))

  if (checkSdundo) {
    thisSd <- -1
    if (appendLeft) {
      leftIndex <- apndIdx
      rightIndex <- apndIdx + 1
    }
    else {
      leftIndex <- apndIdx - 2
      rightIndex <- apndIdx - 1
    }
    ##thisSd <- sd(ratioData[tempShort$seg.start[leftIndex]:tempShort$seg.start[rightIndex], "lowratio"])
    thisSd <- mad(diff(ratioData[, "lowratio"])) / sqrt(2)

    if (abs(tempShort$seg.mean[leftIndex] -
            tempShort$seg.mean[rightIndex]) < (sd.undo * thisSd) ) {

      cat("left",
          tempShort[leftIndex, "chrom"],
          tempShort[leftIndex, "loc.start"],
          tempShort[leftIndex, "loc.end"],
          tempShort[leftIndex, "num.mark"],
          tempShort[leftIndex, "seg.mean"],
          tempShort[leftIndex, "seg.start"],
          tempShort[leftIndex, "seg.end"], "\n")

      cat("right",
          tempShort[rightIndex, "chrom"],
          tempShort[rightIndex, "loc.start"],
          tempShort[rightIndex, "loc.end"],
          tempShort[rightIndex, "num.mark"],
          tempShort[rightIndex, "seg.mean"],
          tempShort[rightIndex, "seg.start"],
          tempShort[rightIndex, "seg.end"], "\n")

      ##  remove breakpoint
      tempShort[leftIndex, "loc.end"] <- tempShort[rightIndex, "loc.end"]
      tempShort[leftIndex, "seg.end"] <- tempShort[rightIndex, "seg.end"]
      tempShort[leftIndex, "num.mark"] <- tempShort[leftIndex, "num.mark"] + tempShort[rightIndex, "num.mark"]
      tempShort[leftIndex, "seg.mean"] <- mean(log(ratioData$lowratio[tempShort[leftIndex, "seg.start"]:tempShort[rightIndex, "seg.end"]], base=2))
      tempShort <- tempShort[-rightIndex, ]
      tempShort$segnum <- seq(1:nrow(tempShort))
    }
  }

  return(tempShort)
}


sdundo.all <- function (sdShort, ratioData, sd.undo) {

  tempShort <- sdShort
  thisSd <- mad(diff(ratioData[, "lowratio"])) / sqrt(2)

  while ( TRUE ) {

    chrom <- tempShort$chrom
    chrom.shift <- c(tempShort$chrom[-1], tempShort$chrom[1])

### RISH: This doesn't seem to be working
    ## breakpoints <- which(chrom == chrom.shift)
### RISH: Fixed by adding as.numeric
    breakpoints <- which(as.numeric(chrom) == chrom.shift)
    cat("sdundo.all intrachrom breakpoints", length(breakpoints), "\n")

    if (length(breakpoints) < 1) {
      break
    }

    breakpoints.shift <- breakpoints + 1

    undo.breakpoints <- breakpoints[which(abs(tempShort$seg.mean[breakpoints] -
                                              tempShort$seg.mean[breakpoints.shift]) < thisSd*sd.undo)]

    cat("sdundo.all undo breakpoints", length(undo.breakpoints), "\n")

    if (length(undo.breakpoints) < 1) {
      break
    }

    undo.breakpoints.shift <- undo.breakpoints + 1

    undo.df <- tempShort[undo.breakpoints, ]
    undo.df$seg.mean.diff <- abs(tempShort$seg.mean[undo.breakpoints] -
                                 tempShort$seg.mean[undo.breakpoints.shift])

    min.index <- which.min(undo.df$seg.mean.diff)

    leftIndex <- undo.df$segnum[min.index]
    rightIndex <- leftIndex + 1

    cat("sdundo.all left",
        tempShort[leftIndex, "chrom"],
        tempShort[leftIndex, "loc.start"],
        tempShort[leftIndex, "loc.end"],
        tempShort[leftIndex, "num.mark"],
        tempShort[leftIndex, "seg.mean"],
        tempShort[leftIndex, "seg.start"],
        tempShort[leftIndex, "seg.end"], "\n")

    cat("sdundo.all right",
        tempShort[rightIndex, "chrom"],
        tempShort[rightIndex, "loc.start"],
        tempShort[rightIndex, "loc.end"],
        tempShort[rightIndex, "num.mark"],
        tempShort[rightIndex, "seg.mean"],
        tempShort[rightIndex, "seg.start"],
        tempShort[rightIndex, "seg.end"], "\n")

    tempShort[leftIndex, "loc.end"] <- tempShort[rightIndex, "loc.end"]
    tempShort[leftIndex, "seg.end"] <- tempShort[rightIndex, "seg.end"]
    tempShort[leftIndex, "num.mark"] <- tempShort[leftIndex, "num.mark"] + tempShort[rightIndex, "num.mark"]
    tempShort[leftIndex, "seg.mean"] <- mean(log(ratioData$lowratio[tempShort[leftIndex, "seg.start"]:tempShort[rightIndex, "seg.end"]], base=2))
    tempShort <- tempShort[-rightIndex, ]
    tempShort$segnum <- seq(1:nrow(tempShort))
  }

  return(tempShort)
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
  thisShort <- segment.smoothed.CNA.object[[2]]

  ## RISH: This can probably be removed by ordering the CNV object
  sortcol <- thisShort$chrom
  sortcol <- gsub("chr", "", sortcol)
  sortcol <- gsub("p", "", sortcol)
  sortcol <- gsub("q", "", sortcol)
  thisShort <- thisShort[order(as.numeric(sortcol)), ]

#####  NEW STUFF  also check min.width=2 above

  workShort <- thisShort
  workShort$segnum <- 0
  workShort$seg.start <- 0
  workShort$seg.end <- 0
  prevEnd <- 0
  for (i in 1:nrow(thisShort)) {
    thisStart <- prevEnd + 1
    thisEnd <- prevEnd + thisShort$num.mark[i]
    workShort$seg.start[i] <- thisStart
    workShort$seg.end[i] <- thisEnd
    workShort$segnum[i] <- i
    prevEnd = thisEnd
  }

  discardSegments <- TRUE
  while (discardSegments) {
    orderShort <- workShort[order(workShort$num.mark, abs(workShort$seg.mean)), ]
    if (orderShort[1, "num.mark"] < min.width) {
      workShort <- remove.segment(workShort, orderShort[1, "segnum"], thisRatioNobig, undo.SD)
    } else {
      discardSegments <- FALSE
    }
  }
  workShort <- sdundo.all(workShort, thisRatioNobig, undo.SD)
  thisShort <- workShort

#####  END NEW STUFF

  m <- matrix(data=0, nrow=nrow(thisRatioNobig), ncol=1)
  prevEnd <- 0
  for (i in 1:nrow(thisShort)) {
    thisStart <- prevEnd + 1
    thisEnd <- prevEnd + thisShort$num.mark[i]
    m[thisStart:thisEnd, 1] <- 2^thisShort$seg.mean[i]
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
  write.table(thisShort, sep="\t",
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
