counts2logR <-function(mat,f=0.2){

    out <- mat[mat$keep==1,]
    out$gcpct <- rep(NA_real_, nrow(out))
    # loop thru chromosomes
    nchr <- max(mat$chrom) # IMPACT doesn't have X so only 22
    for (i in 1:nchr) {
        ii <- out$chrom==i
        jj <- ceiling((out$maploc[ii]-450)/100)
        jj[jj < 1] <- 1 # fix issues (maploc < 500) with chr17
        out$gcpct[ii] <- gcpctdb[[i]][jj]
    }
    ##### log-ratio with gc correction and maf log-odds ratio steps
    chrom <- out$chrom
    maploc <- out$maploc
    rCountN <- out$rCountN
    rCountT <- out$rCountT
    gcpct <- out$gcpct

    # compute gc bias
    ncount <- tapply(rCountN, gcpct, sum)
    tcount <- tapply(rCountT, gcpct, sum)
    pctgc <- as.numeric(names(ncount))
    tscl <- sum(ncount)/sum(tcount)
    gcb <- lowess(pctgc, log2(tcount*tscl)-log2(ncount), f=f)
    jj <- match(gcpct, gcb$x)
    gcbias <- gcb$y[jj]

    # compute cn log-ratio (gc corrected)
    cnlr <- log2(1+rCountT*tscl) - log2(1+rCountN) - gcbias

    out$cnlr <- out$gcbias <- rep(NA_real_, nrow(out))
    out$gcbias <- gcbias
    out$cnlr <- cnlr

    out

}


doCNA<-function(zz,sampleid,minBinCount) {

    # weights
    # zwts <- log2(zz[,"rCountN"]+1-minBinCount)
    zwts <- log2(1/{1/(zz[,"rCountT"]+1) + 1/(zz[,"rCountN"]+1)})
    ii=which(zwts>0)

    zwts=zwts[ii]

    zchr <- zz[ii,"chrom"]
    zpos <- zz[ii,"maploc"]
    # log-ratio
    zlr <- zz[ii,"cnlr"]
	
    zcna <- CNA(zlr, zchr, zpos, sampleid=sampleid)
    zout <- segment(zcna, weights=zwts, undo.splits="sdundo",undo.SD=2)
    #zout <- segment(zcna)
    zout
}

