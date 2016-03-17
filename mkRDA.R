args=commandArgs(trailing=T)
fname=args[1]
print(fname)
tag=gsub("_.*","",fname)
print(tag)

dd=read.delim(fname)

chromosomes=unique(dd$X.1_usercol)

gcpctdb=list()
for(ci in seq(len(chromosomes))){
    chr=chromosomes[ci]
    print(chr)
    ii=which(dd$X.1_usercol==chr)
    ij=which(dd$X10_num_N[ii]!=0 | dd$X12_seq_len[ii]<1000)
    pctgc=dd$X5_pct_gc[ii]
    pctgc[ij]=NA
    gcpctdb[[as.character(ci)]]=pctgc
}

save(gcpctdb,file=paste(tag,"gcpct.rda",sep=""),compress=T)

