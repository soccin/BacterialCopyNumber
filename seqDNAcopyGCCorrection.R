library(DNAcopy)
source("lib.R")
# load gc data
load("mycogcpct.rda")

####
# Edit these
#

COUNTFILE="counts.txt"
BAMFILES="bamList2"
BAMS=scan(BAMFILES,"")
normal="s_Msmeg_mutant1_gris0_S2"

#
# May need to fix this depending on how BAM files are named
#
samples=gsub("___MD.*","",gsub("out.*/s_","s_",BAMS))

######################################################################
######################################################################
######################################################################

counts=read.delim(COUNTFILE,header=F)
colnames(counts)=c("chrom","start","stop",samples)
counts$pos=(counts$start+counts$stop)/2e6
minBinCount=15

for(si in samples){
    if(si != normal){
        print(si)

        mat=data.frame(chrom=as.numeric(factor(counts$chrom)),
                        maploc=floor((counts$start+counts$stop)/2),
                        rCountN=counts[,normal],
                        rCountT=counts[,si])
        mat$keep=rep(1,nrow(mat))
        mat$keep[mat$rCountN<minBinCount]=0

        out=counts2logR(mat)
		
        zout=doCNA(out,paste(si,"vs",normal),minBinCount)
        plot(zout,xmaploc=T,ylim=c(-4,4))
        grid()
		dev.copy2pdf(file=cc("dnaCopy",si,"_vs_",normal,".pdf"),width=11,height=8)
        write.xls(zout$output,file=cc("dnaCopy",si,"_vs_",normal,".seg"))
    }
}

