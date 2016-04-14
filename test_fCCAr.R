# test_fCCAr.R
source("fCCAr.R")

bigwigs <- c("H3K4me3_act_1.bw", "H3K4me3_act_2.bw", "H3K4me3_act_3.bw" )  
labels <- c( "H3K4me3","H3K4me3","H3K4me3" )	
fCCAr(bar=NULL, main="H3K4me3 peaks", peaks="noY_chr_merged_act_K4.bed", bigwigs=bigwigs, labels=labels ,  splines=15, nbins=100, ncan=15)


