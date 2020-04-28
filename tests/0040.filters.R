##Print the FILTERs in the VCF header
# load rbcf
library(rbcf)
# we don't need the index for this file
fp <- bcf.open("./data/gnomad.exomes.r2.0.1.sites.bcf",FALSE)
# error on opening (exit 0 for tests)
if(is.null(fp)) quit(save="no",status=0,runLast=FALSE)
# print FILTERs
bcf.filters(fp)
# dispose the vcf reader
bcf.close(fp)

