##Print the FORMATs in the VCF header
# load rbcf
library(rbcf)
# we don't need the index for this file
fp <- bcf.open("./data/rotavirus_rf.01.vcf",FALSE)
# error on opening (exit 0 for tests)
if(is.null(fp)) quit(save="no",status=0,runLast=FALSE)
# print FORMAT
bcf.formats(fp)
# dispose the vcf reader
bcf.close(fp)
