##Print the INFOs in the VCF header
# load rbcf
library(rbcf)
# we don't need the index for this file
fp <- bcf.open("./data/rotavirus_rf.01.vcf",FALSE)
# error on opening (exit 0 for tests)
if(is.null(fp)) quit(save="no",status=0,runLast=FALSE)
# print INFO
bcf.infos(fp)
# dispose the vcf reader
bcf.close(fp)
# print the table
