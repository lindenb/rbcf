##Open and close a VCF file
# load rbcf
library(rbcf)
# we don't need the index for this file
fp <- bcf.open("./data/rotavirus_rf.01.vcf",FALSE)
# error (exit 0 for tests)
if(is.null(fp)) quit(save="no",status=0,runLast=FALSE)
# dispose the vcf reader
bcf.close(fp)
print("Done.")

