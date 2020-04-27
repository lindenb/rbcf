##Print the Indexed Chromosomes
# load rbcf
library(rbcf)
# Open the indexed VCF
fp <- bcf.open("./data/rotavirus_rf.02.vcf.gz")
# error on opening (exit 0 for tests)
if(is.null(fp)) quit(save="no",status=0,runLast=FALSE)
# get the indexed contigs
bcf.contigs(fp)
# dispose the vcf reader
bcf.close(fp)

