##Print the Indexed Chromosomes
# load rbcf
library(rbcf)
# Open the indexed VCF
fp <- bcf.open("./data/rotavirus_rf.02.vcf.gz")
# get the indexed contigs
bcf.contigs(fp)
# dispose the vcf reader
bcf.close(fp)

