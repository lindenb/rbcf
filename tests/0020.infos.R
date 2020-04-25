##Print the INFOs in the VCF header
# load rbcf
library(rbcf)
# we don't need the index for this file
fp <- bcf.open("./data/rotavirus_rf.01.vcf",FALSE)
bcf.infos(fp)
# dispose the vcf reader
bcf.close(fp)
# print the table
