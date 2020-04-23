##Print the FORMATs in the VCF header
# load rbcf
library(rbcf)
# we don't need the index for this file
fp <- bcf.open("../tests/data/rotavirus_rf.01.vcf",FALSE)
fmts <- bcf.formats(fp)
# dispose the vcf reader
bcf.close(fp)
# print the table
fmts
