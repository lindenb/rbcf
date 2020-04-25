##Print the FILTERs in the VCF header
# load rbcf
library(rbcf)
# we don't need the index for this file
fp <- bcf.open("./data/gnomad.exomes.r2.0.1.sites.bcf",FALSE)
bcf.filters(fp)
# dispose the vcf reader
bcf.close(fp)

