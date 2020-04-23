##Print the FILTERs in the VCF header
# load rbcf
library(rbcf)
# we don't need the index for this file
fp <- bcf.open("../tests/data/gnomad.exomes.r2.0.1.sites.vcf",FALSE)
flt <- bcf.filters(fp)
# dispose the vcf reader
bcf.close(fp)
# print the table
flt
