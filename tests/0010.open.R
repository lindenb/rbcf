##Open and close a VCF file
# load rbcf
library(rbcf)
# we don't need the index for this file
fp <- bcf.open("./data/rotavirus_rf.01.vcf",FALSE)
# dispose the vcf reader
bcf.close(fp)
print("Done.")

