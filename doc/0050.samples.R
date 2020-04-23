##Print the Samples in the VCF header
# load rbcf
library(rbcf)
# we don't need the index for this file
fp <- bcf.open("../tests/data/rotavirus_rf.01.vcf",FALSE)
# print the number of samples
cat(paste("Num. Samples=",bcf.nsamples(fp),".\n"))
# get the name for the 1st sample
cat(paste("First sample is ",bcf.sample1(fp,1),".\n"))
# get the 1-based index for the samples
bcf.sample2index(fp,c("S1","S2","S3","undefined"))
# get the samples
samples <- bcf.samples(fp)
# dispose the vcf reader
bcf.close(fp)
# print the list
samples
