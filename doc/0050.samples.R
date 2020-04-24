##Print the Samples in the VCF header
#' The samples are defined in the '\#CHROM' line of the VCF
# load rbcf
library(rbcf)
# we don't need the index for this file
fp <- bcf.open("../tests/data/rotavirus_rf.01.vcf",FALSE)
# print the number of samples
bcf.nsamples(fp)
# get the name for the 1st sample
bcf.sample.at(fp,1)
# get the 1-based index for the samples
bcf.sample2index(fp,c("S1","S2","S3","missing"))
# get all the samples
bcf.samples(fp)
# dispose the vcf reader
bcf.close(fp)
