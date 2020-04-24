##Print a SNPEFF table for a Variant
# load rbcf
library(rbcf)
# A vcf
filename <- "../tests/data/rotavirus_rf.ann.vcf.gz"
# we don't need the index for this file
fp <- bcf.open(filename,FALSE)
# current variant
vc <- NULL
while(!is.null(vc<-bcf.next(fp))) {
	#find the first variant having an INFO/ANN attribute
	if(variant.has.attribute(vc,"ANN")) break;
	}
if(!is.null(vc)) {
	# get SNPEFF table
	predictions<-variant.snpeff(vc)
	}
# dispose the vcf reader
bcf.close(fp)
# show
predictions
