##Print a SNPEFF table for a Variant
# load rbcf
library(rbcf)
# A vcf
filename <- "./data/rotavirus_rf.ann.vcf.gz"
# we don't need the index for this file
fp <- bcf.open(filename,FALSE)
# error on opening (exit 0 for tests)
if(is.null(fp)) quit(save="no",status=0,runLast=FALSE)
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
