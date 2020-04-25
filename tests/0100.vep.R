##Print a VEP table for a Variant
# load rbcf
library(rbcf)
# A vcf
filename <- "./data/gnomad.exomes.r2.0.1.sites.bcf"
# we don't need the index for this file
fp <- bcf.open(filename,FALSE)
# current variant
vc <- NULL
while(!is.null(vc<-bcf.next(fp))) {
	#find the first variant having an INFO/CSQ attribute
	if(variant.has.attribute(vc,"CSQ")) break;
	}

if(!is.null(vc)) {
	# get the VEP table for the variant
	predictions<-variant.vep(vc)
	}

# dispose the vcf reader
bcf.close(fp)
# show
predictions
