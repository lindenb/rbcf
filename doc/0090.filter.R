##Scanning the variants
# load rbcf
library(rbcf)

# create a function counting variants in a VCF
count.variants<-function(filename,predicate) {
	# we don't need the index for this file
	fp <- bcf.open(filename,FALSE)
	# number of variants
	n<-0
	# loop while we can read a variant
	while(!is.null(vc<-bcf.next(fp))) {
		# test the variant
		if(predicate(vc)) {
			# increment the count
			n<-n+1
			}
	}
	# dispose the vcf reader
	bcf.close(fp)
	# return the number of variant
	n
}

# A vcf
filename <- "../tests/data/gnomad.exomes.r2.0.1.sites.vcf"
# filters
filters<-list(
	list("desc"="accept all","predicate"=function(ctx) {TRUE} ),
	list("desc"="accept none","predicate"=function(ctx) {FALSE} ),
	list("desc"="CHROM is '1'","predicate"=function(ctx) { variant.contig(ctx)=="1"} ),
	list("desc"="POS is even","predicate"=function(ctx) { (variant.pos(ctx)%%2)==1} ),
	list("desc"="PASS filter","predicate"=function(ctx) {!variant.is.filtered(ctx)} ),
	list("desc"="FILTER contains SEGDUP","predicate"=function(ctx) {variant.has.filter(ctx,"SEGDUP")} ),
	list("desc"="SNP","predicate"=function(ctx) {variant.is.snp(ctx)} ),
	list("desc"="not diallelic","predicate"=function(ctx) {variant.nalleles(ctx)!=2} ),
	list("desc"="REF is 'A'","predicate"=function(ctx) {variant.reference(ctx)=="A"} )
	)

# count the variant for each filter
for(flt in filters) {
	cat(paste(flt[["desc"]]," ",count.variants(filename,flt[["predicate"]]),"\n"))
	}

	
