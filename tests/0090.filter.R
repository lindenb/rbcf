##Scanning the variants
# load rbcf
library(rbcf)

# create a function counting variants in a VCF
count.variants<-function(filename,predicate) {
	# we don't need the index for this file
	fp <- bcf.open(filename,FALSE)
	# error on opening
	if(is.null(fp)) return(-1)
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
filename <- "./data/gnomad.exomes.r2.0.1.sites.bcf"
# filters
filters<-list(
	list("desc"="accept all","predicate"=function(ctx) {TRUE} ),
	list("desc"="accept none","predicate"=function(ctx) {FALSE} ),
	list("desc"="CHROM is '1'","predicate"=function(ctx) { variant.contig(ctx)=="1"} ),
	list("desc"="POS is even","predicate"=function(ctx) { (variant.pos(ctx)%%2)==1} ),
	list("desc"="PASS filter","predicate"=function(ctx) {!variant.is.filtered(ctx)} ),
	list("desc"="count(FILTER)>1","predicate"=function(ctx) {length(variant.filters(ctx))>1} ),
	list("desc"="FILTER contains SEGDUP","predicate"=function(ctx) {variant.has.filter(ctx,"SEGDUP")} ),
	list("desc"="SNP","predicate"=function(ctx) {variant.is.snp(ctx)} ),
	list("desc"="POS!=END","predicate"=function(ctx) { variant.pos(ctx)!=variant.end(ctx)} ),
	list("desc"="not diallelic","predicate"=function(ctx) {variant.nalleles(ctx)!=2} ),
	list("desc"="REF is 'A'","predicate"=function(ctx) {variant.reference(ctx)=="A"} ),
	list("desc"="any allele is 'A'","predicate"=function(ctx) {"A"  %in% variant.alleles(ctx)} ),
	list("desc"="any ALT allele is 'A'","predicate"=function(ctx) {"A"  %in% variant.alt.alleles(ctx)} ),
	list("desc"="No QUAL","predicate"=function(ctx) {!variant.has.qual(ctx)} ),
	list("desc"="variant has ID","predicate"=function(ctx) {variant.has.id(ctx)}),
	list("desc"="variant ID match 'rs1*' ","predicate"=function(ctx) {grepl("^rs1",variant.id(ctx))}),
	list("desc"="variant has INFO/AF_NFE","predicate"=function(ctx) {variant.has.attribute(ctx,"AF_NFE")}),
	list("desc"="variant has INFO/AF_NFE > 1E-5","predicate"=function(ctx) {variant.has.attribute(ctx,"AF_NFE") && length(which(variant.float.attribute(ctx,"AF_NFE") > 1E-5))>0}),
	list("desc"="Missense in PLEKHN1 (VEP)","predicate"=function(ctx) {
		# NO VEP annotation ?
		if(!variant.has.attribute(ctx,"CSQ")) return(FALSE);
		# get VEP annotation
		predictions <- variant.vep(ctx)
		# In SCN5A
		predictions <- predictions[which(predictions$SYMBOL=="PLEKHN1"),]
		# Consequence must contain missense
		predictions <- predictions[grep("missense_variant",predictions$Consequence),]
		nrow(predictions)>0
		})
	)

# count the variant for each filter
for(flt in filters) {
	print(paste(basename(filename)," filter:",flt[["desc"]]," count:",count.variants(filename,flt[["predicate"]]),"\n"))
	}

	
