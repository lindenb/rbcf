##Working with Genotypes
# load rbcf
library(rbcf)

# find given variant
find.variant<-function(fp,contig,pos) {
	if(!bcf.query(fp,paste(contig,":",pos,"-",pos,sep=""))) return(NULL)
	# loop while we can read a variant
	while(!is.null(vc<-bcf.next(fp))) {
		return(vc)
	}
	return(NULL)
}
filename<-"./data/1000G.ALL.2of4intersection.20100804.genotypes.bcf"
# open the VCF with index
fp <- bcf.open(filename)
# error on opening (exit 0 for tests)
if(is.null(fp)) quit(save="no",status=0,runLast=FALSE)
# find a variant
ctx <-find.variant(fp,"1",10583)
print(paste("Number of genotypes ",variant.nsamples(ctx)))

dp <- variant.genotypes.int.attribute(ctx, "DP")
stopifnot(length(dp) == variant.nsamples(ctx))
stopifnot(is.integer(dp))
cat("DP: ", paste(head(dp), collapse = ", "), "\n")

ad <- variant.genotypes.int.attribute(ctx, "AD")
stopifnot(length(ad) == 2*variant.nsamples(ctx))
stopifnot(is.integer(ad))
cat("AD: ", paste(ad, collapse = ", "), "\n")

gt <- variant.genotypes.allele.idx0(ctx)
stopifnot(length(gt) == 2*variant.nsamples(ctx))
stopifnot(is.integer(gt))
cat("GT - Integers: ", paste(gt, collapse = ", "), "\n")

gt <- variant.genotypes.allele.strings(ctx)
stopifnot(length(gt) == variant.nsamples(ctx))
stopifnot(is.character(gt))
cat("GT - Strings: ", paste(gt, collapse = ", "), "\n")



# dispose the vcf reader
bcf.close(fp)

