##Attribute in INFO
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
filename<-"./data/gnomad.exomes.r2.0.1.sites.bcf"
# open the VCF with index
fp <- bcf.open(filename)
# error on opening (exit 0 for tests)
if(is.null(fp)) quit(save="no",status=0,runLast=FALSE)

ctx <-find.variant(fp,"1",905608)
stopifnot(variant.has.attribute(ctx,"CSQ"))
print(paste("CSQ(no split) ",variant.string.attribute(ctx,"CSQ",split=FALSE)))
print(paste("CSQ(split) ",variant.string.attribute(ctx,"CSQ")))
stopifnot(variant.has.attribute(ctx,"AN_POPMAX"))
print(paste("AN_POPMAX:",variant.int.attribute(ctx,"AN_POPMAX")))
stopifnot(variant.has.attribute(ctx,"AF_POPMAX"))
print(paste("AF_POPMAX:",variant.float.attribute(ctx,"AF_POPMAX")))
print(paste("flag:VQSR_NEGATIVE_TRAIN_SITE:",variant.flag.attribute(ctx,"VQSR_NEGATIVE_TRAIN_SITE")))
# dispose the vcf reader
bcf.close(fp)

