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
# get 10-th genotype
gt<-variant.genotype(ctx,10)
print(paste("sample ",genotype.sample(gt)))
# get genotype by name
gt<-variant.genotype(ctx,"NA18997")
print(paste("sample ",genotype.sample(gt)))
print(paste("alleles ",genotype.alleles.idx0(gt)))
print(paste("genotype ploidy ? ",genotype.ploidy(gt)))
print(paste("genotype is hom ref ? ",genotype.homref(gt)))
print(paste("genotype is het ? ",genotype.het(gt)))
print(paste("genotype is het-non-ref ? ",genotype.hetnonref(gt)))
print(paste("genotype is phased ? ",genotype.phased(gt)))
print(paste("genotype is no call ? ",genotype.nocall(gt)))
print(paste("genotype FORMAT/OG ? ",genotype.string.attribute(gt,"OG")))
print(paste("genotype FORMAT/GQ ? ",genotype.int.attribute(gt,"GQ")))# hum spec says gt should be integer
print(paste("genotype has GQ ? ",genotype.has.gq(gt)))
print(paste("genotype GQ ",genotype.gq(gt)))
print(paste("genotype has DP ? ",genotype.has.dp(gt)))
print(paste("genotype DP ",genotype.int.attribute(gt,"DP")))
print(paste("genotype DP ",genotype.dp(gt)))
print(paste("genotype has PL ? ",genotype.has.pl(gt)))
print(paste("genotype PL ",genotype.pl(gt)))
print(paste("genotype has AD ? ",genotype.has.ad(gt)))
print(paste("genotype AD ",genotype.ad(gt)))

# dispose the vcf reader
bcf.close(fp)

