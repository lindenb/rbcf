##Scanning the variants
# load rbcf
library(rbcf)

# create a function counting variants in a VCF
count.variants<-function(filename) {
	# we don't need the index for this file
	fp <- bcf.open(filename,FALSE)
	# number of variants
	n<-0
	# loop while we can read a variant
	while(!is.null(vc<-bcf.next(fp))) {
		# increment the count
		n<-n+1
	}
	# dispose the vcf reader
	bcf.close(fp)
	# return the number of variant
	n
}

# filenames
vcfs<-c(
	"./data/gnomad.exomes.r2.0.1.sites.bcf",
	"./data/rotavirus_rf.01.vcf",
	"./data/rotavirus_rf.02.vcf.gz",
	"./data/rotavirus_rf.03.vcf.gz",
	"./data/rotavirus_rf.04.bcf"
	)
# print the number of variants for each vcf
for(f in vcfs) {
	cat(paste(f," ",count.variants(f),"\n"))
	}
