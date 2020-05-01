##Writing variants to a new VCF/BCF file
# load rbcf
library(rbcf)
# vcf input filename
filenamein = "./data/rotavirus_rf.01.vcf"
# output vcf filename. "-" is standard output
filenameout =  "-"

fp <- bcf.open(filenamein,FALSE)
# error on opening (exit 0 for tests)
if(is.null(fp)) quit(save="no",status=0,runLast=FALSE)

# create a new VCF writer using the header from 'fp'
out <- bcf.new.writer(fp,filenameout)
# error on opening (exit 0 for tests)
if(is.null(out)) quit(save="no",status=0,runLast=FALSE)

# loop while we can read a variant
while(!is.null(vc<-bcf.next(fp))) {
	# only write POS%10==0
	if(variant.pos(vc)%%10==0) { 
		# write variant
		bcf.write.variant(out,vc);
		}
	}
# dispose the vcf reader
bcf.close(fp)
# dispose the vcf rwriter
bcf.close(out);
