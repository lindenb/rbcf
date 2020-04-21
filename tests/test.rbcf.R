library(rbcf)

test01<-function(filename,hasIndex) {
	fp <- bcf.open(filename)
	stopifnot(fp!=NULL)
	hdr <- bcf.hdr(fp)
	stopifnot(bcf.hdr.nsamples(hdr)==5)
	stopifnot(bcf.hdr.nsamples(hdr)==length(bcf.hdr.samples(hdr)))
	stopifnot(match("S4",bcf.hdr.samples(hdr))>0)
	bcf.close(fp)
	}
	
test01("data/rotavirus_rf.01.vcf",FALSE)
test01("data/rotavirus_rf.02.vcf.gz",TRUE)
test01("data/rotavirus_rf.03.vcf.gz",TRUE)
test01("data/rotavirus_rf.04.bcf",TRUE)
