library(rbcf)

test01<-function(filename,requireIndex) {
	fp <- bcf.open(filename,requireIndex)
	stopifnot(fp!=NULL)
	hdr <- bcf.hdr(fp)
	stopifnot(bcf.hdr.nsamples(hdr)==5)
	stopifnot(bcf.hdr.nsamples(hdr)==length(bcf.hdr.samples(hdr)))
	cat(bcf.hdr.samples(hdr),file=stderr())
	stopifnot(match("S4",bcf.hdr.samples(hdr))>0)
	stopifnot(bcf.hdr.sample1(hdr,1)=="S1")
	dict<-bcf.hdr.dictionary(hdr)
	head(dict)
	
	if(requireIndex) {
		contigs <-bcf.contigs(fp);
		stopifnot(bcf.query(fp,"RF02:1-1000"))
		}
	while(!is.null(vc<-bcf.next(fp))) {
		cat(variant.tid(hdr,vc),file=stderr())
		cat(" ",file=stderr())
		cat(variant.contig(hdr,vc),file=stderr())
		cat(" ",file=stderr())
		cat(variant.start(hdr,vc),file=stderr())
		cat(" ",file=stderr())
		cat(variant.stop(hdr,vc),file=stderr())
		cat(" ",file=stderr())
		cat(variant.has.id(hdr,vc),file=stderr())
		cat(" ",file=stderr())
		cat(variant.id(hdr,vc),file=stderr())
		cat(" n.alleles:",file=stderr())
		cat(variant.nalleles(hdr,vc),file=stderr())
		cat(" alleles:",file=stderr())
		cat(variant.alleles(hdr,vc),file=stderr())
		cat(" ref:",file=stderr())
		cat(variant.reference(hdr,vc),file=stderr())
		cat(" ",file=stderr())
		cat(variant.has.qual(hdr,vc),file=stderr())
		cat(" ",file=stderr())
		cat(variant.qual(hdr,vc),file=stderr())
		cat(" ",file=stderr())
		cat(variant.is.filtered(hdr,vc),file=stderr())
		cat(" ",file=stderr())
		cat(variant.filters(hdr,vc),file=stderr())
		cat(" ",file=stderr())
		cat(variant.types(hdr,vc),file=stderr())
		cat(" ",file=stderr())
		cat(variant.is.snp(hdr,vc),file=stderr())
		cat("\n",file=stderr())
		cat(variant.max.ploidy(hdr,vc),file=stderr())
		cat("\n",file=stderr())
		}
	bcf.close(fp)
	}
	
test01("data/rotavirus_rf.01.vcf",FALSE)
test01("data/rotavirus_rf.02.vcf.gz",TRUE)
test01("data/rotavirus_rf.03.vcf.gz",TRUE)
test01("data/rotavirus_rf.04.bcf",TRUE)

fp <- bcf.open("data/rotavirus_rf.03.vcf.gz")
stopifnot(fp!=NULL)
hdr <- bcf.hdr(fp)
dict<-bcf.hdr.dictionary(hdr)
print(bcf.contigs(fp))
print(typeof(dict))
print(dict)


bcf.close(fp)

