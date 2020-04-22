library(rbcf)

test01<-function(filename,requireIndex) {
	fp <- bcf.open(filename,requireIndex)
	stopifnot(fp!=NULL)
	hdr <- bcf.hdr(fp)
	stopifnot(bcf.hdr.nsamples(hdr)==length(bcf.hdr.samples(hdr)))
	cat(bcf.hdr.samples(hdr),file=stderr())
	if(bcf.hdr.nsamples(hdr)>0) {
		stopifnot(match(bcf.hdr.sample1(hdr,1),bcf.hdr.samples(hdr))>0)
		}
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
		cat(" filters(",file=stderr())
		cat(variant.filters(hdr,vc),file=stderr())
		cat(") type:",file=stderr())
		cat(variant.types(hdr,vc),file=stderr())
		cat(" snp:",file=stderr())
		cat(variant.is.snp(hdr,vc),file=stderr())
		cat(" ploidy:",file=stderr())
		cat(variant.max.ploidy(hdr,vc),file=stderr())
		cat(" att.str:[",file=stderr())
		cat(variant.string.attributes(hdr,vc,"CSQ"),file=stderr())
		cat("] ",file=stderr())
		cat(" att.str:[",file=stderr())
		cat(variant.int.attributes(hdr,vc,"DP"),file=stderr())
		cat("] ",file=stderr())
		
		
		cat("\n",file=stderr())
		snidx<-1
		while(snidx<=bcf.hdr.nsamples(hdr)) {
			cat(bcf.hdr.sample1(hdr,snidx),file=stderr())
			cat(" alleles.idx[",file=stderr())
			cat(variant.gt.alleles.idx0(hdr,vc,snidx),file=stderr())
			cat("] ploidy:",file=stderr())
			cat(variant.gt.ploidy(hdr,vc,snidx),file=stderr())
			cat(" homref:",file=stderr())
			cat(variant.gt.homref(hdr,vc,snidx),file=stderr())
			cat(" het:",file=stderr())
			cat(variant.gt.het(hdr,vc,snidx),file=stderr())
			cat(" homvar:",file=stderr())
			cat(variant.gt.homvar(hdr,vc,snidx),file=stderr())
			cat(" hetnonref:",file=stderr())
			cat(variant.gt.hetnonref(hdr,vc,snidx),file=stderr())
			cat(" nocall:",file=stderr())
			cat(variant.gt.nocall(hdr,vc,snidx),file=stderr())
			cat(" phased:",file=stderr())
			cat(variant.gt.phased(hdr,vc,snidx),file=stderr())
			cat("\n",file=stderr())
			snidx<- snidx+1
			}
		cat("\n",file=stderr())
		}
	bcf.close(fp)
	}

test01("data/gnomad.exomes.r2.0.1.sites.vcf",FALSE)
test01("data/rotavirus_rf.01.vcf",FALSE)
test01("data/rotavirus_rf.02.vcf.gz",TRUE)
test01("data/rotavirus_rf.03.vcf.gz",TRUE)
test01("data/rotavirus_rf.04.bcf",TRUE)


