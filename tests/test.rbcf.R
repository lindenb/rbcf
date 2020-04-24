library(rbcf)

test01<-function(filename,requireIndex) {
	fp <- bcf.open(filename,requireIndex)
	stopifnot(fp!=NULL)
	stopifnot(bcf.nsamples(fp)==length(bcf.samples(fp)))
	cat(bcf.samples(fp),file=stderr())
	if(bcf.nsamples(fp)>0) {
		stopifnot(match(bcf.sample.at(fp,1),bcf.samples(fp))>0)
		}
	dict<-bcf.dictionary(fp)
	head(dict)
	
	if(requireIndex) {
		contigs <-bcf.contigs(fp);
		stopifnot(bcf.query(fp,"RF02:1-1000"))
		}
	while(!is.null(vc<-bcf.next(fp))) {
		cat(variant.tid(vc),file=stderr())
		cat(" ",file=stderr())
		cat(variant.contig(vc),file=stderr())
		cat(" ",file=stderr())
		cat(variant.start(vc),file=stderr())
		cat(" ",file=stderr())
		cat(variant.stop(vc),file=stderr())
		cat(" ",file=stderr())
		cat(variant.has.id(vc),file=stderr())
		cat(" ",file=stderr())
		cat(variant.id(vc),file=stderr())
		cat(" n.alleles:",file=stderr())
		cat(variant.nalleles(vc),file=stderr())
		cat(" alleles:",file=stderr())
		cat(variant.alleles(vc),file=stderr())
		cat(" ref:",file=stderr())
		cat(variant.reference(vc),file=stderr())
		cat(" ",file=stderr())
		cat(variant.has.qual(vc),file=stderr())
		cat(" ",file=stderr())
		cat(variant.qual(vc),file=stderr())
		cat(" ",file=stderr())
		cat(variant.is.filtered(vc),file=stderr())
		cat(" filters(",file=stderr())
		cat(variant.filters(vc),file=stderr())
		cat(") type:",file=stderr())
		cat(variant.types(vc),file=stderr())
		cat(" snp:",file=stderr())
		cat(variant.is.snp(vc),file=stderr())
		cat(" ploidy:",file=stderr())
		cat(variant.max.ploidy(vc),file=stderr())
		cat(" att.str:[",file=stderr())
		cat(variant.attribute(vc,"CSQ"),file=stderr())
		cat("] ",file=stderr())
		cat(" att.has.attr:",file=stderr())
		cat(variant.has.attribute(vc,"DP"),file=stderr())
		cat(" ",file=stderr())
		cat(" att.str:[",file=stderr())
		cat(variant.attribute(vc,"DP"),file=stderr())
		cat("] ",file=stderr())
		cat(" att.str.float:[",file=stderr())
		cat(variant.attribute(vc,"AF"),file=stderr())
		cat("] ",file=stderr())
		cat(" att.str.flag:[",file=stderr())
		cat(variant.attribute(vc,"X1"),file=stderr())
		cat("] ",file=stderr())
		cat("] info.ids:[",file=stderr())
		cat(variant.info.ids(vc),file=stderr())
		cat("] info.formats:[",file=stderr())
		cat(variant.format.ids(vc),file=stderr())
		cat("] ",file=stderr())
		
		
		
		cat("\n",file=stderr())
		snidx<-1
		while(snidx<=bcf.nsamples(fp)) {
			gt <- variant.genotype(vc,snidx)
			cat(bcf.sample.at(fp,snidx),file=stderr())
			cat(" alleles.idx[",file=stderr())
			cat(genotype.alleles.idx0(gt),file=stderr())
			cat("] ploidy:",file=stderr())
			cat(genotype.ploidy(gt),file=stderr())
			cat(" homref:",file=stderr())
			cat(genotype.homref(gt),file=stderr())
			cat(" het:",file=stderr())
			cat(genotype.het(gt),file=stderr())
			cat(" homvar:",file=stderr())
			cat(genotype.homvar(gt),file=stderr())
			cat(" hetnonref:",file=stderr())
			cat(genotype.hetnonref(gt),file=stderr())
			cat(" nocall:",file=stderr())
			cat(genotype.nocall(gt),file=stderr())
			cat(" phased:",file=stderr())
			cat(genotype.phased(gt),file=stderr())
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

fp <- bcf.open("data/rotavirus_rf.01.vcf",FALSE)
bcf.filters(fp);
bcf.infos(fp);
bcf.formats(fp)
bcf.close(fp);


