#' @return the version of htslib
htslib.version <-function()  {
	.Call("HtslibGetVersion");
	}

#' @return the version of rcbf
rcbf.version <-function()  {
        .Call("RBCFGetVersion");
        }

#' Open a VCF or a BCF file
#' 
#' fp<-bcf.open("my.vcf.gz",TRUE)
#' hts.close(fp)
#' fp<-bcf.open("my.bcf")
#' hts.close(fp)
#' 
#' @param filename the path to the vcf file
#' @param requireIndex load associated vcf index
#' @return the new VCF reader
bcf.open<-function(filename,requireIndex=TRUE)
	{
	stopifnot(is.character(filename))
	stopifnot(is.logical(requireIndex))
	.Call("RBcfFileOpen",filename,requireIndex);
	}

#' Close a VCF reader
#' 
#' fp<-bcf.open("my.vcf.gz",TRUE)
#' bcf.close(fp)
#' 
#' @param fp the vcf reader
#' @return true on success
bcf.close<-function(fp)
	{
	if(!is.null(fp)) .Call("RBcfFileClose",fp);
	}


#' Open a new VCF writer.
#' Must be closed with bcf.close
#' 
#' @param fp the vcf reader
#' @param fname the name of the output vcf file
#' @return the writer or NULL on failure
bcf.new.writer<-function(fp,fname="-")
	{
	.Call("RBcfNewWriter",fp,fname);
	}

#' Save a Variant in a VCF writer
#' 
#' @param fp the vcf reader
#' @param vc the variant
#' @return true on success
bcf.write.variant<-function(fp,vc)
	{
	.Call("RBcfFileWriteCtx",fp,vc);
	}

#' @param fp the vcf reader
#' @return the number of samples
bcf.nsamples<-function(fp)
	{
	.Call("RBcfNSamples",fp);
	}

#' @param fp the vcf reader
#' @return samples
bcf.samples<-function(fp)
	{
	.Call("RBcfSamples",fp);
	}

#'
#' @param fp the vcf reader
#' @param sn the sample name
bcf.sample2index<-function(fp,sn) {
	sapply(sn,function(S) {
		.Call("BcfConvertSampleToIndex0",fp,S) + 1 
		})
	}
#' list the indexed chromosomes 
#' @param fp the vcf reader
#' @return a list of chromosome
bcf.chromosomes<-function(fp)
	{
	.Call("RBcfSeqNames",fp)
	}

#' alias of bcf.chromosomes
#' @param fp the vcf reader
#' @return a list of chromosome
bcf.contigs<-function(fp)
	{
	bcf.chromosomes(fp)
	}

#' @param fp the vcf reader
#' @param idx the 1-based index
#' @return the idx-th sample (1-based)
bcf.sample.at<-function(fp,idx)
	{
	stopifnot(idx>0)
	.Call("RBcfSampleAtIndex0",fp,idx-1)
	}

#' Close a VCF reader
#' 
#' fp<-bcf.open("my.vcf.gz",TRUE)
#' bcf.close(fp)
#'
#' @param fp the vcf reader
#' @return true on success
bcf.close<-function(fp)
	{
	.Call("RBcfFileClose",fp);
	}

#' get dictionary from a VCF reader.
#'
#' 
#' fp<-bcf.open("my.vcf.gz",TRUE)
#' dict<-bcf.dictionary(fp)
#' print(dict)
#'      chrom size
#' RF01  RF01 3302
#' RF02  RF02 2687
#' RF03  RF03 2592
#' RF04  RF04 2362
#' RF05  RF05 1579
#' RF06  RF06 1356
#' RF07  RF07 1074
#' RF08  RF08 1059
#' RF09  RF09 1062
#' RF10  RF10  751
#' RF11  RF11  666
#' 
#' @param fp the vcf reader
#' @return the vcf dictionary as a table
bcf.dictionary<-function(fp)
	{
	.Call("RBcfHeaderDict",fp);
	}


#' prepare the VCF reader for a new vcf iteration over a given interval.
#' VCF reader must be associated and opened with a valid index.
#' 
#' @param fp the vcf reader
#' @param interval the genomic interval to be scanned
#' @param collect if TRUE return a list of variants in the region or NULL on failure
#' @return TRUE on success or a list of variant if collect=TRUE
#' 
bcf.query<-function(fp,interval,collect=FALSE) {
	stopifnot(is.character(interval))
	ret<-.Call("RBcfQueryRegion",fp,interval);
	
	if(collect) {
		if(!ret) return(NULL);
		# initialize
		variants <-list()
		n<- 1
		# loop while we can read a variant
			while(!is.null(vc<-bcf.next(fp))) {
				#append
				variants[[n]] = vc
				n <- n + 1
				}
		return(variants)
		}
	ret
	}
	
#' read the next Variant Context in the VCF reader
#' 
#' 
#' fp <- bcf.open("in.vcf.gz")
#' bcf.query(fp,"RF02:1-1000")
#' while(!is.null(vc<-bcf.next(fp))) {
#'      ## do something with vc
#'      }
#' bcf.close(fp)
#'
#' @param fp the vcf reader
#' @return an opaque vcf context or NULL 
#
bcf.next<-function(fp) {
	.Call("RBcfNextLine",fp);
	}

#' get the numeric index (tid) of the chromosome for this variant
#' 
#' @param vc the variant
#' @return the numeric index for this variant
#
variant.tid<-function(vc) {
	.Call("RBcfCtxRid",vc);
	}
#' get chromosome name for this variant
#' 
#' @param vc the variant
#' @return the chromosome name
#
variant.contig<-function(vc) {
	.Call("RBcfCtxSeqName",vc);
	}
#' alias of variant.contig
#' 
#' @param vc the variant
#' @return the chromosome name
variant.chrom<-function(vc) {
	variant.contig(vc)
	}
#' get chromosome POS for this variant
#' 
#' @param vc the variant
#' @return the starting position for this variant
#
variant.pos<-function(vc) {
	.Call("RBcfCtxPos",vc);
	}
#' alias of variant.pos
#' 
#' @param vc the variant
#' @return the starting position for this variant
variant.start<-function(vc) {
	variant.pos(vc)
	}

#' return whether variant have got an ID
#' 
#' @param vc the variant
#' @return the TRUE if variant have got an ID
variant.has.id<-function(vc) {
	.Call("RBcfCtxHasId",vc);
	}
	
#' return the variant ID
#' 
#' @param vc the variant
#' @return the variant ID or NULL
variant.id<-function(vc) {
	.Call("RBcfCtxId",vc);
	}

#' get chromosome END for this variant
#' 
#' @param vc the variant
#' @return the END position for this variant
#
variant.end<-function(vc) {
	.Call("RBcfCtxEnd",vc);
	}

#' alias if variant.end
#' 
#' @param vc the variant
#' @return the END position for this variant
#
variant.stop<-function(vc) {
	variant.end(vc);
	}


#' get the number of alleles for this variant
#' 
#' @param vc the variant
#' @return the number of alleles for this variant
#
variant.nalleles<-function(vc) {
	.Call("RBcfCtxNAlleles",vc);
	}

#' get the alleles for this variant
#' 
#' @param vc the variant
#' @return the alleles for this variant
#
variant.alleles<-function(vc) {
	.Call("RBcfCtxAlleles",vc);
	}

#' get the reference allele for this variant
#' 
#' @param vc the variant
#' @return the reference allele for this variant
#
variant.reference<-function(vc) {
	.Call("RBcfCtxReference",vc);
	}

#' get the alternate alleles for this variant
#' 
#' @param vc the variant
#' @return the alternate alleles for this variant
#
variant.alt.alleles<-function(vc) {
	.Call("RBcfCtxAlternateAlleles",vc);
	}


#' test if QUAL is available for this variant
#' 
#' @param vc the variant
#' @return QUAL is available for this variant
#
variant.has.qual<-function(vc) {
	.Call("RBcfCtxHasQual",vc);
	}

#' test the QUAL for this variant
#' 
#' @param vc the variant
#' @return QUAL or NULL
#
variant.qual<-function(vc) {
	.Call("RBcfCtxQual",vc);
	}
	

#' test is variant is filtered
#' 
#' @param vc the variant
#' @return  variant is filtered
#
variant.is.filtered<-function(vc) {
	.Call("RBcfCtxFiltered",vc);
	}

#' return the FILTERS for this variant
#' 
#' @param vc the variant
#' @return the filters
#
variant.filters<-function(vc) {
	.Call("RBcfCtxFilters",vc);
	}


#' @param vc the variant
#' @param fn filter name
#' @return true if variant is filtered with 'fn'
#
variant.has.filter<-function(vc,fn) {
	.Call("VariantHasFilter",vc,fn);
	}


#' return the types for this variant (as defined in htslib)
#' 
#' @param vc the variant
#' @return the types for this variant
#
variant.types<-function(vc) {
	.Call("RBcfCtxVariantTypes",vc);
	}


#' return true if variant is a SNP
#' 
#' @param vc the variant
#' @return true if the variant is a SNP
#
variant.is.snp<-function(vc) {
	.Call("RBcfCtxVariantIsSnp",vc);
	}

#' max ploidy for this variant
#' 
#' @param vc the variant
#' @return max ploidy
#
variant.max.ploidy<-function(vc) {
	.Call("RBcfCtxVariantMaxPloidy",vc);
	}

#' @param vc the variant
#' @return the number of samples/genotypes for this variant
variant.nsamples <-function(vc) {
	.Call("VariantNSamples",vc);
	}

#' @param vc the variant
#' @return the number of samples/genotypes for this variant
variant.genotypes <-function(vc) {
	sapply(1:variant.nsamples(vc),function(idx) {
		variant.genotype(vc,idx)
		})
        }


#' return the genotype for a variant
#' 
#' @param vc the variant
#' @param nameOrIdx the sample name (slower) or the 1-based sample index 
#' @return the given genotype
#
variant.genotype<-function(vc,nameOrIdx) {
	if(is.numeric(nameOrIdx)) {
		stopifnot(nameOrIdx>0)
		nameOrIdx<-as.integer(nameOrIdx)-1
		}
	else {
		stopifnot(is.character(nameOrIdx))
		}
	.Call("VariantGetGenotype",vc,nameOrIdx);
	}

#' return the 0-based alleles indexes for the given genotype
#' 
#' @param gt the genotype
#' @return max ploidy
#
genotype.alleles.idx0 <-function(gt) {
	.Call("RBcfCtxVariantGtAllelesIndexes0",gt);
	}

variant.genotypes.alleles.idx0 <-function(vc) {
	.Call("RBcfCtxVariantAllGtAllelesIndexes0", vc);
	}

variant.genotypes.allele.strings <-function(vc) {
	.Call("RBcfCtxVariantAllGtStrings", vc);
	}

#' @param gt the genotype
#' @return the number of alleles for the genotypes
#
genotype.ploidy <-function(gt) {
	alleles<-genotype.alleles.idx0(gt)
	v = 0;
	if(!is.null(alleles)) v= length(alleles)
	v
	}

#' @param gt the g==enotype
#' @return TRUE if genotype is diploid and all alleles are reference
#
genotype.homref <-function(gt) {
	alleles<-genotype.alleles.idx0(gt)
	!is.null(alleles) && length(alleles)==2 && alleles[1]==0 && alleles[2]==0 
	}

#' @param gt the genotype
#' @return TRUE if genotype is diploid and heterozygous
#
genotype.het <-function(gt) {
	alleles<-genotype.alleles.idx0(gt)
	!is.null(alleles) && length(alleles)==2 && alleles[1]!=alleles[2] && alleles[1]>=0 && alleles[2]>=0
	}

#' @param gt the genotype
#' @return TRUE if genotype is diploid and homozygous on alt allele
#
genotype.homvar <-function(gt) {
	alleles<-genotype.alleles.idx0(gt)
	!is.null(alleles) && length(alleles)==2 && alleles[1]>0 && alleles[2]>0 && alleles[1]==alleles[2]
	}
	
#' @param gt the genotype
#' @return TRUE if genotype is diploid and heterozygous and doesn't contains the alt allele
#
genotype.hetnonref  <-function(gt) {
	alleles<-genotype.alleles.idx0(gt)
	!is.null(alleles) && length(alleles)==2 && alleles[1]!=alleles[2] && alleles[1]>0 && alleles[2]>0
	}

#' @param gt the genotype
#' @return TRUE if genotypes contains no allele '' or any is no call '.'
genotype.nocall  <-function(gt) {
	alleles<-genotype.alleles.idx0(gt)
	is.null(alleles) || length(alleles)==0 || -1 %in% alleles
	}

#' @param gt the genotype
#' @return TRUE if genotype is phased
#
genotype.phased  <-function(gt) {
	.Call("RBcfCtxVariantGtPhased",gt);
	}

#' @param gt the genotype
#' @return the name of the sample associated to the genotype
genotype.sample  <-function(gt) {
	.Call("GenotypeSample",gt)
	}


#' @param vc the variant
#' @param att the INFO/Attribute
#' @return true if variant has INFO attribute
variant.has.attribute <-function(vc,att) {
	stopifnot(is.character(att))
	.Call("VariantHasAttribute",vc,att)
	}


#' @param vc the variant
#' @param att the INFO/Attribute
#' @param split split strings using commas
#' @return the the INFO attribute for the given key.
variant.string.attribute <-function(vc,att, split = TRUE) {
	stopifnot(is.character(att))
	s <-.Call("VariantStringAttribute",vc,att)
	if( !is.null(s) && length(s)>0 && split) {
		s <- unlist(strsplit(s, split=","))
		}
	s
	}

#' @param vc the variant
#' @param att the INFO/Attribute
#' @return the the INFO attribute for the given key.
variant.int.attribute <-function(vc,att) {
	.Call("VariantIntAttribute",vc,att)
	}

#' @param vc the variant
#' @param att the INFO/Attribute
#' @return the the INFO attribute for the given key.
variant.float.attribute <-function(vc,att) {
	.Call("VariantFloatAttribute",vc,att)
	}

#' @param vc the variant
#' @param att the INFO/Attribute
#' @return the the INFO attribute for the given key.
variant.flag.attribute <-function(vc,att) {
	.Call("VariantFlagAttribute",vc,att)
	}


#' @param fp the vcf reader
#' @return a table of filters
bcf.filters <-function(fp) {
	 .Call("BcfFilterTable",fp)
}

#' @param fp the vcf reader
#' @return a table of filters
bcf.infos <-function(fp) {
	 .Call("BcfInfoTable",fp)
}

#' @param fp the vcf reader
#' @return a table of filters
bcf.formats <-function(fp) {
	 .Call("BcfFormatTable",fp)
	}

#' @param vc the variant
#' @return the list INFOs for this variant
variant.info.ids <-function(vc) {
	.Call("VariantGetInfoKeySet",vc)
	}

#' @param vc the variant
#' @return the list FORMATs for this variant
variant.format.ids <-function(vc) {
	.Call("VariantGetFormatKeySet",vc)
	}
#' @param vc the variant
#' @return VEP table for this variant or NULL
variant.vep <-function(vc) {
	.Call("VariantVepTable",vc)
	}
#' @param vc the variant
#' @return SNPEFF table for this variant or NULL
variant.snpeff <-function(vc) {
	.Call("VariantSnpEffTable",vc)
	}

#' @param vc the variant
#' @return vector of numerics containing attribute values for all genotypes
variant.genotypes.flag.attribute <- function(vc, att) {
	.Call("VariantGenotypesFlagAttribute", vc, att)
}

#' @param vc the variant
#' @return vector of numerics containing attribute values for all genotypes
variant.genotypes.int.attribute <- function(vc, att) {
	.Call("VariantGenotypesInt32Attribute", vc, att)
}

#' @param vc the variant
#' @return vector of numerics containing attribute values for all genotypes
variant.genotypes.float.attribute <- function(vc, att) {
	.Call("VariantGenotypesFloatAttribute", vc, att)
}

#' @param gt the genotype
#' @param att the key
#' @return the values for this key+genotype
#
genotype.int.attribute <-function(gt,att) {
	.Call("GenotypeInt32Attribute",gt,att)
	}


#' @param gt the genotype
#' @param att the key
#' @return a String for this flag+genotype or NULL
#
genotype.string.attribute <-function(gt,att) {
	.Call("GenotypeStringAttribute",gt,att)
	}

#' @param gt the genotype
#' @return TRUE if genotype is filtered (FORMAT/FT is set)
#
genotype.filtered <-function(gt) {
	flt<- genotype.string.attribute(gt,"FT")
	!(is.null(flt) || length(flt)==0 || flt==".")
	}

#' @param gt the genotype
#' @return the DP or -1
genotype.dp <-function(gt) {
	v<-genotype.int.attribute(gt,"DP")
	if(length(v)!=1) return(-1)
	v[1]	
	}
#' @param gt the genotype
#' @return TRUE if DP is available 
genotype.has.dp <-function(gt) {
	genotype.dp(gt)>=0	
	}
#' @param gt the genotype
#' @return the GQ
genotype.gq <-function(gt) {
	v<-genotype.int.attribute(gt,"GQ")
	if(length(v)!=1) return(-1)
	v[1]	
	}
#' @param gt the genotype
#' @return TRUE if GQ is available 
genotype.has.gq <-function(gt) {
	genotype.gq(gt)>=0	
	}
#' @param gt the genotype
#' @return PL or empty vector
genotype.pl <-function(gt) {
	genotype.int.attribute(gt,"PL")
	}
#' @param gt the genotype
#' @return TRUE if PL is available 
genotype.has.pl<-function(gt) {
	length(genotype.pl(gt))>0	
	}

#' @param gt the genotype
#' @return AD or rempty vector
genotype.ad <-function(gt) {
	genotype.int.attribute(gt,"AD")
	}

#' @param gt the genotype
#' @return TRUE if AD is available 
genotype.has.ad <-function(gt) {
	length(genotype.ad(gt))>0	
	}


#' @param gt the genotype
#' @param att the key
#' @return the values for this key+genotype
#
genotype.float.attribute <-function(gt,att) {
	.Call("GenotypeFloatAttribute",gt,att)
	}
