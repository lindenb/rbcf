#' Open a VCF or a BCF file
#' 
#' @param filename the path to the vcf file
#' @param requireIndex load associated vcf index
#' @return the new VCF reader
#' @examples
#' 
#' bcf.open("my.vcf.gz",T)
#' bcf.open("my.bcf")
bcf.open<-function(filename,requireIndex=TRUE)
	{
	stopifnot(is.character(filename))
	stopifnot(is.logical(requireIndex))
	.Call("RBcfFileOpen",filename,requireIndex);
	}

#' Close a VCF reader
#' 
#' @param fp the vcf reader
#' @return true on success
#' @examples
#' 
#' fp<-bcf.open("my.vcf.gz",T)
#' bcf.close(fp)
bcf.close<-function(fp)
	{
	.Call("RBcfFileClose",fp);
	}

	
bcf.nsamples<-function(fp)
	{
	.Call("RBcfNSamples",fp);
	}

bcf.samples<-function(fp)
	{
	.Call("RBcfSamples",fp);
	}

bcf.chromosomes<-function(fp)
	{
	.Call("RBcfSeqNames",fp)
	}

bcf.contigs<-function(fp)
	{
	bcf.chromosomes(fp)
	}

bcf.sample1<-function(fp,idx)
	{
	stopifnot(idx>0)
	.Call("RBcfSampleAtIndex0",fp,idx-1)
	}

#' Close a VCF reader
#' 
#' @param fp the vcf reader
#' @return true on success
#' @examples
#' 
#' fp<-bcf.open("my.vcf.gz",T)
#' bcf.close(fp)
bcf.close<-function(fp)
	{
	.Call("RBcfFileClose",fp);
	}

#' get dictionary from a VCF reader.
#' 
#' @param fp the vcf reader
#' @return the vcf dictionary as a table
#' @examples
#' 
#' fp<-bcf.open("my.vcf.gz",T)
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
bcf.dictionary<-function(fp)
	{
	.Call("RBcfHeaderDict",fp);
	}

#' prepare the VCF reader for a new vcf iteration over a given interval.
#' VCF reader must be associated and opened with a valid index.
#' 
#' @param fp the vcf reader
#' @param interval the genomic interval to be scanned
#' @return TRUE on success
#' @examples
bcf.query<-function(fp,interval) {
	stopifnot(is.character(interval))
	.Call("RBcfQueryRegion",fp,interval);
	}
	
#' read the next Variant Context in the VCF reader
#' 
#' @param fp the vcf reader
#' @return an opaque vcf context or NULL 
#' @examples
#' fp <- bcf.open("in.vcf.gz")
#' bcf.query(fp,"RF02:1-1000")
#' while(!is.null(vc<-bcf.next(fp))) {
#'      ## do something with vc
#'      }
#' bcf.close(fp)
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


#' return the genotype for a variant
#' 
#' @param vc the variant
#' @param sn the sample name (slower) or the 1-based sample index 
#' @return max ploidy
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


#' @param gt the genotype
#' @return the number of alleles for the genotypes
#
genotype.ploidy <-function(gt) {
	alleles<-genotype.alleles.idx0(gt)
	v = 0;
	if(!is.null(alleles)) v= length(alleles)
	v
	}

#' @param gt the genotype
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
#
genotype.nocall  <-function(gt) {
	alleles<-genotype.alleles.idx0(gt)
	!is.null(alleles) || length(alleles)==0 || match(-1,alleles)<=0
	}

#' @param gt the genotype
#' @return TRUE if genotype is phased
#
genotype.phased  <-function(gt) {
	.Call("RBcfCtxVariantGtPhased",gt);
	}



#' @param vc the variant
#' @param att the INFO/Attribute
#' @return the the INFO attribute for the given key
variant.string.attributes <-function(vc,att) {
	stopifnot(is.character(att))
	s <- .Call("RBcfCtxVariantAttributeAsString",vc,att)
	if( is.null(s) || length(s)==0) {
		c()
	} else {
	     unlist(strsplit(s, split=","))
	     }
	}


#' @param vc the variant
#' @param att the INFO/Attribute
#' @return the the INFO attribute for the given key
variant.int.attributes <-function(vc,att) {
	stopifnot(is.character(att))
	.Call("RBcfCtxVariantAttributeAsInt32",vc,att)
	}

#' @param vc the variant
#' @param att the INFO/Attribute
#' @return the the INFO attribute for the given key
variant.float.attributes <-function(vc,att) {
	stopifnot(is.character(att))
	.Call("RBcfCtxVariantAttributeAsFloat",vc,att)
	}
