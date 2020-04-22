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

#' get the VCF header from a VCF reader.
#' this is a copy of the original header.
#' 
#' @param fp the vcf reader
#' @return the vcf header
#' @examples
#' 
#' fp<-bcf.open("my.vcf.gz",T)
#' hdr<-bcf.hdr(fp)
#' bcf.close(fp)
bcf.hdr<-function(fp)
	{
	.Call("RBcfHeader",fp);
	}
	
bcf.hdr.nsamples<-function(fp)
	{
	.Call("RBcfNSamples",fp);
	}

bcf.hdr.samples<-function(fp)
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

bcf.hdr.sample1<-function(hdr,idx)
	{
	stopifnot(idx>0)
	.Call("RBcfSampleAtIndex0",hdr,idx-1)
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
#' @param hdr the vcf header
#' @return the vcf dictionary as a table
#' @examples
#' 
#' fp<-bcf.open("my.vcf.gz",T)
#' hdr<-bcf.hdr(fp)
#' dict<-bcf.hdr.dictionary(hdr)
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
bcf.hdr.dictionary<-function(hdr)
	{
	.Call("RBcfHeaderDict",hdr);
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
#' @param hdr the vcf header
#' @param vcf the variant
#' @return the numeric index for this variant
#
variant.tid<-function(hdr,vc) {
	.Call("RBcfCtxRid",hdr,vc);
	}
#' get chromosome name for this variant
#' 
#' @param hdr the vcf header
#' @param vcf the variant
#' @return the chromosome name
#
variant.contig<-function(hdr,vc) {
	.Call("RBcfCtxSeqName",hdr,vc);
	}
#' alias of variant.contig
#' 
#' @param hdr the vcf header
#' @param vcf the variant
#' @return the chromosome name
variant.chrom<-function(hdr,vc) {
	variant.contig(vc)
	}
#' get chromosome POS for this variant
#' 
#' @param hdr the vcf header
#' @param vcf the variant
#' @return the starting position for this variant
#
variant.pos<-function(hdr,vc) {
	.Call("RBcfCtxPos",hdr,vc);
	}
#' alias of variant.pos
#' 
#' @param hdr the vcf header
#' @param vcf the variant
#' @return the starting position for this variant
variant.start<-function(hdr,vc) {
	variant.pos(hdr,vc)
	}

#' return whether variant have got an ID
#' 
#' @param hdr the vcf header
#' @param vcf the variant
#' @return the TRUE if variant have got an ID
variant.has.id<-function(hdr,vc) {
	.Call("RBcfCtxHasId",hdr,vc);
	}
	
#' return the variant ID
#' 
#' @param hdr the vcf header
#' @param vcf the variant
#' @return the variant ID or NULL
variant.id<-function(hdr,vc) {
	.Call("RBcfCtxId",hdr,vc);
	}

#' get chromosome END for this variant
#' 
#' @param hdr the vcf header
#' @param vcf the variant
#' @return the END position for this variant
#
variant.end<-function(hdr,vc) {
	.Call("RBcfCtxEnd",hdr,vc);
	}

#' alias if variant.end
#' 
#' @param hdr the vcf header
#' @param vcf the variant
#' @return the END position for this variant
#
variant.stop<-function(hdr,vc) {
	variant.end(hdr,vc);
	}


#' get the number of alleles for this variant
#' 
#' @param hdr the vcf header
#' @param vcf the variant
#' @return the number of alleles for this variant
#
variant.nalleles<-function(hdr,vc) {
	.Call("RBcfCtxNAlleles",hdr,vc);
	}

#' get the alleles for this variant
#' 
#' @param hdr the vcf header
#' @param vcf the variant
#' @return the alleles for this variant
#
variant.alleles<-function(hdr,vc) {
	.Call("RBcfCtxAlleles",hdr,vc);
	}

#' get the reference allele for this variant
#' 
#' @param hdr the vcf header
#' @param vcf the variant
#' @return the reference allele for this variant
#
variant.reference<-function(hdr,vc) {
	.Call("RBcfCtxReference",hdr,vc);
	}

#' get the alternate alleles for this variant
#' 
#' @param hdr the vcf header
#' @param vcf the variant
#' @return the alternate alleles for this variant
#
variant.alt.alleles<-function(hdr,vc) {
	.Call("RBcfCtxAlternateAlleles",hdr,vc);
	}


#' test if QUAL is available for this variant
#' 
#' @param hdr the vcf header
#' @param vcf the variant
#' @return QUAL is available for this variant
#
variant.has.qual<-function(hdr,vc) {
	.Call("RBcfCtxHasQual",hdr,vc);
	}

#' test the QUAL for this variant
#' 
#' @param hdr the vcf header
#' @param vcf the variant
#' @return QUAL or NULL
#
variant.qual<-function(hdr,vc) {
	.Call("RBcfCtxQual",hdr,vc);
	}
	

#' test is variant is filtered
#' 
#' @param hdr the vcf header
#' @param vcf the variant
#' @return  variant is filtered
#
variant.is.filtered<-function(hdr,vc) {
	.Call("RBcfCtxFiltered",hdr,vc);
	}

#' return the FILTERS for this variant
#' 
#' @param hdr the vcf header
#' @param vcf the variant
#' @return the filters
#
variant.filters<-function(hdr,vc) {
	.Call("RBcfCtxFilters",hdr,vc);
	}


#' return the types for this variant (as defined in htslib)
#' 
#' @param hdr the vcf header
#' @param vcf the variant
#' @return the types for this variant
#
variant.types<-function(hdr,vc) {
	.Call("RBcfCtxVariantTypes",hdr,vc);
	}


#' return true if variant is a SNP
#' 
#' @param hdr the vcf header
#' @param vcf the variant
#' @return true if the variant is a SNP
#
variant.is.snp<-function(hdr,vc) {
	.Call("RBcfCtxVariantIsSnp",hdr,vc);
	}

#' max ploidy for this variant
#' 
#' @param hdr the vcf header
#' @param vcf the variant
#' @return max ploidy
#
variant.max.ploidy<-function(hdr,vc) {
	.Call("RBcfCtxVariantMaxPloidy",hdr,vc);
	}

#' return the 0-based alleles indexes for the given genotype
#' 
#' @param hdr the vcf header
#' @param vcf the variant
#' @param sn the sample name (slower) or the 1-based sample index 
#' @return max ploidy
#
variant.gt.alleles.idx0 <-function(hdr,vc,nameOrIdx) {
	if(is.numeric(nameOrIdx)) {
		stopifnot(nameOrIdx>0)
		nameOrIdx<-as.integer(nameOrIdx)-1
		}
	else {
		stopifnot(is.character(nameOrIdx))
		}
	.Call("RBcfCtxVariantGtAllelesIndexes0",hdr,vc,nameOrIdx);
	}


#' @param hdr the vcf header
#' @param vcf the variant
#' @param sn the sample name or the 1-based sample index 
#' @return the number of alleles for the genotypes
#
variant.gt.ploidy <-function(hdr,vc,nameOrIdx) {
	alleles<-variant.gt.alleles.idx0(hdr,vc,nameOrIdx)
	v = 0;
	if(!is.null(alleles)) v= length(alleles)
	v
	}

#' @param hdr the vcf header
#' @param vcf the variant
#' @param sn the sample name or the 1-based sample index 
#' @return TRUE if genotype is diploid and all alleles are reference
#
variant.gt.homref <-function(hdr,vc,nameOrIdx) {
	alleles<-variant.gt.alleles.idx0(hdr,vc,nameOrIdx)
	!is.null(alleles) && length(alleles)==2 && alleles[1]==0 && alleles[2]==0 
	}

#' @param hdr the vcf header
#' @param vcf the variant
#' @param sn the sample name or the 1-based sample index 
#' @return TRUE if genotype is diploid and heterozygous
#
variant.gt.het <-function(hdr,vc,nameOrIdx) {
	alleles<-variant.gt.alleles.idx0(hdr,vc,nameOrIdx)
	!is.null(alleles) && length(alleles)==2 && alleles[1]!=alleles[2] && alleles[1]>=0 && alleles[2]>=0
	}

#' @param hdr the vcf header
#' @param vcf the variant
#' @param sn the sample name or the 1-based sample index 
#' @return TRUE if genotype is diploid and homozygous on alt allele
#
variant.gt.homvar <-function(hdr,vc,nameOrIdx) {
	alleles<-variant.gt.alleles.idx0(hdr,vc,nameOrIdx)
	!is.null(alleles) && length(alleles)==2 && alleles[1]>0 && alleles[2]>0 && alleles[1]==alleles[2]
	}
	
#' @param hdr the vcf header
#' @param vcf the variant
#' @param sn the sample name or the 1-based sample index 
#' @return TRUE if genotype is diploid and heterozygous and doesn't contains the alt allele
#
variant.gt.hetnonref  <-function(hdr,vc,nameOrIdx) {
	alleles<-variant.gt.alleles.idx0(hdr,vc,nameOrIdx)
	!is.null(alleles) && length(alleles)==2 && alleles[1]!=alleles[2] && alleles[1]>0 && alleles[2]>0
	}

#' @param hdr the vcf header
#' @param vcf the variant
#' @param sn the sample name or the 1-based sample index 
#' @return TRUE if genotypes contains no allele '' or any is no call '.'
#
variant.gt.nocall  <-function(hdr,vc,nameOrIdx) {
	alleles<-variant.gt.alleles.idx0(hdr,vc,nameOrIdx)
	!is.null(alleles) || length(alleles)==0 || match(-1,alleles)<=0
	}

#' @param hdr the vcf header
#' @param vcf the variant
#' @param sn the sample name or the 1-based sample index 
#' @return TRUE if genotype is phased
#
variant.gt.phased  <-function(hdr,vc,nameOrIdx) {
	if(is.numeric(nameOrIdx)) {
		stopifnot(nameOrIdx>0)
		nameOrIdx<-as.integer(nameOrIdx)-1
		}
	else {
		stopifnot(is.character(nameOrIdx))
		}
	.Call("RBcfCtxVariantGtPhased",hdr,vc,nameOrIdx);
	}


#' @param hdr the vcf header
#' @param vcf the variant
#' @param att the INFO/Attribute
#' @return the the INFO attribute for the given key
variant.string.attributes <-function(hdr,vc,att) {
	stopifnot(is.character(att))
	s <- .Call("RBcfCtxVariantAttributeAsString",hdr,vc,att)
	v <- NULL
	if( is.null(s) || length(s)==0) {
		NULL
	} else {
	     unlist(strsplit(s, split=","))
	     }
	}

#' @param hdr the vcf header
#' @param vcf the variant
#' @param att the INFO/Attribute
#' @return the the INFO attribute for the given key
variant.int.attributes <-function(hdr,vc,att) {
	stopifnot(is.character(att))
	.Call("RBcfCtxVariantAttributeAsInt32",hdr,vc,att)
	}
