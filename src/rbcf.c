/*
The MIT License (MIT)

Copyright (c) 2020 Pierre Lindenbaum PhD.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.

*/
#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include "htslib/synced_bcf_reader.h"
#include "htslib/khash.h"
KHASH_MAP_INIT_STR(vdict, bcf_idinfo_t)

#define WHERENL REprintf("[%s:%d] ",__FILE__,__LINE__)
#define DIE_FAILURE(FormatLiteral,...) do { WHERENL; REprintf("Exiting: " FormatLiteral "\n", ##__VA_ARGS__); exit(EXIT_FAILURE);} while(0)
#define LOG(FormatLiteral, ...) do { WHERENL;REprintf("" FormatLiteral "\n", ##__VA_ARGS__);} while(0)

// java habits
#define final
#define null NULL

/**
 * Structure holding the VCF file , the header and the index
 */ 
typedef struct rbcffile_t
	{
	bcf_srs_t *files;
	/* header */
	bcf_hdr_t *hdr;
	/** map sample2index */
	void* sample2idx;
	}RBcfFile,*RBcfFilePtr;

static void RBcfFileFree(final RBcfFilePtr ptr) {
	if (ptr==NULL) return;
	if (ptr->files!=NULL) {
		 bcf_sr_destroy(ptr->files);
		}
	if(ptr->sample2idx!=NULL) {
		khash_str2int_destroy(ptr->sample2idx);
		}
	Free(ptr);
	}

/**
 * Close resources associated 
 */
SEXP RBcfFileClose(SEXP handle) {	
	void *p = R_ExternalPtrAddr(handle);
	RBcfFileFree((RBcfFilePtr)p);
	R_ClearExternalPtr(handle);
	return ScalarLogical(1);
	}

/**
 * destructor
 * 
 */
static void RBcfFileFinalizer(SEXP handle)
	{
	RBcfFileClose(handle);
	}


/**
 * Open BCF file
 */
SEXP RBcfFileOpen(SEXP Rfilename)
	{
	int i;
	RBcfFilePtr handler =  NULL;
	
	const char* filename= CHAR(asChar(Rfilename));
	if(filename==NULL) {
		LOG("Filename is null");
		return R_NilValue;
		}
	
	handler = (RBcfFilePtr)(Calloc(1UL,RBcfFile));
	if(handler==NULL) {
		LOG("Out of memory");
		goto die;
		}
	handler->files = bcf_sr_init();
	if (handler->files==NULL) {
		LOG("Cannot init %s.",filename,"");
		goto die;
		}
	if ( !bcf_sr_add_reader(handler->files, filename) ) {
		LOG("Failed to read from %s: %s\n",
			filename,
			bcf_sr_strerror(handler->files->errnum)
			);
		goto die;
		}
	
		
	handler->hdr =  handler->files->readers[0].header;
	if(handler->hdr==NULL) {
		LOG("Cannot read header from %s.",filename);
		goto die;
		}
	/** fill map(sampleid,index) */
	handler->sample2idx = khash_str2int_init();
	ASSERT_NOT_NULL(handler->sample2idx);
    for (i=0; i<bcf_hdr_nsamples(handler->hdr); i++) {
        khash_str2int_inc(handler->sample2idx, bcf_hdr_int2id(handler->hdr,BCF_DT_SAMPLE,i));
		}
	/* wrap pointer in SEXP */
	SEXP ext = PROTECT(R_MakeExternalPtr(handler, R_NilValue, R_NilValue));
	/* register destructor */
	R_RegisterCFinalizerEx(ext,RBcfFileFinalizer, TRUE);
     UNPROTECT(1);
	return ext;
	
	die:
		if(handler!=NULL) {
			RBcfFileFree(handler);
			}
	return R_NilValue;
	}

void RBcfHeaderFinalizer(SEXP handler) {
	void* p=R_ExternalPtrAddr(handler);
	if(p==NULL) return;
	bcf_hdr_destroy((bcf_hdr_t*)p);
	}
	
SEXP RBcfHeader(SEXP handler) {
	void *p = NULL;
	int n=0;
	PROTECT(handler);
	p=R_ExternalPtrAddr(handler);
	if(p==NULL) error("Null object");
	bcf_hdr_t *hdr = ((RBcfFilePtr)p)->hdr;
	if(hdr==NULL) error("Cannot get header");
	hdr = bcf_hdr_dup(hdr);
	if(hdr==NULL) error("Cannot dup header");
	SEXP ext = PROTECT(R_MakeExternalPtr(hdr, R_NilValue, R_NilValue));
	R_RegisterCFinalizerEx(ext,RBcfHeaderFinalizer, TRUE);
	UNPROTECT(2);
	return ext;
	}

SEXP RBcfNSamples(SEXP handler) {
	void *p = NULL;
	int n=0;
	PROTECT(handler);
	p=R_ExternalPtrAddr(handler);
	if(p==NULL) error("Null object");
	n = (int)bcf_hdr_nsamples((bcf_hdr_t*)p);
	UNPROTECT(1);
	return ScalarInteger(n);
	}

SEXP RBcfSamples(SEXP handler) {
	void *p = NULL;
	int i, n=0;
	PROTECT(handler);
	p=R_ExternalPtrAddr(handler);
	if(p==NULL) error("Null object");
	bcf_hdr_t *hdr = (bcf_hdr_t*)p;
	SEXP sNames = PROTECT(allocVector(STRSXP,bcf_hdr_nsamples(hdr)));
	for(i=0;i< bcf_hdr_nsamples(hdr);++i) {
		SET_STRING_ELT(sNames, i , mkChar(hdr->samples[i]));
		}
	UNPROTECT(2);
	return sNames;
	}

SEXP RBcfQueryRegion(SEXP handler,SEXP sexpInterval,SEXP sexpRegionIsFile) {
	RBcfFilePtr ptr = NULL;
	int n=0;
	PROTECT(handler);
	PROTECT(sexpInterval);
	PROTECT(sexpRegionIsFile);
	ptr=(RBcfFilePtr)R_ExternalPtrAddr(handler);
	if(ptr==NULL) error("Null object");
	if(sexpInterval==R_NilValue)
		{
		const char* interval= CHAR(asChar(sexpInterval));
		if(interval==NULL) {
			error("interval is null");
			return R_NilValue;
			}
		bcf_sr_set_regions(ptr->files, interval, asLogical(sexpRegionIsFile) )<0 )
		}
	else
		{
		
		}
	bcf_hdr_t *hdr = ((RBcfFilePtr)p)->hdr;
	if(hdr==NULL) error("Cannot get header");
	hdr = bcf_hdr_dup(hdr);
	if(hdr==NULL) error("Cannot dup header");
	SEXP ext = PROTECT(R_MakeExternalPtr(hdr, R_NilValue, R_NilValue));
	R_RegisterCFinalizerEx(ext,RBcfHeaderFinalizer, TRUE);
	UNPROTECT(3);
	return ScalarLogical(1);
	}


static void RBcf1Finalizer(SEXP handler) {
	void* p=R_ExternalPtrAddr(handler);
	if(p==NULL) return;
	bcf_destroy((bcf1_t*)p);
	}
	

SEXP RBcfNextLine(SEXP handler) {
	RBcfFilePtr ptr = NULL;
	int protect=0;
	SEXP ext;
	PROTECT(handler);protect++;
	
	ptr=(RBcfFilePtr)R_ExternalPtrAddr(handler);
	
	 if ( bcf_sr_next_line(ptr->files) ) {
	 	bcf1_t *line = ptr->files->readers[0].buffer[0];
	 	ASSERT_NOT_NULL(line);
	 	line = bcf_dup(line);
	 	ASSERT_NOT_NULL(line);
	 	ext = PROTECT(R_MakeExternalPtr(line, R_NilValue, R_NilValue));protect++;
	 	R_RegisterCFinalizerEx(ext,RBcf1Finalizer(ext), TRUE);
	 	}
	 else
		 {
		 ext = R_NilValue;
		 }
	UNPROTECT(protect);
	return ext;
	}

SEXP RBcfCtxRid(SEXP sexpheader,SEXP sexpctx) {
	int protect=0;
	int tid = -1;
	bcf1_t *ctx;
	SEXP ext;
	PROTECT(sexpctx);protect++;
	ctx=(bcf1_t*)R_ExternalPtrAddr(sexpctx);
	ASSERT_NOT_NULL(ctx);
	ext = ScalarInteger(ctx->rid);
	UNPROTECT(protect);
	return ext;
	}


SEXP RBcfCtxSeqName(SEXP sexpheader,SEXP sexpctx) {
	int protect=0;
	int tid = -1;
	bcf1_t *ctx = NULL;
	bcf_hdr_t* hdr = NULL;
	SEXP ext;
	PROTECT(sexpctx);protect++;
	PROTECT(sexpheader);protect++;
	ctx=(bcf1_t*)R_ExternalPtrAddr(sexpctx);
	ASSERT_NOT_NULL(ctx);
	hdr=(bcf_hdr_t*)R_ExternalPtrAddr(sexpheader);
	ASSERT_NOT_NULL(ctx);
	ext = mkChar(bcf_seqname(hdr,ctx));
	UNPROTECT(protect);
	return ext;
	}

SEXP RBcfCtxPos(SEXP sexpheader,SEXP sexpctx) {
	int protect=0;
	int tid = -1;
	bcf1_t *ctx;
	SEXP ext;
	PROTECT(sexpctx);protect++;
	ctx=(bcf1_t*)R_ExternalPtrAddr(sexpctx);
	ASSERT_NOT_NULL(ctx);
	ext = ScalarInteger(ctx->pos + 1);
	UNPROTECT(protect);
	return ext;
	}

SEXP RBcfCtxHasId(SEXP sexpheader,SEXP sexpctx) {
	int protect=0;
	int tid = -1;
	bcf1_t *ctx = NULL;
	SEXP ext;
	PROTECT(sexpctx);protect++;
	ctx=(bcf1_t*)R_ExternalPtrAddr(sexpctx);
	ASSERT_NOT_NULL(ctx);
	ext = ScalarLogical(ctx->d.id);
	UNPROTECT(protect);
	return ext;
	}

SEXP RBcfCtxId(SEXP sexpheader,SEXP sexpctx) {
	int protect=0;
	int tid = -1;
	bcf1_t *ctx = NULL;
	SEXP ext;
	PROTECT(sexpctx);protect++;
	ctx=(bcf1_t*)R_ExternalPtrAddr(sexpctx);
	ASSERT_NOT_NULL(ctx);
	if(ctx->d.id) {
		ext = mkChar(ctx->d.id);
		}
	else
		{
		ext  =	R_NilValue;
		}
	UNPROTECT(protect);
	return ext;
	}


SEXP RBcfCtxEnd(SEXP sexpheader,SEXP sexpctx) {
	int protect=0;
	int tid = -1;
	bcf1_t *ctx;
	SEXP ext;
	PROTECT(sexpctx);protect++;
	ctx=(bcf1_t*)R_ExternalPtrAddr(sexpctx);
	ASSERT_NOT_NULL(ctx);
	ext = ScalarInteger(line->pos /* 0 based */ + ctx->rlen);
	UNPROTECT(protect);
	return ext;
	}

SEXP RBcfCtxNAlleles(SEXP sexpheader,SEXP sexpctx) {
	int protect=0;
	int tid = -1;
	bcf1_t *ctx;
	SEXP ext;
	PROTECT(sexpctx);protect++;
	ctx=(bcf1_t*)R_ExternalPtrAddr(sexpctx);
	ASSERT_NOT_NULL(ctx);
	ext = ScalarInteger(ctx->n_allele);
	UNPROTECT(protect);
	return ext;
	}

SEXP RBcfCtxAlleles(SEXP sexpheader,SEXP sexpctx) {
	int protect=0;
	int tid = -1;
	bcf1_t *ctx;
	SEXP ext;
	PROTECT(sexpctx);protect++;
	ctx=(bcf1_t*)R_ExternalPtrAddr(sexpctx);
	ASSERT_NOT_NULL(ctx);
	ext = ScalarInteger(ctx->n_allele);
	UNPROTECT(protect);
	return ext;
	}

SEXP RBcfCtxAlleles(SEXP sexpheader,SEXP sexpctx) {
	int protect=0;
	bcf1_t *ctx;
	SEXP ext;
	PROTECT(sexpctx);protect++;
	ctx=(bcf1_t*)R_ExternalPtrAddr(sexpctx);
	ASSERT_NOT_NULL(ctx);
	ext = PROTECT(allocVector(STRSXP,ctx->n_allele));protect++;
	for(i=0;i< ctx->n_allele;++i) {
		SET_STRING_ELT(ext, i , mkChar(ctx->allele[i]));
		}
	UNPROTECT(protect);
	return ext;
	}

SEXP RBcfCtxReference(SEXP sexpheader,SEXP sexpctx) {
	int protect=0;
	bcf1_t *ctx;
	SEXP ext;
	PROTECT(sexpctx);protect++;
	ctx=(bcf1_t*)R_ExternalPtrAddr(sexpctx);
	ASSERT_NOT_NULL(ctx);
	if(ctx->n_allele>0) {
		ext = mkChar(ctx->allele[0]);
		}
	else
		{
		ext  = R_NilValue;
		}
	UNPROTECT(protect);
	return ext;
	}

SEXP RBcfCtxAlternateAlleles(SEXP sexpheader,SEXP sexpctx) {
	int protect=0;
	bcf1_t *ctx;
	SEXP ext;
	PROTECT(sexpctx);protect++;
	ctx=(bcf1_t*)R_ExternalPtrAddr(sexpctx);
	ASSERT_NOT_NULL(ctx);
	ext = PROTECT(allocVector(STRSXP,(ctx->n_allele==0?0:ctx->n_allele-1));protect++;
	for(i=1;i< ctx->n_allele;++i) {
		SET_STRING_ELT(ext, i-1, mkChar(ctx->allele[i]));
		}
	UNPROTECT(protect);
	return ext;
	}
	
SEXP RBcfCtxHasQual(SEXP sexpheader,SEXP sexpctx) {
	int protect=0;
	bcf1_t *ctx;
	SEXP ext;
	PROTECT(sexpctx);protect++;
	ctx=(bcf1_t*)R_ExternalPtrAddr(sexpctx);
	ASSERT_NOT_NULL(ctx);
	ext = ScalarLogical(!bcf_float_is_missing(ctx->qual));
	UNPROTECT(protect);
	return ext;
	}

SEXP RBcfCtxQual(SEXP sexpheader,SEXP sexpctx) {
	int protect=0;
	bcf1_t *ctx;
	SEXP ext;
	PROTECT(sexpctx);protect++;
	ctx=(bcf1_t*)R_ExternalPtrAddr(sexpctx);
	ASSERT_NOT_NULL(ctx);
	if(bcf_float_is_missing(ctx->qual)) {
		ctx = R_NilValue;
		}
	else
		{
		ext = ext = ScalarNumeric(ctx->qual);
		}
	UNPROTECT(protect);
	return ext;
	}


SEXP RBcfCtxFiltered(SEXP sexpheader,SEXP sexpctx) {
	int protect=0;
	bcf1_t *ctx;
	bcf_hdr_t* hdr;
	SEXP ext;
	PROTECT(sexpctx);protect++;
	PROTECT(sexpheader);protect++;
	ctx=(bcf1_t*)R_ExternalPtrAddr(sexpctx);
	ASSERT_NOT_NULL(ctx);
	hdr=(bcf_hdr_t*)R_ExternalPtrAddr(sexpheader);
	ASSERT_NOT_NULL(hdr);
	ext = ScalerLogical(!bcf_has_filter(hdr,ctx,"PASS"));
	UNPROTECT(protect);
	return ext;
	}

SEXP RBcfCtxFilters(SEXP sexpheader,SEXP sexpctx) {
	int protect=0;
	bcf1_t *ctx;
	bcf_hdr_t* hdr;
	SEXP ext;
	PROTECT(sexpctx);protect++;
	PROTECT(sexpheader);protect++;
	
	ctx=(bcf1_t*)R_ExternalPtrAddr(sexpctx);
	ASSERT_NOT_NULL(ctx);
	hdr=(bcf_hdr_t*)R_ExternalPtrAddr(sexpheader);
	ASSERT_NOT_NULL(hdr);
	
	if(bcf_has_filter(hdr,ctx,"PASS")) {
		ext = PROTECT(allocVector(STRSXP,0));protect++;
		}
	else
		{
		ext = PROTECT(allocVector(STRSXP,ctx->n_flt));protect++;
		for(i=0;i< ctx->n_flt;++i) {
			SET_STRING_ELT(ext, i , mkChar(hdr->id[BCF_DT_ID][ctx->d.flt[i]].key));
			}
		}
	UNPROTECT(protect);
	return ext;
	}

SEXP RBcfCtxVariantTypes(SEXP sexpheader,SEXP sexpctx) {
	int protect=0;
	bcf1_t *ctx;
	SEXP ext;
	PROTECT(sexpctx);protect++;
	ctx=(bcf1_t*)R_ExternalPtrAddr(sexpctx);
	ASSERT_NOT_NULL(ctx);
	ext = ScalarInteger(bcf_get_variant_types(ctx));
	UNPROTECT(protect);
	return ext;
	}
	
SEXP RBcfCtxVariantIsSnp(SEXP sexpheader,SEXP sexpctx) {
	int protect=0;
	bcf1_t *ctx;
	SEXP ext;
	PROTECT(sexpctx);protect++;
	ctx=(bcf1_t*)R_ExternalPtrAddr(sexpctx);
	ASSERT_NOT_NULL(ctx);
	ext = ScalarLogical(bcf_is_snp(ctx));
	UNPROTECT(protect);
	return ext;
	}

SEXP RBcfCtxVariantMaxPloidy(SEXP sexpheader,SEXP sexpctx) {
	int protect=0;
	bcf1_t *ctx;
	bcf_hdr_t* hdr;
	int32_t *gt_arr = NULL, ngt_arr = 0;
	SEXP ext;
	PROTECT(sexpctx);protect++;
	PROTECT(sexpheader);protect++;
	
	ctx=(bcf1_t*)R_ExternalPtrAddr(sexpctx);
	ASSERT_NOT_NULL(ctx);
	hdr=(bcf_hdr_t*)R_ExternalPtrAddr(sexpheader);
	ASSERT_NOT_NULL(hdr);
	
	ngt = bcf_get_genotypes(hdr,ctx, &gt_arr, &ngt_arr);
	if ( ngt<=0 ) {
		ext = ScalarInteger(0);
		}
	else
		{
		int nsmpl = bcf_hdr_nsamples(hdr);
	 	ext = ScalarInteger(ngt/nsmpl);
	 	}
	UNPROTECT(protect);
	return ext;
	}

struct GenotypeShuttle {
	int phased;
	int nocall;
	int* alleles;
	size_t allele_count;
	size_t allele_capacity;
	int error_flag;
	};



static void scanGenotype(SEXP sexpheader,SEXP sexpctx,SEXP sexpgtidx,struct GenotypeShuttle* shuttle) {
	int protect=0;
	bcf1_t *ctx;
	bcf_hdr_t* hdr;
	int sample_idx = 0;
	int32_t *gt_arr = NULL, ngt_arr = 0;
	SEXP ext;
	memset((void*)shuttle,0,sizeof(struct GenotypeShuttle));
	PROTECT(sexpctx);protect++;
	PROTECT(sexpheader);protect++;
	PROTECT(sexpgtidx);protect++;
	
	ctx=(bcf1_t*)R_ExternalPtrAddr(sexpctx);
	ASSERT_NOT_NULL(ctx);
	hdr=(bcf_hdr_t*)R_ExternalPtrAddr(sexpheader);
	ASSERT_NOT_NULL(hdr);
	
	if(IS_CHARACTER(sexpgtidx))
		{
		khint_t k;
		vdict_t* d = NULL;
		const char* sn = CHAR(asChar(sexpgtidx));
		if(sn==NULL) {
			shuttle->error_flag  = 1;
			return;
			}
		sample_idx = bcf_hdr_id2int(hdr,BCF_DT_SAMPLE,sn);
		if(sample_idx<0) {
			warning("unknown sample %s",sn);
			shuttle->error_flag  = 1;
			return;
			}
        }
	else
		{
		sample_idx = asInteger(sexpgtidx);
		}
	
	if(sample_idx<0 || sample_idx>=bcf_hdr_nsamples(hdr)) {
		warning("bad sample index 0<=%d<%d",sample_idx,bcf_hdr_nsamples(hdr));
		shuttle->error_flag  = 1;
		return;
		}
	ngt = bcf_get_genotypes(hdr,ctx, &gt_arr, &ngt_arr);
	 int max_ploidy = ngt/nsmpl;
     
      int32_t *ptr = gt + i*sample_idx;
      for (j=0; j<max_ploidy; j++)
          {
	         // if true, the sample has smaller ploidy
	         if ( ptr[j]==bcf_int32_vector_end ) break;
  
	         // missing allele
	         if ( bcf_gt_is_missing(ptr[j]) ) continue;
  
  			 // is phased?
  			 if(bcf_gt_is_phased(ptr[j]) {
	         	shuttle->phased = 1;
	         	}
  
	         // the VCF 0-based allele index
	         int allele_index = bcf_gt_allele(ptr[j]);
  			 if(shuttle->alleles != NULL) {
  			 	if(shuttle->allele_count + 1 >= shuttle->allele_capacity) {
  			 		shuttle->allele_capacity++;
  			 		shuttle->alleles = (int*)Realloc(shuttle->alleles,shuttle->allele_capacity);
  			 		ASSERT_NOT_NULL(shuttle->alleles);	
  			 		}
  			 	shuttle->alleles[shuttle->allele_count] = allele_index ;
  			 	}
  			shuttle->allele_count++;
  		    }
	
	UNPROTECT(protect);
	}

