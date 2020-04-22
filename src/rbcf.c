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
#ifdef SRS_BCF
#include "htslib/synced_bcf_reader.h"
#else
#include "htslib/vcf.h"
#endif
#include "htslib/tbx.h"
#include "htslib/kseq.h"

#define WHERENL REprintf("[%s:%d] ",__FILE__,__LINE__)
#define LOG(FormatLiteral, ...) do { WHERENL;REprintf("[LOG]" FormatLiteral "\n", ##__VA_ARGS__);} while(0)
#define BCF_WARNING(FormatLiteral, ...) do { WHERENL;REprintf("[WARN]" FormatLiteral "\n", ##__VA_ARGS__);} while(0)
#define BCF_ERROR(FormatLiteral, ...) do { WHERENL;REprintf("[FATAL]" FormatLiteral "\n", ##__VA_ARGS__);error("Fatal Error");} while(0)
#define ASSERT_NOT_NULL(ptr) do {if(ptr==NULL) BCF_ERROR("NULL pointer exception"); } while(0)
// java habits
#define final
#define null NULL


/**
 * Structure holding the VCF file , the header and the index
 */ 
typedef struct rbcffile_t
	{

	htsFile *fp;
	char* filename;
	tbx_t *tbx;
	hts_idx_t *idx;
	hts_itr_t *itr;
	bcf1_t* tmp_ctx; // tmp variant ctx for reading
	kstring_t tmp_line; // for reading line

	/* header */
	bcf_hdr_t *hdr;
	}RBcfFile,*RBcfFilePtr;

static void RBcfFileFree(final RBcfFilePtr ptr) {
	if (ptr==NULL) return;

	free(ptr->filename);
    free(ptr->tmp_line.s);
	if ( ptr->tmp_ctx) bcf_destroy1(ptr->tmp_ctx);
	if ( ptr->itr ) hts_itr_destroy(ptr->itr);
	if ( ptr->tbx ) tbx_destroy(ptr->tbx);
    if ( ptr->idx ) hts_idx_destroy(ptr->idx);
	if ( ptr->hdr ) bcf_hdr_destroy(ptr->hdr);
	if ( ptr->fp ) hts_close(ptr->fp);	
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
SEXP RBcfFileOpen(SEXP Rfilename,SEXP sexpRequireIdx) {
	RBcfFilePtr handler =  NULL;
	int nprotect=0;
	PROTECT(Rfilename);nprotect++;
	PROTECT(sexpRequireIdx);nprotect++;
	
	const char* filename= CHAR(asChar(Rfilename));
	if(filename==NULL) {
		LOG("Filename is null");
		return R_NilValue;
		}
	
	handler = (RBcfFilePtr)(Calloc(1UL,RBcfFile));
	if(handler==NULL) {
		BCF_WARNING("Out of memory");
		goto die;
		}
	
	handler->tmp_ctx = bcf_init1();
	if(handler->tmp_ctx ==NULL) {
		BCF_WARNING("Out of memory");
		goto die;
		}
	
	handler->filename = strdup(filename);
	if(handler->filename==NULL) {
		BCF_WARNING("Out of memory");
		goto die;
		}
	
	handler->fp = hts_open(filename,"rb");
	if (handler->fp==NULL) {
		BCF_WARNING("Failed to read from %s.\n",filename);
		goto die;
		}
	if(asLogical(sexpRequireIdx)) {
	     enum htsExactFormat format = hts_get_format(handler->fp)->format;
	     switch(format) {
	       case vcf:
	         if ( handler->fp->format.compression!=bgzf )
                {
                BCF_WARNING("not bgzf format %s",filename);
				goto die;
                }
			 handler->tbx = tbx_index_load(filename);
			 if(handler->tbx==NULL) {
				  BCF_WARNING("Cannot open tabix index for %s",filename);
				  goto die;
				  }
			 break;
	       case bcf:
	         if ( handler->fp->format.compression!=bgzf )
                {
                BCF_WARNING("not bgzf format %s",filename);
				goto die;
                }
			 handler->idx = bcf_index_load(filename);
			 if(handler->idx==NULL) {
				  BCF_WARNING("Cannot open idx index for %s",filename);
				  goto die;
				  }
			break; 
	     default:
	     	BCF_WARNING("unknown vcf format for %s.",filename);
	       	goto die;
			break;
	    }// end switch
	 }//end if index required

	handler->hdr =	bcf_hdr_read(handler->fp);
	if(handler->hdr==NULL) {
		LOG("Cannot read header from %s.",filename);
		goto die;
		}

	/* wrap pointer in SEXP */
	SEXP ext = PROTECT(R_MakeExternalPtr(handler, R_NilValue, R_NilValue));nprotect++;
	/* register destructor */
	R_RegisterCFinalizerEx(ext,RBcfFileFinalizer, TRUE);
    UNPROTECT(nprotect);
	return ext;
	
	die:
		LOG("error opening file");
		if(handler!=NULL) {
			RBcfFileFree(handler);
			}
	UNPROTECT(nprotect);
	return R_NilValue;
	}

void RBcfHeaderFinalizer(SEXP handler) {
	void* p=R_ExternalPtrAddr(handler);
	if(p==NULL) return;
	bcf_hdr_destroy((bcf_hdr_t*)p);
	}


SEXP RBcfSeqNames(SEXP sexpFile) {
	int i, n=0;
	SEXP ext;
	int nprotect=0;
	PROTECT(sexpFile);nprotect++;
	RBcfFilePtr p=(RBcfFilePtr)R_ExternalPtrAddr(sexpFile);
	ASSERT_NOT_NULL(p);
	const char **names;
	
	if(p->tbx) {
		names =  tbx_seqnames(p->tbx, &n);
		}
	else if( p->idx) {
		names = bcf_hdr_seqnames(p->hdr, &n);
		}
	else
		{
		BCF_WARNING("VCF File \"%s\" is not indexed.",p->filename);
		names = NULL;
		}
	if(names) {
		ext = PROTECT(allocVector(STRSXP,n));nprotect++;
		for(i=0;i< n ;++i) {
			SET_STRING_ELT(ext, i , mkChar(names[i]));
			}
		free(names);
		}
	else
		{
		ext = R_NilValue;
		}
	UNPROTECT(nprotect);
	return ext;
	}
	
SEXP RBcfHeader(SEXP handler) {
	void *p = NULL;
	PROTECT(handler);
	p=R_ExternalPtrAddr(handler);
	if(p==NULL) BCF_ERROR("Null object");
	bcf_hdr_t *hdr = ((RBcfFilePtr)p)->hdr;
	if(hdr==NULL) BCF_ERROR("Cannot get header");
	hdr = bcf_hdr_dup(hdr);
	if(hdr==NULL) BCF_ERROR("Cannot dup header");
	SEXP ext = PROTECT(R_MakeExternalPtr(hdr, R_NilValue, R_NilValue));
	R_RegisterCFinalizerEx(ext,RBcfHeaderFinalizer, TRUE);
	UNPROTECT(2);
	return ext;
	}

SEXP RBcfNSamples(SEXP sexpHeader) {
	int n=0;
	PROTECT(sexpHeader);
	bcf_hdr_t* hdr=(bcf_hdr_t*)R_ExternalPtrAddr(sexpHeader);
	ASSERT_NOT_NULL(hdr);
	n = (int)bcf_hdr_nsamples(hdr);
	UNPROTECT(1);
	return ScalarInteger(n);
	}

SEXP RBcfSamples(SEXP handler) {
	int i;
	PROTECT(handler);
	bcf_hdr_t *hdr = (bcf_hdr_t*)R_ExternalPtrAddr(handler);
	ASSERT_NOT_NULL(hdr);
	SEXP sNames = PROTECT(allocVector(STRSXP,bcf_hdr_nsamples(hdr)));
	for(i=0;i< bcf_hdr_nsamples(hdr);++i) {
		SET_STRING_ELT(sNames, i , mkChar(hdr->samples[i]));
		}
	UNPROTECT(2);
	return sNames;
	}

SEXP RBcfSampleAtIndex0(SEXP handler,SEXP sexpIndex) {
	int nprotect=0;
	SEXP ext;
	PROTECT(handler);++nprotect;
	PROTECT(sexpIndex);++nprotect;
	bcf_hdr_t *hdr = (bcf_hdr_t*)R_ExternalPtrAddr(handler);
	ASSERT_NOT_NULL(hdr);
	int idx = asInteger(sexpIndex);
	if(idx<0 || idx >= bcf_hdr_nsamples(hdr)) {
		BCF_WARNING("sample idx out of range %d/%d\n",bcf_hdr_nsamples(hdr));
		ext = R_NilValue;
		}
	else
		{
		ext =  mkString(hdr->samples[idx]);
		}
	UNPROTECT(nprotect);
	return ext;
	}

SEXP RBcfHeaderDict(SEXP sexpHeader) {
	int i,nprotect=0;
	SEXP res,vChrom,vSize,vRowNames;
	PROTECT(sexpHeader);++nprotect;
	bcf_hdr_t *hdr = (bcf_hdr_t*)R_ExternalPtrAddr(sexpHeader);
	ASSERT_NOT_NULL(hdr);
	int dict_size=  hdr->n[BCF_DT_CTG];
	PROTECT(res = Rf_allocVector(VECSXP,2));nprotect++;
	PROTECT(vChrom = Rf_allocVector(STRSXP,dict_size)); nprotect++;
   	PROTECT(vSize = Rf_allocVector(VECSXP, dict_size)); nprotect++;
	PROTECT(vRowNames = Rf_allocVector(STRSXP, dict_size)); nprotect++;
	
	for (i=0; i< dict_size; i++) {
		if ( !hdr->id[BCF_DT_CTG][i].val ) continue;
		const char* chrom_name = hdr->id[BCF_DT_CTG][i].key;
		SET_STRING_ELT(vChrom, i,  mkChar(chrom_name) );
		SET_STRING_ELT(vRowNames, i,  mkChar(chrom_name) );
		int len = hdr->id[BCF_DT_CTG][i].val->info[0];
		SET_VECTOR_ELT(vSize, i, ScalarInteger(len) );
	    }

	
	/* set the columns of the table */
	SET_VECTOR_ELT(res, 0, vChrom);
  	SET_VECTOR_ELT(res, 1, vSize);

	/* set the columns' name of the table */

	SEXP sNames = PROTECT(allocVector(STRSXP, 2));nprotect++;
	SET_STRING_ELT(sNames, 0, mkChar("chrom"));
	SET_STRING_ELT(sNames, 1, mkChar("size"));
	namesgets(res, sNames);
	
	// https://stackoverflow.com/questions/23547625
   SEXP cls = PROTECT(allocVector(STRSXP, 1));nprotect++;
   SET_STRING_ELT(cls, 0, mkChar("data.frame"));
   classgets(res, cls);
	
	setAttrib(res, R_RowNamesSymbol, vRowNames);
	
	UNPROTECT(nprotect);
	return res;
	}


SEXP RBcfQueryRegion(SEXP sexpFile,SEXP sexpInterval) {
	int nprotect=0;
	SEXP ext;
	PROTECT(sexpFile);nprotect++;
	PROTECT(sexpInterval);nprotect++;
	RBcfFilePtr ptr=(RBcfFilePtr)R_ExternalPtrAddr(sexpFile);
	ASSERT_NOT_NULL(ptr);
	
	
	const char* interval= CHAR(asChar(sexpInterval));
	ASSERT_NOT_NULL(interval);
	if(ptr->itr) hts_itr_destroy(ptr->itr);
	ptr->itr = NULL;
	
	if( ptr->tbx!=NULL) {
		ptr->itr = tbx_itr_querys(ptr->tbx,interval);
		}
	else if(ptr->idx!=NULL) {
		ptr->itr = bcf_itr_querys(ptr->idx,ptr->hdr,interval);
		}
	else
		{
		BCF_ERROR("Cannot query vcf file \"%s\" (no index available)",ptr->filename);
		}
	if(ptr->itr==NULL) {
		BCF_WARNING("Query failed \"%s\" for \"%s\".",interval, ptr->filename);
		}
	ext = ScalarLogical( ptr->itr!=NULL );
	UNPROTECT(nprotect);
	return ext;
	}


static void RBcf1Finalizer(SEXP handler) {
	bcf1_t* p=(bcf1_t*)R_ExternalPtrAddr(handler);
	if(p) bcf_destroy(p);
	}
	

SEXP RBcfNextLine(SEXP handler) {
	int ret=0;
	int nprotect=0;
	SEXP ext;
	PROTECT(handler);nprotect++;
	RBcfFilePtr reader=(RBcfFilePtr)R_ExternalPtrAddr(handler);
	bcf1_t *line = NULL;
	ASSERT_NOT_NULL(reader);
	
	if ( reader->fp->format.format==vcf ) {
  		if( reader->itr) {
  			//LOG("using tbx_itr_next");
  			ret= tbx_itr_next(reader->fp, reader->tbx, reader->itr, &(reader->tmp_line));
  			}
  		else //streaming
  			{
  			//LOG("using getline");
	  		ret= hts_getline(reader->fp, KS_SEP_LINE,&(reader->tmp_line));
	  		}
        if ( ret < 0 ) {
            line = NULL;
        	}
        else
        	{
        	ASSERT_NOT_NULL(reader->tmp_line.s);
        	ASSERT_NOT_NULL(reader->hdr);
        	ASSERT_NOT_NULL(reader->tmp_ctx);
            ret = vcf_parse1(&(reader->tmp_line), reader->hdr, reader->tmp_ctx);
            if ( ret<0 ) BCF_ERROR("error while reading vcf line \"%s\" (error=%d)",reader->tmp_line.s,ret);
            line = reader->tmp_ctx;
            }
    	}
    else if ( reader->fp->format.format==bcf ) {
    	if( reader->itr) {
        	ret = bcf_itr_next(reader->fp, reader->itr, reader->tmp_ctx);
        	}
        else
        	{
        	ret = bcf_read1(reader->fp, reader->hdr, reader->tmp_ctx);
        	}
        if ( ret < -1 ) {
        	BCF_ERROR("error while reading bcf line");
            }
        if ( ret >= 0 ) {
        	line = reader->tmp_ctx;
        	}
    	}
    else
    	{
        BCF_ERROR("illegal state");
    	}
	
	 if ( line !=NULL ) {
	 	line = bcf_dup(line);
	 	ASSERT_NOT_NULL(line);
	 	
	 	ext = PROTECT(R_MakeExternalPtr(line, R_NilValue, R_NilValue));nprotect++;
	 	R_RegisterCFinalizerEx(ext,RBcf1Finalizer, TRUE);
	 	}
	 else
		 {
		 ext = R_NilValue;		 
		 }
	
	UNPROTECT(nprotect);
	return ext;
	}

SEXP RBcfCtxRid(SEXP sexpheader,SEXP sexpctx) {
	int nprotect=0;
	bcf1_t *ctx;
	SEXP ext;
	PROTECT(sexpctx);nprotect++;
	ctx=(bcf1_t*)R_ExternalPtrAddr(sexpctx);
	ASSERT_NOT_NULL(ctx);
	ext = ScalarInteger(ctx->rid);
	UNPROTECT(nprotect);
	return ext;
	}


SEXP RBcfCtxSeqName(SEXP sexpheader,SEXP sexpctx) {
	int nprotect=0;
	bcf1_t *ctx = NULL;
	bcf_hdr_t* hdr = NULL;
	SEXP ext;
	PROTECT(sexpctx);nprotect++;
	PROTECT(sexpheader);nprotect++;
	ctx=(bcf1_t*)R_ExternalPtrAddr(sexpctx);
	ASSERT_NOT_NULL(ctx);
	hdr=(bcf_hdr_t*)R_ExternalPtrAddr(sexpheader);
	ASSERT_NOT_NULL(ctx);
	ext = mkString(bcf_seqname(hdr,ctx));
	UNPROTECT(nprotect);
	return ext;
	}

SEXP RBcfCtxPos(SEXP sexpheader,SEXP sexpctx) {
	int nprotect=0;
	bcf1_t *ctx;
	SEXP ext;
	PROTECT(sexpctx);nprotect++;
	ctx=(bcf1_t*)R_ExternalPtrAddr(sexpctx);
	ASSERT_NOT_NULL(ctx);
	ext = ScalarInteger(ctx->pos + 1);
	UNPROTECT(nprotect);
	return ext;
	}

SEXP RBcfCtxHasId(SEXP sexpheader,SEXP sexpctx) {
	int nprotect=0;
	bcf1_t *ctx = NULL;
	SEXP ext;
	PROTECT(sexpctx);nprotect++;
	ctx=(bcf1_t*)R_ExternalPtrAddr(sexpctx);
	ASSERT_NOT_NULL(ctx);
	ext = ScalarLogical(ctx->d.id!=NULL);
	UNPROTECT(nprotect);
	return ext;
	}

SEXP RBcfCtxId(SEXP sexpheader,SEXP sexpctx) {
	int nprotect=0;
	bcf1_t *ctx = NULL;
	SEXP ext;
	PROTECT(sexpctx);nprotect++;
	ctx=(bcf1_t*)R_ExternalPtrAddr(sexpctx);
	ASSERT_NOT_NULL(ctx);
	if(ctx->d.id) {
		ext = mkString(ctx->d.id);
		}
	else
		{
		ext  =	R_NilValue;
		}
	UNPROTECT(nprotect);
	return ext;
	}


SEXP RBcfCtxEnd(SEXP sexpheader,SEXP sexpctx) {
	int nprotect=0;
	bcf1_t *ctx;
	SEXP ext;
	PROTECT(sexpctx);nprotect++;
	ctx=(bcf1_t*)R_ExternalPtrAddr(sexpctx);
	ASSERT_NOT_NULL(ctx);
	ext = ScalarInteger(ctx->pos /* 0 based */ + ctx->rlen);
	UNPROTECT(nprotect);
	return ext;
	}

SEXP RBcfCtxNAlleles(SEXP sexpheader,SEXP sexpctx) {
	int nprotect=0;
	bcf1_t *ctx;
	SEXP ext;
	PROTECT(sexpctx);nprotect++;
	ctx=(bcf1_t*)R_ExternalPtrAddr(sexpctx);
	ASSERT_NOT_NULL(ctx);
	ext = ScalarInteger(ctx->n_allele);
	UNPROTECT(nprotect);
	return ext;
	}


SEXP RBcfCtxAlleles(SEXP sexpheader,SEXP sexpctx) {
	int i=0;
	int nprotect=0;
	SEXP ext;
	PROTECT(sexpctx);nprotect++;
	bcf1_t * ctx=(bcf1_t*)R_ExternalPtrAddr(sexpctx);
	ASSERT_NOT_NULL(ctx);
	bcf_unpack(ctx,BCF_UN_STR);
	ext = PROTECT(allocVector(STRSXP,ctx->n_allele));nprotect++;
	for(i=0;i< ctx->n_allele;++i) {
		ASSERT_NOT_NULL(ctx->d.allele);
		ASSERT_NOT_NULL(ctx->d.allele[i]);
		SET_STRING_ELT(ext, i , mkChar(ctx->d.allele[i]));
		}
	UNPROTECT(nprotect);
	return ext;
	}

SEXP RBcfCtxReference(SEXP sexpheader,SEXP sexpctx) {
	int nprotect=0;
	bcf1_t *ctx;
	SEXP ext;
	PROTECT(sexpctx);nprotect++;
	ctx=(bcf1_t*)R_ExternalPtrAddr(sexpctx);
	ASSERT_NOT_NULL(ctx);
	bcf_unpack(ctx,BCF_UN_STR);
	if(ctx->n_allele>0) {
		ext = mkString(ctx->d.allele[0]);
		}
	else
		{
		ext  = R_NilValue;
		}
	UNPROTECT(nprotect);
	return ext;
	}

SEXP RBcfCtxAlternateAlleles(SEXP sexpheader,SEXP sexpctx) {
	int nprotect=0;
	int i;
	bcf1_t *ctx;
	SEXP ext;
	PROTECT(sexpctx);nprotect++;
	ctx=(bcf1_t*)R_ExternalPtrAddr(sexpctx);
	bcf_unpack(ctx,BCF_UN_STR);
	ASSERT_NOT_NULL(ctx);
	ext = PROTECT(allocVector(STRSXP,(ctx->n_allele<1?0:ctx->n_allele-1)));nprotect++;
	for(i=1;i< ctx->n_allele;++i) {
		SET_STRING_ELT(ext, i-1, mkChar(ctx->d.allele[i]));
		}
	UNPROTECT(nprotect);
	return ext;
	}
	
SEXP RBcfCtxHasQual(SEXP sexpheader,SEXP sexpctx) {
	int nprotect=0;
	bcf1_t *ctx;
	SEXP ext;
	PROTECT(sexpctx);nprotect++;
	ctx=(bcf1_t*)R_ExternalPtrAddr(sexpctx);
	ASSERT_NOT_NULL(ctx);
	ext = ScalarLogical(!bcf_float_is_missing(ctx->qual));
	UNPROTECT(nprotect);
	return ext;
	}

SEXP RBcfCtxQual(SEXP sexpheader,SEXP sexpctx) {
	int nprotect=0;
	bcf1_t *ctx;
	SEXP ext;
	PROTECT(sexpctx);nprotect++;
	ctx=(bcf1_t*)R_ExternalPtrAddr(sexpctx);
	ASSERT_NOT_NULL(ctx);
	if(bcf_float_is_missing(ctx->qual)) {
		ext = R_NilValue;
		}
	else
		{
		ext = ScalarReal(ctx->qual);
		}
	UNPROTECT(nprotect);
	return ext;
	}


SEXP RBcfCtxFiltered(SEXP sexpheader,SEXP sexpctx) {
	int nprotect=0;
	bcf1_t *ctx;
	bcf_hdr_t* hdr;
	SEXP ext;
	PROTECT(sexpctx);nprotect++;
	PROTECT(sexpheader);nprotect++;
	ctx=(bcf1_t*)R_ExternalPtrAddr(sexpctx);
	ASSERT_NOT_NULL(ctx);
	bcf_unpack(ctx,BCF_UN_FLT);
	hdr=(bcf_hdr_t*)R_ExternalPtrAddr(sexpheader);
	ASSERT_NOT_NULL(hdr);
	ext = ScalarLogical(!bcf_has_filter(hdr,ctx,"PASS"));
	UNPROTECT(nprotect);
	return ext;
	}

SEXP RBcfCtxFilters(SEXP sexpheader,SEXP sexpctx) {
	int nprotect=0;
	int i=0;
	bcf1_t *ctx;
	bcf_hdr_t* hdr;
	SEXP ext;
	PROTECT(sexpctx);nprotect++;
	PROTECT(sexpheader);nprotect++;
	
	ctx=(bcf1_t*)R_ExternalPtrAddr(sexpctx);
	ASSERT_NOT_NULL(ctx);
	hdr=(bcf_hdr_t*)R_ExternalPtrAddr(sexpheader);
	ASSERT_NOT_NULL(hdr);
	bcf_unpack(ctx,BCF_UN_FLT);
	
	if(bcf_has_filter(hdr,ctx,"PASS")) {
		ext = PROTECT(allocVector(STRSXP,0));nprotect++;
		}
	else
		{
		ext = PROTECT(allocVector(STRSXP,ctx->d.n_flt));nprotect++;
		for(i=0;i< ctx->d.n_flt;++i) {
			SET_STRING_ELT(ext, i , mkChar(hdr->id[BCF_DT_ID][ctx->d.flt[i]].key));
			}
		}
	UNPROTECT(nprotect);
	return ext;
	}

SEXP RBcfCtxVariantTypes(SEXP sexpheader,SEXP sexpctx) {
	int nprotect=0;
	bcf1_t *ctx;
	SEXP ext;
	PROTECT(sexpctx);nprotect++;
	ctx=(bcf1_t*)R_ExternalPtrAddr(sexpctx);
	ASSERT_NOT_NULL(ctx);
	ext = ScalarInteger(bcf_get_variant_types(ctx));
	UNPROTECT(nprotect);
	return ext;
	}
	
SEXP RBcfCtxVariantIsSnp(SEXP sexpheader,SEXP sexpctx) {
	int nprotect=0;
	bcf1_t *ctx;
	SEXP ext;
	PROTECT(sexpctx);nprotect++;
	ctx=(bcf1_t*)R_ExternalPtrAddr(sexpctx);
	ASSERT_NOT_NULL(ctx);
	ext = ScalarLogical(bcf_is_snp(ctx));
	UNPROTECT(nprotect);
	return ext;
	}

SEXP RBcfCtxVariantMaxPloidy(SEXP sexpheader,SEXP sexpctx) {
	int nprotect=0;
	bcf1_t *ctx;
	bcf_hdr_t* hdr;
	int32_t *gt_arr = NULL, ngt_arr = 0;
	SEXP ext;
	PROTECT(sexpctx);nprotect++;
	PROTECT(sexpheader);nprotect++;
	
	ctx=(bcf1_t*)R_ExternalPtrAddr(sexpctx);
	ASSERT_NOT_NULL(ctx);
	bcf_unpack(ctx,BCF_UN_IND);
	
	hdr=(bcf_hdr_t*)R_ExternalPtrAddr(sexpheader);
	ASSERT_NOT_NULL(hdr);
	
	int ngt = bcf_get_genotypes(hdr,ctx, &gt_arr, &ngt_arr);
	if ( ngt<=0 ) {
		ext = ScalarInteger(0);
		}
	else
		{
		int nsmpl = bcf_hdr_nsamples(hdr);
	 	ext = ScalarInteger(ngt/nsmpl);
	 	}
	UNPROTECT(nprotect);
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
	int nprotect=0;
	bcf1_t *ctx;
	bcf_hdr_t* hdr;
	int sample_idx = 0;
	int32_t *gt_arr = NULL, ngt_arr = 0;

	PROTECT(sexpctx);nprotect++;
	PROTECT(sexpheader);nprotect++;
	PROTECT(sexpgtidx);nprotect++;
	
	ctx=(bcf1_t*)R_ExternalPtrAddr(sexpctx);
	ASSERT_NOT_NULL(ctx);
	bcf_unpack(ctx,BCF_UN_IND);
	
	hdr=(bcf_hdr_t*)R_ExternalPtrAddr(sexpheader);
	ASSERT_NOT_NULL(hdr);
	
	if(IS_CHARACTER(sexpgtidx))
		{
		const char* sn = CHAR(asChar(sexpgtidx));
		if(sn==NULL) {
			shuttle->error_flag  = 1;
			return;
			}
		sample_idx = bcf_hdr_id2int(hdr,BCF_DT_SAMPLE,sn);
		if(sample_idx<0) {
			BCF_WARNING("unknown sample %s",sn);
			shuttle->error_flag  = 1;
			return;
			}
        }
	else
		{
		sample_idx = asInteger(sexpgtidx);
		}
	
	if(sample_idx<0 || sample_idx>=bcf_hdr_nsamples(hdr)) {
		BCF_WARNING("bad sample index 0<=%d<%d",sample_idx,bcf_hdr_nsamples(hdr));
		shuttle->error_flag  = 1;
		return;
		}
	int j;
	int ngt = bcf_get_genotypes(hdr,ctx, &gt_arr, &ngt_arr);
	int nsmpl = bcf_hdr_nsamples(hdr);
	int max_ploidy = ngt/nsmpl;
     
      int32_t *ptr = gt_arr + max_ploidy*sample_idx;
      for (j=0; j<max_ploidy; j++)
          {
	         // if true, the sample has smaller ploidy
	         if ( ptr[j]==bcf_int32_vector_end ) break;
  
	         
  
  			 // is phased?
  			 if(bcf_gt_is_phased(ptr[j])) {
	         	shuttle->phased = 1;
	         	}
  
	         // the VCF 0-based allele index
	         int allele_index ;
	         
	         // missing allele
	         if ( bcf_gt_is_missing(ptr[j]) ) {
	         	allele_index = -1;
	         	}
	         else
	         	{
		        allele_index= bcf_gt_allele(ptr[j]);
		        }
		        
  			 if(shuttle->alleles != NULL) {
  			 	if(shuttle->allele_count + 1 >= shuttle->allele_capacity) {
  			 		shuttle->allele_capacity++;
  			 		shuttle->alleles =Realloc(shuttle->alleles,shuttle->allele_capacity,int);
  			 		ASSERT_NOT_NULL(shuttle->alleles);	
  			 		}
  			 	shuttle->alleles[shuttle->allele_count] = allele_index ;
  			 	}
  			shuttle->allele_count++;
  		    }
	
	UNPROTECT(nprotect);
	}

