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
#include "htslib/version.h"
#include "rbcf_version.h"


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
	tbx_t *tbx;
	hts_idx_t *idx;
	hts_itr_t *itr;
	bcf1_t* tmp_ctx; // tmp variant ctx for reading
	kstring_t tmp_line; // for reading line
	int query_failed;
	}RBcfFile,*RBcfFilePtr;



#ifndef RBCF_VERSION
	#define RBCF_VERSION "undefined"
#endif

SEXP RBCFGetVersion() {
	return mkString(RBCF_VERSION);
	}

/* version of htslib */
SEXP HtslibGetVersion() {
	return mkString(HTS_VERSION_TEXT);
	}

static void RBcfFileFree(final RBcfFilePtr ptr) {
	if (ptr==NULL) return;

	free(ptr->tmp_line.s);
	if ( ptr->tmp_ctx) bcf_destroy1(ptr->tmp_ctx);
	if ( ptr->itr ) hts_itr_destroy(ptr->itr);
	if ( ptr->tbx ) tbx_destroy(ptr->tbx);
    if ( ptr->idx ) hts_idx_destroy(ptr->idx);
	if ( ptr->fp ) hts_close(ptr->fp);	
	Free(ptr);
	}

/**
 * Close resources associated 
 */
SEXP RBcfFileClose(SEXP sexpFile) {
	PROTECT(sexpFile);
	SEXP fp = VECTOR_ELT(sexpFile,0);
	RBcfFilePtr p = (RBcfFilePtr)R_ExternalPtrAddr(fp);
	RBcfFileFree(p);
	R_ClearExternalPtr(fp);
	SEXP ext= ScalarLogical(1);
	UNPROTECT(1);
	return ext;
	}

void RBcfHeaderFinalizer(SEXP handler) {
	void* p=R_ExternalPtrAddr(handler);
	if(p==NULL) return;
	bcf_hdr_destroy((bcf_hdr_t*)p);
	}

static void RBcfFileFinalizer(SEXP handle)
	{
	RBcfFilePtr p = (RBcfFilePtr)R_ExternalPtrAddr(handle);
	RBcfFileFree(p);
	}




/**
 * Open BCF file
 */
SEXP RBcfFileOpen(SEXP Rfilename,SEXP sexpRequireIdx) {
	bcf_hdr_t* hdr = NULL;
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
	handler->query_failed = 0;
	
	handler->tmp_ctx = bcf_init1();
	if(handler->tmp_ctx ==NULL) {
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

	hdr =	bcf_hdr_read(handler->fp);
	if(hdr==NULL) {
		LOG("Cannot read header from %s.",filename);
		goto die;
		}

	SEXP sexpheader = PROTECT(R_MakeExternalPtr(hdr, R_NilValue, R_NilValue));nprotect++;
	R_RegisterCFinalizerEx(sexpheader,RBcfHeaderFinalizer, TRUE);

	SEXP sexpFile = PROTECT(R_MakeExternalPtr(handler, R_NilValue, R_NilValue));nprotect++;
	R_RegisterCFinalizerEx(sexpFile,RBcfFileFinalizer, TRUE);
	
	SEXP ext = PROTECT(allocVector(VECSXP, 3));nprotect++;
	SET_VECTOR_ELT(ext, 0, sexpFile);
	SET_VECTOR_ELT(ext, 1, sexpheader);
	SET_VECTOR_ELT(ext, 2, Rfilename);
    UNPROTECT(nprotect);
	return ext;
	
	die:
		LOG("error opening file");
		if(handler!=NULL) {
			RBcfFileFree(handler);
			}
		if(hdr!=NULL) bcf_hdr_destroy(hdr);
	UNPROTECT(nprotect);
	return R_NilValue;
	}




SEXP RBcfSeqNames(SEXP sexpFile) {
	int i, n=0;
	SEXP ext;
	int nprotect=0;
	PROTECT(sexpFile);nprotect++;
	
	RBcfFilePtr p = (RBcfFilePtr)R_ExternalPtrAddr(VECTOR_ELT(sexpFile,0));
	ASSERT_NOT_NULL(p);
	bcf_hdr_t* hdr =(bcf_hdr_t*)R_ExternalPtrAddr(VECTOR_ELT(sexpFile,1));
	ASSERT_NOT_NULL(hdr);
	
	const char **names;
	
	if(p->tbx) {
		names =  tbx_seqnames(p->tbx, &n);
		}
	else if( p->idx) {
		names = bcf_hdr_seqnames(hdr, &n);
		}
	else
		{
		BCF_WARNING("VCF File \"%s\" is not indexed.",CHAR(asChar(VECTOR_ELT(sexpFile,2))));
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

SEXP RBcfNSamples(SEXP sexpFile) {
	int n=0;
	PROTECT(sexpFile);
	bcf_hdr_t* hdr =(bcf_hdr_t*)R_ExternalPtrAddr(VECTOR_ELT(sexpFile,1));
	ASSERT_NOT_NULL(hdr);	
	n = (int)bcf_hdr_nsamples(hdr);
	UNPROTECT(1);
	return ScalarInteger(n);
	}

SEXP RBcfSamples(SEXP sexpFile) {
	int i;
	int nprotect=0;
	PROTECT(sexpFile);nprotect++;
	bcf_hdr_t* hdr =(bcf_hdr_t*)R_ExternalPtrAddr(VECTOR_ELT(sexpFile,1));
	ASSERT_NOT_NULL(hdr);
	
	SEXP sNames = PROTECT(allocVector(STRSXP,bcf_hdr_nsamples(hdr)));nprotect++;
	for(i=0;i< bcf_hdr_nsamples(hdr);++i) {
		SET_STRING_ELT(sNames, i , mkChar(hdr->samples[i]));
		}
	UNPROTECT(nprotect);
	return sNames;
	}

SEXP RBcfSampleAtIndex0(SEXP sexpFile,SEXP sexpIndex) {
	int nprotect=0;
	SEXP ext;
	PROTECT(sexpFile);nprotect++;
	bcf_hdr_t* hdr =(bcf_hdr_t*)R_ExternalPtrAddr(VECTOR_ELT(sexpFile,1));
	ASSERT_NOT_NULL(hdr);
	
	PROTECT(sexpIndex);++nprotect;
	
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

/** number of vcfheaderline with this type */
static int my_count_hdr_type(bcf_hdr_t* hdr,int type) {
	int count  = 0;
	int i = 0;
	for (i=0; i<hdr->nhrec; i++) {
        if ( hdr->hrec[i]->type!=type ) continue;
        count++;
        }
    return count;
	}

SEXP BcfFilterTable(SEXP sexpFile) {
	int  row_idx=0,hidx=0,nprotect=0;
	SEXP res,vID,vDesc,vRowNames;
	PROTECT(sexpFile);nprotect++;
	bcf_hdr_t* hdr =(bcf_hdr_t*)R_ExternalPtrAddr(VECTOR_ELT(sexpFile,1));
	ASSERT_NOT_NULL(hdr);
	
	final int header_type = BCF_HL_FLT;
	int nrows = my_count_hdr_type(hdr,header_type);
	
	PROTECT(res = Rf_allocVector(VECSXP,2));nprotect++;
	PROTECT(vID = Rf_allocVector(STRSXP,nrows)); nprotect++;
	PROTECT(vDesc = Rf_allocVector(STRSXP,nrows)); nprotect++;
	PROTECT(vRowNames = Rf_allocVector(STRSXP,nrows)); nprotect++;
	
	for(hidx=0; hidx <hdr->nhrec; hidx++) {
		bcf_hrec_t *hrec = hdr->hrec[hidx];
        if ( hrec->type!= header_type ) continue;
		/* ID */
			{
			int i = bcf_hrec_find_key(hrec,"ID");
			SEXP sexpCell = mkChar(i<0?".":hrec->vals[i]);
			SET_STRING_ELT(vID, row_idx,  sexpCell);
			SET_STRING_ELT(vRowNames, row_idx, sexpCell );
			}
		/* Description */
			{
			int i = bcf_hrec_find_key(hrec,"Description");
			SEXP sexpCell = mkChar(i<0?".":hrec->vals[i]);
			SET_STRING_ELT(vDesc,row_idx,sexpCell);
			}
		row_idx++;
		}

	/* set the columns of the table */
	SET_VECTOR_ELT(res, 0, vID);
  	SET_VECTOR_ELT(res, 1, vDesc);

	/* set the columns' name of the table */

	SEXP sNames = PROTECT(allocVector(STRSXP, 2));nprotect++;
	SET_STRING_ELT(sNames, 0, mkChar("ID"));
	SET_STRING_ELT(sNames, 1, mkChar("Description"));
	namesgets(res, sNames);
	
	// https://stackoverflow.com/questions/23547625
   SEXP cls = PROTECT(allocVector(STRSXP, 1));nprotect++;
   SET_STRING_ELT(cls, 0, mkChar("data.frame"));
   classgets(res, cls);
	
	setAttrib(res, R_RowNamesSymbol, vRowNames);
	
	UNPROTECT(nprotect);
	return res;
	}
	
static SEXP BcfInfoOrFormatTable(SEXP sexpFile,int header_type) {
	int  row_idx=0,hidx=0,nprotect=0;
	SEXP res,vID,vNumber,vType,vDesc,vRowNames;
	PROTECT(sexpFile);nprotect++;
	bcf_hdr_t* hdr =(bcf_hdr_t*)R_ExternalPtrAddr(VECTOR_ELT(sexpFile,1));
	ASSERT_NOT_NULL(hdr);
	
	int nrows = my_count_hdr_type(hdr,header_type);
	
	PROTECT(res = Rf_allocVector(VECSXP,4));nprotect++;
	PROTECT(vID = Rf_allocVector(STRSXP,nrows)); nprotect++;
	PROTECT(vDesc = Rf_allocVector(STRSXP,nrows)); nprotect++;
	PROTECT(vNumber = Rf_allocVector(STRSXP,nrows)); nprotect++;
	PROTECT(vType = Rf_allocVector(STRSXP,nrows)); nprotect++;
	PROTECT(vRowNames = Rf_allocVector(STRSXP,nrows)); nprotect++;
	
	for(hidx=0; hidx <hdr->nhrec; hidx++) {
		bcf_hrec_t *hrec = hdr->hrec[hidx];
        if ( hrec->type!= header_type ) continue;
		/* ID */
			{
			int i = bcf_hrec_find_key(hrec,"ID");
			SEXP sexpCell = mkChar(i<0?".":hrec->vals[i]);
			SET_STRING_ELT(vID, row_idx,  sexpCell);
			SET_STRING_ELT(vRowNames, row_idx, sexpCell );
			}
		/* Number */
			{
			int i = bcf_hrec_find_key(hrec,"Number");
			SEXP sexpCell = mkChar(i<0?".":hrec->vals[i]);
			SET_STRING_ELT(vNumber, row_idx,  sexpCell);
			}	
		/* Type */
			{
			int i = bcf_hrec_find_key(hrec,"Type");
			SEXP sexpCell = mkChar(i<0?".":hrec->vals[i]);
			SET_STRING_ELT(vType, row_idx,  sexpCell);
			}
		/* Description */
			{
			int i = bcf_hrec_find_key(hrec,"Description");
			SEXP sexpCell = mkChar(i<0?".":hrec->vals[i]);
			SET_STRING_ELT(vDesc,row_idx,sexpCell);
			}
		row_idx++;
		}

	/* set the columns of the table */
	SET_VECTOR_ELT(res, 0, vID);
	SET_VECTOR_ELT(res, 1, vNumber);
	SET_VECTOR_ELT(res, 2, vType);
  	SET_VECTOR_ELT(res, 3, vDesc);

	/* set the columns' name of the table */
	SEXP sNames = PROTECT(allocVector(STRSXP, 4));nprotect++;
	SET_STRING_ELT(sNames, 0, mkChar("ID"));
	SET_STRING_ELT(sNames, 1, mkChar("Number"));
	SET_STRING_ELT(sNames, 2, mkChar("Type"));
	SET_STRING_ELT(sNames, 3, mkChar("Description"));
	namesgets(res, sNames);
	
	// https://stackoverflow.com/questions/23547625
   SEXP cls = PROTECT(allocVector(STRSXP, 1));nprotect++;
   SET_STRING_ELT(cls, 0, mkChar("data.frame"));
   classgets(res, cls);
	
	setAttrib(res, R_RowNamesSymbol, vRowNames);
	
	UNPROTECT(nprotect);
	return res;
	}

SEXP BcfInfoTable(SEXP sexpFile) {
	return BcfInfoOrFormatTable(sexpFile,BCF_HL_INFO);
	}

SEXP BcfFormatTable(SEXP sexpFile) {
	return BcfInfoOrFormatTable(sexpFile,BCF_HL_FMT);
	}

SEXP RBcfHeaderDict(SEXP sexpFile) {
	int i,nprotect=0;
	SEXP res,vChrom,vSize,vRowNames;
	PROTECT(sexpFile);nprotect++;
	bcf_hdr_t* hdr =(bcf_hdr_t*)R_ExternalPtrAddr(VECTOR_ELT(sexpFile,1));
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


SEXP BcfConvertSampleToIndex0(SEXP sexpFile,SEXP sexpSample) {
	int nprotect=0;
	PROTECT(sexpFile);nprotect++;
	bcf_hdr_t* hdr =(bcf_hdr_t*)R_ExternalPtrAddr(VECTOR_ELT(sexpFile,1));
	ASSERT_NOT_NULL(hdr);
	const char* sn = CHAR(asChar(sexpSample));
	int sample_idx=-1;

	if(sn!=NULL) {
		sample_idx = bcf_hdr_id2int(hdr,BCF_DT_SAMPLE,sn);
       		}
	UNPROTECT(nprotect);
	return ScalarInteger(sample_idx);
	}


SEXP RBcfQueryRegion(SEXP sexpFile,SEXP sexpInterval) {
	int nprotect=0;
	SEXP ext;
	PROTECT(sexpFile);nprotect++;
	PROTECT(sexpInterval);nprotect++;
	
	RBcfFilePtr ptr=(RBcfFilePtr)R_ExternalPtrAddr(VECTOR_ELT(sexpFile,0));
	ASSERT_NOT_NULL(ptr);
	bcf_hdr_t* hdr =(bcf_hdr_t*)R_ExternalPtrAddr(VECTOR_ELT(sexpFile,1));
	ASSERT_NOT_NULL(hdr);
	const char* filename  = CHAR(asChar(VECTOR_ELT(sexpFile,2)));
	
	const char* interval= CHAR(asChar(sexpInterval));
	
	
	ASSERT_NOT_NULL(interval);
	if(ptr->itr) hts_itr_destroy(ptr->itr);
	ptr->itr = NULL;
	ptr->query_failed=0;
	
	if( ptr->tbx!=NULL) {
		ptr->itr = tbx_itr_querys(ptr->tbx,interval);
		}
	else if(ptr->idx!=NULL) {
		ptr->itr = bcf_itr_querys(ptr->idx,hdr,interval);
		}
	else
		{
		BCF_ERROR("Cannot query vcf file \"%s\" (no index available)",filename);
		}
	if(ptr->itr==NULL) {
		BCF_WARNING("Query failed \"%s\" for \"%s\".",interval, filename);
		ptr->query_failed = 1;
		}
	ext = ScalarLogical( ptr->query_failed==0 );
	UNPROTECT(nprotect);
	return ext;
	}


static void RBcf1Finalizer(SEXP handler) {
	bcf1_t* p=(bcf1_t*)R_ExternalPtrAddr(handler);
	if(p) bcf_destroy(p);
	}


SEXP RBcfNextLine(SEXP sexpFile) {
	int ret=0;
	int nprotect=0;
	bcf1_t *line = NULL;
	
	SEXP ext;
	PROTECT(sexpFile);nprotect++;
	RBcfFilePtr reader=(RBcfFilePtr)R_ExternalPtrAddr(VECTOR_ELT(sexpFile,0));
	ASSERT_NOT_NULL(reader);
	bcf_hdr_t* hdr =(bcf_hdr_t*)R_ExternalPtrAddr(VECTOR_ELT(sexpFile,1));
	ASSERT_NOT_NULL(hdr);
	
	if( reader->query_failed!=0) {
		line=NULL;
		}
	else if ( reader->fp->format.format==vcf ) {
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
        	ASSERT_NOT_NULL(hdr);
        	ASSERT_NOT_NULL(reader->tmp_ctx);
            ret = vcf_parse1(&(reader->tmp_line), hdr, reader->tmp_ctx);
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
        	ret = bcf_read1(reader->fp,hdr, reader->tmp_ctx);
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
	 	
	 	SEXP sexpctx = PROTECT(R_MakeExternalPtr(line, R_NilValue, R_NilValue));nprotect++;
	 	R_RegisterCFinalizerEx(sexpctx,RBcf1Finalizer, TRUE);
	 	
	 	ext = PROTECT(allocVector(VECSXP, 2));nprotect++;
		SET_VECTOR_ELT(ext, 0, VECTOR_ELT(sexpFile,1));
		SET_VECTOR_ELT(ext, 1, sexpctx);
	 	}
	else
		 {
		 ext = R_NilValue;		 
		 }
	
	UNPROTECT(nprotect);
	return ext;
	}

SEXP RBcfCtxRid(SEXP sexpCtx) {
	int nprotect=0;
	bcf1_t *ctx;
	SEXP ext;
	PROTECT(sexpCtx);nprotect++;
	ctx=(bcf1_t*)R_ExternalPtrAddr(VECTOR_ELT(sexpCtx,1));
	ASSERT_NOT_NULL(ctx);
	ext = ScalarInteger(ctx->rid);
	UNPROTECT(nprotect);
	return ext;
	}


SEXP RBcfCtxSeqName(SEXP sexpCtx) {
	int nprotect=0;
	bcf1_t *ctx = NULL;
	bcf_hdr_t* hdr = NULL;
	SEXP ext;
	PROTECT(sexpCtx);nprotect++;
	hdr = (bcf_hdr_t*)R_ExternalPtrAddr(VECTOR_ELT(sexpCtx,0));
	ctx = (bcf1_t*)R_ExternalPtrAddr(VECTOR_ELT(sexpCtx,1));
	ext = mkString(bcf_seqname(hdr,ctx));
	UNPROTECT(nprotect);
	return ext;
	}

SEXP RBcfCtxPos(SEXP sexpCtx) {
	int nprotect=0;
	bcf1_t *ctx;
	SEXP ext;
	PROTECT(sexpCtx);nprotect++;
	ctx = (bcf1_t*)R_ExternalPtrAddr(VECTOR_ELT(sexpCtx,1));
	ASSERT_NOT_NULL(ctx);
	ext = ScalarInteger(ctx->pos + 1);
	UNPROTECT(nprotect);
	return ext;
	}

SEXP RBcfCtxHasId(SEXP sexpCtx) {
	int nprotect=0;
	bcf1_t *ctx = NULL;
	SEXP ext;
	PROTECT(sexpCtx);nprotect++;
	ctx = (bcf1_t*)R_ExternalPtrAddr(VECTOR_ELT(sexpCtx,1));
	ASSERT_NOT_NULL(ctx);
	ext = ScalarLogical(ctx->d.id!=NULL);
	UNPROTECT(nprotect);
	return ext;
	}

SEXP RBcfCtxId(SEXP sexpCtx) {
	int nprotect=0;
	bcf1_t *ctx = NULL;
	SEXP ext;
	PROTECT(sexpCtx);nprotect++;
	ctx = (bcf1_t*)R_ExternalPtrAddr(VECTOR_ELT(sexpCtx,1));
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


SEXP RBcfCtxEnd(SEXP sexpCtx) {
	int nprotect=0;
	bcf1_t *ctx;
	SEXP ext;
	PROTECT(sexpCtx);nprotect++;
	ctx = (bcf1_t*)R_ExternalPtrAddr(VECTOR_ELT(sexpCtx,1));
	ASSERT_NOT_NULL(ctx);
	ext = ScalarInteger(ctx->pos /* 0 based */ + ctx->rlen);
	UNPROTECT(nprotect);
	return ext;
	}

SEXP RBcfCtxNAlleles(SEXP sexpCtx) {
	int nprotect=0;
	bcf1_t *ctx;
	SEXP ext;
	PROTECT(sexpCtx);nprotect++;
	ctx = (bcf1_t*)R_ExternalPtrAddr(VECTOR_ELT(sexpCtx,1));
	ASSERT_NOT_NULL(ctx);
	ext = ScalarInteger(ctx->n_allele);
	UNPROTECT(nprotect);
	return ext;
	}


SEXP RBcfCtxAlleles(SEXP sexpCtx) {
	int i=0;
	int nprotect=0;
	SEXP ext;
	PROTECT(sexpCtx);nprotect++;
	bcf1_t* ctx = (bcf1_t*)R_ExternalPtrAddr(VECTOR_ELT(sexpCtx,1));
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

SEXP RBcfCtxReference(SEXP sexpCtx) {
	int nprotect=0;
	bcf1_t *ctx;
	SEXP ext;
	PROTECT(sexpCtx);nprotect++;
	ctx = (bcf1_t*)R_ExternalPtrAddr(VECTOR_ELT(sexpCtx,1));
	bcf_unpack(ctx,BCF_UN_STR);
	ASSERT_NOT_NULL(ctx);
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

SEXP RBcfCtxAlternateAlleles(SEXP sexpCtx) {
	int nprotect=0;
	int i;
	bcf1_t *ctx;
	SEXP ext;
	PROTECT(sexpCtx);nprotect++;
	ctx = (bcf1_t*)R_ExternalPtrAddr(VECTOR_ELT(sexpCtx,1));
	ASSERT_NOT_NULL(ctx);
	bcf_unpack(ctx,BCF_UN_STR);
	
	ext = PROTECT(allocVector(STRSXP,(ctx->n_allele<1?0:ctx->n_allele-1)));nprotect++;
	for(i=1;i< ctx->n_allele;++i) {
		SET_STRING_ELT(ext, i-1, mkChar(ctx->d.allele[i]));
		}
	UNPROTECT(nprotect);
	return ext;
	}
	
SEXP RBcfCtxHasQual(SEXP sexpCtx) {
	int nprotect=0;
	bcf1_t *ctx;
	SEXP ext;
	PROTECT(sexpCtx);nprotect++;
	ctx = (bcf1_t*)R_ExternalPtrAddr(VECTOR_ELT(sexpCtx,1));
	ASSERT_NOT_NULL(ctx);
	ext = ScalarLogical(!bcf_float_is_missing(ctx->qual));
	UNPROTECT(nprotect);
	return ext;
	}

SEXP RBcfCtxQual(SEXP sexpCtx) {
	int nprotect=0;
	bcf1_t *ctx;
	SEXP ext;
	PROTECT(sexpCtx);nprotect++;
	ctx = (bcf1_t*)R_ExternalPtrAddr(VECTOR_ELT(sexpCtx,1));
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


SEXP RBcfCtxFiltered(SEXP sexpCtx) {
	int nprotect=0;
	bcf1_t *ctx;
	bcf_hdr_t* hdr;
	SEXP ext;
	PROTECT(sexpCtx);nprotect++;
	hdr = (bcf_hdr_t*)R_ExternalPtrAddr(VECTOR_ELT(sexpCtx,0));
	ctx = (bcf1_t*)R_ExternalPtrAddr(VECTOR_ELT(sexpCtx,1));
	ASSERT_NOT_NULL(ctx);
	ASSERT_NOT_NULL(hdr);
	bcf_unpack(ctx,BCF_UN_FLT);
	ext = ScalarLogical(!bcf_has_filter(hdr,ctx,"PASS"));
	UNPROTECT(nprotect);
	return ext;
	}

SEXP RBcfCtxFilters(SEXP sexpCtx) {
	int nprotect=0;
	int i=0;
	bcf1_t *ctx;
	bcf_hdr_t* hdr;
	SEXP ext;
	PROTECT(sexpCtx);nprotect++;
	
	hdr = (bcf_hdr_t*)R_ExternalPtrAddr(VECTOR_ELT(sexpCtx,0));
	ctx = (bcf1_t*)R_ExternalPtrAddr(VECTOR_ELT(sexpCtx,1));
	ASSERT_NOT_NULL(ctx);
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

SEXP RBcfCtxVariantTypes(SEXP sexpCtx) {
	int nprotect=0;
	bcf1_t *ctx;
	SEXP ext;
	PROTECT(sexpCtx);nprotect++;
	ctx = (bcf1_t*)R_ExternalPtrAddr(VECTOR_ELT(sexpCtx,1));
	ext = ScalarInteger(bcf_get_variant_types(ctx));
	UNPROTECT(nprotect);
	return ext;
	}
	
SEXP RBcfCtxVariantIsSnp(SEXP sexpCtx) {
	int nprotect=0;
	bcf1_t *ctx;
	SEXP ext;
	PROTECT(sexpCtx);nprotect++;
	ctx = (bcf1_t*)R_ExternalPtrAddr(VECTOR_ELT(sexpCtx,1));
	ASSERT_NOT_NULL(ctx);
	ext = ScalarLogical(bcf_is_snp(ctx));
	UNPROTECT(nprotect);
	return ext;
	}

SEXP RBcfCtxVariantMaxPloidy(SEXP sexpCtx) {
	int nprotect=0;
	bcf1_t *ctx;
	bcf_hdr_t* hdr;
	int32_t *gt_arr = NULL, ngt_arr = 0;
	SEXP ext;
	PROTECT(sexpCtx);nprotect++;

	hdr = (bcf_hdr_t*)R_ExternalPtrAddr(VECTOR_ELT(sexpCtx,0));
	ctx = (bcf1_t*)R_ExternalPtrAddr(VECTOR_ELT(sexpCtx,1));
	ASSERT_NOT_NULL(ctx);
	bcf_unpack(ctx,BCF_UN_IND);

	
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



SEXP VariantGetGenotype(SEXP sexpCtx,SEXP sexpgtidx) {
	int nprotect=0;
	SEXP ext;
	PROTECT(sexpCtx);nprotect++;
	
	bcf_hdr_t* hdr = (bcf_hdr_t*)R_ExternalPtrAddr(VECTOR_ELT(sexpCtx,0));
	int sample_idx = -1;
	
	if(IS_CHARACTER(sexpgtidx))
		{
		const char* sn = CHAR(asChar(sexpgtidx));
		if(sn==NULL) {
			UNPROTECT(nprotect);
			return R_NilValue;
			}
		sample_idx = bcf_hdr_id2int(hdr,BCF_DT_SAMPLE,sn);
		if(sample_idx<0) {
			BCF_WARNING("unknown sample %s",sn);
			UNPROTECT(nprotect);
			return R_NilValue;
			}
        }
	else
		{
		sample_idx = asInteger(sexpgtidx);
		}
	
	if(sample_idx<0 || sample_idx>=bcf_hdr_nsamples(hdr)) {
		BCF_WARNING("bad sample index 0<=%d<%d",sample_idx,bcf_hdr_nsamples(hdr));
		ext= R_NilValue;
		}
	else
		{
		ext = PROTECT(allocVector(VECSXP, 3));nprotect++;
		SET_VECTOR_ELT(ext, 0, VECTOR_ELT(sexpCtx,0));//hdr
		SET_VECTOR_ELT(ext, 1, VECTOR_ELT(sexpCtx,1));//variant
		SET_VECTOR_ELT(ext, 2, ScalarInteger(sample_idx));
		}
	UNPROTECT(nprotect);
	return ext;
	}

static void scanGenotype(bcf_hdr_t* hdr,bcf1_t *ctx,int sample_idx, struct GenotypeShuttle* shuttle) {
	int32_t *gt_arr = NULL, ngt_arr = 0;

	bcf_unpack(ctx,BCF_UN_IND);
	
	
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
	
	}

SEXP RBcfCtxVariantGtAllelesIndexes0(SEXP sexpGt) {
	int nprotect=0;
	struct GenotypeShuttle shuttle;
	SEXP ext;
	PROTECT(sexpGt);nprotect++;
	
	bcf_hdr_t*	hdr = (bcf_hdr_t*)R_ExternalPtrAddr(VECTOR_ELT(sexpGt,0));
	bcf1_t* ctx = (bcf1_t*)R_ExternalPtrAddr(VECTOR_ELT(sexpGt,1));
	int sample_index = asInteger(VECTOR_ELT(sexpGt,2));
	

	memset((void*)&shuttle,0,sizeof(struct GenotypeShuttle));
	shuttle.allele_capacity=10;
	shuttle.alleles =Calloc(shuttle.allele_capacity,int);
	scanGenotype(hdr,ctx,sample_index,&shuttle);
	if(shuttle.error_flag)
		{
		ext = R_NilValue;
		}
	else
		{
		PROTECT(ext = Rf_allocVector(INTSXP,shuttle.allele_count)); nprotect++;
		for(int i=0;i< shuttle.allele_count;i++) {
			INTEGER(ext)[i] = shuttle.alleles[i];
			}
		}
	Free(shuttle.alleles);
	UNPROTECT(nprotect);
	return ext;
	}

SEXP RBcfCtxVariantGtPhased(SEXP sexpGt) {
	int nprotect=0;
	
	struct GenotypeShuttle shuttle;
	PROTECT(sexpGt);nprotect++;
	bcf_hdr_t*	hdr = (bcf_hdr_t*)R_ExternalPtrAddr(VECTOR_ELT(sexpGt,0));
	bcf1_t* ctx = (bcf1_t*)R_ExternalPtrAddr(VECTOR_ELT(sexpGt,1));
	int sample_index = asInteger(VECTOR_ELT(sexpGt,2));
	
	memset((void*)&shuttle,0,sizeof(struct GenotypeShuttle));
	scanGenotype(hdr,ctx,sample_index,&shuttle);
	UNPROTECT(nprotect);
	return ScalarLogical(shuttle.phased==1);
	}


SEXP VariantHasAttribute(SEXP sexpCtx,SEXP sexpatt) {
	int nprotect=0;
	SEXP ext;
	PROTECT(sexpCtx);nprotect++;
	bcf_hdr_t*	hdr = (bcf_hdr_t*)R_ExternalPtrAddr(VECTOR_ELT(sexpCtx,0));
	bcf1_t* ctx = (bcf1_t*)R_ExternalPtrAddr(VECTOR_ELT(sexpCtx,1));
	const char* att = CHAR(asChar(sexpatt));
	ASSERT_NOT_NULL(ctx);
	ASSERT_NOT_NULL(hdr);
	ASSERT_NOT_NULL(att);
	int tag_id = bcf_hdr_id2int(hdr,BCF_DT_ID,att);
	
	if(!(bcf_hdr_idinfo_exists(hdr,BCF_HL_INFO,tag_id))) {
		ext = ScalarLogical(0);
		}
	else 
		{
		bcf_info_t* info = bcf_get_info_id(ctx,  tag_id);
		ext = ScalarLogical(info!=NULL);
		}
	UNPROTECT(nprotect);
	return ext;
	}
	
SEXP VariantGetInfoKeySet(SEXP sexpCtx) {
	int i,nprotect=0;
	PROTECT(sexpCtx);nprotect++;
	bcf_hdr_t*	hdr = (bcf_hdr_t*)R_ExternalPtrAddr(VECTOR_ELT(sexpCtx,0));
	bcf1_t* ctx = (bcf1_t*)R_ExternalPtrAddr(VECTOR_ELT(sexpCtx,1));

	ASSERT_NOT_NULL(ctx);
	ASSERT_NOT_NULL(hdr);
	
	bcf_unpack(ctx,BCF_UN_INFO);
	int count=0;
	for(i=0;i< ctx->n_info ;++i) {
		bcf_info_t *z = &(ctx->d.info[i]);
        if ( !z->vptr ) continue;
		count++;
		}
	
	
	SEXP ext = PROTECT(allocVector(STRSXP,count));nprotect++;
	count=0;//reset
	for(i=0;i< ctx->n_info ;++i) {
		bcf_info_t *z = &(ctx->d.info[i]);
        if ( !z->vptr ) continue;
        const char * id = hdr->id[BCF_DT_ID][z->key].key;
		SET_STRING_ELT(ext, count , mkChar(id));
		count++;
		}
	
	UNPROTECT(nprotect);
	return ext;
	}



SEXP VariantGetFormatKeySet(SEXP sexpCtx) {
	int i,nprotect=0;
	PROTECT(sexpCtx);nprotect++;
	bcf_hdr_t*	hdr = (bcf_hdr_t*)R_ExternalPtrAddr(VECTOR_ELT(sexpCtx,0));
	bcf1_t* ctx = (bcf1_t*)R_ExternalPtrAddr(VECTOR_ELT(sexpCtx,1));

	ASSERT_NOT_NULL(ctx);
	ASSERT_NOT_NULL(hdr);
	bcf_unpack(ctx,BCF_UN_FMT);
	int count=0;
	bcf_fmt_t *fmt = ctx->d.fmt;
	for(i=0;i< ctx->n_fmt ;++i) {
        if ( !fmt[i].p || fmt[i].id<0) continue;
		count++;
		}
	
	
	SEXP ext = PROTECT(allocVector(STRSXP,count));nprotect++;
	count=0;//reset
	for(i=0;i< ctx->n_fmt ;++i) {
        if ( !fmt[i].p || fmt[i].id<0 ) continue;
        const char * id = hdr->id[BCF_DT_ID][fmt[i].id].key;
		SET_STRING_ELT(ext, count , mkChar(id));
		count++;
		}
	
	UNPROTECT(nprotect);
	return ext;
	}

SEXP VariantGetAttribute(SEXP sexpCtx,SEXP sexpatt) {
	int nprotect=0;
	
	PROTECT(sexpCtx);nprotect++;
	bcf_hdr_t*	hdr = (bcf_hdr_t*)R_ExternalPtrAddr(VECTOR_ELT(sexpCtx,0));
	bcf1_t* ctx = (bcf1_t*)R_ExternalPtrAddr(VECTOR_ELT(sexpCtx,1));
	const char* att = CHAR(asChar(sexpatt));
	ASSERT_NOT_NULL(ctx);
	ASSERT_NOT_NULL(hdr);
	ASSERT_NOT_NULL(att);
	
	int tag_id = bcf_hdr_id2int(hdr, BCF_DT_ID, att);
	

	if(!(bcf_hdr_idinfo_exists(hdr,BCF_HL_INFO,tag_id))) {
		SEXP ext = PROTECT(allocVector(STRSXP,0));nprotect++;
		UNPROTECT(nprotect);
		return ext;
		}
	
	if ( !(ctx->unpacked & BCF_UN_INFO) ) bcf_unpack(ctx, BCF_UN_INFO);
	
	
	uint32_t type = bcf_hdr_id2type(hdr,BCF_HL_INFO, tag_id);
	
	bcf_info_t* info = bcf_get_info_id(ctx,  tag_id);
	if(info==NULL) {
		SEXP ext = PROTECT(allocVector(STRSXP,0));nprotect++;
		UNPROTECT(nprotect);
		return ext;
		}
	switch(type) {
		case BCF_HT_STR: {
			SEXP ext;
			char* dst=NULL;
			int ndst=0;
			int ret=bcf_get_info_string(hdr,ctx,att,(void**)&dst,&ndst);
			if(ret<0 || dst==NULL) {
				ext = R_NilValue;
				}
			else
				{
				ext= mkString(dst);
				}
			free(dst);
			UNPROTECT(nprotect);
			return ext;
			break;
			}
		case BCF_HT_INT: {
			SEXP ext;
			int ndst=0;
			int32_t* dst=NULL;
			int ret=bcf_get_info_int32(hdr,ctx,att,(void**)&dst,&ndst);
			if(ret<0 || dst==NULL ) {
				PROTECT(ext = Rf_allocVector(INTSXP,0)); nprotect++;
				}
			else
				{
				PROTECT(ext = Rf_allocVector(INTSXP,ndst)); nprotect++;
				        for(int i=0;i< ndst;i++) {
				               	INTEGER(ext)[i] = dst[i];
				                }
		
				}
			free(dst);
			UNPROTECT(nprotect);
			return ext;
			break;
			}
		case BCF_HT_REAL: {
			SEXP ext;
			int ndst=0;
			float* dst=NULL;
			int ret=bcf_get_info_float(hdr,ctx,att,(void**)&dst,&ndst);
			if(ret<0 || dst==NULL ) {
				PROTECT(ext = Rf_allocVector(REALSXP,0)); nprotect++;
				}
			else
				{
				PROTECT(ext = Rf_allocVector(REALSXP,ndst)); nprotect++;
				        for(int i=0;i< ndst;i++) {
				               	REAL(ext)[i] = dst[i];
				                }
		
				}
			free(dst);
			UNPROTECT(nprotect);
			return ext;
			break;
			}
		case BCF_HT_FLAG: {
			int is_set = bcf_get_info_flag(hdr,ctx,att,NULL,NULL)==1;
			SEXP ext = ScalarLogical(is_set==1);
			UNPROTECT(nprotect);
			return ext;
			break;
			}
		default:
			{
			BCF_WARNING("unknown INFO type");
			SEXP ext = PROTECT(allocVector(STRSXP,0));nprotect++;
			UNPROTECT(nprotect);
			return ext;
			break;
			}
		}
	return R_NilValue;
	}

#ifdef XXX
SEXP GenotypeAttributes(SEXP sexpGt,SEXP sexpatt) {
	int nprotect=0;
	
	PROTECT(sexpGt);nprotect++;
	bcf_hdr_t*	hdr = (bcf_hdr_t*)R_ExternalPtrAddr(VECTOR_ELT(sexpGt,0));
	bcf1_t* ctx = (bcf1_t*)R_ExternalPtrAddr(VECTOR_ELT(sexpGt,1));
	int sample_index = asInteger(VECTOR_ELT(sexpGt,2));
	const char* att = CHAR(asChar(sexpatt));
	ASSERT_NOT_NULL(ctx);
	ASSERT_NOT_NULL(hdr);
	ASSERT_NOT_NULL(att);
	
	
	int tag_id = bcf_hdr_id2int(hdr, BCF_DT_ID, att);
	

	if(!(bcf_hdr_idinfo_exists(hdr,BCF_HL_FMT,tag_id))) {
		SEXP ext = PROTECT(allocVector(STRSXP,0));nprotect++;
		UNPROTECT(nprotect);
		return ext;
		}
	
	int nsmpl = bcf_hdr_nsamples(args->hdr);
    int len = bcf_hdr_id2length(args->hdr,BCF_HL_FMT,fmt->id);
	
	
	if ( !(ctx->unpacked & BCF_UN_FMT) ) bcf_unpack(ctx, BCF_UN_FMT);
	
	
	uint32_t type = bcf_hdr_id2type(hdr,BCF_HL_FMT, tag_id);
	bcf_fmt_t * fmt = bcf_get_fmt_id(ctx,  tag_id);
	if(fmt==NULL) {
		SEXP ext = PROTECT(allocVector(STRSXP,0));nprotect++;
		UNPROTECT(nprotect);
		return ext;
		}
	
	LOG("format/%s has type = %d",att,type);
	
	switch(type) {
		case BCF_HT_STR: {
			SEXP ext;
			char* dst=NULL;
			int ndst=0;
			int nret = bcf_get_format_char(hdr,ctx,att,dst,&ndst);
			if(nret<0 || dst==NULL) {
				ext = R_NilValue;
				}
			else
				{
				ext= mkString(dst[sample_index]);
				}
			free(dst);
			UNPROTECT(nprotect);
			return ext;
			break;
			}
		case BCF_HT_INT: {
			BCF_ERROR("HT_INT %s",att);
			SEXP ext;
			int ndst=0;
			int32_t* dst=NULL;
			int ret=bcf_get_format_int32(hdr,ctx,att,(void**)&dst,&ndst);
			if(ret<0 || dst==NULL ) {
				PROTECT(ext = Rf_allocVector(INTSXP,0)); nprotect++;
				}
			else
				{
				PROTECT(ext = Rf_allocVector(INTSXP,ndst)); nprotect++;
				        for(int i=0;i< ndst;i++) {
				               	INTEGER(ext)[i] = dst[i];
				                }
		
				}
			free(dst);
			UNPROTECT(nprotect);
			return ext;
			break;
			}
		case BCF_HT_REAL: {
			BCF_ERROR("HT_REAL %s",att);
			SEXP ext;
			int ndst=0;
			float* dst=NULL;
			int ret=bcf_get_format_float(hdr,ctx,att,(void**)&dst,&ndst);
			if(ret<0 || dst==NULL ) {
				PROTECT(ext = Rf_allocVector(REALSXP,0)); nprotect++;
				}
			else
				{
				PROTECT(ext = Rf_allocVector(REALSXP,ndst)); nprotect++;
				        for(int i=0;i< ndst;i++) {
				               	REAL(ext)[i] = dst[i];
				                }
		
				}
			free(dst);
			UNPROTECT(nprotect);
			return ext;
			break;
			}
		case BCF_HT_FLAG: {
			BCF_ERROR("HT_FLAG %s",att);
			int is_set = bcf_get_format_flag(hdr,ctx,att,NULL,NULL)==1;
			SEXP ext = ScalarLogical(is_set==1);
			UNPROTECT(nprotect);
			return ext;
			break;
			}
		default:
			{
			BCF_WARNING("unknown FORMAT type");
			SEXP ext = PROTECT(allocVector(STRSXP,0));nprotect++;
			UNPROTECT(nprotect);
			return ext;
			break;
			}
		}
	return R_NilValue;
	}
#endif	

