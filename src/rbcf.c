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
#include <math.h>
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

#define VEP_CSQ_KEY "CSQ"
#define VEP_FORMAT "Format: "
#define WHERENL REprintf("[%s:%d] ",__FILE__,__LINE__)
#define LOG(FormatLiteral, ...) do { WHERENL;REprintf("[LOG]" FormatLiteral "\n", ##__VA_ARGS__);} while(0)
#define BCF_WARNING(FormatLiteral, ...) do { WHERENL;REprintf("[WARN]" FormatLiteral "\n", ##__VA_ARGS__);} while(0)
#define BCF_ERROR(FormatLiteral, ...) do { WHERENL;REprintf("[FATAL]" FormatLiteral "\n", ##__VA_ARGS__);error("Fatal Error");} while(0)
#define ASSERT_NOT_NULL(ptr) do {if(ptr==NULL) BCF_ERROR("NULL pointer exception"); } while(0)
#define IF_NULL_UNPROTECT_AND_RETURN_NULL(P) do {if(isNull(P)) {BCF_WARNING("variable <"  #P  "> is NULL.");UNPROTECT(nprotect); return R_NilValue;} } while(0)
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
	int is_writer;
	}RBcfFile,*RBcfFilePtr;

typedef struct array_of_strings {
	char* dup;
	char** tokens;
	int count;
	}Tokens,*TokensPtr;
	
static void trim(char *s) {
	char * p = s;
	size_t l = strlen(p);

	while(l>0 && isspace(p[l - 1])) p[--l] = 0;
	while(l>0 && *p && isspace(*p)) ++p, --l;

	memmove(s, p, l + 1);
	}

static int endsWith(const char *str, const char *suffix) {
    if (str==NULL || suffix==NULL) return 0;
    size_t lenstr = strlen(str);
    size_t lensuffix = strlen(suffix);
    if (lensuffix >  lenstr) return 0;
    return strncmp(str + lenstr - lensuffix, suffix, lensuffix) == 0;
    }

#define MAKE_ARRAY_API(NAME,TYPE) \
typedef struct NAME##array_t {\
	TYPE* data;\
	size_t capacity;\
	size_t size;\
	} NAME##Array,* NAME##ArrayPtr;\
NAME##ArrayPtr NAME##ArrayNew() {\
	NAME##ArrayPtr ptr = (NAME##ArrayPtr)calloc(1,sizeof(NAME##Array));\
	ptr->capacity = 10L;\
	ptr->data = (TYPE*)calloc(ptr->capacity,sizeof(TYPE));\
	return ptr;\
	}\
void NAME##ArrayFree(NAME##ArrayPtr ptr) {\
	if(ptr==NULL) return;\
	free(ptr->data);\
	free(ptr);\
	}\
NAME##ArrayPtr NAME##ArrayPush(NAME##ArrayPtr ptr,TYPE v) {\
	if(ptr->size>=ptr->capacity) {\
		ptr->capacity+=10;\
		ptr->data = (TYPE*)realloc(ptr->data,sizeof(TYPE)*(ptr->capacity));\
		if(ptr->data==NULL) {BCF_WARNING("Out of memory");return NULL;}\
		}\
	ptr->data[ptr->size]=v;\
	ptr->size++;\
	return ptr;\
	}


MAKE_ARRAY_API(Int32,int32_t)
MAKE_ARRAY_API(Float,float)

static TokensPtr NewTokensPtr(const char* str,char delim) {
	TokensPtr dest = (TokensPtr)malloc(sizeof(Tokens));
	if(dest==NULL) {
		BCF_WARNING("Out of memory");
		return NULL;
		}
	memset((void*)dest,0,sizeof(Tokens));
	dest->dup = strdup(str);
	char* p = dest->dup;
	char* prev = p;
	for(;;) {
		if(*p==delim || *p==0) {
			dest->tokens = (char**)realloc(dest->tokens,sizeof(char*)*(1+dest->count));
			if(dest->tokens==NULL) {
				BCF_WARNING("Out of memory");
				return NULL;
				}
			dest->tokens[dest->count] = prev;
			dest->count++;
			if(*p==0) break;
			*p=0;
			prev=p+1;
			}
		++p;
		}
	return dest;
	}

void FreeTokensPtr(TokensPtr dest) {
	if(dest==NULL) return;
	free(dest->dup);
	free(dest->tokens);
	free(dest);
	}
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
	if(sexpFile!=R_NilValue) {
		SEXP fp = VECTOR_ELT(sexpFile,0);
		RBcfFilePtr p = (RBcfFilePtr)R_ExternalPtrAddr(fp);
		RBcfFileFree(p);
		R_ClearExternalPtr(fp);
	}
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
	
	IF_NULL_UNPROTECT_AND_RETURN_NULL(Rfilename);
	IF_NULL_UNPROTECT_AND_RETURN_NULL(sexpRequireIdx);
	
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
	handler->is_writer = 0;
	
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
/**
 * Open BCF file
 */
SEXP RBcfNewWriter(SEXP sexpIn,SEXP Rfilename) {
	int nprotect=0;
	RBcfFilePtr handler =  NULL;
	
	


	PROTECT(sexpIn);nprotect++;
	PROTECT(Rfilename);nprotect++;
	IF_NULL_UNPROTECT_AND_RETURN_NULL(sexpIn);
	IF_NULL_UNPROTECT_AND_RETURN_NULL(Rfilename);	
	
	RBcfFilePtr p = (RBcfFilePtr)R_ExternalPtrAddr(VECTOR_ELT(sexpIn,0));
	ASSERT_NOT_NULL(p);
	SEXP sexpheader = VECTOR_ELT(sexpIn,1);
	bcf_hdr_t* hdr =(bcf_hdr_t*)R_ExternalPtrAddr(sexpheader);
	ASSERT_NOT_NULL(hdr);

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
	handler->is_writer = 1;
	
	char modew[8];
	strcpy(modew, "w");

	if(endsWith(filename,".bcf")) strcat(modew, "b"); // uncompressed BCF
	else if(endsWith(filename,".gz") || endsWith(filename,".vcfz")) strcat(modew, "z");// compressed VCF
    	
	handler->fp = hts_open(filename,modew);
	if (handler->fp==NULL) {
		BCF_WARNING("Failed to open writer %s.\n",filename);
		goto die;
		}

	if ( bcf_hdr_write(handler->fp,hdr)!=0 ) {
		BCF_WARNING("Failed to write header for writer %s.\n",filename);
		goto die;
		}

	SEXP sexpFile = PROTECT(R_MakeExternalPtr(handler, R_NilValue, R_NilValue));nprotect++;
	R_RegisterCFinalizerEx(sexpFile,RBcfFileFinalizer, TRUE);
	
	SEXP ext = PROTECT(allocVector(VECSXP, 3));nprotect++;
	SET_VECTOR_ELT(ext, 0, sexpFile);
	SET_VECTOR_ELT(ext, 1, sexpheader);
	SET_VECTOR_ELT(ext, 2, Rfilename);
    	UNPROTECT(nprotect);
	return ext;


	die:
		BCF_WARNING("error opening writer");
		if(handler!=NULL) {
			RBcfFileFree(handler);
			}
		if(hdr!=NULL) bcf_hdr_destroy(hdr);
	UNPROTECT(nprotect);
	return R_NilValue;
	}

/**
 * Write variant to output vcf
 */
SEXP RBcfFileWriteCtx(SEXP sexpOut,SEXP sexpCtx) {
	int nprotect=0;
	PROTECT(sexpOut);nprotect++;
	PROTECT(sexpCtx);nprotect++;
	IF_NULL_UNPROTECT_AND_RETURN_NULL(sexpOut);
	IF_NULL_UNPROTECT_AND_RETURN_NULL(sexpCtx);
	
	RBcfFilePtr p = (RBcfFilePtr)R_ExternalPtrAddr(VECTOR_ELT(sexpOut,0));
	ASSERT_NOT_NULL(p);
	bcf_hdr_t* hdr1 = (bcf_hdr_t*)R_ExternalPtrAddr(VECTOR_ELT(sexpOut,1));


	if(p->is_writer!=1) {
		BCF_WARNING("Cannot save to a reader");
		UNPROTECT(nprotect);
		return ScalarLogical(0);
		}

	bcf_hdr_t* hdr2 = (bcf_hdr_t*)R_ExternalPtrAddr(VECTOR_ELT(sexpCtx,0));

	if(hdr1!=hdr2) {
		BCF_WARNING("header from reader is not the same as writer.");
		UNPROTECT(nprotect);
		return ScalarLogical(0);
		}

	bcf1_t* ctx = (bcf1_t*)R_ExternalPtrAddr(VECTOR_ELT(sexpCtx,1));
	int ret = bcf_write1(p->fp, hdr2, ctx);
	UNPROTECT(nprotect);
	return ScalarLogical(ret==0);
	}

SEXP RBcfSeqNames(SEXP sexpFile) {
	int i, n=0;
	SEXP ext;
	int nprotect=0;
	PROTECT(sexpFile);nprotect++;
	IF_NULL_UNPROTECT_AND_RETURN_NULL(sexpFile);
	
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
	int nprotect=0;
	PROTECT(sexpFile);nprotect++;
	IF_NULL_UNPROTECT_AND_RETURN_NULL(sexpFile);
	
	bcf_hdr_t* hdr =(bcf_hdr_t*)R_ExternalPtrAddr(VECTOR_ELT(sexpFile,1));
	ASSERT_NOT_NULL(hdr);	
	n = (int)bcf_hdr_nsamples(hdr);
	UNPROTECT(nprotect);
	return ScalarInteger(n);
	}

SEXP RBcfSamples(SEXP sexpFile) {
	int i;
	int nprotect=0;
	PROTECT(sexpFile);nprotect++;
	IF_NULL_UNPROTECT_AND_RETURN_NULL(sexpFile);
	
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
	IF_NULL_UNPROTECT_AND_RETURN_NULL(sexpFile);
	
	bcf_hdr_t* hdr =(bcf_hdr_t*)R_ExternalPtrAddr(VECTOR_ELT(sexpFile,1));
	ASSERT_NOT_NULL(hdr);
	
	PROTECT(sexpIndex);++nprotect;
	
	int idx = asInteger(sexpIndex);
	
	if(idx<0 || idx >= bcf_hdr_nsamples(hdr)) {
		BCF_WARNING("sample idx out of range %d/%d\n",idx,bcf_hdr_nsamples(hdr));
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
	IF_NULL_UNPROTECT_AND_RETURN_NULL(sexpFile);
	
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
	IF_NULL_UNPROTECT_AND_RETURN_NULL(sexpFile);
	
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
	IF_NULL_UNPROTECT_AND_RETURN_NULL(sexpFile);
	
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
	PROTECT(sexpSample);nprotect++;
	if(isNull(sexpFile) || isNull(sexpSample)) {
		UNPROTECT(nprotect);
		return ScalarInteger(-1);
		}
	
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
	if(isNull(sexpFile) || isNull(sexpInterval)) {
		BCF_WARNING("null parameter");
		UNPROTECT(nprotect);
		return ScalarLogical(0);
		}
	
	
	RBcfFilePtr ptr=(RBcfFilePtr)R_ExternalPtrAddr(VECTOR_ELT(sexpFile,0));
	ASSERT_NOT_NULL(ptr);
	bcf_hdr_t* hdr =(bcf_hdr_t*)R_ExternalPtrAddr(VECTOR_ELT(sexpFile,1));
	ASSERT_NOT_NULL(hdr);
	const char* filename  = CHAR(asChar(VECTOR_ELT(sexpFile,2)));
	
	const char* interval= CHAR(asChar(sexpInterval));
	ASSERT_NOT_NULL(interval);

	if( ptr->is_writer == 1 ) {
		BCF_WARNING("cannot query a writer");
		UNPROTECT(nprotect);
		return ScalarLogical(0);
		}	


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
	IF_NULL_UNPROTECT_AND_RETURN_NULL(sexpFile);
	
	RBcfFilePtr reader=(RBcfFilePtr)R_ExternalPtrAddr(VECTOR_ELT(sexpFile,0));
	ASSERT_NOT_NULL(reader);
	bcf_hdr_t* hdr =(bcf_hdr_t*)R_ExternalPtrAddr(VECTOR_ELT(sexpFile,1));
	ASSERT_NOT_NULL(hdr);
	
	if( reader->is_writer == 1 ) {
		BCF_WARNING("try to read from a writer");
		UNPROTECT(nprotect);
		return R_NilValue;
	}

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
	IF_NULL_UNPROTECT_AND_RETURN_NULL(sexpCtx);
	
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
	IF_NULL_UNPROTECT_AND_RETURN_NULL(sexpCtx);
	
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
	IF_NULL_UNPROTECT_AND_RETURN_NULL(sexpCtx);
	
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
	IF_NULL_UNPROTECT_AND_RETURN_NULL(sexpCtx);
	
	ctx = (bcf1_t*)R_ExternalPtrAddr(VECTOR_ELT(sexpCtx,1));
	ASSERT_NOT_NULL(ctx);
	bcf_unpack(ctx,BCF_UN_STR);
	ext = ScalarLogical(ctx->d.id!=NULL && strcmp(ctx->d.id,".")!=0 && strcmp(ctx->d.id,"")!=0);
	UNPROTECT(nprotect);
	return ext;
	}

SEXP RBcfCtxId(SEXP sexpCtx) {
	int nprotect=0;
	bcf1_t *ctx = NULL;
	SEXP ext;
	PROTECT(sexpCtx);nprotect++;
	IF_NULL_UNPROTECT_AND_RETURN_NULL(sexpCtx);
	
	ctx = (bcf1_t*)R_ExternalPtrAddr(VECTOR_ELT(sexpCtx,1));
	ASSERT_NOT_NULL(ctx);
	bcf_unpack(ctx,BCF_UN_STR);
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
	IF_NULL_UNPROTECT_AND_RETURN_NULL(sexpCtx);
	
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
	IF_NULL_UNPROTECT_AND_RETURN_NULL(sexpCtx);
	
	ctx = (bcf1_t*)R_ExternalPtrAddr(VECTOR_ELT(sexpCtx,1));
	ASSERT_NOT_NULL(ctx);
	bcf_unpack(ctx,BCF_UN_STR);
	ext = ScalarInteger(ctx->n_allele);
	UNPROTECT(nprotect);
	return ext;
	}


SEXP RBcfCtxAlleles(SEXP sexpCtx) {
	int i=0;
	int nprotect=0;
	SEXP ext;
	PROTECT(sexpCtx);nprotect++;
	IF_NULL_UNPROTECT_AND_RETURN_NULL(sexpCtx);
	
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
	IF_NULL_UNPROTECT_AND_RETURN_NULL(sexpCtx);
	
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
	IF_NULL_UNPROTECT_AND_RETURN_NULL(sexpCtx);
	
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
	IF_NULL_UNPROTECT_AND_RETURN_NULL(sexpCtx);
	
	ctx = (bcf1_t*)R_ExternalPtrAddr(VECTOR_ELT(sexpCtx,1));
	ASSERT_NOT_NULL(ctx);
	ext = ScalarLogical(bcf_float_is_missing(ctx->qual)==0 && !isnan(ctx->qual));
	UNPROTECT(nprotect);
	return ext;
	}

SEXP RBcfCtxQual(SEXP sexpCtx) {
	int nprotect=0;
	bcf1_t *ctx;
	SEXP ext;
	PROTECT(sexpCtx);nprotect++;
	IF_NULL_UNPROTECT_AND_RETURN_NULL(sexpCtx);
	
	ctx = (bcf1_t*)R_ExternalPtrAddr(VECTOR_ELT(sexpCtx,1));
	ASSERT_NOT_NULL(ctx);
	ext = ScalarReal(ctx->qual);
	UNPROTECT(nprotect);
	return ext;
	}


SEXP RBcfCtxFiltered(SEXP sexpCtx) {
	int nprotect=0;
	bcf1_t *ctx;
	bcf_hdr_t* hdr;
	SEXP ext;
	PROTECT(sexpCtx);nprotect++;
	IF_NULL_UNPROTECT_AND_RETURN_NULL(sexpCtx);
	
	hdr = (bcf_hdr_t*)R_ExternalPtrAddr(VECTOR_ELT(sexpCtx,0));
	ctx = (bcf1_t*)R_ExternalPtrAddr(VECTOR_ELT(sexpCtx,1));
	ASSERT_NOT_NULL(ctx);
	ASSERT_NOT_NULL(hdr);
	bcf_unpack(ctx,BCF_UN_FLT);
	ext = ScalarLogical(bcf_has_filter(hdr,ctx,"PASS")==0);
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
	IF_NULL_UNPROTECT_AND_RETURN_NULL(sexpCtx);
	
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

SEXP VariantHasFilter(SEXP sexpCtx,SEXP sexpFilterName) {
	int nprotect=0;
	bcf1_t *ctx;
	bcf_hdr_t* hdr;
	SEXP ext;
	PROTECT(sexpCtx);nprotect++;
	PROTECT(sexpFilterName);nprotect++;
	if(isNull(sexpCtx) || isNull(sexpFilterName)) {
		UNPROTECT(nprotect);
		BCF_WARNING("NULL parameters");
		return ScalarLogical(0);
		}
	
	hdr = (bcf_hdr_t*)R_ExternalPtrAddr(VECTOR_ELT(sexpCtx,0));
	ctx = (bcf1_t*)R_ExternalPtrAddr(VECTOR_ELT(sexpCtx,1));
	const char* filterName= CHAR(asChar(sexpFilterName));
	
	ASSERT_NOT_NULL(ctx);
	ASSERT_NOT_NULL(hdr);
	
	bcf_unpack(ctx,BCF_UN_FLT);
	
	if(filterName!=NULL && bcf_has_filter(hdr,ctx,(char*)filterName)) {
		ext = ScalarLogical(1);
		}
	else
		{
		ext = ScalarLogical(0);
		}
	UNPROTECT(nprotect);
	return ext;
	}
	
SEXP RBcfCtxVariantTypes(SEXP sexpCtx) {
	int nprotect=0;
	bcf1_t *ctx;
	SEXP ext;
	PROTECT(sexpCtx);nprotect++;
	if(isNull(sexpCtx)) {
		UNPROTECT(nprotect);
		BCF_WARNING("NULL parameters");
		return ScalarLogical(0);
		}
	
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
	if(isNull(sexpCtx)) {
		UNPROTECT(nprotect);
		BCF_WARNING("NULL parameters");
		return ScalarLogical(0);
		}
	
	ctx = (bcf1_t*)R_ExternalPtrAddr(VECTOR_ELT(sexpCtx,1));
	ASSERT_NOT_NULL(ctx);
	ext = ScalarLogical(bcf_is_snp(ctx)!=0);
	UNPROTECT(nprotect);
	return ext;
	}

SEXP VariantNSamples(SEXP sexpCtx) {
	int nprotect=0;
	PROTECT(sexpCtx);nprotect++;
	if(isNull(sexpCtx)) {
		UNPROTECT(nprotect);
		BCF_WARNING("NULL parameters");
		return ScalarInteger(0);
		}
	
	bcf_hdr_t* hdr = (bcf_hdr_t*)R_ExternalPtrAddr(VECTOR_ELT(sexpCtx,0));
   	SEXP ext = ScalarInteger(bcf_hdr_nsamples(hdr));
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
	if(isNull(sexpCtx)) {
		UNPROTECT(nprotect);
		BCF_WARNING("NULL parameters");
		return ScalarInteger(0);
		}
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
	PROTECT(sexpgtidx);nprotect++;
	IF_NULL_UNPROTECT_AND_RETURN_NULL(sexpCtx);
	IF_NULL_UNPROTECT_AND_RETURN_NULL(sexpgtidx);
	
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
	IF_NULL_UNPROTECT_AND_RETURN_NULL(sexpGt);
	
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

SEXP RBcfCtxVariantAllGtAllelesIndexes0(SEXP sexpCtx) {
	int nprotect=0;
	int i,j,k;
	PROTECT(sexpCtx);nprotect++;
	
	bcf_hdr_t*	hdr = (bcf_hdr_t*)R_ExternalPtrAddr(VECTOR_ELT(sexpCtx,0));
	bcf1_t* ctx = (bcf1_t*)R_ExternalPtrAddr(VECTOR_ELT(sexpCtx,1));
	
	// Total number of samples
	int nsmpl = bcf_hdr_nsamples(hdr);

	// Identify the max-ploidy
	int32_t *gt_arr = NULL, ngt_arr = 0;
	int ngt = bcf_get_genotypes(hdr,ctx, &gt_arr, &ngt_arr);
	int max_ploidy = 0;
	if ( ngt > 0 ) {
		max_ploidy = ngt / nsmpl;
	}
  
	SEXP ext = PROTECT(allocVector(INTSXP,ngt));nprotect++;	

	for (i=0; i < nsmpl; i ++) {
		for (j=0; j<max_ploidy; j++) {
		  k = i * max_ploidy + j;
		  // if true, the sample has smaller ploidy
		  if (gt_arr[k]==bcf_int32_vector_end ) {
		     INTEGER(ext)[k] = -2;
		  }
		  else if ( bcf_gt_is_missing(gt_arr[k]) ) {
		    INTEGER(ext)[k] = -1;
		  } else {
		    INTEGER(ext)[k] = bcf_gt_allele(gt_arr[k]);
		  }
		}
   }

	UNPROTECT(nprotect);
	return ext;
	}

SEXP RBcfCtxVariantAllGtStrings(SEXP sexpCtx) {
	int nprotect=0;
	int i,j,k;
	PROTECT(sexpCtx);nprotect++;
	
	bcf_hdr_t*	hdr = (bcf_hdr_t*)R_ExternalPtrAddr(VECTOR_ELT(sexpCtx,0));
	bcf1_t* ctx = (bcf1_t*)R_ExternalPtrAddr(VECTOR_ELT(sexpCtx,1));
	
	// Total number of samples
	int nsmpl = bcf_hdr_nsamples(hdr);

	// Identify the max-ploidy
	int32_t *gt_arr = NULL, ngt_arr = 0;
	int ngt = bcf_get_genotypes(hdr,ctx, &gt_arr, &ngt_arr);
	int max_ploidy = 0;
	if ( ngt > 0 ) {
		max_ploidy = ngt / nsmpl;
	}
	char *buf = malloc(max_ploidy * 2 * sizeof(char));
  
	SEXP ext = PROTECT(allocVector(STRSXP,nsmpl));nprotect++;

	for (i=0; i < nsmpl; i ++) {
		for (j=0; j<max_ploidy; j++) {
		  k = i * max_ploidy + j;
		  // if true, the sample has smaller ploidy
		  if (gt_arr[k]==bcf_int32_vector_end ) {
		     buf[j*2] = '.';
		  }
		  else if ( bcf_gt_is_missing(gt_arr[k]) ) {
		     buf[j*2] = '.';
		  } else {
				sprintf(buf + (j*2), "%d", bcf_gt_allele(gt_arr[k]));
		  }
		  if (j + 1 >= max_ploidy) {
					 buf[j*2 + 1] = '\0';
		  } else if (bcf_gt_is_phased(gt_arr[k]) ){
					 buf[j*2 + 1] = '|';
			} else {
					 buf[j*2 + 1] = '/';
		  }
		}
		SET_STRING_ELT(ext, i, mkChar(buf));
   }
	free(buf);

	UNPROTECT(nprotect);
	return ext;
	}

SEXP GenotypeSample(SEXP sexpGt) {
        int nprotect=0;
		SEXP ext;
        PROTECT(sexpGt);nprotect++;
        IF_NULL_UNPROTECT_AND_RETURN_NULL(sexpGt);
        
        bcf_hdr_t*	hdr = (bcf_hdr_t*)R_ExternalPtrAddr(VECTOR_ELT(sexpGt,0));
        int sample_index = asInteger(VECTOR_ELT(sexpGt,2));

        if( sample_index <0 || sample_index >= bcf_hdr_nsamples(hdr)) {
                BCF_WARNING("sample idx out of range %d/%d\n",sample_index,bcf_hdr_nsamples(hdr));
                ext = R_NilValue;
                }
        else
                {
                ext =  mkString(hdr->samples[sample_index]);
	        }

	UNPROTECT(nprotect);
        return ext;
	}

SEXP RBcfCtxVariantGtPhased(SEXP sexpGt) {
	int nprotect=0;
	
	struct GenotypeShuttle shuttle;
	PROTECT(sexpGt);nprotect++;
	if(isNull(sexpGt)) {
		BCF_WARNING("parameter is null");
		UNPROTECT(nprotect);
		return ScalarLogical(0);
	}
	
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
	PROTECT(sexpatt);nprotect++;
	
	if(isNull(sexpCtx) || isNull(sexpatt)) {
			BCF_WARNING("parameter is null");
			UNPROTECT(nprotect);
			return ScalarLogical(0);
		}
	
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
	IF_NULL_UNPROTECT_AND_RETURN_NULL(sexpCtx);
	
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
	IF_NULL_UNPROTECT_AND_RETURN_NULL(sexpCtx);
	
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

SEXP VariantVepTable(SEXP sexpCtx) {
	int nprotect=0;
	SEXP ext = R_NilValue;
	char* csq_str =NULL;
	PROTECT(sexpCtx);nprotect++;
	IF_NULL_UNPROTECT_AND_RETURN_NULL(sexpCtx);
	
	bcf_hdr_t*	hdr = (bcf_hdr_t*)R_ExternalPtrAddr(VECTOR_ELT(sexpCtx,0));
	bcf1_t* ctx = (bcf1_t*)R_ExternalPtrAddr(VECTOR_ELT(sexpCtx,1));

	/** check variant has VEP tag */
	int tag_id = bcf_hdr_id2int(hdr, BCF_DT_ID, VEP_CSQ_KEY);
	if(!(bcf_hdr_idinfo_exists(hdr,BCF_HL_INFO,tag_id))) goto theend;

	/** search VEP_CSQ_KEY in header */
	char* format = NULL;
	
	for(int hidx=0; hidx <hdr->nhrec; hidx++) {
		bcf_hrec_t *hrec = hdr->hrec[hidx];
        if ( hrec->type!= BCF_HL_INFO) continue;
        int i = bcf_hrec_find_key(hrec,"ID");
		if(i<0 || strcmp(hrec->vals[i],VEP_CSQ_KEY)!=0) continue;
        i = bcf_hrec_find_key(hrec,"Description");
		if(i<0)  {
			BCF_WARNING("Cannot find Description, in INFO/"VEP_CSQ_KEY);
			break;
			}
		char* format_p = strstr(hrec->vals[i],VEP_FORMAT);
		if(format_p==NULL) {
			BCF_WARNING("Cannot find " VEP_FORMAT " in INFO/"VEP_CSQ_KEY);
			break;
			}
		format_p += strlen(VEP_FORMAT);
		format = format_p;
		break;
		}
	if(format==NULL) goto theend;

	
	bcf_unpack(ctx, BCF_UN_INFO);
	bcf_info_t* info = bcf_get_info_id(ctx,  tag_id);
	if(info==NULL) goto theend;
	
	
	int ndst=0;
	int ret=bcf_get_info_string(hdr,ctx,VEP_CSQ_KEY,(void**)& csq_str,&ndst);
	if(ret<0 ||  csq_str==NULL) goto theend;

	
	TokensPtr fmtTokens = NewTokensPtr(format,'|');
	TokensPtr transcriptTokens = NewTokensPtr(csq_str,',');
	
	/** alloc table */
	ext = PROTECT(Rf_allocVector(VECSXP,fmtTokens->count));nprotect++;
	SEXP colNames = PROTECT(allocVector(STRSXP, fmtTokens->count));nprotect++;
	SEXP rowNames = PROTECT(Rf_allocVector(STRSXP,transcriptTokens->count)); nprotect++;
	/* create columns */
	SEXP* columns=(SEXP*)malloc(sizeof(SEXP)*fmtTokens->count);
	ASSERT_NOT_NULL(columns);
	for(int x=0;x<fmtTokens->count;++x)
		{
		columns[x]=PROTECT(Rf_allocVector(STRSXP,transcriptTokens->count));nprotect++;
		SET_VECTOR_ELT(ext, x, columns[x]);
		// set header
		SET_STRING_ELT(colNames, x, mkChar(fmtTokens->tokens[x]));
		}
	

	/* scan each transcript */
	for(int y=0;y< transcriptTokens->count;++y) {
		const char* transcript = transcriptTokens->tokens[y];
		TokensPtr components = NewTokensPtr(transcript,'|');
		for(int x=0;x< fmtTokens->count;++x)
			{
			SET_STRING_ELT(columns[x],y,mkChar(x<components->count?components->tokens[x]:""));
			}
		FreeTokensPtr(components);
		//set name for this row
		char tmp[20];
		sprintf(tmp,"%d",(y+1));
		SET_STRING_ELT(rowNames, y, mkChar(tmp));

		}
	/* set the columns' name of the table */	
	namesgets(ext, colNames);
   // https://stackoverflow.com/questions/23547625
   SEXP cls = PROTECT(allocVector(STRSXP, 1));nprotect++;
   SET_STRING_ELT(cls, 0, mkChar("data.frame"));
   classgets(ext, cls);
   setAttrib(ext, R_RowNamesSymbol, rowNames);
	
	free(columns);
	FreeTokensPtr(fmtTokens);
	FreeTokensPtr(transcriptTokens);
	
	theend:
		UNPROTECT(nprotect);
		return ext;
	}

#define SNPEFF_ANN_KEY "ANN"

SEXP VariantSnpEffTable(SEXP sexpCtx) {
	int nprotect=0;
	SEXP ext = R_NilValue;
	char* csq_str =NULL;
	PROTECT(sexpCtx);nprotect++;
	IF_NULL_UNPROTECT_AND_RETURN_NULL(sexpCtx);
	
	bcf_hdr_t*	hdr = (bcf_hdr_t*)R_ExternalPtrAddr(VECTOR_ELT(sexpCtx,0));
	bcf1_t* ctx = (bcf1_t*)R_ExternalPtrAddr(VECTOR_ELT(sexpCtx,1));

	/** check variant has ANN tag */
	int tag_id = bcf_hdr_id2int(hdr, BCF_DT_ID, SNPEFF_ANN_KEY);
	if(!(bcf_hdr_idinfo_exists(hdr,BCF_HL_INFO,tag_id))) goto theend;

	/** search SNPEFF_ANN_KEY in header */
	char* format = NULL;
	
	for(int hidx=0; hidx <hdr->nhrec; hidx++) {
		bcf_hrec_t *hrec = hdr->hrec[hidx];
        if ( hrec->type!= BCF_HL_INFO) continue;
        int i = bcf_hrec_find_key(hrec,"ID");
		if(i<0 || strcmp(hrec->vals[i],SNPEFF_ANN_KEY)!=0) continue;
        i = bcf_hrec_find_key(hrec,"Description");
		if(i<0)  {
			BCF_WARNING("Cannot find Description, in INFO/"SNPEFF_ANN_KEY);
			break;
			}
		char* format_p = strchr(hrec->vals[i],'\'');
		if(format_p==NULL) {
			BCF_WARNING("Cannot find \' in INFO/"SNPEFF_ANN_KEY);
			break;
			}
		format = format_p+1;
		break;
		}
	if(format==NULL) goto theend;

	
	bcf_unpack(ctx, BCF_UN_INFO);
	bcf_info_t* info = bcf_get_info_id(ctx,  tag_id);
	if(info==NULL) goto theend;
	
	
	int ndst=0;
	int ret=bcf_get_info_string(hdr,ctx,SNPEFF_ANN_KEY,(void**)& csq_str,&ndst);
	if(ret<0 ||  csq_str==NULL) goto theend;

	
	TokensPtr fmtTokens = NewTokensPtr(format,'|');
	for(int i=0;i< fmtTokens->count;i++) {
		trim(fmtTokens->tokens[i]);
		char* p= strchr(fmtTokens->tokens[i],'\''); //last contains quote : ARNINGS / INFO'">
		if(p!=NULL && *(p+1)==0) *p=0;
		}
	
	TokensPtr transcriptTokens = NewTokensPtr(csq_str,',');
	
	/** alloc table */
	ext = PROTECT(Rf_allocVector(VECSXP,fmtTokens->count));nprotect++;
	SEXP colNames = PROTECT(allocVector(STRSXP, fmtTokens->count));nprotect++;
	SEXP rowNames = PROTECT(Rf_allocVector(STRSXP,transcriptTokens->count)); nprotect++;
	/* create columns */
	SEXP* columns=(SEXP*)malloc(sizeof(SEXP)*fmtTokens->count);
	ASSERT_NOT_NULL(columns);
	for(int x=0;x<fmtTokens->count;++x)
		{
		columns[x]=PROTECT(Rf_allocVector(STRSXP,transcriptTokens->count));nprotect++;
		SET_VECTOR_ELT(ext, x, columns[x]);
		// set header
		SET_STRING_ELT(colNames, x, mkChar(fmtTokens->tokens[x]));
		}
	

	/* scan each transcript */
	for(int y=0;y< transcriptTokens->count;++y) {
		const char* transcript = transcriptTokens->tokens[y];
		TokensPtr components = NewTokensPtr(transcript,'|');
		for(int x=0;x< fmtTokens->count;++x)
			{
			SET_STRING_ELT(columns[x],y,mkChar(x<components->count?components->tokens[x]:""));
			}
		FreeTokensPtr(components);
		//set name for this row
		char tmp[20];
		sprintf(tmp,"%d",(y+1));
		SET_STRING_ELT(rowNames, y, mkChar(tmp));

		}
	/* set the columns' name of the table */	
	namesgets(ext, colNames);
   // https://stackoverflow.com/questions/23547625
   SEXP cls = PROTECT(allocVector(STRSXP, 1));nprotect++;
   SET_STRING_ELT(cls, 0, mkChar("data.frame"));
   classgets(ext, cls);
   setAttrib(ext, R_RowNamesSymbol, rowNames);
	
	free(columns);
	FreeTokensPtr(fmtTokens);
	FreeTokensPtr(transcriptTokens);
	
	theend:
		UNPROTECT(nprotect);
		return ext;
	}

SEXP VariantStringAttribute(SEXP sexpCtx,SEXP sexpatt) {
	int nprotect=0;
	SEXP ext= R_NilValue;
	PROTECT(sexpCtx);nprotect++;
	PROTECT(sexpatt);nprotect++;
	IF_NULL_UNPROTECT_AND_RETURN_NULL(sexpCtx);
	IF_NULL_UNPROTECT_AND_RETURN_NULL(sexpatt);
	
	bcf_hdr_t*	hdr = (bcf_hdr_t*)R_ExternalPtrAddr(VECTOR_ELT(sexpCtx,0));
	bcf1_t* ctx = (bcf1_t*)R_ExternalPtrAddr(VECTOR_ELT(sexpCtx,1));
	const char* att = CHAR(asChar(sexpatt));
	ASSERT_NOT_NULL(ctx);
	ASSERT_NOT_NULL(hdr);
	ASSERT_NOT_NULL(att);
	
	bcf_unpack(ctx, BCF_UN_INFO);
	
	char* dst=NULL;
	int ndst=0;
	int ret=bcf_get_info_string(hdr,ctx,att,(void**)&dst,&ndst);
	if(!(ret<0 || dst==NULL)) {
		ext= mkString(dst);
		}
	free(dst);
	UNPROTECT(nprotect);
	return ext;
	}

SEXP VariantIntAttribute(SEXP sexpCtx,SEXP sexpatt) {
	int nprotect=0;
	PROTECT(sexpCtx);nprotect++;
	PROTECT(sexpatt);nprotect++;
	IF_NULL_UNPROTECT_AND_RETURN_NULL(sexpCtx);
	IF_NULL_UNPROTECT_AND_RETURN_NULL(sexpatt);
	
	bcf_hdr_t*	hdr = (bcf_hdr_t*)R_ExternalPtrAddr(VECTOR_ELT(sexpCtx,0));
	bcf1_t* ctx = (bcf1_t*)R_ExternalPtrAddr(VECTOR_ELT(sexpCtx,1));
	const char* att = CHAR(asChar(sexpatt));
	ASSERT_NOT_NULL(ctx);
	ASSERT_NOT_NULL(hdr);
	ASSERT_NOT_NULL(att);
	
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
	}

SEXP VariantFloatAttribute(SEXP sexpCtx,SEXP sexpatt) {
	int nprotect=0;
	PROTECT(sexpCtx);nprotect++;
	PROTECT(sexpatt);nprotect++;
	IF_NULL_UNPROTECT_AND_RETURN_NULL(sexpCtx);
	IF_NULL_UNPROTECT_AND_RETURN_NULL(sexpatt);
	
	bcf_hdr_t*	hdr = (bcf_hdr_t*)R_ExternalPtrAddr(VECTOR_ELT(sexpCtx,0));
	bcf1_t* ctx = (bcf1_t*)R_ExternalPtrAddr(VECTOR_ELT(sexpCtx,1));
	const char* att = CHAR(asChar(sexpatt));
	ASSERT_NOT_NULL(ctx);
	ASSERT_NOT_NULL(hdr);
	ASSERT_NOT_NULL(att);
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
	}
		
SEXP VariantFlagAttribute(SEXP sexpCtx,SEXP sexpatt) {
	int nprotect=0;
	PROTECT(sexpCtx);nprotect++;
	PROTECT(sexpatt);nprotect++;
	IF_NULL_UNPROTECT_AND_RETURN_NULL(sexpCtx);
	IF_NULL_UNPROTECT_AND_RETURN_NULL(sexpatt);
	
	bcf_hdr_t*	hdr = (bcf_hdr_t*)R_ExternalPtrAddr(VECTOR_ELT(sexpCtx,0));
	bcf1_t* ctx = (bcf1_t*)R_ExternalPtrAddr(VECTOR_ELT(sexpCtx,1));
	const char* att = CHAR(asChar(sexpatt));
	ASSERT_NOT_NULL(ctx);
	ASSERT_NOT_NULL(hdr);
	ASSERT_NOT_NULL(att);
	
	int is_set = bcf_get_info_flag(hdr,ctx,att,NULL,NULL)==1;
	SEXP ext = ScalarLogical(is_set==1);
	UNPROTECT(nprotect);
	return ext;		
	}




SEXP GenotypeStringAttribute(SEXP sexpGt,SEXP sexpatt) {
	int nprotect=0;	
	SEXP ext= R_NilValue;
	PROTECT(sexpGt);nprotect++;
	PROTECT(sexpatt);nprotect++;
	IF_NULL_UNPROTECT_AND_RETURN_NULL(sexpGt);
	IF_NULL_UNPROTECT_AND_RETURN_NULL(sexpatt);
	
	bcf_hdr_t*	hdr = (bcf_hdr_t*)R_ExternalPtrAddr(VECTOR_ELT(sexpGt,0));
	bcf1_t* ctx = (bcf1_t*)R_ExternalPtrAddr(VECTOR_ELT(sexpGt,1));
	int sample_index = asInteger(VECTOR_ELT(sexpGt,2));
	const char* att = CHAR(asChar(sexpatt));
	ASSERT_NOT_NULL(ctx);
	ASSERT_NOT_NULL(hdr);
	ASSERT_NOT_NULL(att);
	bcf_unpack(ctx, BCF_UN_FMT);
	int ndst=0;
	char** dst=NULL;
	int ret= bcf_get_format_string(hdr,ctx,att,&dst,&ndst);
	if(ret>=0 ) {
		if(dst!=NULL && sample_index < ndst && dst[sample_index]!=NULL)
			{
			ext = mkString(dst[sample_index]);
			}
		free(dst[0]);
		}
	free(dst);
	UNPROTECT(nprotect);
	return ext;
	}

SEXP GenotypeInt32Attribute(SEXP sexpGt,SEXP sexpatt) {
	int nprotect=0;
	
	PROTECT(sexpGt);nprotect++;
	PROTECT(sexpatt);nprotect++;
	IF_NULL_UNPROTECT_AND_RETURN_NULL(sexpGt);
	IF_NULL_UNPROTECT_AND_RETURN_NULL(sexpatt);
	
	bcf_hdr_t*	hdr = (bcf_hdr_t*)R_ExternalPtrAddr(VECTOR_ELT(sexpGt,0));
	bcf1_t* ctx = (bcf1_t*)R_ExternalPtrAddr(VECTOR_ELT(sexpGt,1));
	int sample_index = asInteger(VECTOR_ELT(sexpGt,2));
	const char* att = CHAR(asChar(sexpatt));
	ASSERT_NOT_NULL(ctx);
	ASSERT_NOT_NULL(hdr);
	ASSERT_NOT_NULL(att);
	Int32ArrayPtr array = Int32ArrayNew();


	bcf_unpack(ctx, BCF_UN_FMT);
	int tag_id = bcf_hdr_id2int(hdr, BCF_DT_ID, att);
	

	if(bcf_hdr_idinfo_exists(hdr,BCF_HL_FMT,tag_id)) {
		int ndst=0;
		int32_t* dst=NULL;
		int ret=bcf_get_format_int32(hdr,ctx,att,(void**)&dst,&ndst);
		if(ret>=0 && dst!=NULL ) {
			int per_sample = ndst/bcf_hdr_nsamples(hdr);
			int32_t *ptr = dst + sample_index*per_sample;
			for (int j=0; j< per_sample ; j++) {
			    int32_t val=ptr[j];
			    if(val==bcf_int32_vector_end ) break;
			    //IGNORE if(val==bcf_int32_missing ) call->qsum[j] = 0;
			    Int32ArrayPush(array,val);
			    }
			}
		free(dst);
		}

	SEXP ext = PROTECT(allocVector(INTSXP,array->size));nprotect++;	
 	for(int i=0;i< array->size;i++) {
	       	INTEGER(ext)[i] = array->data[i];
		}

	Int32ArrayFree(array);
	
	UNPROTECT(nprotect);
	return ext;
	}



SEXP GenotypeFloatAttribute(SEXP sexpGt,SEXP sexpatt) {
	int nprotect=0;
	
	PROTECT(sexpGt);nprotect++;
	PROTECT(sexpatt);nprotect++;
	IF_NULL_UNPROTECT_AND_RETURN_NULL(sexpGt);
	IF_NULL_UNPROTECT_AND_RETURN_NULL(sexpatt);
	
	bcf_hdr_t*	hdr = (bcf_hdr_t*)R_ExternalPtrAddr(VECTOR_ELT(sexpGt,0));
	bcf1_t* ctx = (bcf1_t*)R_ExternalPtrAddr(VECTOR_ELT(sexpGt,1));
	int sample_index = asInteger(VECTOR_ELT(sexpGt,2));
	const char* att = CHAR(asChar(sexpatt));
	ASSERT_NOT_NULL(ctx);
	ASSERT_NOT_NULL(hdr);
	ASSERT_NOT_NULL(att);
	FloatArrayPtr array = FloatArrayNew();


	bcf_unpack(ctx, BCF_UN_FMT);
	int tag_id = bcf_hdr_id2int(hdr, BCF_DT_ID, att);
	

	if(bcf_hdr_idinfo_exists(hdr,BCF_HL_FMT,tag_id)) {
		int ndst=0;
		float* dst=NULL;
		int ret=bcf_get_format_float(hdr,ctx,att,(void**)&dst,&ndst);
		if(ret>=0 && dst!=NULL ) {
			int per_sample = ndst/bcf_hdr_nsamples(hdr);
			float *ptr = dst + sample_index*per_sample;
			for (int j=0; j< per_sample ; j++) {
			    float val=ptr[j];
			    if(val==bcf_float_vector_end ) break;
			    //IGNORE if(val==bcf_float_missing ) call->qsum[j] = 0;
			    FloatArrayPush(array,val);
			    }
			}
		free(dst);
		}

	SEXP ext = PROTECT(allocVector(REALSXP,array->size));nprotect++;	
 	for(int i=0;i< array->size;i++) {
	       	REAL(ext)[i] = array->data[i];
		}

	FloatArrayFree(array);
	
	UNPROTECT(nprotect);
	return ext;
	}


SEXP VariantGenotypesFlagAttribute(SEXP sexpCtx,SEXP sexpatt) {
	int nprotect=0;
	PROTECT(sexpCtx);nprotect++;
	bcf_hdr_t*	hdr = (bcf_hdr_t*)R_ExternalPtrAddr(VECTOR_ELT(sexpCtx,0));
	bcf1_t* ctx = (bcf1_t*)R_ExternalPtrAddr(VECTOR_ELT(sexpCtx,1));
	const char* att = CHAR(asChar(sexpatt));
	ASSERT_NOT_NULL(ctx);
	ASSERT_NOT_NULL(hdr);
	ASSERT_NOT_NULL(att);
	Int32ArrayPtr array = Int32ArrayNew();


	bcf_unpack(ctx, BCF_UN_FMT);
	int tag_id = bcf_hdr_id2int(hdr, BCF_DT_ID, att);

	if(bcf_hdr_idinfo_exists(hdr,BCF_HL_FMT,tag_id)) {
		int ndst=0;
		int32_t* dst=NULL;
		int ret=bcf_get_format_int32(hdr,ctx,att,(void**)&dst,&ndst);
		if(ret>=0 && dst!=NULL ) {
			for (int j=0; j< ndst; j++) {
			    int32_t val=dst[j];
			    if(val==bcf_int32_vector_end ) break;
			    //IGNORE if(val==bcf_float_missing ) call->qsum[j] = 0;
			    Int32ArrayPush(array,val == 0 ? 0 : 1);
			}
		}
		free(dst);
	}

	SEXP ext = PROTECT(allocVector(INTSXP,array->size));nprotect++;	
 	for(int i=0;i< array->size;i++) {
	       	INTEGER(ext)[i] = array->data[i];
		}

	Int32ArrayFree(array);
	
	UNPROTECT(nprotect);
	return ext;
	}


SEXP VariantGenotypesInt32Attribute(SEXP sexpCtx,SEXP sexpatt) {
	int nprotect=0;
	PROTECT(sexpCtx);nprotect++;
	bcf_hdr_t*	hdr = (bcf_hdr_t*)R_ExternalPtrAddr(VECTOR_ELT(sexpCtx,0));
	bcf1_t* ctx = (bcf1_t*)R_ExternalPtrAddr(VECTOR_ELT(sexpCtx,1));
	const char* att = CHAR(asChar(sexpatt));
	ASSERT_NOT_NULL(ctx);
	ASSERT_NOT_NULL(hdr);
	ASSERT_NOT_NULL(att);
	Int32ArrayPtr array = Int32ArrayNew();


	bcf_unpack(ctx, BCF_UN_FMT);
	int tag_id = bcf_hdr_id2int(hdr, BCF_DT_ID, att);
	

	if(bcf_hdr_idinfo_exists(hdr,BCF_HL_FMT,tag_id)) {
		int ndst=0;
		int32_t* dst=NULL;
		int ret=bcf_get_format_int32(hdr,ctx,att,(void**)&dst,&ndst);
		if(ret>=0 && dst!=NULL ) {
			for (int j=0; j< ndst; j++) {
			    int32_t val=dst[j];
			    if(val==bcf_int32_vector_end ) break;
			    //IGNORE if(val==bcf_float_missing ) call->qsum[j] = 0;
			    Int32ArrayPush(array,val);
			}
		}
		free(dst);
	}

	SEXP ext = PROTECT(allocVector(INTSXP,array->size));nprotect++;	
 	for(int i=0;i< array->size;i++) {
	       	INTEGER(ext)[i] = array->data[i];
		}

	Int32ArrayFree(array);
	
	UNPROTECT(nprotect);
	return ext;
	}


SEXP VariantGenotypesFloatAttribute(SEXP sexpCtx,SEXP sexpatt) {
	int nprotect=0;
	PROTECT(sexpCtx);nprotect++;
	bcf_hdr_t*	hdr = (bcf_hdr_t*)R_ExternalPtrAddr(VECTOR_ELT(sexpCtx,0));
	bcf1_t* ctx = (bcf1_t*)R_ExternalPtrAddr(VECTOR_ELT(sexpCtx,1));
	const char* att = CHAR(asChar(sexpatt));
	ASSERT_NOT_NULL(ctx);
	ASSERT_NOT_NULL(hdr);
	ASSERT_NOT_NULL(att);
	FloatArrayPtr array = FloatArrayNew();


	bcf_unpack(ctx, BCF_UN_FMT);
	int tag_id = bcf_hdr_id2int(hdr, BCF_DT_ID, att);
	

	if(bcf_hdr_idinfo_exists(hdr,BCF_HL_FMT,tag_id)) {
		int ndst=0;
		float* dst=NULL;
		int ret=bcf_get_format_float(hdr,ctx,att,(void**)&dst,&ndst);
		if(ret>=0 && dst!=NULL ) {
			for (int j=0; j< ndst; j++) {
			    float val=dst[j];
			    if(val==bcf_float_vector_end ) break;
			    //IGNORE if(val==bcf_float_missing ) call->qsum[j] = 0;
			    FloatArrayPush(array,val);
			}
		}
		free(dst);
	}

	SEXP ext = PROTECT(allocVector(REALSXP,array->size));nprotect++;	
 	for(int i=0;i< array->size;i++) {
	       	REAL(ext)[i] = array->data[i];
		}

	FloatArrayFree(array);
	
	UNPROTECT(nprotect);
	return ext;
	}


