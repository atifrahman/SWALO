//
//  main.cpp
//  align
//
//  Created by Atif Rahman on 4/17/13.
//  Copyright (c) 2013 Atif Rahman. All rights reserved.
//
#include <iostream>
#include <string>
#include <limits>
#include <vector>
#include <algorithm>

using namespace std;
#include <time.h>
#include <pthread.h>
#include "cgal.h"
#include <emmintrin.h>

#include "swsse2.h"
#include "swstriped.h"

#define MAX_REC_LEN 1024

int noContigs=0;
long int noReads=0;
long int contigLength=0;

FILE *contigFile;
FILE *readFile;
FILE *summaryFile;
FILE *outFile;
char *contigFileName;
char *readFileName;

vector<char*> contigs;
vector<char*> contigNames;
vector<int> contigLengths;

long int totalReads;


long int totalReadCount;
int maxFragmentLength, maxReadLength;
int toAlign=100;

pthread_mutex_t readFile_mutex = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t outFile_mutex = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t totalReads_mutex = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t toAlign_mutex = PTHREAD_MUTEX_INITIALIZER;


//--------------------------------------------------


#define ALPHA_SIZE 6

const char AMINO_ACIDS[ALPHA_SIZE] = {
    'A', 'C', 'G', 'T', 'N', '_'
};

const int AMINO_ACID_VALUE[256] = {
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1,  0,  5,  1,  5,  -1,  -1,  2,  5, -1, -1,  5, -1, 5, 4, -1,
    -1, -1, 5, 5, 3, 5, 5, 5, 5, 5, -1, -1, -1, -1, -1, -1,
    -1,  0,  5,  1,  5,  -1,  -1,  2,  5,  -1, -1,  5, -1, 5, 4, -1,
    -1, -1, 5, 5, 3, 5, 5, 5, 5, 5, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1
};

typedef struct {
    __m128i        *pvbQueryProf;
    __m128i        *pvsQueryProf;
    __m128i        *pvH1;
    __m128i        *pvH2;
    __m128i        *pvE;
    unsigned char  *pData;
    unsigned short  bias;
} SwStripedData;


void printByte(__m128i v, int i)
{
	union u
    {
        __m128i m;
        char  c[16];
    } x;
    
    x.m = v;
	
	printf("%d",x.c[i]);
    
}

signed char matrix[36]={1,-1,-1,-1,-1,-1,
    -1,1,-1,-1,-1,-1,
    -1,-1,1,-1,-1,-1,
    -1,-1,-1,1,-1,-1,
    -1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1};

void *
swStripedInit(unsigned char   *querySeq2,
              int              queryLength,
              signed char     *matrix)
{
    
	char *querySeq=new char[queryLength+1];
    
	for(int i=0;i<queryLength;i++)
	{
		querySeq[i]=AMINO_ACID_VALUE[querySeq2[i]];
	}
	querySeq[queryLength]='\0';
    
    
    int i, j, k;
    
    int segSize;
    int nCount;
    
    int bias;
    
    int lenQryByte;
    int lenQryShort;
    
    int weight;
    
    short *ps;
    char *pc;
    
    signed char *matrixRow;
    
    size_t aligned;
    
    SwStripedData *pSwData;
    
    lenQryByte = (queryLength + 15) / 16;
    lenQryShort = (queryLength + 7) / 8;
    
    pSwData = (SwStripedData *) malloc (sizeof (SwStripedData));
    if (!pSwData) {
        fprintf (stderr, "Unable to allocate memory for SW data\n");
        exit (-1);
    }
    
    nCount = 64 +                             /* slack bytes */
    lenQryByte * ALPHA_SIZE +        /* query profile byte */
    lenQryShort * ALPHA_SIZE +       /* query profile short */
    (lenQryShort * 3);               /* vH1, vH2 and vE */
    
    pSwData->pData = (unsigned char *) calloc (nCount, sizeof (__m128i));
    if (!pSwData->pData) {
        fprintf (stderr, "Unable to allocate memory for SW data buffers\n");
        exit (-1);
    }
    
    /* since we might port this to another platform, lets align the data */
    /* to 16 byte boundries ourselves */
    aligned = ((size_t) pSwData->pData + 15) & ~(0x0f);
    
    pSwData->pvbQueryProf = (__m128i *) aligned;
    pSwData->pvsQueryProf = pSwData->pvbQueryProf + lenQryByte * ALPHA_SIZE;
    
    pSwData->pvH1 = pSwData->pvsQueryProf + lenQryShort * ALPHA_SIZE;
    pSwData->pvH2 = pSwData->pvH1 + lenQryShort;
    pSwData->pvE  = pSwData->pvH2 + lenQryShort;
    
    /* Use a scoring profile for the SSE2 implementation, but the layout
     * is a bit strange.  The scoring profile is parallel to the query, but is
     * accessed in a stripped pattern.  The query is divided into equal length
     * segments.  The number of segments is equal to the number of elements
     * processed in the SSE2 register.  For 8-bit calculations, the query will
     * be divided into 16 equal length parts.  If the query is not long enough
     * to fill the last segment, it will be filled with neutral weights.  The
     * first element in the SSE register will hold a value from the first segment,
     * the second element of the SSE register will hold a value from the
     * second segment and so on.  So if the query length is 288, then each
     * segment will have a length of 18.  So the first 16 bytes will  have
     * the following weights: Q1, Q19, Q37, ... Q271; the next 16 bytes will
     * have the following weights: Q2, Q20, Q38, ... Q272; and so on until
     * all parts of all segments have been written.  The last seqment will
     * have the following weights: Q18, Q36, Q54, ... Q288.  This will be
     * done for the entire alphabet.
     */
    
    /* Find the bias to use in the substitution matrix */
    bias = 127;
    for (i = 0; i < ALPHA_SIZE * ALPHA_SIZE; i++) {
        if (matrix[i] < bias) {
            bias = matrix[i];
        }
    }
    if (bias > 0) {
        bias = 0;
    }
    
    /* Fill in the byte query profile */
    pc = (char *) pSwData->pvbQueryProf;
    segSize = (queryLength + 15) / 16;
    nCount = segSize * 16;
    for (i = 0; i < ALPHA_SIZE; ++i) {
        matrixRow = matrix + i * ALPHA_SIZE;
        for (j = 0; j < segSize; ++j) {
            for (k = j; k < nCount; k += segSize) {
                if (k >= queryLength) {
                    weight = 0;
                } else {
                    weight = matrixRow[*(querySeq + k)];
                }
                *pc++ = (char) (weight - bias);
            }
        }
    }
	
    
    //	printByte(pSwData->pvbQueryProf[1],0);
    //	cout<<endl;
    
    /* Fill in the short query profile */
    ps = (short *) pSwData->pvsQueryProf;
    segSize = (queryLength + 7) / 8;
    nCount = segSize * 8;
    for (i = 0; i < ALPHA_SIZE; ++i) {
        matrixRow = matrix + i * ALPHA_SIZE;
        for (j = 0; j < segSize; ++j) {
            for (k = j; k < nCount; k += segSize) {
                if (k >= queryLength) {
                    weight = 0;
                } else {
                    weight = matrixRow[*(querySeq + k)];
                }
                *ps++ = (unsigned short) weight;
            }
        }
    }
    
    pSwData->bias = (unsigned short) -bias;
    
	delete []querySeq;
    
    return pSwData;
}

int
swStripedWord(unsigned char   *querySeq2,
              int              queryLength,
              unsigned char   *dbSeq2,
              int              dbLength,
              unsigned short   gapOpen,
              unsigned short   gapExtend,
              __m128i         *pvQueryProf,
              __m128i         *pvHLoad,
              __m128i         *pvHStore,
              __m128i         *pvE)
{
	char *querySeq=new char[queryLength+1];
    
	for(int i=0;i<queryLength;i++)
		querySeq[i]=AMINO_ACID_VALUE[querySeq2[i]];
    
	querySeq[queryLength]='\0';
    
	char *dbSeq=new char[dbLength+1];
    
	for(int i=0;i<dbLength;i++)
		dbSeq[i]=AMINO_ACID_VALUE[dbSeq2[i]];
    
	dbSeq[dbLength]='\0';
    
    
    int     i, j;
    int     score;
    
    int     cmp;
    int     iter = (queryLength + 7) / 8;
    
    __m128i *pv;
    
    __m128i vE=_mm_setzero_si128();
	__m128i vF=_mm_setzero_si128();
	__m128i vH=_mm_setzero_si128();
    
    __m128i vMaxScore=_mm_setzero_si128();
    __m128i vBias=_mm_setzero_si128();
    __m128i vGapOpen=_mm_setzero_si128();
    __m128i vGapExtend=_mm_setzero_si128();
    
    __m128i vTemp=_mm_setzero_si128();
    __m128i vZero=_mm_setzero_si128();
    
    __m128i vMin=_mm_setzero_si128();
    __m128i vMinimums=_mm_setzero_si128();
    
    __m128i *pvScore;
    
    
    
    /* remove unreferenced warning */
    querySeq;
    
    /* Load gap opening penalty to all elements of a constant */
    vGapOpen = _mm_insert_epi16 (vGapOpen, gapOpen, 0);
    vGapOpen = _mm_shufflelo_epi16 (vGapOpen, 0);
    vGapOpen = _mm_shuffle_epi32 (vGapOpen, 0);
    
    /* Load gap extension penalty to all elements of a constant */
    vGapExtend = _mm_insert_epi16 (vGapExtend, gapExtend, 0);
    vGapExtend = _mm_shufflelo_epi16 (vGapExtend, 0);
    vGapExtend = _mm_shuffle_epi32 (vGapExtend, 0);
    
    /*  load vMaxScore with the zeros.  since we are using signed */
    /*  math, we will bias the maxscore to -32768 so we have the */
    /*  full range of the short. */
    vMaxScore = _mm_cmpeq_epi16 (vMaxScore, vMaxScore);
    vMaxScore = _mm_slli_epi16 (vMaxScore, 15);
    
    vMinimums = _mm_shuffle_epi32 (vMaxScore, 0);
    
    vMin = _mm_shuffle_epi32 (vMaxScore, 0);
    vMin = _mm_srli_si128 (vMin, 14);
    
    /* Zero out the storage vector */
    for (i = 0; i < iter; i++)
    {
        _mm_store_si128 (pvE + i, vMaxScore);
        _mm_store_si128 (pvHStore + i, vMaxScore);
    }
    
    for (i = 0; i < dbLength; ++i)
    {
        /* fetch first data asap. */
        pvScore = pvQueryProf + dbSeq[i] * iter;
        
        /* zero out F. */
        vF = _mm_cmpeq_epi16 (vF, vF);
        vF = _mm_slli_epi16 (vF, 15);
        
        /* load the next h value */
        vH = _mm_load_si128 (pvHStore + iter - 1);
        vH = _mm_slli_si128 (vH, 2);
        vH = _mm_or_si128 (vH, vMin);
        
        pv = pvHLoad;
        pvHLoad = pvHStore;
        pvHStore = pv;
        
        for (j = 0; j < iter; j++)
        {
            /* load values of vF and vH from previous row (one unit up) */
            vE = _mm_load_si128 (pvE + j);
            
            /* add score to vH */
            vH = _mm_adds_epi16 (vH, *pvScore++);
            
            /* Update highest score encountered this far */
            vMaxScore = _mm_max_epi16 (vMaxScore, vH);
            
            /* get max from vH, vE and vF */
            vH = _mm_max_epi16 (vH, vE);
            vH = _mm_max_epi16 (vH, vF);
            
            /* save vH values */
            _mm_store_si128 (pvHStore + j, vH);
            
            /* update vE value */
            vH = _mm_subs_epi16 (vH, vGapOpen);
            vE = _mm_subs_epi16 (vE, vGapExtend);
            vE = _mm_max_epi16 (vE, vH);
            
            /* update vF value */
            vF = _mm_subs_epi16 (vF, vGapExtend);
            vF = _mm_max_epi16 (vF, vH);
            
            /* save vE values */
            _mm_store_si128 (pvE + j, vE);
            
            /* load the next h value */
            vH = _mm_load_si128 (pvHLoad + j);
        }
        
        /* reset pointers to the start of the saved data */
        j = 0;
        vH = _mm_load_si128 (pvHStore + j);
        
        /*  the computed vF value is for the given column.  since */
        /*  we are at the end, we need to shift the vF value over */
        /*  to the next column. */
        vF = _mm_slli_si128 (vF, 2);
        vF = _mm_or_si128 (vF, vMin);
        vTemp = _mm_subs_epi16 (vH, vGapOpen);
        vTemp = _mm_cmpgt_epi16 (vF, vTemp);
        cmp  = _mm_movemask_epi8 (vTemp);
        while (cmp != 0x0000)
        {
            vE = _mm_load_si128 (pvE + j);
            
            vH = _mm_max_epi16 (vH, vF);
            
            /* save vH values */
            _mm_store_si128 (pvHStore + j, vH);
            
            /*  update vE incase the new vH value would change it */
            vH = _mm_subs_epi16 (vH, vGapOpen);
            vE = _mm_max_epi16 (vE, vH);
            _mm_store_si128 (pvE + j, vE);
            
            /* update vF value */
            vF = _mm_subs_epi16 (vF, vGapExtend);
            
            j++;
            if (j >= iter)
            {
                j = 0;
                vF = _mm_slli_si128 (vF, 2);
                vF = _mm_or_si128 (vF, vMin);
            }
            
            vH = _mm_load_si128 (pvHStore + j);
            
            vTemp = _mm_subs_epi16 (vH, vGapOpen);
            vTemp = _mm_cmpgt_epi16 (vF, vTemp);
            cmp  = _mm_movemask_epi8 (vTemp);
        }
    }
    
    /* find largest score in the vMaxScore vector */
    vTemp = _mm_srli_si128 (vMaxScore, 8);
    vMaxScore = _mm_max_epi16 (vMaxScore, vTemp);
    vTemp = _mm_srli_si128 (vMaxScore, 4);
    vMaxScore = _mm_max_epi16 (vMaxScore, vTemp);
    vTemp = _mm_srli_si128 (vMaxScore, 2);
    vMaxScore = _mm_max_epi16 (vMaxScore, vTemp);
    
    /* store in temporary variable */
    score = (short) _mm_extract_epi16 (vMaxScore, 0);
    
    
	delete []querySeq;
	delete []dbSeq;
    
    /* return largest score */
    return score + SHORT_BIAS;
}

int
swStripedByte(unsigned char   *querySeq2,
              int              queryLength,
              unsigned char   *dbSeq2,
              int              dbLength,
              unsigned short   gapOpen,
              unsigned short   gapExtend,
              __m128i         *pvQueryProf,
              __m128i         *pvHLoad,
              __m128i         *pvHStore,
              __m128i         *pvE,
              unsigned short   bias)
{
	char *querySeq=new char[queryLength+1];
    
	for(int i=0;i<queryLength;i++)
		querySeq[i]=AMINO_ACID_VALUE[querySeq2[i]];
    
	querySeq[queryLength]='\0';
    
	char *dbSeq=new char[dbLength+1];
    
	for(int i=0;i<dbLength;i++)
		dbSeq[i]=AMINO_ACID_VALUE[dbSeq2[i]];
    
	dbSeq[dbLength]='\0';
    
    int     i, j;
    int     score;
    
    int     dup;
    int     cmp;
    int     iter = (queryLength + 15) / 16;
    
    
    __m128i *pv;
    
    __m128i vE=_mm_setzero_si128();
	__m128i vF=_mm_setzero_si128();
	__m128i vH=_mm_setzero_si128();
    
    __m128i vMaxScore=_mm_setzero_si128();
    __m128i vBias=_mm_setzero_si128();
    __m128i vGapOpen=_mm_setzero_si128();
    __m128i vGapExtend=_mm_setzero_si128();
    
    __m128i vTemp=_mm_setzero_si128();
    __m128i vZero=_mm_setzero_si128();
    
    __m128i *pvScore;
    
    
    /* remove unreferenced warning */
    querySeq;
    
    /* Load the bias to all elements of a constant */
    dup    = (bias << 8) | (bias & 0x00ff);
    vBias = _mm_insert_epi16 (vBias, dup, 0);
    vBias = _mm_shufflelo_epi16 (vBias, 0);
    vBias = _mm_shuffle_epi32 (vBias, 0);
    
    /* Load gap opening penalty to all elements of a constant */
    dup    = (gapOpen << 8) | (gapOpen & 0x00ff);
    vGapOpen = _mm_insert_epi16 (vGapOpen, dup, 0);
    vGapOpen = _mm_shufflelo_epi16 (vGapOpen, 0);
    vGapOpen = _mm_shuffle_epi32 (vGapOpen, 0);
    
    /* Load gap extension penalty to all elements of a constant */
    dup    = (gapExtend << 8) | (gapExtend & 0x00ff);
    vGapExtend = _mm_insert_epi16 (vGapExtend, dup, 0);
    vGapExtend = _mm_shufflelo_epi16 (vGapExtend, 0);
    vGapExtend = _mm_shuffle_epi32 (vGapExtend, 0);
    
    vMaxScore = _mm_xor_si128 (vMaxScore, vMaxScore);
    
    vZero = _mm_xor_si128 (vZero, vZero);
    
    /* Zero out the storage vector */
    for (i = 0; i < iter; i++)
    {
        _mm_store_si128 (pvE + i, vMaxScore);
        _mm_store_si128 (pvHStore + i, vMaxScore);
    }
    
    
    
    for (i = 0; i < dbLength; ++i)
    {
        /* fetch first data asap. */
        pvScore = pvQueryProf + dbSeq[i] * iter;
        
		
        
        /* zero out F. */
        vF = _mm_xor_si128 (vF, vF);
        
        /* load the next h value */
        vH = _mm_load_si128 (pvHStore + iter - 1);
        vH = _mm_slli_si128 (vH, 1);
        
        pv = pvHLoad;
        pvHLoad = pvHStore;
        pvHStore = pv;
        
        
        for (j = 0; j < iter; j++)
        {
            /* load values of vF and vH from previous row (one unit up) */
            vE = _mm_load_si128 (pvE + j);
            
            /* add score to vH */
            vH = _mm_adds_epu8 (vH, *(pvScore + j));
            vH = _mm_subs_epu8 (vH, vBias);
            
            
            /* Update highest score encountered this far */
            vMaxScore = _mm_max_epu8 (vMaxScore, vH);
            
            /* get max from vH, vE and vF */
            vH = _mm_max_epu8 (vH, vE);
            vH = _mm_max_epu8 (vH, vF);
            
            /* save vH values */
            _mm_store_si128 (pvHStore + j, vH);
            
            /* update vE value */
            vH = _mm_subs_epu8 (vH, vGapOpen);
            vE = _mm_subs_epu8 (vE, vGapExtend);
            vE = _mm_max_epu8 (vE, vH);
            
            /* update vF value */
            vF = _mm_subs_epu8 (vF, vGapExtend);
            vF = _mm_max_epu8 (vF, vH);
            
            /* save vE values */
            _mm_store_si128 (pvE + j, vE);
            
            /* load the next h value */
            vH = _mm_load_si128 (pvHLoad + j);
        }
        
        /* reset pointers to the start of the saved data */
        j = 0;
        vH = _mm_load_si128 (pvHStore + j);
        
        /*  the computed vF value is for the given column.  since */
        /*  we are at the end, we need to shift the vF value over */
        /*  to the next column. */
        vF = _mm_slli_si128 (vF, 1);
        vTemp = _mm_subs_epu8 (vH, vGapOpen);
        vTemp = _mm_subs_epu8 (vF, vTemp);
        vTemp = _mm_cmpeq_epi8 (vTemp, vZero);
        cmp  = _mm_movemask_epi8 (vTemp);
        
        while (cmp != 0xffff)
        {
            vE = _mm_load_si128 (pvE + j);
            
            vH = _mm_max_epu8 (vH, vF);
            
            /* save vH values */
            _mm_store_si128 (pvHStore + j, vH);
            
            /*  update vE incase the new vH value would change it */
            vH = _mm_subs_epu8 (vH, vGapOpen);
            vE = _mm_max_epu8 (vE, vH);
            _mm_store_si128 (pvE + j, vE);
            
            /* update vF value */
            vF = _mm_subs_epu8 (vF, vGapExtend);
            
            j++;
            if (j >= iter)
            {
                j = 0;
                vF = _mm_slli_si128 (vF, 1);
            }
            
            vH = _mm_load_si128 (pvHStore + j);
            
            vTemp = _mm_subs_epu8 (vH, vGapOpen);
            vTemp = _mm_subs_epu8 (vF, vTemp);
            vTemp = _mm_cmpeq_epi8 (vTemp, vZero);
            cmp  = _mm_movemask_epi8 (vTemp);
        }
		
    }
    
    
    
    /* find largest score in the vMaxScore vector */
    vTemp = _mm_srli_si128 (vMaxScore, 8);
    vMaxScore = _mm_max_epu8 (vMaxScore, vTemp);
    vTemp = _mm_srli_si128 (vMaxScore, 4);
    vMaxScore = _mm_max_epu8 (vMaxScore, vTemp);
    vTemp = _mm_srli_si128 (vMaxScore, 2);
    vMaxScore = _mm_max_epu8 (vMaxScore, vTemp);
    vTemp = _mm_srli_si128 (vMaxScore, 1);
    vMaxScore = _mm_max_epu8 (vMaxScore, vTemp);
    
    /* store in temporary variable */
    score = _mm_extract_epi16 (vMaxScore, 0);
    
	
	score = score & 0x00ff;
    
    /*  check if we might have overflowed */
    
	
	if (score + bias >= 255)
    {
        score = 255;
    }
    
    
	delete []querySeq;
	delete []dbSeq;
    
    /* return largest score */
    return score;
}

//----------------------------------------------------


void reverse(char *reverse, char *read)
{
    
	char ch='A';
	int readLength=strlen(read);
	for(int i=1; i<=readLength;i++)
	{
		ch=read[i-1];
		if(ch=='A')
			reverse[readLength-i]='T';
		else if(ch=='C')
			reverse[readLength-i]='G';
		else if(ch=='G')
			reverse[readLength-i]='C';
		else if(ch=='T')
			reverse[readLength-i]='A';
		else
			reverse[readLength-i]='N';
	}
	reverse[readLength]='\0';
    
}



void allocDPTable(int **dis, int m,int n)
{
	dis = new int *[m+1] ;
    
	for( int i = 0 ; i <= m ; i++ )
		dis[i] = new int[n+1];
	
}

void freeDPTable(int **dis, int m, int n)
{
	for( int i = 0 ; i <= m ; i++ )
        delete [] dis[i] ;
	delete [] dis ;
}


int getMatch(int **dis, char *s, char * t, int n)
{
	int m=strlen(s);
	
	int maxi=0;
	int maxj=0;
	int maxVal=0;
    
	int val1;
	int val2;
	int val3;
	
	int i,j;
    
	for(i=0;i<=m;i++)
        dis[i][0]=0;
    
    
	for(j=0;j<=n;j++)
        dis[0][j]=0;
    
	for(i=1;i<=m;i++)
	{
		for(j=1;j<=n;j++)
		{
			val1=dis[i-1][j-1];
			val2=dis[i-1][j]-1;
			val3=dis[i][j-1]-1;
			val2=val2>val3?val2:val3;
			val2=val2>0?val2:0;
			if(s[i-1]==t[j-1])
			{
				val1++;
                
				val3=val1>val2?val1:val2;
			}
			else
			{
				val1--;
                
				val3=val1>val2?val1:val2;
			}
			dis[i][j]=val3;
			if(val3>maxVal)
			{
                
				maxVal=val3;
			}
		}
		
	}
    
	return 	maxVal;
    
}

int align(int **dis, char *read, char *read2, char *qname1, char *qname2, char *quality1, char *quality2)
{
	int maxMatch=0;
	int tempMatch;
    
	int maxMatch1;
	int maxMatch2;
    
	int readLength1=strlen(read);
	int readLength2=strlen(read2);
    
	char *read2Reverse=new char[readLength2+1];
    
	reverse(read2Reverse,read2);
    
	int isReverse=0;
	int maxContig=0;
	int maxStart1=0;
	int maxStart2=0;
	int tempMaxStart1=0;
	int tempMaxStart2=0;
    
	void *swData=swStripedInit((unsigned char*)read,readLength1,matrix);
	SwStripedData *stripedData1 = (SwStripedData *) swData;
    
	swData=swStripedInit((unsigned char*)read2Reverse,readLength2,matrix);
	SwStripedData *stripedData2 = (SwStripedData *) swData;
    
	
	for(int i=0;i<contigs.size();i++)
	{
		int length=contigLengths[i];
		int start=0;
		
		maxMatch1=0;
        
		while(start<length-1)
		{
            
            //		tempMatch=getMatch(dis,read,contigs[i]+start,min(5000,length-start-1));
			tempMatch = swStripedByte ((unsigned char*)read,readLength1,
                                       (unsigned char*)(contigs[i]+start),min(5000,length-start-1),
                                       1, 1,
                                       stripedData1->pvbQueryProf,
                                       stripedData1->pvH1,
                                       stripedData1->pvH2,
                                       stripedData1->pvE,
                                       stripedData1->bias);
            
            
			if(tempMatch>maxMatch1)
			{
				maxMatch1=tempMatch;
				tempMaxStart1=start;
			}
			if(maxMatch1==readLength1)
				break;
			start+=4800;
		}
        
		
		start=0;
		maxMatch2=0;
		while(start<length-1)
		{
            
            //		tempMatch=getMatch(dis,read2Reverse,contigs[i]+start,min(5000,length-start-1));
			tempMatch = swStripedByte ((unsigned char*)read2Reverse,readLength2,
                                       (unsigned char*)(contigs[i]+start),min(5000,length-start-1),
                                       1, 1,
                                       stripedData2->pvbQueryProf,
                                       stripedData2->pvH1,
                                       stripedData2->pvH2,
                                       stripedData2->pvE,
                                       stripedData2->bias);
            
		    if(tempMatch>maxMatch2)
			{
				maxMatch2=tempMatch;
				tempMaxStart2=start;
			}
			if(maxMatch2==readLength2)
				break;
			start+=4800;
		}
        
		if(maxMatch1+maxMatch2>maxMatch)
		{
			maxMatch=maxMatch1+maxMatch2;
			maxContig=i;
			maxStart1=tempMaxStart1;
			maxStart2=tempMaxStart2;
		}
		
		if(maxMatch==readLength1+readLength2)
			break;
	}
	
    
	char *read1Reverse=new char[readLength1+1];
	reverse(read1Reverse,read);
	
	swData=swStripedInit((unsigned char*)read1Reverse,readLength1,matrix);
	stripedData1 = (SwStripedData *) swData;
    
	swData=swStripedInit((unsigned char*)read2,readLength2,matrix);
	stripedData2 = (SwStripedData *) swData;
    
    
	if (maxMatch<readLength1+readLength2)
	{
		
		int reverseMaxMatch=0;
        
		for(int i=0;i<contigs.size();i++)
		{
            
			int length=contigLengths[i];
			int start=0;
            
			maxMatch1=0;
            
			while(start<length-1)
			{
                
                //		tempMatch=getMatch(dis,read1Reverse,contigs[i]+start,min(5000,length-start-1));
				tempMatch = swStripedByte ((unsigned char*)read1Reverse,readLength1,
                                           (unsigned char*)contigs[i]+start,min(5000,length-start-1),
                                           1, 1,
                                           stripedData1->pvbQueryProf,
                                           stripedData1->pvH1,
                                           stripedData1->pvH2,
                                           stripedData1->pvE,
                                           stripedData1->bias);
                
                if(tempMatch>maxMatch1)
				{
					maxMatch1=tempMatch;
					tempMaxStart1=start;
				}
				if(maxMatch1==readLength1)
					break;
				start+=4800;
			}
            
            
			start=0;
			maxMatch2=0;
			while(start<length-1)
			{
                //		tempMatch=getMatch(dis,read2,contigs[i]+start,min(5000,length-start-1));
				tempMatch = swStripedByte ((unsigned char*)read2,readLength2,
                                           (unsigned char*)contigs[i]+start,min(5000,length-start-1),
                                           1, 1,
                                           stripedData2->pvbQueryProf,
                                           stripedData2->pvH1,
                                           stripedData2->pvH2,
                                           stripedData2->pvE,
                                           stripedData2->bias);
                
				if(tempMatch>maxMatch2)
				{
					maxMatch2=tempMatch;
					tempMaxStart2=start;
				}
				if(maxMatch2==readLength2)
					break;
				start+=4800;
			}
            
			if(maxMatch1+maxMatch2>maxMatch)
			{
				maxMatch=maxMatch1+maxMatch2;
				isReverse=1;
				maxContig=i;
				maxStart1=tempMaxStart1;
				maxStart2=tempMaxStart2;
			}
			if(maxMatch==readLength1+readLength2)
				break;
		}
        
	}
	if(isReverse==0)
	{
		strcpy(read2,read2Reverse);
	}
	else
	{
		strcpy(read,read1Reverse);
	}
    
	delete [] read2Reverse;
	delete [] read1Reverse;
    
    
	
	int length=contigLengths[maxContig];
    
	char *s=read;
	char t[10000];
    
    
	int n=min(5000,length-maxStart1-1);
    
	strncpy(t,contigs[maxContig]+maxStart1,n);
	int m=readLength1;
	
    
	char * contig=t;
    
	int maxi=0;
	int maxj=0;
	int maxVal=0;
    
	int val1;
	int val2;
	int val3;
    
    
	for(int i=0;i<=m;i++)
		dis[i][0]=0;
    
	for(int j=0;j<=n;j++)
		dis[0][j]=0;
    
	for(int i=1;i<=m;i++)
	{
		for(int j=1;j<=n;j++)
		{
			val1=dis[i-1][j-1];
			val2=dis[i-1][j]-1;
			val3=dis[i][j-1]-1;
			val2=val2>val3?val2:val3;
			val2=val2>0?val2:0;
			if(s[i-1]==t[j-1])
			{
                
				val1++;
				dis[i][j]=val1>val2?val1:val2;
			}
			else
			{
                
				val1--;
				dis[i][j]=val1>val2?val1:val2;
			}
			if(dis[i][j]>maxVal)
			{
				maxi=i;
				maxj=j;
				maxVal=dis[i][j];
			}
		}
	}
    
    
	
    
	int nm=(readLength1-maxMatch1)/2;
    
	char md[500];
	
	char cigar[500];
	char temp[500];
	char current[2];
	char dummy[2];
	int count=0;
	int matchCount=0;
	
	current[0]='M';
	current[1]='\0';
	dummy[0]='N';
	dummy[1]='\0';
	
	md[0]='\0';
	cigar[0]='\0';
	temp[0]='\0';
    
    
	int i=maxi;
	int j=maxj;
	int	value=dis[i][j];
	
	
	if(i<readLength1)
	{
		
        //	cigar=(readLength-i).toString+"I"
		itoa(readLength1-i,cigar,10);
		strcat(cigar,"I");
        
	}
    
	while(value>0)
	{
		if(read[i-1]==contig[j-1])
		{
			if(current[0]!='M')
			{
				if(current[0]=='D')
				{
                    //	md="^"+md
					strcpy(temp,md);
					strcpy(md,"^");
					strcat(md,temp);
				}
                //		cigar=count+current+cigar
                
				strcpy(temp,cigar);
				itoa(count,cigar,10);
				strcat(cigar,current);
				strcat(cigar,temp);
				
				current[0]='M';
				count=1;
			}
			else
			{
				count+=1;
			}
			matchCount+=1;
			i=i-1;
			j=j-1;
			value=dis[i][j];
		}
		else if(dis[i][j]==dis[i-1][j-1]-1)
		{
			if(current[0]!='M')
			{
                //	cigar=count+current+cigar
                
				strcpy(temp,cigar);
				itoa(count,cigar,10);
				strcat(cigar,current);
				strcat(cigar,temp);
				
				current[0]='M';
				count=1;
			}
			else
			{
				count+=1;
			}
			if(matchCount!=0)
			{
                //	md=matchCount+md
				strcpy(temp,md);
				itoa(matchCount,md,10);
				strcat(md,temp);
				
				matchCount=0;
			}
            //	md=contig.charAt(j-1)+md
			strcpy(temp,md);
			dummy[0]=contig[j-1];
			strcpy(md,dummy);
			strcat(md,temp);
			
			i=i-1;
			j=j-1;
			value=dis[i][j];
		}
		else if(dis[i][j]==dis[i][j-1]-1)
		{
			if(current[0]!='D')
			{
                //	cigar=count+current+cigar
                
				strcpy(temp,cigar);
				itoa(count,cigar,10);
				strcat(cigar,current);
				strcat(cigar,temp);
				
				current[0]='D';
				count=1;
			}
			else
			{
				count+=1;
			}
			if(matchCount!=0)
			{
                //	md=matchCount+md
				strcpy(temp,md);
				itoa(matchCount,md,10);
				strcat(md,temp);
				
				matchCount=0;
			}
            //	md=contig.charAt(j)+md
			strcpy(temp,md);
			dummy[0]=contig[j];
			strcpy(md,dummy);
			strcat(md,temp);
            
			i=i;
			j=j-1;
			value=dis[i][j];
		}
		else if(dis[i][j]==dis[i-1][j]-1)
		{
			if(current[0]!='I')
			{
                //	cigar=count+current+cigar
                
				strcpy(temp,cigar);
				itoa(count,cigar,10);
				strcat(cigar,current);
				strcat(cigar,temp);
				
				current[0]='I';
				count=1;
			}
			else
			{
				count+=1;
			}
			i=i-1;
			j=j;
			value=dis[i][j];
		}
		
	}
    
	
    //	cigar=count+current+cigar
	
	strcpy(temp,cigar);
	itoa(count,cigar,10);
	strcat(cigar,current);
	strcat(cigar,temp);
	
    
    //	md="MD:Z:"+matchCount+md
	strcpy(temp,md);
	itoa(matchCount,md,10);
	strcat(md,temp);
    
	int	pos=maxStart1+j+1;
    
	int	flag=isReverse<<4;
    
	if(i>0)
	{
        
		strcpy(temp,cigar);
		itoa(i,cigar,10);
		current[0]='I';
		strcat(cigar,current);
		strcat(cigar,temp);
        
        //	cigar=i+"I"+cigar
	}
    
	s=read2;
    
	n=min(5000,length-maxStart2-1);
	strncpy(t,contigs[maxContig]+maxStart2,n);
	m=readLength2;
    
	contig=t;
    
	maxi=0;
	maxj=0;
	maxVal=0;
    
	for(int i=0;i<= m;i++)
        dis[i][0]=0;
    
	for(int j=0;j<=n;j++)
        dis[0][j]=0;
    
	for(int i=1;i<=m;i++)
	{
		for(int j=1;j<=n;j++)
		{
			val1=dis[i-1][j-1];
			val2=dis[i-1][j]-1;
			val3=dis[i][j-1]-1;
			val2=val2>val3?val2:val3;
			val2=val2>0?val2:0;
			if(s[i-1]==t[j-1])
			{
                
				val1++;
				dis[i][j]=val1>val2?val1:val2;
			}
			else
			{
                
				val1--;
				dis[i][j]=val1>val2?val1:val2;
			}
			if(dis[i][j]>maxVal)
			{
				maxi=i;
				maxj=j;
				maxVal=dis[i][j];
			}
		}
	}
    
    
	
    
	int nm2=(readLength2-maxMatch2)/2;
    
	char md2[500];
	char cigar2[500];
	count=0;
	matchCount=0;
	
	current[0]='M';
	current[1]='\0';
	dummy[0]='N';
	dummy[1]='\0';
	
	md2[0]='\0';
	cigar2[0]='\0';
	temp[0]='\0';
    
    
	i=maxi;
	j=maxj;
	value=dis[i][j];
	
	
	if(i<readLength2)
	{
        //	cigar=(readLength-i).toString+"I"
		itoa(readLength2-i,cigar2,10);
		strcat(cigar2,"I");
	}
    
	while(value>0)
	{
        //		cout<<strlen(temp)<<endl;
		if(read2[i-1]==contig[j-1])
		{
			if(current[0]!='M')
			{
				if(current[0]=='D')
				{
                    //	md="^"+md
					strcpy(temp,md2);
					strcpy(md2,"^");
					strcat(md2,temp);
				}
                //		cigar=count+current+cigar
				strcpy(temp,cigar2);
				itoa(count,cigar2,10);
				strcat(cigar2,current);
				strcat(cigar2,temp);
				
				current[0]='M';
				count=1;
			}
			else
			{
				count+=1;
			}
			matchCount+=1;
			i=i-1;
			j=j-1;
			value=dis[i][j];
		}
		else if(dis[i][j]==dis[i-1][j-1]-1)
		{
			if(current[0]!='M')
			{
                //	cigar=count+current+cigar
				strcpy(temp,cigar2);
				itoa(count,cigar2,10);
				strcat(cigar2,current);
				strcat(cigar2,temp);
				
				current[0]='M';
				count=1;
			}
			else
			{
				count+=1;
			}
			if(matchCount!=0)
			{
                //	md=matchCount+md
				strcpy(temp,md2);
				itoa(matchCount,md2,10);
				strcat(md2,temp);
				
				matchCount=0;
			}
            //	md=contig.charAt(j-1)+md
			strcpy(temp,md2);
			dummy[0]=contig[j-1];
			strcpy(md2,dummy);
			strcat(md2,temp);
			
			i=i-1;
			j=j-1;
			value=dis[i][j];
		}
		else if(dis[i][j]==dis[i][j-1]-1)
		{
			if(current[0]!='D')
			{
                //	cigar=count+current+cigar
				strcpy(temp,cigar2);
				itoa(count,cigar2,10);
				strcat(cigar2,current);
				strcat(cigar2,temp);
				
				current[0]='D';
				count=1;
			}
			else
			{
				count+=1;
			}
			if(matchCount!=0)
			{
                //	md=matchCount+md
				strcpy(temp,md2);
				itoa(matchCount,md2,10);
				strcat(md2,temp);
				
				matchCount=0;
			}
            //	md=contig.charAt(j)+md
			strcpy(temp,md2);
			dummy[0]=contig[j];
			strcpy(md2,dummy);
			strcat(md2,temp);
            
			i=i;
			j=j-1;
			value=dis[i][j];
		}
		else if(dis[i][j]==dis[i-1][j]-1)
		{
			if(current[0]!='I')
			{
                //	cigar=count+current+cigar
				strcpy(temp,cigar2);
				itoa(count,cigar2,10);
				strcat(cigar2,current);
				strcat(cigar2,temp);
				
				current[0]='I';
				count=1;
			}
			else
			{
				count+=1;
			}
			i=i-1;
			j=j;
			value=dis[i][j];
		}
		
	}
    
	
    //	cigar=count+current+cigar
	strcpy(temp,cigar2);
	itoa(count,cigar2,10);
	strcat(cigar2,current);
	strcat(cigar2,temp);
    
    //	md="MD:Z:"+matchCount+md
	strcpy(temp,md2);
	itoa(matchCount,md2,10);
	strcat(md2,temp);
    
	int pos2=maxStart2+j+1;
    
	int flag2=(!isReverse)<<4;
    
	if(i>0)
	{
        //		cigar=i+"I"+cigar
		strcpy(temp,cigar2);
		itoa(i,cigar2,10);
		current[0]='I';
		strcat(cigar2,current);
		strcat(cigar2,temp);
        
	}
    
	int tlen=0;
	
	if(pos2>=pos)
		tlen=pos2-pos+readLength2;
	else
		tlen=pos2-pos-readLength1;
    
	pthread_mutex_lock(&outFile_mutex);
    
    //	fprintf(outFile,"%s\t%d\t%s\t%d\t%d\t%s\t%s\t%d\t%d\t%s\t%s\tNM:i:%d\tNH:i:1\tIH:i:1\tMD:Z:%s\n",qname1,flag,contigNames[maxContig],pos,0,cigar,contigNames[maxContig],pos2,pos2-pos,read,quality1,nm,md);
    //	fprintf(outFile,"%s\t%d\t%s\t%d\t%d\t%s\t%s\t%d\t%d\t%s\t%s\tNM:i:%d\tNH:i:1\tIH:i:1\tMD:Z:%s\n",qname2,flag2,contigNames[maxContig],pos2,0,cigar2,contigNames[maxContig],pos,pos-pos2,read2,quality2,nm2,md2);
    
	fprintf(outFile,"%s\t%d\t%s\t%d\t%s\t%d\t%s\tIH:i:1\tMD:Z:%s\n",qname1,flag,"dummy",pos,cigar,tlen,read,md);
	fprintf(outFile,"%s\t%d\t%s\t%d\t%s\t%d\t%s\tIH:i:1\tMD:Z:%s\n",qname2,flag2,"dummy",pos2,cigar2,-tlen,read2,md2);
	
	pthread_mutex_unlock(&outFile_mutex);
	
	fflush(outFile);
    
    
	return maxMatch;
}

void * map(void *threadid)
{
	char qname1[500],qname2[500],read1[500],read2[500],quality1[500],quality2[500];
	char *line= new char[MAX_REC_LEN];
	
	int MAX_FILE_READ=MAX_REC_LEN/sizeof(line[0]);
    
	
	int m=maxReadLength;
	int n=5000;
	
	
	int **dis;
	
	dis = new int *[m+1] ;
    
	for( int i = 0 ; i <= m ; i++ )
		dis[i] = new int[n+1];
    
	
	for(int i=0;i<=m;i++)
        dis[i][0]=0;
    
	for(int j=0;j<=n;j++)
        dis[0][j]=0;
	
    
    
    
	pthread_mutex_lock(&readFile_mutex);
	
	while(fgets(line, MAX_FILE_READ, readFile)!=NULL && toAlign>0)
	{
		strcpy(qname1,line+1);
		fgets(read1, MAX_FILE_READ, readFile);
		fgets(line, MAX_FILE_READ, readFile);
		fgets(quality1, MAX_FILE_READ, readFile);
        
		fgets(line, MAX_FILE_READ, readFile);
		strcpy(qname2,line+1);
		fgets(read2, MAX_FILE_READ, readFile);
		fgets(line, MAX_FILE_READ, readFile);
		fgets(quality2, MAX_FILE_READ, readFile);
        
		pthread_mutex_unlock(&readFile_mutex);
        
		
		read1[strlen(read1)-1]='\0';
		read2[strlen(read2)-1]='\0';
        
		qname1[strlen(qname1)-1]='\0';
		qname2[strlen(qname2)-1]='\0';
        
		quality1[strlen(quality1)-1]='\0';
		quality2[strlen(quality2)-1]='\0';
        
		double rnd=rand()/(double(RAND_MAX)+1);
        
		if(rnd<toAlign/(double)totalReads)
		{
			int readLength1=strlen(read1);
			int readLength2=strlen(read2);
			if(readLength1>maxReadLength || readLength2>maxReadLength)
			{
				maxReadLength=max(readLength1,readLength2);
				freeDPTable(dis,m,n);
				m=maxReadLength;
				allocDPTable(dis,m,n);
			}
			pthread_mutex_lock(&toAlign_mutex);
			toAlign--;
			pthread_mutex_unlock(&toAlign_mutex);
			align(dis,read1,read2,qname1,qname2,quality1,quality2);
			
			
            
		}
		pthread_mutex_lock(&totalReads_mutex);
		totalReads--;
		pthread_mutex_unlock(&totalReads_mutex);
        
		pthread_mutex_lock(&readFile_mutex);
        
	}
    
	pthread_mutex_unlock(&readFile_mutex);
	
	
	freeDPTable(dis,m,n);
    
	pthread_exit(NULL);
    
}

void printHelp()
{
    
    
	cout<<"swalo-"<<VERSION<<endl;
	cout<<"----------------"<<endl;
	cout<<endl;
	cout<<"align - maps random subset of reads not mapped by alignment tool using Smith-Waterman algorithm"<<endl;
	cout<<"Usage:"<<endl;
	cout<<"align  <contigFile> "<<endl;
	cout<<endl;
	cout<<"Required arguments:"<<endl;
	cout<<"<contigFile>\t Contig file in fasta format"<<endl;
	cout<<"Options:"<<endl;
	cout<<"-h [--help]\t\t Prints this message"<<endl;
	cout<<endl;
	exit(1);
    
}



int main(int argc, char *argv[])
{
	/*	input contig file name, read file name
     contig file - fasta format
     read file - fastq format
     */
    
    

	if(argc<2)
		printHelp();
    
	if(strcmp(argv[1],"--help")==0 || strcmp(argv[1],"-h")==0)
		printHelp();
    
    
	contigFileName=argv[1];
    
	toAlign=200;//atoi(argv[2]);
    
	srand(SEED);


//	if(argc==4)
//		NUM_THREADS=atoi(argv[3]);
//	else
//		NUM_THREADS=2;
    
    
//	readFileName="unmapped.txt";
    
	int unmapped=toAlign;
	
	contigFile=fopen(contigFileName, "r");
    
    
	if (contigFile == NULL)
	{
		printf("Can't open contig file\n");
		exit(1);
	}
    
	readFile=fopen("unmapped.txt", "r");
    
	if (readFile == NULL)
	{
		printf("Can't open read file\n");
		exit(1);
	}
    
	summaryFile=fopen("stat.txt", "r");
    
	if (summaryFile == NULL)
	{
		printf("Can't open stat file\n");
		exit(1);
	}
    
	outFile=fopen("unmappedOut.sam", "w");
    
	char *line= new char[MAX_REC_LEN];
    
	char *line1= new char[MAX_REC_LEN];
	char *line2= new char[MAX_REC_LEN];
	
	int read;
	int MAX_FILE_READ=MAX_REC_LEN/sizeof(line[0]);
    
	long int bufferLength=1024;
    
	char *contig=new char[bufferLength];
	contig[0]='\0';
	char *newcontig;
	char *contigName;
	contigLength=0;
    
    
	long int tempContigLength=0;
	
    
	while(fgets(line, MAX_FILE_READ, contigFile)!=NULL)
	{
		if(line[0]==';')
		{
			continue;
		}
		else if(line[0]=='>')
		{
			contigName=new char[strlen(line)];
			strcpy(contigName,line+1);
			contigName[strlen(contigName)-1]='\0';
			contigNames.push_back(strtok(contigName," \t\n"));
			if(contigLength>0)
			{
				noContigs++;
				contigs.push_back(toupper(contig));
				contigLengths.push_back(contigLength);
				contigLength=0;
				bufferLength=1024;
				contig=new char[bufferLength];
				contig[0]='\0';
			}
		}
		else
		{
			read=strlen(line);
			tempContigLength=contigLength;
			if(read<MAX_FILE_READ-1)
			{	
				contigLength+=(read-1);
			}
			else
			{
				contigLength+=MAX_FILE_READ-1;
				read++;
				
			}
			if(contigLength>bufferLength)
			{
				bufferLength=max(bufferLength*2,contigLength+1);
				newcontig=new char[bufferLength];
				strcpy(newcontig,contig);
				line[read-1]='\0';
				strcat(newcontig, line);
				delete []contig;
				contig=newcontig;
			}
			else
			{
				line[read-1]='\0';
				strcpy(contig+tempContigLength, line);	
			}
            
		}
		
	}
	noContigs++;
	contigs.push_back(toupper(contig));
	contigLengths.push_back(contigLength);
	
    
    
    
	fclose(contigFile);
    
    
    
    
    
	fscanf(summaryFile,"%ld %ld %d %d",&totalReadCount,&totalReads,&maxReadLength,&maxFragmentLength);
    
	long int totalUnmapped=totalReads;
    
	
	pthread_t threads[NUM_THREADS];
	int rc;
	long t;
	void *status;
	for(t=0; t<NUM_THREADS; t++)
	{
        printf("In main: creating thread %ld\n", t);
        rc = pthread_create(&threads[t], NULL, map, (void *)t);
        if (rc){
            printf("ERROR; return code from pthread_create() is %d\n", rc);
            exit(-1);
        }
	}
	
	for(t=0; t<NUM_THREADS; t++) 
	{
		rc = pthread_join(threads[t], &status);
		if (rc) 
	  	{
            printf("ERROR; return code from pthread_join() is %d\n", rc);
            exit(-1);
        }
        printf("Main: completed join with thread %ld having a status of %ld\n",t,(long)status);
    }
	
	fclose(readFile);
	fclose(outFile);
	fclose(summaryFile);
    
	summaryFile=fopen("stat.txt","w");
    
	fprintf(summaryFile,"%ld %ld %d %d %d",totalReadCount,totalUnmapped,unmapped, maxReadLength,maxFragmentLength);
    
	fclose(summaryFile);	
    
    
	pthread_exit(NULL);
    
	
    
	return 0;
}

