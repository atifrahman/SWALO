//
//  bowtieconvert.cpp
//  swalo
//

#include <iostream>
#include <string>
#include <limits>
#include <vector>
using namespace std;

#include "cgal.h"
#include <math.h>
#include <limits.h>


#define MAX_REC_LEN 4096
#define MAX_READLENGTH 1000
#define MAX_NAMELENGTH 1000

struct SAM
{
	char qname[MAX_NAMELENGTH];
	int flag;
	char rname[MAX_NAMELENGTH];
	int pos;
	int mapq;
	char cigar[MAX_READLENGTH];
	char rnext[MAX_NAMELENGTH];
	int pnext;
	int tlen;
	char seq[MAX_NAMELENGTH];
	char qual[MAX_NAMELENGTH];
	char md[MAX_NAMELENGTH];
	unsigned long ih;
	int nm;
	long int contigNo;
};


FILE *outFile;
FILE *linkingFile;
FILE *singletonFile;
FILE *mapFile1;
FILE *mapFile2;
FILE *contigFile;
FILE *unFile;
FILE *statFile;
FILE *scaffoldFile;



vector <SAM *> mixedReads1;
vector <SAM *> mixedReads2;



char line1[MAX_REC_LEN];
char line2[MAX_REC_LEN];

vector<char*> contigs;
vector<char*> contigNames;
vector<unsigned long> contigLengths;
long int noContigs;

double *contigReadCounts;



//int **readCounts;

long int unCount=0;
long int totalCount=0;
unsigned long maxReadLength=0;
double insertSizeMean=0;
double insertSizeVar=0;
double squaredError=0;

long int MAX_FRAGMENT_SIZE=5000;

struct Contig
{
	char contigName[1000];
	long int contigNo;
};


#define HASH_TABLE_SIZE 10001
vector<Contig *> hashTable[HASH_TABLE_SIZE];



int isJump=0;
long int X;

void reverse(char *reverse, char *read)
{
    
	char ch='A';
	unsigned long readLength=strlen(read);
	for(unsigned long i=1; i<=readLength;i++)
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

unsigned long getHash(char *str)
{
        unsigned long hash = 5381;
        int c;

        while (c = *str++)
            hash = ((hash << 5) + hash) + c; /* hash * 33 + c */

        return hash;
}

void initHashTable()
{

	unsigned long index;
	for(long int i=0;i<noContigs;i++)
	{
		index=getHash(contigNames[i]) % HASH_TABLE_SIZE;
		Contig *c=new Contig();
		strcpy(c->contigName,contigNames[i]);
		c->contigNo=i;
		hashTable[index].push_back(c);
	}
}

int getContigNo(char *contigName)
{
	unsigned long index=getHash(contigName) % HASH_TABLE_SIZE;
    
	for(int i=0;i<hashTable[index].size();i++)
	{
		if(strcmp(hashTable[index][i]->contigName,contigName)==0)
			return hashTable[index][i]->contigNo;
	}
	return -1;
    
}

/*
int getContigNo(char *contigName)
{
    
	for(int i=0;i<contigNames.size();i++)
	{
		if(strcmp(contigNames[i],contigName)==0)
			return i;
	}
	return -1;
    
}
*/

void writeSam(SAM* read, FILE *out)
{
    /*	fprintf(out,"%s\t%d\t%s\t%d\t%d\t%s\t%s\t%d\t%d\t%s\t%s\t\%s\tIH:i:%d\n",read->qname,read->flag,read->rname,read->pos,
     read->mapq,read->cigar,read->rnext,read->pnext,read->tlen,read->seq,read->qual,read->md,read->ih);
     */
    
	fprintf(out,"%s\t%d\t%ld\t%d\t%s\t%d\t%s\t\%s\tIH:i:%ld\n",read->qname,read->flag,read->contigNo,read->pos,
            read->cigar,read->tlen,read->seq,read->md,read->ih);
    
}

void printSam(SAM* read)
{
	printf("%s\t%d\t%s\t%d\t%d\t%s\t%s\t%d\t%d\t%s\t%s\t\%s\tIH:i:%ld\n",read->qname,read->flag,read->rname,read->pos,
           read->mapq,read->cigar,read->rnext,read->pnext,read->tlen,read->seq,read->qual,read->md,read->ih);
    
    
    /*	printf("%s\t%d\t%d\t%s\t%d\t%s\t\%s\tIH:i:%d\n",read->qname,read->flag,read->pos,
     read->cigar,read->tlen,read->seq,read->md,read->ih);
     */
}
/*
void printVectors(FILE * out, FILE * u)
{
    
	SAM *read1, *read2;
	unsigned long readLength1,readLength2;
	int pos1, pos2;
    //	int insertSize;
    
	char temp[MAX_READLENGTH];
    
	unsigned long ih=reads1.size();
    
	for(unsigned long i=0;i<ih;i++)
	{
		read1=reads1[i];
		read2=reads2[i];
        
		
		if(strcmp(read1->rname,"*")==0 || strcmp(read2->rname,"*")==0)
		{
            //		int i=0;
            //		int j=0;
			int ncount1=0;
			int ncount2=0;
            
			readLength1=strlen(read1->seq);
			if(readLength1>maxReadLength)
				maxReadLength=readLength1;
            
			readLength2=strlen(read2->seq);
			if(readLength2>maxReadLength)
				maxReadLength=readLength2;
            
            
			for(int i=0;i<readLength1;i++)
			{
				if(read1->seq[i]=='N'||read1->seq[i]=='n')
				{
					ncount1++;
				}
			}
            
			for(int i=0;i<readLength2;i++)
			{
                
				if(read2->seq[i]=='N'||read2->seq[i]=='n')
				{
					ncount2++;
				}
			}
            
			if(ncount1/(double)readLength1 < 0.8 && ncount2/(double)readLength2 < 0.8)
			{
				unCount++;
				totalCount++;
                
                
				fputs("@",u);
				fputs(read1->qname,u);
				fputs("\n",u);
                
                
				int strandNo=(read1->flag&16)>>4;
				if(strandNo==1)
				{
					reverse(temp,read1->seq);
					fputs(temp,u);
				}
				else
				{
					fputs(read1->seq,u);
				}
				fputs("\n",u);
                
				fputs("+",u);
				fputs(read1->qname,u);
				fputs("\n",u);
                
				fputs(read1->qual,u);
				fputs("\n",u);
                
				fputs("@",u);
				fputs(read2->qname,u);
				fputs("\n",u);
                
				strandNo=(read2->flag&16)>>4;
				if(strandNo==1)
				{
					reverse(temp,read2->seq);
					fputs(temp,u);
				}
				else
				{
					fputs(read2->seq,u);
                    
				}
				fputs("\n",u);
                
				fputs("+",u);
				fputs(read2->qname,u);
				fputs("\n",u);
                
				fputs(read2->qual,u);
				fputs("\n",u);
				
				for(int i=0;i<reads1.size();i++)
					delete reads1[i];
                
				for(int i=0;i<reads2.size();i++)
					delete reads2[i];
                
                
				reads1.clear();
				reads2.clear();
				
				return;
			}
            
		}
		else if(strcmp(read1->rname,read2->rname)!=0)
		{
            
//			readCounts[read1->contigNo][read2->contigNo]++;
			printSam(read1);
			printSam(read2);
//			getchar();
            
		}
		else
		{
//			readCounts[read1->contigNo][read2->contigNo]++;
            
			pos1=read1->pos;
			readLength1=strlen(read1->seq);
			if(readLength1>maxReadLength)
				maxReadLength=readLength1;
			read1->ih=ih;
            
            
			pos2=read2->pos;
			readLength2=strlen(read2->seq);
			if(readLength2>maxReadLength)
				maxReadLength=readLength2;
			read2->ih=ih;
            
			writeSam(read1,out);
			writeSam(read2,out);
            
		}
	}
	totalCount++;
    
	
    
	for(int i=0;i<reads1.size();i++)
		delete reads1[i];
    
	for(int i=0;i<reads2.size();i++)
		delete reads2[i];
	
    
	reads1.clear();
	reads2.clear();
    
}
*/

void printMixedVectors(FILE * linking, FILE * single, FILE * u)
{
    
    
	SAM *read1, *read2;
	unsigned long readLength1,readLength2;
    //	int pos1, pos2;
    //	int insertSize;
    
    //	char temp[MAX_READLENGTH];
    
    
    char temp[MAX_READLENGTH];
    
	int ih1=mixedReads1.size();
    int ih2=mixedReads2.size();
    
    if(ih1>1 || ih2>1)
    {
        for(int i=0;i<mixedReads1.size();i++)
        {
            read1=mixedReads1[i];
//            cout<<read1->qname<<" ";
        }
        for(int j=0;j<mixedReads2.size();j++)
        {
            read1=mixedReads2[j];
//            cout<<read1->qname<<" ";
        }
//        cout<<endl;
    }
    
    
	for(int i=0;i<mixedReads1.size();i++)
	{
		read1=mixedReads1[i];
        
		for(int j=0;j<mixedReads2.size();j++)
		{
			read2=mixedReads2[j];
            
            
            if(((read1->flag) & 4) != 0 || ((read2->flag)& 4) != 0)
			{
                for(int i=0;i<mixedReads1.size();i++)
                    delete mixedReads1[i];
                for(int i=0;i<mixedReads2.size();i++)
                    delete mixedReads2[i];
                
                
				mixedReads1.clear();
				mixedReads2.clear();
				
				return;
			}
            
            if(i==0 && j==0)
            {
                readLength1=strlen(read1->seq);
                if(readLength1>maxReadLength)
                    maxReadLength=readLength1;
                
                readLength2=strlen(read2->seq);
                if(readLength2>maxReadLength)
                    maxReadLength=readLength2;
                
                
                int ncount1=0;
                int ncount2=0;
                
                
                for(int i=0;i<readLength1;i++)
                {
                    if(read1->seq[i]=='N'||read1->seq[i]=='n')
                    {
                        ncount1++;
                    }
                }
                
                for(int i=0;i<readLength2;i++)
                {
                    
                    if(read2->seq[i]=='N'||read2->seq[i]=='n')
                    {
                        ncount2++;
                    }
                }
                
                if(ncount1/(double)readLength1 < 0.8 && ncount2/(double)readLength2 < 0.8)
                {
                    unCount++;
                    totalCount++;
                    
                    
                    fputs("@",u);
                    fputs(read1->qname,u);
                    fputs("\n",u);
                    
                    
                    int strandNo=(read1->flag&16)>>4;
                    if(strandNo==1)
                    {
                        reverse(temp,read1->seq);
                        fputs(temp,u);
                    }
                    else
                    {
                        fputs(read1->seq,u);
                    }
                    fputs("\n",u);
                    
                    fputs("+",u);
                    fputs(read1->qname,u);
                    fputs("\n",u);
                    
                    fputs(read1->qual,u);
                    fputs("\n",u);
                    
                    fputs("@",u);
                    fputs(read2->qname,u);
                    fputs("\n",u);
                    
                    strandNo=(read2->flag&16)>>4;
                    if(strandNo==1)
                    {
                        reverse(temp,read2->seq);
                        fputs(temp,u);
                    }
                    else
                    {
                        fputs(read2->seq,u);
                        
                    }
                    fputs("\n",u);
                    
                    fputs("+",u);
                    fputs(read2->qname,u);
                    fputs("\n",u);
                    
                    fputs(read2->qual,u);
                    fputs("\n",u);
                }
                else
                {
                    for(int i=0;i<mixedReads1.size();i++)
                        delete mixedReads1[i];
                    
                    for(int i=0;i<mixedReads2.size();i++)
                        delete mixedReads2[i];
                    
                    
                    mixedReads1.clear();
                    mixedReads2.clear();
                    
                    return;
                    
                }
            }
            /*            if(mixedReads1.size()>1 || mixedReads2.size()>1)
             {
             for(int i=0;i<mixedReads1.size();i++)
             delete mixedReads1[i];
             
             for(int i=0;i<mixedReads2.size();i++)
             delete mixedReads2[i];
             
             
             mixedReads1.clear();
             mixedReads2.clear();
             
             return;
             
             }
             */
            
            if(((read1->flag) & 4) != 0 && ((read2->flag)& 4) != 0)
			{
                for(int i=0;i<mixedReads1.size();i++)
                    delete mixedReads1[i];
                for(int i=0;i<mixedReads2.size();i++)
                    delete mixedReads2[i];
                
                
				mixedReads1.clear();
				mixedReads2.clear();
				
				return;
			}
            else if(((read1->flag) & 4) == 0 && ((read2->flag)& 4) != 0)
            {
                for(int i=0;i<mixedReads1.size();i++)
                {
                    writeSam(mixedReads1[i], single);
                }
                for(int i=0;i<mixedReads2.size();i++)
                {
                    writeSam(mixedReads2[i], single);
                }
                
                for(int i=0;i<mixedReads1.size();i++)
                    delete mixedReads1[i];
                for(int i=0;i<mixedReads2.size();i++)
                    delete mixedReads2[i];
                
                
                mixedReads1.clear();
                mixedReads2.clear();
                
                return;
                
            }
            else if(((read1->flag) & 4) != 0 && ((read2->flag)& 4) == 0)
            {
                for(int i=0;i<mixedReads1.size();i++)
                {
                    writeSam(mixedReads1[i], single);
                }
                for(int i=0;i<mixedReads2.size();i++)
                {
                    writeSam(mixedReads2[i], single);
                }
                for(int i=0;i<mixedReads1.size();i++)
                    delete mixedReads1[i];
                for(int i=0;i<mixedReads2.size();i++)
                    delete mixedReads2[i];
                
                mixedReads1.clear();
                mixedReads2.clear();
                return;
            }
            else if(strcmp(read1->rname,read2->rname)!=0)
			{
                read1->ih=ih1;
                read2->ih=ih2;
//                readCounts[read1->contigNo][read2->contigNo]++;
                writeSam(read1,linking);
                writeSam(read2, linking);
                //	printSam(read1);
                //	printSam(read2);
                //	getchar();
                
			}
		}
	}
    
	for(int i=0;i<mixedReads1.size();i++)
		delete mixedReads1[i];
	for(int i=0;i<mixedReads2.size();i++)
		delete mixedReads2[i];
	
    
	mixedReads1.clear();
	mixedReads2.clear();
    
}

void testMixedVectors_jump(FILE* out,FILE * linking, FILE * single, FILE * u)
{
    
    SAM *read1, *read2;

    unsigned long readLength1,readLength2;
	int pos1, pos2;
    //	int insertSize;
    
	char temp[MAX_READLENGTH];
    

	int ih1=mixedReads1.size();
	int ih2=mixedReads2.size();

	int ih1_0=0;
	int ih2_0=0;  


	for(int i=0;i<ih1;i++)
	{
		if(strcmp(mixedReads1[i]->rname,"*")!=0 && mixedReads1[i]->nm==0)
		{
			ih1_0++;
		}
	}

	for(int i=0;i<ih2;i++)
	{
		if(strcmp(mixedReads2[i]->rname,"*")!=0  && mixedReads2[i]->nm==0)
		{
			ih2_0++;
		}
	}
	
	if(ih1_0>0)
	{
	for(int i=0;i<ih1;i++)
	{
		if(strcmp(mixedReads1[i]->rname,"*")!=0  && mixedReads1[i]->nm==0)
		{
			contigReadCounts[mixedReads1[i]->contigNo]+=1/(double)ih1_0;
		}
	}
	}
	if(ih2_0>0)
	{
	for(int i=0;i<ih2;i++)
	{
		if(strcmp(mixedReads2[i]->rname,"*")!=0   && mixedReads2[i]->nm==0)
		{
			contigReadCounts[mixedReads2[i]->contigNo]+=1/(double)ih2_0;
		}
	}
	}
     
    
    int pairedCount=0;

    for(int i=0;i<mixedReads1.size();i++)
	{
		read1=mixedReads1[i];
        readLength1=strlen(read1->seq);
        if(readLength1>maxReadLength)
            maxReadLength=readLength1;
        

        
		for(int j=0;j<mixedReads2.size();j++)
		{
			read2=mixedReads2[j];
            readLength2=strlen(read2->seq);

            if(readLength2>maxReadLength)
                maxReadLength=readLength2;
            
            
            if(strcmp(read1->rname,"*")==0 || strcmp(read2->rname,"*")==0)
            {
                //		int i=0;
                //		int j=0;
                int ncount1=0;
                int ncount2=0;
                
                
                
                
                for(int i=0;i<readLength1;i++)
                {
                    if(read1->seq[i]=='N'||read1->seq[i]=='n')
                    {
                        ncount1++;
                    }
                }
                
                for(int i=0;i<readLength2;i++)
                {
                    
                    if(read2->seq[i]=='N'||read2->seq[i]=='n')
                    {
                        ncount2++;
                    }
                }
                
                if(ncount1/(double)readLength1 < 0.8 && ncount2/(double)readLength2 < 0.8)
                {
                    unCount++;
                    totalCount++;
                    
                    
                    fputs("@",u);
                    fputs(read1->qname,u);
                    fputs("\n",u);
                    
                    
                    int strandNo=(read1->flag&16)>>4;
                    if(strandNo==1)
                    {
                        reverse(temp,read1->seq);
                        fputs(temp,u);
                    }
                    else
                    {
                        fputs(read1->seq,u);
                    }
                    fputs("\n",u);
                    
                    fputs("+",u);
                    fputs(read1->qname,u);
                    fputs("\n",u);
                    
                    fputs(read1->qual,u);
                    fputs("\n",u);
                    
                    fputs("@",u);
                    fputs(read2->qname,u);
                    fputs("\n",u);
                    
                    strandNo=(read2->flag&16)>>4;
                    if(strandNo==1)
                    {
                        reverse(temp,read2->seq);
                        fputs(temp,u);
                    }
                    else
                    {
                        fputs(read2->seq,u);
                        
                    }
                    fputs("\n",u);
                    
                    fputs("+",u);
                    fputs(read2->qname,u);
                    fputs("\n",u);
                    
                    fputs(read2->qual,u);
                    fputs("\n",u);
                    
                    for(int i=0;i<mixedReads1.size();i++)
                        delete mixedReads1[i];
                    
                    for(int i=0;i<mixedReads2.size();i++)
                        delete mixedReads2[i];
                    
                    
                    mixedReads1.clear();
                    mixedReads2.clear();
                    
                    return;
                }
            }
            else if(strcmp(read1->rname,read2->rname)==0)
            {
                if(read1->pos-read2->pos+readLength1<X && read1->pos>read2->pos && (read1->flag&16)>>4 == 0 && (read2->flag&16)>>4 == 1)
                {
                    read1->tlen=read1->pos-read2->pos+readLength1;
                    read2->tlen=-read1->tlen;
                    pairedCount++;
                    
                }
                else if(read2->pos-read1->pos+readLength2<X && read2->pos>read1->pos && (read2->flag&16)>>4 == 0 && (read1->flag&16)>>4 == 1)
                {
                    read2->tlen=read2->pos-read1->pos+readLength2;
                    read1->tlen=-read2->tlen;

                    pairedCount++;
                }
                
            }
            
        }
    }
    if(pairedCount==0)
    {
        printMixedVectors(linking, single, u);
        return;
    }
    int ih=pairedCount;
    for(int i=0;i<mixedReads1.size();i++)
	{
		read1=mixedReads1[i];
        readLength1=strlen(read1->seq);
        
        
		for(int j=0;j<mixedReads2.size();j++)
		{
			read2=mixedReads2[j];
            readLength2=strlen(read2->seq);
            
            if(strcmp(read1->rname,read2->rname)==0)
            {
                if(read1->pos-read2->pos+readLength1<X && read1->pos>read2->pos && (read1->flag&16)>>4 == 0 && (read2->flag&16)>>4 == 1)
                {
                    read1->ih=ih;
                    read2->ih=ih;
                    writeSam(read1,out);
                    writeSam(read2,out);
                    
                                       
                }
                else if(read2->pos-read1->pos+readLength2<X && read2->pos>read1->pos && (read2->flag&16)>>4 == 0 && (read1->flag&16)>>4 == 1)
                {
                    read1->ih=ih;
                    read2->ih=ih;
                    writeSam(read1,out);
                    writeSam(read2,out);
                    
                    
                   
                }
                
            }
            
        }
    }
    totalCount++;
    
	for(int i=0;i<mixedReads1.size();i++)
		delete mixedReads1[i];
	for(int i=0;i<mixedReads2.size();i++)
		delete mixedReads2[i];
	
    
	mixedReads1.clear();
	mixedReads2.clear();
    
}


void testMixedVectors(FILE* out,FILE * linking, FILE * single, FILE * u)
{
    
    SAM *read1, *read2;

    unsigned long readLength1,readLength2;
	int pos1, pos2;
    //	int insertSize;
    
	char temp[MAX_READLENGTH];
    
     
    
    int pairedCount=0;

	

	int ih1=mixedReads1.size();
	int ih2=mixedReads2.size();


	int ih1_0=0;
	int ih2_0=0;  


	for(int i=0;i<ih1;i++)
	{
		if(strcmp(mixedReads1[i]->rname,"*")!=0 && mixedReads1[i]->nm==0)
		{
			ih1_0++;
		}
	}

	for(int i=0;i<ih2;i++)
	{
		if(strcmp(mixedReads2[i]->rname,"*")!=0  && mixedReads2[i]->nm==0)
		{
			ih2_0++;
		}
	}
	
	if(ih1_0>0)
	{
	for(int i=0;i<ih1;i++)
	{
		if(strcmp(mixedReads1[i]->rname,"*")!=0  && mixedReads1[i]->nm==0)
		{
			contigReadCounts[mixedReads1[i]->contigNo]+=1/(double)ih1_0;
		}
	}
	}
	if(ih2_0>0)
	{
	for(int i=0;i<ih2;i++)
	{
		if(strcmp(mixedReads2[i]->rname,"*")!=0   && mixedReads2[i]->nm==0)
		{
			contigReadCounts[mixedReads2[i]->contigNo]+=1/(double)ih2_0;
		}
	}
	}


    for(int i=0;i<mixedReads1.size();i++)
	{
		read1=mixedReads1[i];
        readLength1=strlen(read1->seq);
        if(readLength1>maxReadLength)
            maxReadLength=readLength1;
        

        
		for(int j=0;j<mixedReads2.size();j++)
		{
			read2=mixedReads2[j];
            readLength2=strlen(read2->seq);

            if(readLength2>maxReadLength)
                maxReadLength=readLength2;
            
            
            if(strcmp(read1->rname,"*")==0 || strcmp(read2->rname,"*")==0)
            {
                //		int i=0;
                //		int j=0;
                int ncount1=0;
                int ncount2=0;
                
                
                
                
                for(int i=0;i<readLength1;i++)
                {
                    if(read1->seq[i]=='N'||read1->seq[i]=='n')
                    {
                        ncount1++;
                    }
                }
                
                for(int i=0;i<readLength2;i++)
                {
                    
                    if(read2->seq[i]=='N'||read2->seq[i]=='n')
                    {
                        ncount2++;
                    }
                }
                
                if(ncount1/(double)readLength1 < 0.8 && ncount2/(double)readLength2 < 0.8)
                {
                    unCount++;
                    totalCount++;
                    
                    
                    fputs("@",u);
                    fputs(read1->qname,u);
                    fputs("\n",u);
                    
                    
                    int strandNo=(read1->flag&16)>>4;
                    if(strandNo==1)
                    {
                        reverse(temp,read1->seq);
                        fputs(temp,u);
                    }
                    else
                    {
                        fputs(read1->seq,u);
                    }
                    fputs("\n",u);
                    
                    fputs("+",u);
                    fputs(read1->qname,u);
                    fputs("\n",u);
                    
                    fputs(read1->qual,u);
                    fputs("\n",u);
                    
                    fputs("@",u);
                    fputs(read2->qname,u);
                    fputs("\n",u);
                    
                    strandNo=(read2->flag&16)>>4;
                    if(strandNo==1)
                    {
                        reverse(temp,read2->seq);
                        fputs(temp,u);
                    }
                    else
                    {
                        fputs(read2->seq,u);
                        
                    }
                    fputs("\n",u);
                    
                    fputs("+",u);
                    fputs(read2->qname,u);
                    fputs("\n",u);
                    
                    fputs(read2->qual,u);
                    fputs("\n",u);
                    
                    for(int i=0;i<mixedReads1.size();i++)
                        delete mixedReads1[i];
                    
                    for(int i=0;i<mixedReads2.size();i++)
                        delete mixedReads2[i];
                    
                    
                    mixedReads1.clear();
                    mixedReads2.clear();
                    
                    return;
                }
            }
            else if(strcmp(read1->rname,read2->rname)==0)
            {
                if(read1->pos-read2->pos+readLength1<X && read1->pos>read2->pos && (read1->flag&16)>>4 == 1 && (read2->flag&16)>>4 == 0)
                {
                    read1->tlen=read1->pos-read2->pos+readLength1;
                    read2->tlen=-read1->tlen;
                    pairedCount++;
                    
                }
                else if(read2->pos-read1->pos+readLength2<X && read2->pos>read1->pos && (read2->flag&16)>>4 == 1 && (read1->flag&16)>>4 == 0)
                {
                    read2->tlen=read2->pos-read1->pos+readLength2;
                    read1->tlen=-read2->tlen;

                    pairedCount++;
                }
                
            }
            
        }
    }
    if(pairedCount==0)
    {
        printMixedVectors(linking, single, u);
        return;
    }
    int ih=pairedCount;
    for(int i=0;i<mixedReads1.size();i++)
	{
		read1=mixedReads1[i];
        readLength1=strlen(read1->seq);
        
        
		for(int j=0;j<mixedReads2.size();j++)
		{
			read2=mixedReads2[j];
            readLength2=strlen(read2->seq);
            
            if(strcmp(read1->rname,read2->rname)==0)
            {
                if(read1->pos-read2->pos+readLength1<X && read1->pos>read2->pos && (read1->flag&16)>>4 == 1 && (read2->flag&16)>>4 == 0)
                {
                    read1->ih=ih;
                    read2->ih=ih;
                    writeSam(read1,out);
                    writeSam(read2,out);
                    
                                       
                }
                else if(read2->pos-read1->pos+readLength2<X && read2->pos>read1->pos && (read2->flag&16)>>4 == 1 && (read1->flag&16)>>4 == 0)
                {
                    read1->ih=ih;
                    read2->ih=ih;
                    writeSam(read1,out);
                    writeSam(read2,out);
                    
                    
                   
                }
                
            }
            
        }
    }
    totalCount++;
    
	for(int i=0;i<mixedReads1.size();i++)
		delete mixedReads1[i];
	for(int i=0;i<mixedReads2.size();i++)
		delete mixedReads2[i];
	
    
	mixedReads1.clear();
	mixedReads2.clear();
    
}


SAM *getSAM(char *line)
{
  	SAM *sam=new SAM;
	char *temp;
    
	temp=strtok(line,"\t\n ");
	strcpy(sam->qname,temp);
    
	temp=strtok(NULL,"\t\n ");
	sam->flag=atoi(temp);
    
	temp=strtok(NULL,"\t\n ");
	strcpy(sam->rname,temp);
    
	temp=strtok(NULL,"\t\n ");
	sam->pos=atoi(temp);
    
	temp=strtok(NULL,"\t\n ");
	sam->mapq=atoi(temp);
    
	temp=strtok(NULL,"\t\n ");
	strcpy(sam->cigar,temp);
	
	temp=strtok(NULL,"\t\n ");
	strcpy(sam->rnext,temp);
    
	temp=strtok(NULL,"\t\n ");
	sam->pnext=atoi(temp);
    
	temp=strtok(NULL,"\t\n ");
	sam->tlen=atoi(temp);
    
	temp=strtok(NULL,"\t\n ");
	strcpy(sam->seq,temp);
	
	temp=strtok(NULL,"\t\n ");
	strcpy(sam->qual,temp);
    
	sam->nm=-1;

	while((temp=strtok(NULL,"\t\n "))!=NULL)
	{
		if(temp[0]=='M' && temp[1]=='D')
		{
			strcpy(sam->md,temp);
		}
		if(temp[0]=='N' && temp[1]=='M')
		{
			sam->nm=atoi(&temp[5]);
		}

	}
    	
	sam->contigNo=getContigNo(sam->rname);
    
	return sam;
}

void printHelp()
{
    
	cout<<"swalo-"<<VERSION<<endl;
	cout<<"----------------"<<endl;
	cout<<endl;
	cout<<"bowtieconvert - to preprocess alignments generated with Bowtie"<<endl;
	cout<<"Note: Reads need to be mapped with Bowtie first. Map each end separately. For non repeat rich genomes -a can be used otherwise use -k 5"<<endl; 
	cout<<"Usage:"<<endl;
	cout<<"bowtieconvert <mapFile1> <mapFile2> <contigFile> <maxInsertSize> <insertSizeLimit(>=maxInsertSize)>"<<endl;
	cout<<endl;
	cout<<"Required arguments:"<<endl;
	cout<<"<mapFile1>\t Map file 1 in SAM format"<<endl;
	cout<<"<mapFile2>\t Map file 2 in SAM format"<<endl;
	cout<<"<contigFile>\t Contig file in FASTA format"<<endl;
	cout<<"<maxInsertSize>\t An estimate of maximum insert size. Can be set something large. Does not need to be exact."<<endl;
	cout<<"<insertSizeLimit(>=maxInsertSize)>\t Set equal to <maxInsertSize>. Only used to allow consistency with an evaluation paper."<<endl;
	cout<<endl;
	cout<<"Optional arguments:"<<endl;
	cout<<"--jump\t If the reads are mate-pair reads"<<endl;
	cout<<endl;
	exit(1);
    
}

int main(int argc, char *argv[])
{
	/*	input contig file name, read file name
     contig file - fasta format
     read file - fastq format
     */
    

     if(argc<6)
     printHelp();
     
     X=atoi(argv[4]);

	MAX_FRAGMENT_SIZE=atoi(argv[5]);
	
	for(int i=6;i<argc;i++)
	{
		if(strcmp(argv[i],"-j")==0 || strcmp(argv[i],"--jump")==0)
		{
			isJump=1;
		}
	}     
	

    
	char *line= new char[MAX_REC_LEN];
    char *line1= new char[MAX_REC_LEN];
	char *line2= new char[MAX_REC_LEN];

    //	char *templine= new char[MAX_REC_LEN];
	
	int MAX_FILE_READ=MAX_REC_LEN/sizeof(line[0]);
	
	char * mapFileName1=argv[1];	//"/Users/atif/Studies/mlga/BowtieConvert/BowtieConvert/map_bowtie_v3_1.sam";
    char * mapFileName2=argv[2];	//"/Users/atif/Studies/mlga/BowtieConvert/BowtieConvert/map_bowtie_v3_2.sam";
	char *	contigFileName=argv[3];	//"/Users/atif/Studies/mlga/swalo/swalo/contigs.fa";
    
    
	contigFile=fopen(contigFileName, "r");
    
    
	if (contigFile == NULL)
	{
		printf("Can't open contig file\n");
		exit(1);
	}
    
	noContigs=0;
	
	long int contigLength=0;
	
	unsigned long read;
	char *contig;
	char *newcontig;
	char *contigName;
	long int bufferLength=1024;
    
	contigLength=0;
    
	long int tempContigLength=0;
    
	contig=new char[bufferLength];
	contig[0]='\0';
	
    // read contig file
    
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
				contig=(char *)realloc(contig, (contigLength+1)*sizeof(char));
				contigs.push_back(contig);
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
				strcpy(newcontig+tempContigLength, line);
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
	contig=(char *)realloc(contig, (contigLength+1)*sizeof(char));
	contigs.push_back(contig);
	contigLengths.push_back(contigLength);
	
    
	cout<<noContigs<<endl;
    //	cout<<contig<<endl;
    //	cout<<contigLength<<endl;
    

    
	fclose(contigFile);

	initHashTable();

	contigReadCounts=new double[noContigs];

	for(int i=0;i<noContigs;i++)
	{
		contigReadCounts[i]=0;
	}
    
	statFile=fopen("stat.txt","w");
    
 //   scaffoldFile=fopen("scaffold.txt","w");
    
/*    
	readCounts=new int*[noContigs];
	
	for (int i = 0; i < noContigs; i++)
		readCounts[i] = new int[noContigs];
    
    
    
	for(int i=0;i<noContigs;i++)
	{
		for(int j=0;j<noContigs;j++)
		{
			readCounts[i][j]=0;
		}
	}
 */   
	
	mapFile1=fopen(mapFileName1, "r");
    
	if (mapFile1 == NULL)
	{
		printf("Can't open map file\n");
		exit(1);
	}

    mapFile2=fopen(mapFileName2, "r");
    
	if (mapFile2 == NULL)
	{
		printf("Can't open map file\n");
		exit(1);
	}

    
	outFile=fopen("myout.sam","w");
    linkingFile=fopen("linking.sam","w");
    singletonFile=fopen("singletons.sam","w");
    
    
	unFile=fopen("unmapped.txt","w");
    
    //	char *temp,nhstring[200];
    
    
	int it=0;
    
    int end=0;
    
    
	char preqname1[1000];
	char preqname2[1000];
	
	strcpy(preqname1,"*");
	strcpy(preqname2,"*");
    
    char qname[1000];
    int flag_segment;
    
    // read map file
    
	while(fgets(line1, MAX_FILE_READ, mapFile1)!=NULL)
	{
		
		
		if(line1[0]=='@')
			continue;
        
        fgets(line2, MAX_FILE_READ, mapFile2);
        while(line2[0]=='@')
        {
            fgets(line2, MAX_FILE_READ, mapFile2);
        }
        
		it++;
        
        
        // check for linking reads
		
		SAM *read1=getSAM(line1);
        SAM *read2=getSAM(line2);

        
        while(strcmp(preqname1,read1->qname)==0)
        {
            mixedReads1.push_back(read1);
            strcpy(preqname1,read1->qname);
            
            if(fgets(line1, MAX_FILE_READ, mapFile1)!=NULL)
            {
                
                read1=getSAM(line1);
                
            }
            else
            {
                end=1;
            }
        }
        while(strcmp(preqname2,read2->qname)==0)
        {
            mixedReads2.push_back(read2);
            strcpy(preqname2,read2->qname);
            
            if(fgets(line2, MAX_FILE_READ, mapFile2)!=NULL)
            {
                
                read2=getSAM(line2);
                
            }
            else
            {
                end=1;
            }
        }

 //       printMixedVectors(linkingFile,singletonFile,unFile);

	 if(isJump==1)
        	testMixedVectors_jump(outFile,linkingFile,singletonFile,unFile);
        else
		testMixedVectors(outFile,linkingFile,singletonFile,unFile);

        if(end==0)
        {
            mixedReads1.push_back(read1);
            strcpy(preqname1,read1->qname);
        
            mixedReads2.push_back(read2);
            strcpy(preqname2,read2->qname);
        }
		if(it%10000==0)
			cout<<it<<endl;
      	}
    
    
	if(isJump==1)
     	 	testMixedVectors_jump(outFile,linkingFile,singletonFile,unFile);
       else
		testMixedVectors(outFile,linkingFile,singletonFile,unFile);
    
/*    
    for(int i=0;i<noContigs;i++)
        fprintf(scaffoldFile,";%s",contigNames[i]);
    fprintf(scaffoldFile,";\n");
    
    
	for(int i=0;i<noContigs;i++)
	{
        fprintf(scaffoldFile,"%s;",contigNames[i]);
		for(int j=0;j<noContigs;j++)
		{
            if(i==j)
                fprintf(scaffoldFile,"%d;",0);
			else if(readCounts[i][j]>10)
                fprintf(scaffoldFile,"%d;",readCounts[i][j]);
            else
                fprintf(scaffoldFile,"%d;",0);
			if(i!=j && readCounts[i][j]>0)
				cout<<contigNames[i]<<" "<<contigNames[j]<<" "<<readCounts[i][j]<<endl;
		}
		fprintf(scaffoldFile,"\n");
	}
*/    
	fclose(mapFile1);
    fclose(mapFile2);

	fclose(outFile);
	fclose(unFile);
//	fclose(scaffoldFile);
	
    //    cout<<totalCount<<"\t"<<unCount<<"\t"<<maxReadLength<<"\t"<<MAX_FRAGMENT_SIZE<<endl;
    
    fprintf(statFile,"%ld %ld %ld %ld",totalCount, unCount, maxReadLength, MAX_FRAGMENT_SIZE);
    
    fclose(statFile);


	FILE *countFile=fopen("counts.txt","w");
    
	for(int i=0;i<noContigs;i++)
	{
		fprintf(countFile,"%d %ld %lf\n",i, contigLengths[i], contigReadCounts[i]);
    		cout<<contigNames[i]<<"\t"<<contigLengths[i]<<"\t"<<contigReadCounts[i]<<endl;
	}
	fclose(countFile);    

	return 0;
}


