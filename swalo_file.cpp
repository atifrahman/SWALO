//
//  swalo_file.cpp
//  swalo
//

#include <iostream>
#include <string>
#include <limits>
#include <vector>
#include <cmath>
#include <cfloat>
using namespace std;

#define _USE_MATH_DEFINES
#include <pthread.h>
#include <math.h>
#include "cgal.h"

#define MAX_REC_LEN 1024

int noContigs=0;
int noReads=0;
long int contigLength=0;
int maxReadLength=0;

long int totalContigLength=0;
vector<char*> contigs;
vector<char*> contigNames;
vector<long int> contigLengths;
vector<long int> contigLengthsOriginals;
double *nc;

FILE *contigFile;
FILE *mapFile;
FILE *summaryFile;
FILE *outFile;

char *contigFileName;
char const *mapFileName;

long int *insertCounts;
int maxInsertSize=0;
int MAX_INSERT_SIZE;
double insertSizeMean;
double insertSizeVar;
double insertSizeSD;
int insertSizeMode;
int insertCutoffMax=0;
int insertCutoffMin=0;
int insertThresholdMax=0;
int insertThresholdMin=0;
int tolerance;
int insertCountMax;

int priorSize=0;

long int *insertCountsMapped;
long int *insertCountsUnmapped;

long int errorTypes[5][5];
long int baseCounts[5];
long int *errorPos;
long int *inPos;
long int *inLengths;
long int *delPos;
long int *delLengths;
long int *readLengths;

long int *effectiveLengths;

double errorTypeProbs[5][5];
double baseErrorRates[5];
double *errorPosDist;
double *inPosDist;
double *inLengthDist;
double *delPosDist;
double *delLengthDist;
double *insertLengthDist;
double *insertLengthDistSmoothed;

int windowSize=12;

double *noErrorProbs;
int tmpCount=0;
int toAlign;

long int erroredReads=0;
long int uniqueMappedReads=0;
long int discardedReads=0;
long int totalCount,unCount;

char tempCigar[500], tempMD[500];
char noErrorCigar[500], noErrorMD[500];


struct MAP
{
	double errorProb;
	int insertSize;
    long int pos1;
    long int pos2;
    long int contigNo1;
    long int contigNo2;
    int readLength1;
    int readLength2;
    bool isSameStrand;
};


struct InsertTable
{
	int insertSize;
    int count;
};

struct Join
{
    double likelihood;
    long int contigNo1;
    long int contigNo2;
    int joinType;
    int gap;
    int isExtra;
	int count;
};

vector<Join *> joins;

vector<Join *> intraJoins1;
vector<Join *> intraJoins2;
vector<Join *> extraJoins1;
vector<Join *> extraJoins2;

vector<Join *> failedJoins;


struct Interval
{
    long int contigNo;
    int contigOrientation;
    long int contigStart;
    long int contigEnd;
    double likelihood;
    int gap;
    int end;
    int selected;
    int realGap;
	int count;
};

vector<Interval *> intervals;

struct QEntry
{
    long int contigNo;
    double likelihood;
};

vector<QEntry *> quarantine;

InsertTable *insertTableMapped;
int insertTableMappedSize;
long int insertTableMappedCount;

InsertTable *insertTableUnmapped;
int insertTableUnmappedSize;
long int insertTableUnmappedCount;

class LinkingMaps
{
	public:
	int size;
	vector<int> indexes;
	vector< vector<MAP *> > vals;

	LinkingMaps()
	{
		size=0;
	}

	inline vector<MAP *> &operator[](int index)
	{

		for(int i=0;i<size;i++)
		{
			if(indexes[i]==index)
			{
		
				return vals[i];
			}
		}
		size++;
		indexes.push_back(index);
		vector<MAP *> * vp=new vector<MAP *>;
		vals.push_back(*vp);
		return vals[size-1];
 	}

};

class Overlaps
{
	public:
	int size;
	vector<int> indexes;
	vector< vector<int> > vals;

	Overlaps()
	{
		size=0;
	}

	inline vector<int> &operator[](int index)
	{

		for(int i=0;i<size;i++)
		{
			if(indexes[i]==index)
			{
		
				return vals[i];
			}
		}
		size++;
		indexes.push_back(index);
		vector<int> * vp=new vector<int>;
		vals.push_back(*vp);
		return vals[size-1];
 	}

};

class AdjListInt
{
	public:
	int size;
	vector<int> indexes;
	vector<int> vals;

	AdjListInt()
	{
		size=0;
	}

	inline void set(int index,int val)
	{
		for(int i=0;i<size;i++)
		{
			if(indexes[i]==index)
			{
		
				vals[i]=val;
				return;
			}
		}
		size++;
		indexes.push_back(index);
		vals.push_back(val);

	}

	void reset()
	{
		for(int i=0;i<size;i++)
		{
			vals[i]=0;
		}
	}

	inline int operator[](int index)
	{

		for(int i=0;i<size;i++)
		{
			if(indexes[i]==index)
			{
		
				return vals[i];
			}
		}
		return 0;
 	}

};

class AdjListDouble
{
	public:
	int size;
	vector<int> indexes;
	vector<double> vals;

	AdjListDouble()
	{
		size=0;
	}

	inline void set(int index,double val)
	{
		for(int i=0;i<size;i++)
		{
			if(indexes[i]==index)
			{
		
				vals[i]=val;
				return;
			}
		}
		size++;
		indexes.push_back(index);
		vals.push_back(val);

	}

	void reset()
	{
		for(int i=0;i<size;i++)
		{
			vals[i]=0;
		}
	}



	inline double operator[](int index)
	{

		for(int i=0;i<size;i++)
		{
			if(indexes[i]==index)
			{
		
				return vals[i];
			}
		}
		return 0;
 	}

};


class Likelihoods
{
	public:
	int size;
	long int contigLength;
	double valUnmapped;
	vector<int> indexes;
	vector<double> vals;

	Likelihoods()
	{
		size=0;
	}

	void init(long int length, double valUn)
	{
		contigLength=length;
		valUnmapped=valUn;
	}

	inline void set(int index,double val)
	{
		for(int i=0;i<size;i++)
		{
			if(indexes[i]==index)
			{
		
				vals[i]=val;
				return;
			}
		}
		size++;
		indexes.push_back(index);
		vals.push_back(val);

	}

	void reset()
	{
		for(int i=0;i<size;i++)
		{
			vals[i]=0;
		}
	}

	inline double operator[](int index)
	{
		for(int i=0;i<size;i++)
		{
			if(indexes[i]==index)
			{
				return vals[i];
			}
		}
		return 0;
 	}


	inline double get(int index);
};


//vector<MAP *> ** linkingMaps;
LinkingMaps *linkingMaps;
//vector<MAP *> ** linkingMapsReversed;
LinkingMaps *linkingMapsReversed;

//vector<int> **overlaps;
Overlaps *overlaps;
//vector<int> **overlapsReversed;
Overlaps *overlapsReversed;

pthread_mutex_t overlaps_mutex = PTHREAD_MUTEX_INITIALIZER;

int isJump=0;
int isConservative=0;

int MIN_GAP=-300;
int MAX_GAP=5000;
int MIN_READ=0;
int MIN_READ_JOIN=1;
int NO_INTERVALS=0;
double *scoreDiffMapped;
double *oldScoreUnmapped;
double *newScoreUnmapped;

double *scoreDiffMapped_2;
double *oldScoreUnmapped_2;
double *newScoreUnmapped_2;

double *scoreDiffMapped_3;
double *oldScoreUnmapped_3;
double *newScoreUnmapped_3;


double *insertProbsSums;

//int **gaps;
AdjListInt *gaps;
//double **likelihoods;
Likelihoods *likelihoods;
//double **counts;
AdjListDouble *counts;


//int **gapsReversed;
AdjListInt *gapsReversed;
//double **likelihoodsReversed;
Likelihoods *likelihoodsReversed;
//double **countsReversed;
AdjListDouble *countsReversed;

int *contigOrientations;
long int *contigSuccessors;
long int *contigSuccessorsReversed;
long int *contigRoots;
long int *contigStarts;
long int *contigEnds;

FILE *scaffoldFile;
FILE *joinsFile;


vector<double> likelihoodsDist;
double likelihoodMean;
double likelihoodSD;

double stretchFactor=1.1;

struct ThreadArg {
    double valMapped;
    double valUnmapped;
    int threadID;
};


void initInsertCounts(int max)
{
	maxInsertSize=max;
	insertCounts=new long int[maxInsertSize];
	for(int i=0;i<maxInsertSize;i++)
	{
		insertCounts[i]=1;
	}
}

void updateInsertCounts(int index)
{
    
	if(index<=0)
		return;
	if(index<maxInsertSize)
	{
		insertCounts[index]++;
	}
	else
	{
		
		if(index>MAX_INSERT_SIZE)
		{
			discardedReads++;
			return;
		}
		int tempInsertSize=max(maxInsertSize*2,index);
		long int *tempCounts=new long int[maxInsertSize];
		for(int i=0;i<maxInsertSize;i++)
		{
			tempCounts[i]=insertCounts[i];
		}
		insertCounts=new long int[tempInsertSize];
		for(int i=0;i<maxInsertSize;i++)
		{
			insertCounts[i]=tempCounts[i];
		}
		for(int i=maxInsertSize;i<tempInsertSize;i++)
		{
			insertCounts[i]=1;
		}
        
		insertCounts[index]++;
		maxInsertSize=tempInsertSize;
		delete []tempCounts;
		
	}
    
}

void initErrorTypes(int readLength)
{
	for(int i=0;i<5;i++)
		for(int j=0;j<5;j++)
			errorTypes[i][j]=1;
    
	for(int i=0;i<5;i++)
		baseCounts[i]=1;
    
	errorPos=new long int[readLength];
	inPos=new long int[readLength];
	inLengths=new long int[readLength];
	delPos=new long int[readLength];
	delLengths=new long int[readLength];
	readLengths=new long int[readLength];
	
	for(int i=0;i<readLength;i++)
	{
		errorPos[i]=1;
		inPos[i]=1;
		inLengths[i]=1;
		delPos[i]=1;
		delLengths[i]=1;
		readLengths[i]=0;
	}
}


int getLength(char *read)
{
	
	int i=0;
	while(read[i])
	{
		if(read[i]=='A')
			baseCounts[0]++;
		else if(read[i]=='C')
			baseCounts[1]++;
		else if(read[i]=='G')
			baseCounts[2]++;
		else if(read[i]=='T')
			baseCounts[3]++;
		else
			baseCounts[4]++;
        
		i++;
	}
	
	return i;
}

long int getContigNo(char *contigName)
{
    for(long int i=0;i<contigNames.size();i++)
    {
		if(strcmp(contigNames[i],contigName)==0)
			return i;
	}
	return -1;
    
}



void processErrorTypes(char *cigar, char *md, char *read, int strandNo)
{
    
	int readLength=getLength(read);
	readLengths[readLength-1]++;
    
	if(strcmp(md,noErrorCigar)!=0)
		erroredReads++;
	else
		return;
    
    
	unsigned long mdLength=strlen(md)-5;
	unsigned long tempLength=0;
    
	char *temp;
	int index=0,totalLength=0;
	
    
	int curIndex=0;
	int *inserts=new int[readLength];
	
	for(int i=0;i<readLength;i++)
	{
		inserts[i]=0;
	}
    
	unsigned long cigarLength=strlen(cigar);
    //	char *tempCigar=new char[cigarLength];
	char cigarChar;
    
	strcpy(tempCigar,cigar);
    
    
	temp=strtok(tempCigar,"IDMS^\t\n ");
    
	while(temp!=NULL)
	{
        
		tempLength=atoi(temp);
		totalLength+=strlen(temp);
		cigarChar=cigar[totalLength];
        
		if(cigarChar=='M')
		{
			index+=tempLength;
			curIndex+=tempLength;
		}
		else if(cigarChar=='I' || cigarChar=='S')
		{
			if(strandNo==0)
			{
				inPos[index]++;
				inLengths[tempLength-1]++;
                
			}
			else
			{
                
				inPos[readLength-index-1]++;
				inLengths[tempLength-1]++;
			}
            
			inserts[curIndex]=tempLength;
            
			index+=tempLength;
		}
		else if(cigarChar=='D' )
		{
			if(strandNo==0)
			{
				delPos[index]++;
				delLengths[tempLength-1]++;
                
			}
			else
			{
                
				delPos[readLength-index-1]++;
				delLengths[tempLength-1]++;
			}
		}
		totalLength++;
		temp=strtok(NULL,"IDMS^\t\n ");
	}
    
	
	strcpy(tempMD,md);
    
	strtok(tempMD,":");
	strtok(NULL,":");
    
	
	index=0,totalLength=0,tempLength=0;
    
	int f, t;
    
	while((temp=strtok(NULL,"ACGTN^\t\n "))!=NULL)
	{
		tempLength=strlen(temp);
        
        
		totalLength+=tempLength;
		
        
		if(totalLength<mdLength)
		{
			char from=md[5+totalLength];
			
			
			if(from=='^')
			{
				totalLength++;
				index+=atoi(temp);
				for(int i=totalLength;i<mdLength;i++)
				{
					from=md[5+totalLength];
					if(from=='A' || from=='C' || from=='G' || from=='T'|| from=='N')
						totalLength++;
					else
						break;
                    
				}
			}
			else if(from=='A' || from=='C' || from=='G' || from=='T'|| from=='N')
			{
				totalLength++;
				index+=atoi(temp)+1;
				
				
				
				curIndex=0;
				for(int i=0;i<index;i++)
				{
					curIndex+=inserts[i];
				}
				char to=read[index-1+curIndex];
                
				if(strandNo==0)
					errorPos[index-1+curIndex]++;
				else
					errorPos[readLength-index-curIndex]++;
                
                
				switch(from)
				{
					case 'A':
						f=0;
						break;
					case 'C':
						f=1;
						break;
					case 'G':
						f=2;
						break;
					case 'T':
						f=3;
						break;
					default:
						f=4;
				}
                
				switch(to)
				{
					case 'A':
						t=0;
						break;
					case 'C':
						t=1;
						break;
					case 'G':
						t=2;
						break;
					case 'T':
						t=3;
						break;
					default:
						t=4;
				}
                
				if(f==t)
				{
                    
				}
				else
                    
					errorTypes[f][t]++;
                
			}
			else
				break;
		}
		
	}
	delete []inserts;
    
}

void computeProbabilites()
{
	
	int errorCount=0;
    
	for(int i=0;i<5;i++)
	{
		errorCount=0;
		for(int j=0;j<5;j++)
		{
			errorCount+=errorTypes[i][j];
		}
		for(int j=0;j<5;j++)
		{
			errorTypeProbs[i][j]=(double)errorTypes[i][j]/errorCount;
		}
		
		baseErrorRates[i]=errorCount/(double)baseCounts[i];
	}
    
	double sum=0;
	for(int i=0;i<4;i++)
		sum+=baseErrorRates[i];
    
	for(int i=0;i<4;i++)
	{
		baseErrorRates[i]=4*baseErrorRates[i]/sum;
	}
    
	baseErrorRates[4]=1;
    
	for(int i=maxReadLength-1;i>0;i--)
	{
		readLengths[i-1]=readLengths[i]+readLengths[i-1];
	}
    
	errorPosDist=new double[maxReadLength];
    
	for(int i=0;i<maxReadLength;i++)
	{
		errorPosDist[i]=(double)errorPos[i]/readLengths[i];
	}
    
	inPosDist=new double[maxReadLength];
    
	for(int i=0;i<maxReadLength;i++)
	{
		inPosDist[i]=(double)inPos[i]/readLengths[i];
	}
    
	inLengthDist=new double[maxReadLength];
    
	int inCount=0;
    
	for(int i=0;i<maxReadLength;i++)
	{
		inCount+=inLengths[i];
	}
	
	for(int i=0;i<maxReadLength;i++)
	{
		inLengthDist[i]=(double)inLengths[i]/inCount;
	}
    
	delPosDist=new double[maxReadLength];
    
	for(int i=0;i<maxReadLength;i++)
	{
		delPosDist[i]=(double)delPos[i]/readLengths[i];
	}
    
	delLengthDist=new double[maxReadLength];
    
	int delCount=0;
    
	for(int i=0;i<maxReadLength;i++)
	{
		delCount+=delLengths[i];
	}
	
	for(int i=0;i<maxReadLength;i++)
	{
		delLengthDist[i]=(double)delLengths[i]/delCount;
	}
	
    
	insertLengthDist=new double[maxInsertSize];
    
	long int insCount=discardedReads;
    
	sum=0;
    
	for(int i=0;i<maxInsertSize;i++)
	{
		insCount+=(insertCounts[i]-1);
		sum+=i*(insertCounts[i]-1);
        
	}
	insertSizeMean=sum/insCount;
    
   
    
	sum=0;
    
	for(int i=0;i<maxInsertSize;i++)
	{
		insertLengthDist[i]=(double)insertCounts[i]/insCount;
        
		sum+=(insertCounts[i]-1)*(insertSizeMean-i)*(insertSizeMean-i);
	}
    
	insertSizeVar=sum/insCount;
    
    insertSizeSD=sqrt(insertSizeVar);
    
	noErrorProbs=new double[maxReadLength];
    
	double noErrorProb=1.0;
    
    
	for(int i=0;i<maxReadLength;i++)
	{
		noErrorProb*=(1-errorPosDist[i]-inPosDist[i]-delPosDist[i]);
		noErrorProbs[i]=noErrorProb;
  	}
    
	effectiveLengths=new long int[maxInsertSize];
    
	for(int i=0;i<maxInsertSize;i++)
	{
		effectiveLengths[i]=-1;
	}
    
	long int totalContigLength=0;
	for(int i=0;i<contigLengths.size();i++)
	{
		totalContigLength+=contigLengths[i];
	}
	effectiveLengths[0]=totalContigLength;
	
    
    insertCountsMapped=new long int[maxInsertSize];
    insertCountsUnmapped=new long int[maxInsertSize];
    
    
    for(long int i=0;i<maxInsertSize;i++)
    {
        insertCountsMapped[i]=0;
        insertCountsUnmapped[i]=0;
    }
    
    insertLengthDistSmoothed=new double[maxInsertSize];
    double windowSum=0;
    
    for(int i=0;i<windowSize;i++)
    {
        insertLengthDistSmoothed[i]=insertLengthDist[i];
        
    }
    
    for(int i=0;i<2*windowSize+1;i++)
    {
        windowSum+=insertLengthDist[i];
        
    }
    insertLengthDistSmoothed[windowSize]=windowSum/(2*windowSize+1);
    
    for(int i=windowSize+1;i<maxInsertSize-windowSize;i++)
    {
        windowSum-=insertLengthDist[i-windowSize-1];
        windowSum+=insertLengthDist[i+windowSize];
        insertLengthDistSmoothed[i]=windowSum/(2*windowSize+1);
    }
    for(int i=maxInsertSize-windowSize;i<maxInsertSize;i++)
    {
        insertLengthDistSmoothed[i]=insertLengthDist[i];
    }

	for(int i=0;i<maxInsertSize;i++)
	{
		insertLengthDistSmoothed[i]=insertLengthDistSmoothed[i]-1/(double)(insCount)+(1/(double)maxInsertSize)/(double)(insCount+1);

	}
    
    int count=0;
    
    for(int i=insertSizeMean;i<maxInsertSize;i++)
    {
        if(insertCounts[i]<=1)
        {
            count++;
            if(count==10)
            {
                insertCutoffMax=i;
                break;
            }
        }
        else
        {
            count=0;
        }
        
    }
    
    count=0;
    
    for(int i=insertSizeMean;i>=0;i--)
    {
        if(insertCounts[i]<=1)
        {
            count++;
            if(count==10)
            {
                insertCutoffMin=i;
                break;
            }
        }
        else
        {
            count=0;
        }
        
    }
    
    insertCountMax=0;
    for(int i=insertCutoffMin;i<insertCutoffMax;i++)
    {
        if(insertCounts[i]>insertCountMax)
        {
            insertCountMax=insertCounts[i];
            insertSizeMode=i;
            
        }
    }

    
    count=0;

    for(int i=insertSizeMean;i<maxInsertSize;i++)
    {
        if(insertCounts[i]<=max(insertCountMax/100,2))
        {
            count++;
            if(count==2)
            {
                insertThresholdMax=i;
                break;
            }
        }
        else
        {
            count=0;
        }
        
    }
    
    count=0;
    
    for(int i=insertSizeMean;i>=0;i--)
    {
        if(insertCounts[i]<=max(insertCountMax/100,2))
        {
            count++;
            if(count==2)
            {
                insertThresholdMin=i;
                break;
            }
        }
        else
        {
            count=0;
        }
        
    }
    
    //mean used instead of mode
    
    double insertSum=0;
    double insertCount=0;
    for(int i=insertCutoffMin;i<insertCutoffMax;i++)
    {
        insertCount+=insertCounts[i]-1;
        insertSum+=(insertCounts[i]-1)*i;
    }
    insertSizeMode=insertSum/insertCount;

    insertCutoffMax=insertThresholdMax;
    insertCutoffMin=insertThresholdMin;

}

double dnorm(double x,double mean, double variance)
{
	double val=1/sqrt(M_PI*2*variance);
	val*=exp(-((x-mean)*(x-mean))/(2*variance));
	return val;
}

void processMapping(char *line)
{
	
	char * temp;
	char *qname, *rname, *mapq;
	int	pos,flag;
	char * cigar, * readString; // * md, *nhstring;
	long int contigNo;    
		

	char md[500];
	char nhstring[500];
    
	int nh;
    
	int strandNo=0;
    
    
	qname=strtok(line,"\t");
	
    
	temp=strtok(NULL,"\t");
	flag=atoi(temp);
	
    
	strandNo=(flag&16)>>4;
    
    rname=strtok(NULL,"\t");
	
    
	temp=strtok(NULL,"\t");
	pos=atoi(temp);
	
	
	cigar=strtok(NULL,"\t");
	
	
	temp=strtok(NULL,"\t");
	
	readString=strtok(NULL,"\t");
    
	int insertSize=atoi(temp);
    
	while((temp=strtok(NULL,"\t\n"))!=NULL)
	{
		if(temp[0]=='M' && temp[1]=='D')
		{
			strcpy(md,temp);
		}
		else if(temp[0]=='I' && temp[1]=='H')
		{
			strcpy(nhstring,(temp+5));
			nh=atoi(nhstring) ;
		}
        
	}
    
	
	if(nh==1 && md[5]!='^')
	{
		contigNo=getContigNo(rname);
		if(contigLengths[contigNo]>1200)
			updateInsertCounts(insertSize);
		processErrorTypes(cigar,md,readString,strandNo);
		uniqueMappedReads++;
	}
	
    
}

long int getEffectiveLength(int insertSize)
{
	if(insertSize<0)
		return effectiveLengths[0];
    
	if(insertSize>=maxInsertSize)
	{
		long int effectiveLength=0;
		for(int i=0;i<contigLengths.size();i++)
		{
			if(contigLengths[i]>=insertSize)
				effectiveLength+=(contigLengths[i]-insertSize+1);
		}
		return effectiveLength;
        
	}
	if(effectiveLengths[insertSize]==-1)
	{
		long int effectiveLength=0;
		for(int i=0;i<contigLengths.size();i++)
		{
			if(contigLengths[i]>=insertSize)
				effectiveLength+=(contigLengths[i]-insertSize+1);
		}
		effectiveLengths[insertSize]=effectiveLength;
	}
	return effectiveLengths[insertSize];
}

long double computeErrorProb(char *cigar, char *md, char *read, int strandNo)
{
    
    unsigned long readLength=strlen(read);
	
    
	long double errorProb=noErrorProbs[readLength-1];


    
	if(md[5]=='^')
		return errorProb;
	
	
	char tempMD[1000], tempCigar[1000];
    
	unsigned long mdLength=strlen(md)-5;
	unsigned long tempLength=0;
	
	char *temp;
	int index=0,totalLength=0;
	
    
	int curIndex=0;
	int *inserts=new int[readLength];
	
	for(int i=0;i<readLength;i++)
	{
		inserts[i]=0;
	}
    
    //	int cigarLength=strlen(cigar);
	char cigarChar;
    
	strcpy(tempCigar,cigar);
    
	temp=strtok(tempCigar,"IDM^\t\n ");
    
	while(temp!=NULL)
	{
        
		tempLength=atoi(temp);
		totalLength+=strlen(temp);
		cigarChar=cigar[totalLength];
        
		if(cigarChar=='M')
		{
			index+=tempLength;
			curIndex+=tempLength;
		}
		else if(cigarChar=='I')
		{
			unsigned long i;
			if(strandNo==0)
			{
				//look up insert probs
				i=index;
				
			}
			else
			{
				i=readLength-index-1;
			}
            
			errorProb=errorProb*inPosDist[i]*inLengthDist[tempLength-1]/(1-errorPosDist[i]-inPosDist[i]-delPosDist[i]);
            
			inserts[curIndex]=tempLength;
            
			index+=tempLength;
		}
		else if(cigarChar=='D')
		{
			unsigned long i;
			if(strandNo==0)
			{
				i=index;
                //	look up delete probs
			}
			else
			{
				i=readLength-index-1;
			}
            
			errorProb=errorProb*delPosDist[i]*delLengthDist[tempLength-1]/(1-errorPosDist[i]-inPosDist[i]-delPosDist[i]);
		}
		totalLength++;
		temp=strtok(NULL,"IDM^\t\n ");
	}
    
    
	strcpy(tempMD,md);
    
	strtok(tempMD,":");
	strtok(NULL,":");
    
	index=0,totalLength=0,tempLength=0;
    
	int f, t;
    
	while((temp=strtok(NULL,"ACGTN^\t\n "))!=NULL)
	{
		tempLength=strlen(temp);
        
		totalLength+=tempLength;
		
		if(totalLength<mdLength)
		{
			char from=md[5+totalLength];
            
			if(from=='^')
			{
				totalLength++;
				index+=atoi(temp);
				for(int i=totalLength;i<mdLength;i++)
				{
					from=md[5+totalLength];
					if(from=='A' || from=='C' || from=='G' || from=='T'|| from=='N')
						totalLength++;
					else
						break;
				}
			}
			else if(from=='A' || from=='C' || from=='G' || from=='T'|| from=='N')
			{
				totalLength++;
				index+=atoi(temp)+1;
				
                
				curIndex=0;
				for(int i=0;i<index;i++)
				{
					curIndex+=inserts[i];
				}
				char to=read[index-1+curIndex];
                
				int i;
				if(strandNo==0)
					i=index-1+curIndex;
				else
					i=readLength-index-curIndex;
                
                
				errorProb=errorProb*errorPosDist[i]/(1-errorPosDist[i]-inPosDist[i]-delPosDist[i]);
                
                
				switch(from)
				{
					case 'A':
						f=0;
						break;
					case 'C':
						f=1;
						break;
					case 'G':
						f=2;
						break;
					case 'T':
						f=3;
						break;
					default:
						f=4;
				}
                
				switch(to)
				{
					case 'A':
						t=0;
						break;
					case 'C':
						t=1;
						break;
					case 'G':
						t=2;
						break;
					case 'T':
						t=3;
						break;
					default:
						t=4;
				}
                
				if(f==t)
				{
					
                    
				}
				else
				{
					//errorTypeProb
					errorProb*=baseErrorRates[f]*errorTypeProbs[f][t];
				}
			}
			else
				break;
            
		}
		
	}
	
	delete []inserts;
	return errorProb;
}


double computeLikelihood(char const *file)
{
    
	mapFile=fopen(file, "r");
	char *line1= new char[MAX_REC_LEN];
	char *line2= new char[MAX_REC_LEN];
    
	char *qname1,*qname2,preqname1[500],preqname2[500];
    
	int MAX_FILE_READ=MAX_REC_LEN/sizeof(line1[0]);
    
	long double sum=0.0;
	long double logsum=0.0;
    
	char * temp;
	char *rname1, *rname2;
	int	pos1,pos2,flag,strandNo1, strandNo2, insertSize1, insertSize2;
	char *cigar1, *cigar2, *readString1, *readString2, md1[1000], md2[1000];
    
    
	long double insertSizeProb;
    
	long double errorProb1, errorProb2;
    
    int tempInsertSize=0;
    long double tempProb=0;
    
	preqname1[0]=preqname2[0]='*';
	preqname1[1]=preqname2[1]=0;
    
    
	int it=0;
    
	while(fgets(line1, MAX_FILE_READ, mapFile)!=NULL)
	{
		if(line1[0]=='@')
			continue;
        //????
		if(fgets(line2, MAX_FILE_READ, mapFile)==NULL)
			break;
        
		qname1=strtok(line1,"\t");
		
        temp=strtok(NULL,"\t");
        flag=atoi(temp);
		
        
		strandNo1=(flag&16)>>4;
        
        temp=strtok(NULL,"\t");
		
		temp=strtok(NULL,"\t");
		pos1=atoi(temp);
        
		
		cigar1=strtok(NULL,"\t");
        
		
		temp=strtok(NULL,"\t");
		
		
		
		insertSize1=atoi(temp);
        
        
		readString1=strtok(NULL,"\t");
        
        
		while((temp=strtok(NULL,"\t\n"))!=NULL)
		{
			if(temp[0]=='M' && temp[1]=='D')
			{
				strcpy(md1,temp);
			}
		}
        
//        		cout<<insertSize1<<" "<<cigar1<<" "<<md1<<endl;
		
        //second of the pair
        
		qname2=strtok(line2,"\t");
		temp=strtok(NULL,"\t");
		flag=atoi(temp);
		
        
		strandNo2=(flag&16)>>4;
        
        temp=strtok(NULL,"\t");
		
		temp=strtok(NULL,"\t");
		pos2=atoi(temp);
        
		
        
		cigar2=strtok(NULL,"\t");
        
		
		temp=strtok(NULL,"\t");
		
		insertSize2=atoi(temp);
        
		readString2=strtok(NULL,"\t");
        
        
		while((temp=strtok(NULL,"\t\n"))!=NULL)
		{
			if(temp[0]=='M' && temp[1]=='D')
			{
				strcpy(md2,temp);
			}
		}
        
//        		cout<<insertSize2<<" "<<cigar2<<" "<<md2<<endl;
        
		
		
		int insertSize=max(insertSize1, insertSize2);
		
        
		insertSizeProb=0;
        
		if(insertSize>=0 && insertSize<maxInsertSize)
		{
			insertSizeProb=insertLengthDist[insertSize];
		}
        
		if(insertSizeProb==0)
		{
			insertSizeProb=1/(double)uniqueMappedReads;
		}
		
        
		errorProb1=computeErrorProb(cigar1,md1,readString1,strandNo1);
        
        
		errorProb2=computeErrorProb(cigar2,md2,readString2,strandNo2);
        


        
		long int totalEffectiveLength=getEffectiveLength(insertSize);
        
        
		long double prob=(1/(long double)(totalEffectiveLength))*insertSizeProb*errorProb1*errorProb2;
        
        
        
//        cout<<errorProb1<<" "<<errorProb2<<" "<<insertSizeProb<<" "<<prob<<endl;
        
        
        
        
		if(strcmp(qname1,preqname1)==0 && strcmp(qname2,preqname2)==0)
		{
            if(tempProb<prob)
            {
                tempProb=prob;
                tempInsertSize=insertSize;
		  tempInsertSize=tempInsertSize<0?0:tempInsertSize;

            }
            
			sum+=prob;
            
		}
		else if(strcmp("*",preqname1)!=0 && strcmp("*",preqname2)!=0)
		{
			if(sum<1e-320 || isnan(sum))
			{
				sum=1e-320;
			}
			logsum+=log(sum);
            
            
            
            if(tempInsertSize>=maxInsertSize)
                insertCountsMapped[maxInsertSize-1]++;
            else
                insertCountsMapped[tempInsertSize]++;
			
			sum=prob;
            
            tempProb=prob;
            tempInsertSize=insertSize;
		  tempInsertSize=tempInsertSize<0?0:tempInsertSize;


		}
		else
		{
			sum=prob;
			
            tempProb=prob;
            tempInsertSize=insertSize;
		  tempInsertSize=tempInsertSize<0?0:tempInsertSize;

			
		}
        
		strcpy(preqname1,qname1);
		strcpy(preqname2,qname2);
		it++;
        
		
		if(isinf( logsum ))
		{
			cout<<it<<endl;
			exit(1);
		}
        
        
	}
	if(sum!=0)
    {
        if(sum<1e-320 || isnan(sum))
        {
            sum=1e-320;
        }
		logsum+=log(sum);
    }
    
	fclose(mapFile);
	
	return logsum;
}

void printHelp()
{
    
	cout<<"swalo-"<<VERSION<<endl;
	cout<<"----------------"<<endl;
	cout<<endl;
	cout<<"swaloFile - scaffold contigs with precomputed scaffold graph (Use for multiple read libraries)"<<endl;
	cout<<"Note: For this first run SWALO on each library separately in different directories using swalo. This will generate three files with names staring with the same prefix and a text file prefixes.txt containing the prefixes. Now create a new directory, copy into the new directory the three files starting with prefix from each directory and create a new text file prefixes.txt in the new directory containing one prefix in each line. Then in the new directory run swaloFile"<<endl;
	cout<<"Usage:"<<endl;
	cout<<"swalo <contigFile> [options]"<<endl;
	cout<<endl;
	cout<<"Required arguments:"<<endl;
	cout<<"<contigFile>\t Contig file in FASTA format"<<endl;
	cout<<endl;
	cout<<"Options:"<<endl;
	cout<<"-h [--help]\t\t Prints this message"<<endl;
	cout<<"--conservative or -c\t\t Use a conservative mode (please see paper). Recommended if standard deviation of insert library is high (>1000)"<<endl;

	cout<<endl;
	cout<<"Output: "<<endl;
	cout<<"scaffolds.fa will contain the scaffolds in FASTA format."<<endl;
	cout<<endl;
	exit(1);
    
}

void reverse(char *reverse, char *read)
{
    
	char ch='A';
	int readLength=strlen(read);
	for(int i=1; i<=readLength;i++)
	{
		ch=read[i-1];
		if(ch=='A' || ch=='a')
			reverse[readLength-i]='T';
		else if(ch=='C' || ch=='c')
			reverse[readLength-i]='G';
		else if(ch=='G' || ch=='g')
			reverse[readLength-i]='C';
		else if(ch=='T' || ch=='t')
			reverse[readLength-i]='A';
		else
			reverse[readLength-i]='N';
	}
	reverse[readLength]='\0';
    
}


int getDistance(char *s, char * t, int sStart, int sEnd, int tStart, int tEnd, int ** dis)
{
	int m=sEnd-sStart+1;
    int n=tEnd-tStart+1;
	
    
	int val1;
	int val2;
	int val3;

	
    int i,j;
    
	for(i=0;i<=m;i++)
        dis[i][0]=i;
    
	for(j=0;j<=n;j++)
        dis[0][j]=j;
    
	for(i=1;i<=m;i++)
	{
		for(j=1;j<=n;j++)
		{
			val2=dis[i-1][j]+1;
			val3=dis[i][j-1]+1;
			val2=val2<val3?val2:val3;
			
            val1=dis[i-1][j-1];
            if(s[i-1+sStart]!=t[j-1]+tStart)
			{
				val1++;
                
			}
			val3=val1<val2?val1:val2;
			
            dis[i][j]=val3;
			
		}
		
	}
    
    
    
    return 	dis[m][n];
    
}

long int getContigRoot(long int contigNo)
{
    while(contigNo!=contigRoots[contigNo])
    {
        contigNo=contigRoots[contigNo];
        
    }
    return contigNo;
    
}

int getOverlap(char *s, char *t, double error)
{
    int sLength=strlen(s);
    int tLength=strlen(t);
    int distance=0;
    
    int **dis;
    
    dis=new int*[sLength+1];
    
    for(int i=0;i<=sLength;i++)
    {
        dis[i]=new int[tLength+1];
    }
    
    
    for(int i=sLength>tLength?sLength-tLength:0;i<sLength;i++)
    {
        distance=getDistance(s,t , i, sLength-1, 0, sLength-i-1,dis);
        
        if(distance<=error*(sLength-i))
        {
            for(int j=0;j<=sLength;j++)
            {
                delete [] dis[j];
            }
            
            delete [] dis;
            
            return (sLength-i);
        }
    }
    
    for(int i=0;i<=sLength;i++)
    {
        delete [] dis[i];
    }
    
    delete [] dis;
    
    return 0;
}

void computeOverlaps(long int contigNo1, long int contigNo2, double error, int maxOverlap, int isReverse)
{
    char *s;
    char *t;
    long int contigRoot1=getContigRoot(contigNo1);
    long int contigRoot2=getContigRoot(contigNo2);
    
    if(isReverse==1)
    {
        if(contigNo1<contigNo2)
        {
            s=contigs[contigRoot1];
            t=new char[contigLengths[contigRoot2]+1];
            reverse(t,contigs[contigRoot2]);
        }
        else
        {
            s=new char[contigLengths[contigRoot1]+1];
            reverse(s,contigs[contigRoot1]);
            t=contigs[contigRoot2];
        }
        
    }
    else
    {
        s=contigs[contigRoot1];
        t=contigs[contigRoot2];
    }
    int m=contigLengths[contigRoot1]<maxOverlap?contigLengths[contigRoot1]:maxOverlap;
    int n=contigLengths[contigRoot2]<maxOverlap?contigLengths[contigRoot2]:maxOverlap;
    int distance=0;
    
    int **dis;
    
    dis=new int*[m+1];
    
    for(int i=0;i<=m;i++)
    {
        dis[i]=new int[n+1];
    }
    
    long int minLength=m<n?m:n;
    
    for(int i=0;i<minLength;i++)
    {
        distance=getDistance(s,t , contigLengths[contigRoot1]-minLength+i, contigLengths[contigRoot1]-1, 0, minLength-1-i, dis);
        
        if(distance<=error*(minLength-i))
        {
	//     pthread_mutex_lock(&overlaps_mutex);		
            if(isReverse==0)
                overlaps[contigNo1][contigNo2].push_back(minLength-i);
            else
                overlapsReversed[contigNo1][contigNo2].push_back(minLength-i);
	//     pthread_mutex_unlock(&overlaps_mutex);		
	     
        }
    }
//    pthread_mutex_lock(&overlaps_mutex);		
    if(isReverse==0)
        overlaps[contigNo1][contigNo2].push_back(0);
    else
        overlapsReversed[contigNo1][contigNo2].push_back(0);
//    pthread_mutex_unlock(&overlaps_mutex);		
    

    for(int i=0;i<=m;i++)
    {
        delete [] dis[i];
    }
    
    delete [] dis;
    
    if(isReverse==1)
    {
        if(contigNo1<contigNo2)
        {
            delete []t;
        }
        else
        {
            delete []s;
        }
    }
    
}

int getOverlap(long int contigNo1, long int contigNo2, int joinType)
{
    if(contigNo1==contigNo2)
        return 0;
    
    if(joinType==0)
    {
        if(overlaps[contigNo1][contigNo2].size()==0)
        {
            computeOverlaps(contigNo1, contigNo2, 0, -MIN_GAP, 0);
        }
        return overlaps[contigNo1][contigNo2][0];
    }
    else if(joinType==1)
    {
        long int tempContigNo1=contigNo1<contigNo2?contigNo1:contigNo2;
        long int tempContigNo2=contigNo1>contigNo2?contigNo1:contigNo2;
        
        if(overlapsReversed[tempContigNo1][tempContigNo2].size()==0)
        {
            computeOverlaps(tempContigNo1, tempContigNo2, 0, -MIN_GAP, 1);
        }
        return overlapsReversed[tempContigNo1][tempContigNo2][0];
    }
    else if(joinType==2)
    {
        long int tempContigNo1=contigNo1>contigNo2?contigNo1:contigNo2;
        long int tempContigNo2=contigNo1<contigNo2?contigNo1:contigNo2;
        
        if(overlapsReversed[tempContigNo1][tempContigNo2].size()==0)
        {
            computeOverlaps(tempContigNo1, tempContigNo2, 0, -MIN_GAP, 1);
        }
        return overlapsReversed[tempContigNo1][tempContigNo2][0];
    }

    
    return 0;
}


int getOverlap(long int contigNo1, long int contigNo2, int orientation1, int orientation2)
{
    if(orientation1==0 && orientation2==0)
    {
        return getOverlap(contigNo1, contigNo2, 0);
    }
    else if(orientation1==1 && orientation2==1)
    {
        return getOverlap(contigNo2, contigNo1, 0);
    }
    else if(orientation1==0 && orientation2==1)
    {
        return getOverlap(contigNo1, contigNo2, 1);
    }
    else if(orientation1==1 && orientation2==0)
    {
        return getOverlap(contigNo1, contigNo2, 2);
    }

    return 0;
}

void computeScores()
{
    
    scoreDiffMapped=new double[MAX_GAP-MIN_GAP+1];
    oldScoreUnmapped=new double[MAX_GAP-MIN_GAP+1];
    newScoreUnmapped=new double[MAX_GAP-MIN_GAP+1];
    
    scoreDiffMapped_2=new double[insertCutoffMax+1];
    oldScoreUnmapped_2=new double[insertCutoffMax+1];
    newScoreUnmapped_2=new double[insertCutoffMax+1];
    
    scoreDiffMapped_3=new double[insertCutoffMax+1];
    oldScoreUnmapped_3=new double[insertCutoffMax+1];
    newScoreUnmapped_3=new double[insertCutoffMax+1];

    
    
    double oldSum;
    double newSum;
    long int effectiveLength;
    
    int gap=MIN_GAP;
    for(int i=0;i<=MAX_GAP-MIN_GAP;i++)
    {
        oldSum=0;
        newSum=0;
        effectiveLength=0;
        
        for(int j=0;j<insertTableMappedSize;j++)
        {
            effectiveLength=getEffectiveLength(insertTableMapped[j].insertSize);
            
            oldSum+=insertTableMapped[j].count*log(1/(double)effectiveLength);
            
            newSum+=insertTableMapped[j].count*log(1/(double)(effectiveLength+insertTableMapped[j].insertSize+gap));
        }
        
        scoreDiffMapped[i]=oldSum-newSum;
        
        oldSum=0;
        newSum=0;
        effectiveLength=0;
        
        for(int j=0;j<insertTableUnmappedSize;j++)
        {
            effectiveLength=getEffectiveLength(insertTableUnmapped[j].insertSize);
            
            oldSum+=insertTableUnmapped[j].count*log(1/(double)effectiveLength);
            
            newSum+=insertTableUnmapped[j].count*log(1/(double)(effectiveLength+insertTableUnmapped[j].insertSize+gap));
        }
        
        oldScoreUnmapped[i]=oldSum;
        newScoreUnmapped[i]=newSum;
        
        gap++;
    }
    
    for(int i=0;i<=insertCutoffMax;i++)
    {
        oldSum=0;
        newSum=0;
        effectiveLength=0;
        
        for(int j=0;j<insertTableMappedSize;j++)
        {
            effectiveLength=getEffectiveLength(insertTableMapped[j].insertSize);
            
            oldSum+=insertTableMapped[j].count*log(1/(double)effectiveLength);
            
            newSum+=insertTableMapped[j].count*log(1/(double)(effectiveLength+i));
        }
        
        scoreDiffMapped_2[i]=oldSum-newSum;
        
        oldSum=0;
        newSum=0;
        effectiveLength=0;
        
        for(int j=0;j<insertTableUnmappedSize;j++)
        {
            effectiveLength=getEffectiveLength(insertTableUnmapped[j].insertSize);
            
            oldSum+=insertTableUnmapped[j].count*log(1/(double)effectiveLength);
            
            newSum+=insertTableUnmapped[j].count*log(1/(double)(effectiveLength+i));
        }
        
        oldScoreUnmapped_2[i]=oldSum;
        newScoreUnmapped_2[i]=newSum;
        
    }
    
    for(int i=0;i<=insertCutoffMax;i++)
    {
        oldSum=0;
        newSum=0;
        effectiveLength=0;
        
        for(int j=0;j<insertTableMappedSize;j++)
        {
            if(insertTableMapped[j].insertSize<i)
            {
                effectiveLength=getEffectiveLength(insertTableMapped[j].insertSize);
                
                oldSum+=insertTableMapped[j].count*log(1/(double)effectiveLength);
                
                newSum+=insertTableMapped[j].count*log(1/(double)(effectiveLength+i-insertTableMapped[j].insertSize));
            }

        }
        
        scoreDiffMapped_3[i]=oldSum-newSum;
        
        oldSum=0;
        newSum=0;
        effectiveLength=0;
        
        for(int j=0;j<insertTableUnmappedSize;j++)
        {
            if(insertTableUnmapped[j].insertSize<i)
            {
                effectiveLength=getEffectiveLength(insertTableUnmapped[j].insertSize);
                
                oldSum+=insertTableUnmapped[j].count*log(1/(double)effectiveLength);
                
                newSum+=insertTableUnmapped[j].count*log(1/(double)(effectiveLength+i-insertTableUnmapped[j].insertSize));
            }
        }
        
        oldScoreUnmapped_3[i]=oldSum;
        newScoreUnmapped_3[i]=newSum;
        
    }
    

	double sum=0;
    
    insertProbsSums=new double[MAX_INSERT_SIZE];
    
    sum=0;
    for(int i=0;i<MAX_INSERT_SIZE;i++)
    {
        
        sum+=insertLengthDistSmoothed[i];
    }
    
    for(int i=0;i<MAX_INSERT_SIZE;i++)
    {
        insertProbsSums[i]=sum;
        sum-=insertLengthDistSmoothed[i];
        
    }

	insertTableMappedCount=0;
	for(int j=0;j<insertTableMappedSize;j++)
       {
		insertTableMappedCount+=insertTableMapped[j].count;
	}

	insertTableUnmappedCount=0;
	for(int j=0;j<insertTableUnmappedSize;j++)
       {
		insertTableUnmappedCount+=insertTableUnmapped[j].count;
	}
    
    
}

double getScore(int gap,int linkingCount, double valUnmapped, long int contigLength1, long contigLength2)
{
    
    int i=gap-MIN_GAP;
    long int newContigLength=gap+contigLength1+contigLength2;
    long int smallContigLength=contigLength1<contigLength2?contigLength1:contigLength2;
    
    if(newContigLength<insertCutoffMax)
    {
        return -scoreDiffMapped_3[newContigLength]-(oldScoreUnmapped_3[newContigLength]-newScoreUnmapped_3[newContigLength])*(unCount-linkingCount)/toAlign-valUnmapped*linkingCount/unCount;
    }
    else if(smallContigLength+gap<insertCutoffMax)
    {
        return -scoreDiffMapped_2[smallContigLength+gap]-(oldScoreUnmapped_2[smallContigLength+gap]-newScoreUnmapped_2[smallContigLength+gap])*(unCount-linkingCount)/toAlign-valUnmapped*linkingCount/unCount;
    }
    else
    {
        return -scoreDiffMapped[i]-(oldScoreUnmapped[i]-newScoreUnmapped[i])*(unCount-linkingCount)/toAlign-valUnmapped*linkingCount/unCount;
    }
    
}

/*
double getScore(int gap,int linkingCount, double valUnmapped, long int contigLength1, long contigLength2)
{
    
    int i=gap-MIN_GAP;
    long int newContigLength=gap+contigLength1+contigLength2;
    long int smallContigLength=contigLength1<contigLength2?contigLength1:contigLength2;
    
    if(newContigLength<insertCutoffMax)
    {
        return -scoreDiffMapped_3[newContigLength]-(scoreDiffMapped_3[newContigLength])*(unCount-linkingCount)/insertTableMappedCount-valUnmapped*linkingCount/unCount;
    }
    else if(smallContigLength+gap<insertCutoffMax)
    {
        return -scoreDiffMapped_2[smallContigLength+gap]-(scoreDiffMapped_2[smallContigLength+gap])*(unCount-linkingCount)/insertTableMappedCount-valUnmapped*linkingCount/unCount;
    }
    else
    {
        return -scoreDiffMapped[i]-(scoreDiffMapped[i])*(unCount-linkingCount)/insertTableMappedCount-valUnmapped*linkingCount/unCount;
    }
    
}
*/
double getInsertProbsSum(long int minInsertSize,long int maxInsertSize)
{
    maxInsertSize++;
    
    if(maxInsertSize<MAX_INSERT_SIZE)
        return insertProbsSums[minInsertSize]-insertProbsSums[maxInsertSize];//+(maxInsertSize-minInsertSize)/(double)uniqueMappedReads;
    else
        return insertProbsSums[minInsertSize];//+(MAX_INSERT_SIZE-minInsertSize)/(double)uniqueMappedReads;
    
}

double getProb(MAP *map)
{
    
    long int contigNo1=map->contigNo1;
    long int contigNo2=map->contigNo2;
    
    int gap;
    
    long int insertMin1, insertMax1, insertMin2, insertMax2;
    
    
    if(map->isSameStrand==false)
    {
        gap=gaps[contigNo1][contigNo2];
        insertMin1=contigLengths[map->contigNo1]-map->pos1+1+gap+map->readLength2;
        
        insertMax1=contigLengths[map->contigNo1]-map->pos1+1+gap+contigLengths[map->contigNo2];
        
        insertMin2=map->readLength1+gap+map->pos2+map->readLength2-1;
        
        insertMax2=contigLengths[map->contigNo1]+gap+map->pos2+map->readLength2-1;
    }
    else
    {
        gap=gapsReversed[contigNo1][contigNo2];
        insertMin1=map->readLength1+map->pos1-1+gap+map->readLength2;
        
        insertMax1=map->readLength1+map->pos1-1+gap+contigLengths[map->contigNo2];
        
        insertMin2=map->readLength1+gap+map->pos2-1+map->readLength2;
        
        insertMax2=contigLengths[map->contigNo1]+gap+map->pos2-1+map->readLength2;
    }
    
    long int insertSize=map->insertSize+gap;
    
    insertMin1=insertSize<insertMin1?insertSize:insertMin1;
    insertMax1=insertSize>insertMax1?insertSize:insertMax1;
    
    insertMin2=insertSize<insertMin2?insertSize:insertMin2;
    insertMax2=insertSize>insertMax2?insertSize:insertMax2;
    
    //    long int totalEffectiveLength=getEffectiveLength(insertSize);
    
    //    long int totalEffectiveLength=contigLengths[contigNo1]+contigLengths[contigNo2]+gap-insertSize;
    
    //    totalEffectiveLength=totalEffectiveLength>1?totalEffectiveLength:1;
    
    long int totalEffectiveLength1=getEffectiveLength(map->readLength1);
    long int totalEffectiveLength2=getEffectiveLength(map->readLength2);
    
    double insertSizeProb=0;
    
    if(insertSize>=0 && insertSize<maxInsertSize)
    {
        insertSizeProb=insertLengthDistSmoothed[insertSize];
    }
    
    if(insertSizeProb==0)
    {
        insertSizeProb=1/(double)uniqueMappedReads;
    }
    
    
    //    double prob=(1/(double)(totalEffectiveLength))*insertSizeProb*map->errorProb;
    
    double prob=((1/(double)(totalEffectiveLength1))*insertSizeProb/getInsertProbsSum(insertMin1, insertMax1)*(1/(double)(totalEffectiveLength2))*insertSizeProb/getInsertProbsSum(insertMin2, insertMax2))*map->errorProb;
    
    if(prob<0 || prob>1)
    {
        return 0;
    }
    return prob;
    
}






void readFromFile()
{
    
    gaps=new AdjListInt[noContigs];
    gapsReversed=new AdjListInt[noContigs];
    
    
    likelihoods=new Likelihoods[noContigs];
    likelihoodsReversed=new Likelihoods[noContigs];
	
	counts=new AdjListDouble[noContigs];
    countsReversed=new AdjListDouble[noContigs];
    
    linkingMaps=new LinkingMaps[noContigs];
    linkingMapsReversed=new LinkingMaps[noContigs];
    
    overlaps=new Overlaps[noContigs];
    overlapsReversed=new Overlaps[noContigs];

   
    contigOrientations=new int[noContigs];
    contigRoots=new long int[noContigs];
    contigSuccessors=new long int[noContigs];
    contigSuccessorsReversed=new long int[noContigs];
    contigStarts=new long int[noContigs];
    contigEnds=new long int[noContigs];

	nc=new double[noContigs];

    for (int i = 0; i < noContigs; i++)
    {
        contigOrientations[i]=0;
        contigSuccessors[i]=-1;
        contigSuccessorsReversed[i]=-1;
        contigRoots[i]=i;
        contigStarts[i]=i;
        contigEnds[i]=i;
	 nc[i]=0;
    }

    
	for (int i = 0; i < noContigs; i++)
    {
//        gaps[i] = new int[noContigs];

	likelihoods[i].init(contigLengths[i],0);
	contigLengthsOriginals.push_back(contigLengths[i]);

//        likelihoods[i]= new double[noContigs];
//        counts[i]= new double[noContigs];
//        gapsReversed[i] = new int[noContigs];

	likelihoodsReversed[i].init(contigLengths[i],0);
//        likelihoodsReversed[i]= new double[noContigs];
//        countsReversed[i]= new double[noContigs];
//        linkingMaps[i] = new vector<MAP *>[noContigs];
//        linkingMapsReversed[i] = new vector<MAP *>[noContigs];
//        overlaps[i]=new vector<int>[noContigs];
//        overlapsReversed[i]=new vector<int>[noContigs];
      
    }
    
	for(int i=0;i<noContigs;i++)
	{

		gaps[i].reset();
		gapsReversed[i].reset();

		counts[i].reset();
		countsReversed[i].reset();

		likelihoods[i].reset();
		likelihoodsReversed[i].reset();

		for(int j=0;j<noContigs;j++)
		{
		//    	gaps[i][j]=0;
            	//	likelihoods[i][j]=0;
            	//	counts[i][j]=0;
            	//	gapsReversed[i][j]=0;
            	//	likelihoodsReversed[i][j]=0;
            	//	countsReversed[i][j]=0;
            
	       //     	linkingMaps[i][j].clear();
            	//	linkingMapsReversed[i][j].clear();
            
	       //     	overlaps[i][j].clear();
            	//	overlapsReversed[i][j].clear();
            
		}
	}
    
	FILE *prefixFile=fopen("prefixes.txt","r");
	
	
	char *line= new char[MAX_REC_LEN];
	char *line1= new char[MAX_REC_LEN];
	
	int MAX_FILE_READ=MAX_REC_LEN/sizeof(line[0]);


	char graphFileName[100];
	
	FILE *graphFile;

	int prefix,gap,contigNo1,contigNo2;
	double likelihood,count;
	char *temp;    

	while(fgets(line, MAX_FILE_READ, contigFile)!=NULL)
	{

		prefix=atoi(line);
		sprintf(graphFileName,"%d_graph.txt",prefix);
		
		graphFile=fopen(graphFileName,"r");
		
//		while(fgets(line1, MAX_FILE_READ, graphFile)!=NULL)
		while(fscanf(graphFile,"%d %d %lf %d %lf\n",&contigNo1,&contigNo2,&likelihood,&gap,&count)==5)
		{

/*			temp=strtok(line1," ");
			contigNo1=atoi(temp);

			temp=strtok(NULL," ");
			contigNo2=atoi(temp);

			temp=strtok(NULL," ");
			likelihood=atof(temp);
			
			temp=strtok(NULL," ");
			gap=atoi(temp);
			
			temp=strtok(NULL," \n");
			count=atof(temp);

*/

			if(gaps[contigNo1][contigNo2]<0 && gap<0)
			{
				if(likelihood>likelihoods[contigNo1][contigNo2])
				{
					gaps[contigNo1].set(contigNo2,gap);
				}
			}
			else if(gap<0)
			{
				gaps[contigNo1].set(contigNo2,gap);
			}
			else if(gaps[contigNo1][contigNo2]>=0 && gap>=0)
			{
				gaps[contigNo1].set(contigNo2,gaps[contigNo1][contigNo2]*likelihoods[contigNo1][contigNo2]/(likelihoods[contigNo1][contigNo2]+likelihood)+gap*likelihood/(likelihoods[contigNo1][contigNo2]+likelihood));
			} 
			likelihoods[contigNo1].set(contigNo2,likelihoods[contigNo1][contigNo2]+likelihood);
			counts[contigNo1].set(contigNo2,counts[contigNo1][contigNo2]+count);

		}

		sprintf(graphFileName,"%d_graph_rev.txt",prefix);
		
		graphFile=fopen(graphFileName,"r");
		
//		while(fgets(line1, MAX_FILE_READ, graphFile)!=NULL)
		while(fscanf(graphFile,"%d %d %lf %d %lf\n",&contigNo1,&contigNo2,&likelihood,&gap,&count)==5)
		{
/*			temp=strtok(line1," ");
			contigNo1=atoi(temp);

			temp=strtok(NULL," ");
			contigNo2=atoi(temp);

			temp=strtok(NULL," ");
			likelihood=atof(temp);
			
			temp=strtok(NULL," ");
			gap=atoi(temp);
			
			temp=strtok(NULL," \n");
			count=atof(temp);
*/
			if(gapsReversed[contigNo1][contigNo2]<0 && gap<0)
			{
				if(likelihood>likelihoodsReversed[contigNo1][contigNo2])
				{
					gapsReversed[contigNo1].set(contigNo2,gap);
				}
			}
			else if(gap<0)
			{
				gapsReversed[contigNo1].set(contigNo2,gap);
			}
			else if(gapsReversed[contigNo1][contigNo2]>=0 && gap>=0)
			{
				gapsReversed[contigNo1].set(contigNo2,gapsReversed[contigNo1][contigNo2]*likelihoodsReversed[contigNo1][contigNo2]/(likelihoodsReversed[contigNo1][contigNo2]+likelihood)+gap*likelihood/(likelihoodsReversed[contigNo1][contigNo2]+likelihood));
			} 
			likelihoodsReversed[contigNo1].set(contigNo2,likelihoodsReversed[contigNo1][contigNo2]+likelihood);
			countsReversed[contigNo1].set(contigNo2,countsReversed[contigNo1][contigNo2]+count);

		}
		sprintf(graphFileName,"%d_counts.txt",prefix);
		graphFile=fopen(graphFileName,"r");
		
		long int contigNo;
		long int contigLength;
		double readCount;

//		while(fgets(line1, MAX_FILE_READ, graphFile)!=NULL)
		while(fscanf(graphFile,"%ld %ld %lf\n",&contigNo,&contigLength,&readCount)==3)
		{
/*			temp=strtok(line1," ");
			contigNo=atol(temp);

			temp=strtok(NULL," ");
			contigLength=atol(temp);

			temp=strtok(NULL," \n");
			readCount=atof(temp);
*/
			nc[contigNo]+=readCount;
		}
	
	}

	double ncMean=0;
	double ncSD=0;
	double sum=0;
	long int hcCount=0;
	for(long int i=0;i<noContigs;i++)
	{
		sum+=nc[i];
	}
	ncMean=sum/noContigs;
	sum=0;
	
	for(long int i=0;i<noContigs;i++)
	{
		if(nc[i]>=ncMean)
		{
			sum+=(nc[i]-ncMean)*(nc[i]-ncMean);
			hcCount++;
		}
	}
	ncSD=sqrt(sum/hcCount);
	for(long int i=0;i<noContigs;i++)
	{
		if(nc[i]>ncMean+2*ncSD  && nc[i]>1.5*ncMean)
		{
			QEntry *qe=new QEntry();
       	       qe->contigNo=i;
                	qe->likelihood=DBL_MAX;
                	quarantine.push_back(qe);
                
		}
	
	}

}

char *merge(char *s, char *t, int gap)
{
    long int sLen=strlen(s);
    long int tLen=strlen(t);
    
    char *dest=new char[sLen+tLen+gap+1];
    
    long int i;
    long int j;
    
    for(i=0;i<sLen;i++)
    {
        dest[i]=s[i];
    }
    if(gap<0)
    {
        for(j=-gap;j<tLen;j++)
        {
            dest[i++]=t[j];
            
        }
        
    }
    else
    {
        for(j=0;j<gap;j++)
        {
            dest[i++]='N';
            
        }
        
        for(j=0;j<tLen;j++)
        {
            dest[i++]=t[j];
            
        }
    
        
    }
    dest[i]='\0';
    return dest;
}




void reverseContig(long int contigNo)
{
    long int contigRoot=getContigRoot(contigNo);
    char *contigReverse=new char[contigLengths[contigRoot]+1];
    reverse(contigReverse,contigs[contigRoot]);
    delete [] contigs[contigRoot];
    contigs[contigRoot]=contigReverse;
    contigOrientations[contigStarts[contigRoot]]^=1;
    if(contigStarts[contigRoot]!=contigEnds[contigRoot])
    {
        contigOrientations[contigEnds[contigRoot]]^=1;
        long int temp=contigStarts[contigRoot];
        contigStarts[contigRoot]=contigEnds[contigRoot];
        contigEnds[contigRoot]=temp;
    }
}

int checkLikelihood(double likelihood1, double likelihood2)
{
    	double ratio,lambda;
	lambda=-log(2.0/5.0)/(5*likelihoodSD);
	ratio=5*exp(-lambda*likelihood1);
	ratio=ratio>2?ratio:2;
	
    if((likelihood1-likelihood2)/likelihoodSD>2.5 || likelihood1/likelihood2>ratio)
        return 1;
    else
        return 0;
}


int checkQuarantine(long int contigRoot1,long int contigRoot2, double likelihood)
{
    for(int i=0;i<quarantine.size();i++)
    {
        if((contigRoot1==quarantine[i]->contigNo|| contigRoot2==quarantine[i]->contigNo) && checkLikelihood(likelihood, quarantine[i]->likelihood)==0)
        {
            return 0;
        }
        
    }
    return 1;
}


int joinContigs(long int contigNo1, long int contigNo2, int joinType, int gap, int isDummy, double likelihoodIn)
{
    
    long int contigRoot1=getContigRoot(contigNo1);
    long int contigRoot2=getContigRoot(contigNo2);
    
    if(contigRoot1==contigRoot2)
        return 0;
    
    
    char *mergedContig;
    long int mergedLength;

    double likelihood;
    
    if(joinType==0)
    {

        if(gap<MIN_GAP)
        {
            gap=gaps[contigNo1][contigNo2];
        }
        likelihood=likelihoods[contigNo1][contigNo2];
        
        if(likelihoodIn>0)
            likelihood=likelihoodIn;
        
        if(checkQuarantine(contigRoot1, contigRoot2, likelihood)==0)
            return 0;

        

        
        //reverse a contig if opposite orientation, swap if both are reversed
        if(contigOrientations[contigNo1]==0 && contigOrientations[contigNo2]==0)
        {
            if(contigEnds[contigRoot1]!=contigNo1 || contigStarts[contigRoot2]!=contigNo2)
                return 0;
            
        }
        else if(contigOrientations[contigNo1]==1 && contigOrientations[contigNo2]==1)
        {
            if(contigStarts[contigRoot1]!=contigNo1 || contigEnds[contigRoot2]!=contigNo2)
                return 0;
            else
            {
                long int temp=contigNo1;
                contigNo1=contigNo2;
                contigNo2=temp;
                
                temp=contigRoot1;
                contigRoot1=contigRoot2;
                contigRoot2=temp;
            }
        }
        else if(contigOrientations[contigNo1]==0 && contigOrientations[contigNo2]==1)
        {
            if(contigEnds[contigRoot1]!=contigNo1 || contigEnds[contigRoot2]!=contigNo2)
                return 0;
            else
            {
                reverseContig(contigNo2);
            }
        }
        else if(contigOrientations[contigNo1]==1 && contigOrientations[contigNo2]==0)
        {
            if(contigStarts[contigRoot1]!=contigNo1 || contigStarts[contigRoot2]!=contigNo2)
                return 0;
            else
            {
                reverseContig(contigNo1);
            }
        }
        cout<<"joining 0 "<<contigNames[contigNo1]<<" "<<contigNames[contigNo2]<<endl;
        
        fprintf(joinsFile, "%d %s %s %lf %ld %ld\n",gap,contigNames[contigNo1],contigNames[contigNo2],likelihoods[contigNo1][contigNo2],contigLengths[contigRoot1],contigLengths[contigRoot2]);

        if(isDummy==1)
        {
            return 1;
        }
        
        mergedContig=merge(contigs[contigRoot1],contigs[contigRoot2],gap);
        mergedLength=strlen(mergedContig);
        
    }
    else if(joinType==1)
    {
        if(gap<MIN_GAP)
        {
            if(contigNo1<contigNo2)
                gap=gapsReversed[contigNo1][contigNo2];
            else
                gap=gapsReversed[contigNo2][contigNo1];
        }
        if(contigNo1<contigNo2)
            likelihood=likelihoodsReversed[contigNo1][contigNo2];
        else
            likelihood=likelihoodsReversed[contigNo2][contigNo1];
        
        if(likelihoodIn>0)
            likelihood=likelihoodIn;

        
        if(checkQuarantine(contigRoot1, contigRoot2, likelihood)==0)
            return 0;

        
        //reverse a contig if same orientation
        
        if(contigOrientations[contigNo1]==0 && contigOrientations[contigNo2]==0)
        {
            
            if(contigEnds[contigRoot1]!=contigNo1 || contigEnds[contigRoot2]!=contigNo2)
                return 0;
            else
            {
                reverseContig(contigNo2);
            }
            
        }
        else if(contigOrientations[contigNo1]==1 && contigOrientations[contigNo2]==1)
        {
            if(contigStarts[contigRoot1]!=contigNo1 || contigStarts[contigRoot2]!=contigNo2)
                return 0;
            else
            {
                reverseContig(contigNo1);
            }
        }
        else if(contigOrientations[contigNo1]==0 && contigOrientations[contigNo2]==1)
        {
            if(contigEnds[contigRoot1]!=contigNo1 || contigStarts[contigRoot2]!=contigNo2)
                return 0;
            
        }
        else if(contigOrientations[contigNo1]==1 && contigOrientations[contigNo2]==0)
        {
            if(contigStarts[contigRoot1]!=contigNo1 || contigEnds[contigRoot2]!=contigNo2)
                return 0;
            else
            {
                long int temp=contigNo1;
                contigNo1=contigNo2;
                contigNo2=temp;
                
                temp=contigRoot1;
                contigRoot1=contigRoot2;
                contigRoot2=temp;
            }
        }
        
        cout<<"joining 1 "<<contigNames[contigNo1]<<" "<<contigNames[contigNo2]<<endl;
        
        fprintf(joinsFile, "%d %s %s %lf %ld %ld\n",gap,contigNames[contigNo1],contigNames[contigNo2],likelihood,contigLengths[contigRoot1],contigLengths[contigRoot2]);
        
        if(isDummy==1)
        {
            return 1;
        }

        
        mergedContig=merge(contigs[contigRoot1],contigs[contigRoot2],gap);
        mergedLength=strlen(mergedContig);
    }
    else if(joinType==2)
    {
        if(gap<MIN_GAP)
        {
            if(contigNo1>contigNo2)
                gap=gapsReversed[contigNo1][contigNo2];
            else
                gap=gapsReversed[contigNo2][contigNo1];
        }
        if(contigNo1>contigNo2)
            likelihood=likelihoodsReversed[contigNo1][contigNo2];
        else
            likelihood=likelihoodsReversed[contigNo2][contigNo1];
        
        if(likelihoodIn>0)
            likelihood=likelihoodIn;

        
        if(checkQuarantine(contigRoot1, contigRoot2, likelihood)==0)
            return 0;

        
        //reverse a contig if same orientation
        
        if(contigOrientations[contigNo1]==0 && contigOrientations[contigNo2]==0)
        {
            
            if(contigStarts[contigRoot1]!=contigNo1 || contigStarts[contigRoot2]!=contigNo2)
                return 0;
            else
            {
                reverseContig(contigNo1);
            }
            
        }
        else if(contigOrientations[contigNo1]==1 && contigOrientations[contigNo2]==1)
        {
            if(contigEnds[contigRoot1]!=contigNo1 || contigEnds[contigRoot2]!=contigNo2)
                return 0;
            else
            {
                reverseContig(contigNo2);
            }
        }
        else if(contigOrientations[contigNo1]==1 && contigOrientations[contigNo2]==0)
        {
            if(contigEnds[contigRoot1]!=contigNo1 || contigStarts[contigRoot2]!=contigNo2)
                return 0;
            
        }
        else if(contigOrientations[contigNo1]==0 && contigOrientations[contigNo2]==1)
        {
            if(contigStarts[contigRoot1]!=contigNo1 || contigEnds[contigRoot2]!=contigNo2)
                return 0;
            else
            {
                long int temp=contigNo1;
                contigNo1=contigNo2;
                contigNo2=temp;
                
                temp=contigRoot1;
                contigRoot1=contigRoot2;
                contigRoot2=temp;
            }
        }
        
        cout<<"joining 2 "<<contigNames[contigNo1]<<" "<<contigNames[contigNo2]<<endl;
        
        fprintf(joinsFile, "%d %s %s %lf %ld %ld\n",gap,contigNames[contigNo1],contigNames[contigNo2],likelihood,contigLengths[contigRoot1],contigLengths[contigRoot2]);
        
        if(isDummy==1)
        {
            return 1;
        }

        
        mergedContig=merge(contigs[contigRoot1],contigs[contigRoot2],gap);
        mergedLength=strlen(mergedContig);
    }


    delete [] contigs[contigRoot1];
    delete [] contigs[contigRoot2];
    
    
    if(contigRoot1<contigRoot2)
    {
        contigRoots[contigRoot2]=contigRoot1;
        contigEnds[contigRoot1]=contigEnds[contigRoot2];
        contigs[contigRoot1]=mergedContig;
        contigLengths[contigRoot1]=mergedLength;
        contigLengths[contigRoot2]=-1;
    }
    else
    {
        contigRoots[contigRoot1]=contigRoot2;
        contigStarts[contigRoot2]=contigStarts[contigRoot1];
        contigs[contigRoot2]=mergedContig;
        contigLengths[contigRoot2]=mergedLength;
        contigLengths[contigRoot1]=-1;
    }
    return 1;
}

void findUniqueJoins()
{
    
    
      
    
    int *indegrees=new int[noContigs];
    int *outdegrees=new int[noContigs];
    int *indegreesReversed=new int[noContigs];
    int *outdegreesReversed=new int[noContigs];
  

	for(int i=0;i<noContigs;i++)
    	{
  		indegrees[i]=0;
		outdegrees[i]=0;
  		indegreesReversed[i]=0;
		outdegreesReversed[i]=0;		
  	}
    

	for(int i=0;i<noContigs;i++)
    	{
		for(int index=0;index<likelihoods[i].size;index++)
		{
			int j=likelihoods[i].indexes[index];

			if(i==j)
                		continue;

			if(likelihoods[i][j]>0)
            		{
                		outdegrees[i]++;
				indegrees[j]++;
                		contigSuccessors[i]=j;
            		}
			
		}


		for(int index=0;index<likelihoodsReversed[i].size;index++)
		{
			int j=likelihoodsReversed[i].indexes[index];

			if(likelihoodsReversed[i][j]>0)
            		{

				if(i<j)
				{
					outdegreesReversed[i]++;
					outdegreesReversed[j]++;
					contigSuccessors[i]=j;
				//	contigSuccessorsReversed[j]=i;
		
				}
				else if(i>j)
				{
					indegreesReversed[i]++;
					indegreesReversed[j]++;
					contigSuccessorsReversed[i]=j;
				//	contigSuccessors[j]=i;

				}
            		}
			
		}


	}

/*     
    int indegree, outdegree, indegreeReversed, outdegreeReversed;
    
    for(int i=0; i < noContigs; i++)
    {
        indegree=0;
        outdegree=0;
        indegreeReversed=0;
        outdegreeReversed=0;
        
        for(int j=0;j<noContigs;j++)
        {
            if(i==j)
                continue;
            
            if(likelihoods[i][j]>0)
            {
                outdegree++;
                contigSuccessors[i]=j;
            }
            
            if(likelihoods[j][i]>0)
                indegree++;
            
        }
        for(int j=0;j<i;j++)
        {
            if(likelihoodsReversed[j][i]>0)
            {
                outdegreeReversed++;
        //        contigSuccessors[i]=j;
            }
        }
        for(int j=i+1;j<noContigs;j++)
        {
            if(likelihoodsReversed[i][j]>0)
            {
                outdegreeReversed++;
                contigSuccessors[i]=j;
            }
        }
        for(int j=0;j<i;j++)
        {
            if(likelihoodsReversed[i][j]>0)
            {
                indegreeReversed++;
                contigSuccessorsReversed[i]=j;
            }
        }
        for(int j=i+1;j<noContigs;j++)
        {
            if(likelihoodsReversed[j][i]>0)
            {
                indegreeReversed++;
       //         contigSuccessors[i]=j;
            }
        }
        
        outdegrees[i]=outdegree;
        indegrees[i]=indegree;
        outdegreesReversed[i]=outdegreeReversed;
        indegreesReversed[i]=indegreeReversed;
    
    }
*/ 
    int isDummy=0;
    
    for(int i=0;i<noContigs;i++)
    {
        if(outdegrees[i]==1 && outdegreesReversed[i]==0 && contigSuccessors[i]!=-1 && indegrees[contigSuccessors[i]]==1 && indegreesReversed[contigSuccessors[i]]==0)
        {
            int joined;
		if(counts[i][contigSuccessors[i]]>MIN_READ_JOIN)
		{	
			joined=joinContigs(i,contigSuccessors[i],0,MIN_GAP-1,isDummy,0);
            		if(joined==1)
                		likelihoodsDist.push_back(likelihoods[i][contigSuccessors[i]]);
            	}
      //      cout<<likelihoods[i][contigSuccessors[i]]<<endl;
                
        }
        else if(outdegreesReversed[i]==1 && outdegrees[i]==0 && contigSuccessors[i]!=-1 && outdegreesReversed[contigSuccessors[i]]==1 && outdegrees[contigSuccessors[i]]==0)
        {
            int joined;
		if((i<contigSuccessors[i] && countsReversed[i][contigSuccessors[i]]>MIN_READ_JOIN) || (i>contigSuccessors[i] && countsReversed[contigSuccessors[i]][i]>MIN_READ_JOIN)) 
		{
			joined=joinContigs(i,contigSuccessors[i],1,MIN_GAP-1,isDummy,0);
            		if(joined==1)
            		{
                		if(i<contigSuccessors[i])
                    			likelihoodsDist.push_back(likelihoodsReversed[i][contigSuccessors[i]]);
                		else	
                    			likelihoodsDist.push_back(likelihoodsReversed[contigSuccessors[i]][i]);
            		}
		}
      //      cout<<likelihoodsReversed[i][contigSuccessors[i]]<<endl;

        }
        else
            contigSuccessors[i]=-1;
        
        if(indegreesReversed[i]==1 && indegrees[i]==0 && contigSuccessorsReversed[i]!=-1 && indegrees[contigSuccessorsReversed[i]]==0 && indegreesReversed[contigSuccessorsReversed[i]]==1)
        {
            int joined;

	 	if((i>contigSuccessorsReversed[i] && countsReversed[i][contigSuccessorsReversed[i]]>MIN_READ_JOIN) || (i<contigSuccessorsReversed[i] && countsReversed[contigSuccessorsReversed[i]][i]>MIN_READ_JOIN))
		{
			joined=joinContigs(i,contigSuccessorsReversed[i],2,MIN_GAP-1,isDummy,0);
            		if(joined==1)
            		{
                		if(i>contigSuccessorsReversed[i])
                    			likelihoodsDist.push_back(likelihoodsReversed[i][contigSuccessorsReversed[i]]);
                		else
                    			likelihoodsDist.push_back(likelihoodsReversed[contigSuccessorsReversed[i]][i]);
            		}
		}
      //      cout<<likelihoodsReversed[i][contigSuccessors[i]]<<endl;
        }
        else
            contigSuccessorsReversed[i]=-1;
    
    }
/*  
    for(int i=0;i<noContigs;i++)
    {
        
        cout<<i<<" "<<contigSuccessors[i]<<" "<<outdegrees[i]<<" "<<indegrees[i]<<" "<<outdegreesReversed[i]<<" "<<indegreesReversed[i]<<" "<<contigOrientations[i]<<endl;
    }
*/  
   
}

void computeLikelihoodDist()
{
    long int n=likelihoodsDist.size();
    likelihoodMean=0;
    for(int i=0;i<n;i++)
    {
        
        likelihoodMean+=likelihoodsDist[i];
    }
    likelihoodMean/=n;
    likelihoodSD=0;
    for(int i=0;i<n;i++)
    {
        likelihoodSD+=(likelihoodsDist[i]-likelihoodMean)*(likelihoodsDist[i]-likelihoodMean);
    }
    likelihoodSD=sqrt(likelihoodSD/n);
    
}


void insertJoin(Join * newJoin)
{
    for(int i=0;i<joins.size();i++)
    {
        if(joins[i]->contigNo1==newJoin->contigNo1 && joins[i]->contigNo2==newJoin->contigNo2 && joins[i]->joinType==newJoin->joinType)
        {
//            cout<<"deleting join "<<contigNames[joins[i]->contigNo1]<<" "<<contigNames[joins[i]->contigNo2]<<" "<<joins[i]->joinType<<endl;
            delete joins[i];
            joins.erase(joins.begin()+i);
            break;
        }
    }
    
    for(int i=0;i<joins.size();i++)
    {
        if(newJoin->likelihood>=joins[i]->likelihood)
        {
            joins.insert(joins.begin()+i,newJoin);
            return;
            
        }
    }
    joins.push_back(newJoin);
}

int checkFailedJoins(Join * j)
{
    Join* fj;
    
    for(int i=0;i<failedJoins.size();i++)
    {
        fj=failedJoins[i];
        
        if(j->likelihood>fj->likelihood)
            continue;
        
        if(j->joinType==0)
        {
            if((fj->joinType==0 && (fj->contigNo1==j->contigNo1 || fj->contigNo2==j->contigNo2)) || (fj->joinType==1 && (fj->contigNo1==j->contigNo1 || fj->contigNo2==j->contigNo1)) || (fj->joinType==2 && (fj->contigNo1==j->contigNo2 || fj->contigNo2==j->contigNo2)))
            {
                return 0;
                
            }
            
        }
        else if(j->joinType==1)
        {
            if((fj->joinType==0 && (fj->contigNo1==j->contigNo1 || fj->contigNo1==j->contigNo2)) || (fj->joinType==1 && (fj->contigNo1==j->contigNo1 || fj->contigNo2==j->contigNo1 || fj->contigNo1==j->contigNo2 || fj->contigNo2==j->contigNo2)))
            {
                return 0;
                
            }

            
        }
        else if(j->joinType==2)
        {
            if((fj->joinType==0 && (fj->contigNo2==j->contigNo1 || fj->contigNo2==j->contigNo2)) || (fj->joinType==2 && (fj->contigNo1==j->contigNo1 || fj->contigNo2==j->contigNo1 || fj->contigNo1==j->contigNo2 || fj->contigNo2==j->contigNo2)))
            {
                return 0;
                
            }
            
        }

        
    }
    return 1;
}

int isMember(vector<long int> set, long int element)
{
    for(int j=0;j<set.size();j++)
    {
        if(set[j]==element)
        {
            return 1;
        }
        
    }
    return 0;
    
}

void insertInterval(Interval * newInterval)
{
    for(int i=0;i<intervals.size();i++)
    {
        if(newInterval->end<=intervals[i]->end)
        {
            intervals.insert(intervals.begin()+i,newInterval);
            return;
            
        }
    }
    intervals.push_back(newInterval);
}

void selectIntraContigs()
{
    
    //vals indexed from 1 and intervals from 0
    
    double *vals=new double[intervals.size()+1];
    int pre;
    
    vals[0]=0;
    
    for(int i=0;i<intervals.size();i++)
    {
        pre=0;
        for(int j=i-1;j>=0;j--)
        {
            if(intervals[j]->end-getOverlap(intervals[j]->contigNo, intervals[i]->contigNo, intervals[j]->contigOrientation, intervals[i]->contigOrientation)<intervals[i]->gap*stretchFactor)
            {
                pre=j+1;
                break;
            }
        }
        
        if(vals[pre]+intervals[i]->likelihood>=vals[i])
        {
            vals[i+1]=vals[pre]+intervals[i]->likelihood;
        }
        else
        {
            vals[i+1]=vals[i];
        }
    }
    
    for(int i=intervals.size()-1;i>=0;)
    {
        if(vals[i+1]==vals[i])
        {
            intervals[i]->selected=0;
            i--;
        }
        else
        {
            intervals[i]->selected=1;
            pre=0;
            for(int j=i-1;j>=0;j--)
            {
                if(intervals[j]->end-getOverlap(intervals[j]->contigNo, intervals[i]->contigNo, intervals[j]->contigOrientation, intervals[i]->contigOrientation)<intervals[i]->gap*stretchFactor)
                {
                    pre=j+1;
                    break;
                }
            }
            i=pre-1;
        }
    }
    
    
    for(int i=0;i<intervals.size();i++)
    {
        for(int j=i-1;j>=0;j--)
        {
            if(getContigRoot(intervals[i]->contigNo)==getContigRoot(intervals[j]->contigNo))
            {
                if(intervals[i]->selected==1 && intervals[j]->selected==1)
                {
                    int toDelete=intervals[i]->likelihood<intervals[j]->likelihood?i:j;
                    intervals[toDelete]->selected=0;
                }
                
            }
        }
    }
    
    for(int i=0;i<intervals.size();i++)
    {
        if(intervals[i]->selected!=1)
            continue;
        
        for(int j=0;j<intervals.size();j++)
        {
            if(i==j)
                continue;
            
            if(getContigRoot(intervals[i]->contigNo)==getContigRoot(intervals[j]->contigNo))
            {
                intervals[j]->selected=-1;
                
            }
        }
    }

    for(int i=0;i<intervals.size();i++)
    {
        if(intervals[i]->selected!=1)
            continue;

	 if(intervals[i]->count<=MIN_READ_JOIN)
	 {
                intervals[i]->selected=0;	
	 }
    }
    
    for(int i=0;i<intervals.size();i++)
    {
        if(intervals[i]->selected!=1)
            continue;
        
        for(int j=0;j<intervals.size();j++)
        {
            if(intervals[j]->selected!=0)
                continue;
            
            if((intervals[j]->end-getOverlap(intervals[j]->contigNo, intervals[i]->contigNo, intervals[j]->contigOrientation, intervals[i]->contigOrientation)>intervals[i]->gap*stretchFactor && intervals[j]->gap*stretchFactor<intervals[i]->end-getOverlap(intervals[i]->contigNo, intervals[j]->contigNo, intervals[i]->contigOrientation, intervals[j]->contigOrientation)) && checkLikelihood(intervals[i]->likelihood, intervals[j]->likelihood)==0)
            {
                
//                cout<<"**** conflict **** "<<contigNames[intervals[i]->contigNo]<<" "<<contigNames[intervals[j]->contigNo]<<endl;
                intervals[i]->selected=0;
                break;
            }
        }
    }
    
    delete [] vals;
    
    
}

void addExtraJoins()
{
    
    for(int i=0;i<extraJoins1.size();i++)
    {
        Join *newJoin=new Join();
        Join *extraJoin=extraJoins1[i];
        int tempGap;
        // get join type
        
        newJoin->joinType=extraJoin->joinType;
        newJoin->contigNo1=extraJoin->contigNo1;
        newJoin->contigNo2=extraJoin->contigNo2;
        newJoin->isExtra=1;
        
        // check if edge already exists
        // if not set likelihood and gap
        // otherwise add likelihood
        
//        if(extraJoin->gap<-100)
//            continue;
        
        if(extraJoin->gap<10)
        {
            tempGap=-getOverlap(extraJoin->contigNo1, extraJoin->contigNo2, extraJoin->joinType);
            if(tempGap<0)
            {
                extraJoin->gap=tempGap;
            }
            else
            {
                extraJoin->gap=10;
            }
        }
        
        if(newJoin->joinType==0)
        {
            if(likelihoods[newJoin->contigNo1][newJoin->contigNo2]>0)
            {
                newJoin->likelihood=likelihoods[newJoin->contigNo1][newJoin->contigNo2]+extraJoin->likelihood;
			newJoin->count=counts[newJoin->contigNo1][newJoin->contigNo2]+extraJoin->count;
                if(gaps[newJoin->contigNo1][newJoin->contigNo2]>0)
                {
                    gaps[newJoin->contigNo1].set(newJoin->contigNo2,gaps[newJoin->contigNo1][newJoin->contigNo2]*likelihoods[newJoin->contigNo1][newJoin->contigNo2]/newJoin->likelihood + extraJoin->gap*extraJoin->likelihood/newJoin->likelihood);
                    
                }
                
                likelihoods[newJoin->contigNo1].set(newJoin->contigNo2,newJoin->likelihood);
	         counts[newJoin->contigNo1].set(newJoin->contigNo2,newJoin->count);
                
            }
            else
            {
                newJoin->likelihood=extraJoin->likelihood;
                likelihoods[newJoin->contigNo1].set(newJoin->contigNo2,newJoin->likelihood);
                newJoin->count=extraJoin->count;
                counts[newJoin->contigNo1].set(newJoin->contigNo2,newJoin->count);
		
                //  set gap
                gaps[newJoin->contigNo1].set(newJoin->contigNo2,extraJoin->gap);
            }
        }
        else if(newJoin->joinType==1)
        {
            if(likelihoodsReversed[newJoin->contigNo1][newJoin->contigNo2]>0)
            {
                newJoin->likelihood=likelihoodsReversed[newJoin->contigNo1][newJoin->contigNo2]+extraJoin->likelihood;
                newJoin->count=countsReversed[newJoin->contigNo1][newJoin->contigNo2]+extraJoin->count;
                if(gapsReversed[newJoin->contigNo1][newJoin->contigNo2]>0)
                {
                    gapsReversed[newJoin->contigNo1].set(newJoin->contigNo2,gapsReversed[newJoin->contigNo1][newJoin->contigNo2]*likelihoodsReversed[newJoin->contigNo1][newJoin->contigNo2]/newJoin->likelihood + extraJoin->gap*extraJoin->likelihood/newJoin->likelihood);
                    
                }
                likelihoodsReversed[newJoin->contigNo1].set(newJoin->contigNo2,newJoin->likelihood);
                countsReversed[newJoin->contigNo1].set(newJoin->contigNo2,newJoin->count);

            }
            else
            {
                newJoin->likelihood=extraJoin->likelihood;
                likelihoodsReversed[newJoin->contigNo1].set(newJoin->contigNo2,newJoin->likelihood);
                newJoin->count=extraJoin->count;
                countsReversed[newJoin->contigNo1].set(newJoin->contigNo2,newJoin->count);

                //  set gap
                gapsReversed[newJoin->contigNo1].set(newJoin->contigNo2,extraJoin->gap);
            }
        }
        
        else if(newJoin->joinType==2)
        {
            if(likelihoodsReversed[newJoin->contigNo1][newJoin->contigNo2]>0)
            {
                newJoin->likelihood=likelihoodsReversed[newJoin->contigNo1][newJoin->contigNo2]+extraJoin->likelihood;
                newJoin->count=countsReversed[newJoin->contigNo1][newJoin->contigNo2]+extraJoin->count;
                if(gapsReversed[newJoin->contigNo1][newJoin->contigNo2]>0)
                {
                    gapsReversed[newJoin->contigNo1].set(newJoin->contigNo2,gapsReversed[newJoin->contigNo1][newJoin->contigNo2]*likelihoodsReversed[newJoin->contigNo1][newJoin->contigNo2]/newJoin->likelihood + extraJoin->gap*extraJoin->likelihood/newJoin->likelihood);
                    
                }
                likelihoodsReversed[newJoin->contigNo1].set(newJoin->contigNo2,newJoin->likelihood);
                countsReversed[newJoin->contigNo1].set(newJoin->contigNo2,newJoin->count);
            }
            else
            {
                newJoin->likelihood=extraJoin->likelihood;
                likelihoodsReversed[newJoin->contigNo1].set(newJoin->contigNo2,newJoin->likelihood);
                newJoin->count=extraJoin->count;
                countsReversed[newJoin->contigNo1].set(newJoin->contigNo2,newJoin->count);

                //  set gap
                gapsReversed[newJoin->contigNo1].set(newJoin->contigNo2,extraJoin->gap);
            }
        }
        // insert into joins
        
        
//        cout<<"adding to joins "<<contigNames[newJoin->contigNo1]<<" "<<contigNames[newJoin->contigNo2]<<" "<<newJoin->joinType<<" "<<newJoin->likelihood<<" "<<newJoin->gap<<endl;
	if(newJoin->count > MIN_READ_JOIN)        
 		insertJoin(newJoin);
        
        
    }
    
    
    for(int i=0;i<extraJoins2.size();i++)
    {
        Join *newJoin=new Join();
        Join *extraJoin=extraJoins2[i];
        int tempGap;
        
        // get join type
        
        newJoin->joinType=extraJoin->joinType;
        newJoin->contigNo1=extraJoin->contigNo1;
        newJoin->contigNo2=extraJoin->contigNo2;
        newJoin->isExtra=1;
        
    //    if(extraJoin->gap<-100)
    //        continue;
        
        if(extraJoin->gap<10)
        {
            tempGap=-getOverlap(extraJoin->contigNo1, extraJoin->contigNo2, extraJoin->joinType);
            if(tempGap<0)
            {
                extraJoin->gap=tempGap;
            }
            else
            {
                extraJoin->gap=10;
            }
        }
        
        // check if edge already exists
        // if not set likelihood and gap
        // otherwise add likelihood
        if(newJoin->joinType==0)
        {
            if(likelihoods[newJoin->contigNo1][newJoin->contigNo2]>0)
            {
                newJoin->likelihood=likelihoods[newJoin->contigNo1][newJoin->contigNo2]+extraJoin->likelihood;
                newJoin->count=counts[newJoin->contigNo1][newJoin->contigNo2]+extraJoin->count;
                if(gaps[newJoin->contigNo1][newJoin->contigNo2]>0)
                {
                    gaps[newJoin->contigNo1].set(newJoin->contigNo2,gaps[newJoin->contigNo1][newJoin->contigNo2]*likelihoods[newJoin->contigNo1][newJoin->contigNo2]/newJoin->likelihood + extraJoin->gap*extraJoin->likelihood/newJoin->likelihood);
                    
                }
                likelihoods[newJoin->contigNo1].set(newJoin->contigNo2,newJoin->likelihood);
                counts[newJoin->contigNo1].set(newJoin->contigNo2,newJoin->count);

            }
            else
            {
                newJoin->likelihood=extraJoin->likelihood;
                likelihoods[newJoin->contigNo1].set(newJoin->contigNo2,newJoin->likelihood);
                newJoin->count=extraJoin->count;
                counts[newJoin->contigNo1].set(newJoin->contigNo2,newJoin->count);

                //  set gap
                gaps[newJoin->contigNo1].set(newJoin->contigNo2,extraJoin->gap);
            }
        }
        else if(newJoin->joinType==2)
        {
            if(likelihoodsReversed[newJoin->contigNo1][newJoin->contigNo2]>0)
            {
                newJoin->likelihood=likelihoodsReversed[newJoin->contigNo1][newJoin->contigNo2]+extraJoin->likelihood;
                newJoin->count=countsReversed[newJoin->contigNo1][newJoin->contigNo2]+extraJoin->count;
                if(gapsReversed[newJoin->contigNo1][newJoin->contigNo2]>0)
                {
                    gapsReversed[newJoin->contigNo1].set(newJoin->contigNo2,gapsReversed[newJoin->contigNo1][newJoin->contigNo2]*likelihoodsReversed[newJoin->contigNo1][newJoin->contigNo2]/newJoin->likelihood + extraJoin->gap*extraJoin->likelihood/newJoin->likelihood);
                    
                }
                likelihoodsReversed[newJoin->contigNo1].set(newJoin->contigNo2,newJoin->likelihood);
                countsReversed[newJoin->contigNo1].set(newJoin->contigNo2,newJoin->count);

            }
            else
            {
                newJoin->likelihood=extraJoin->likelihood;
                likelihoodsReversed[newJoin->contigNo1].set(newJoin->contigNo2,newJoin->likelihood);
                newJoin->count=extraJoin->count;
                countsReversed[newJoin->contigNo1].set(newJoin->contigNo2,newJoin->count);

                //  set gap
                gapsReversed[newJoin->contigNo1].set(newJoin->contigNo2,extraJoin->gap);
            }
        }
        else if(newJoin->joinType==1)
        {
            if(likelihoodsReversed[newJoin->contigNo1][newJoin->contigNo2]>0)
            {
                newJoin->likelihood=likelihoodsReversed[newJoin->contigNo1][newJoin->contigNo2]+extraJoin->likelihood;
                newJoin->count=countsReversed[newJoin->contigNo1][newJoin->contigNo2]+extraJoin->count;
                if(gapsReversed[newJoin->contigNo1][newJoin->contigNo2]>0)
                {
                    gapsReversed[newJoin->contigNo1].set(newJoin->contigNo2,gapsReversed[newJoin->contigNo1][newJoin->contigNo2]*likelihoodsReversed[newJoin->contigNo1][newJoin->contigNo2]/newJoin->likelihood + extraJoin->gap*extraJoin->likelihood/newJoin->likelihood);
                    
                }
                likelihoodsReversed[newJoin->contigNo1].set(newJoin->contigNo2,newJoin->likelihood);
                countsReversed[newJoin->contigNo1].set(newJoin->contigNo2,newJoin->count);
            }
            else
            {
                newJoin->likelihood=extraJoin->likelihood;
                likelihoodsReversed[newJoin->contigNo1].set(newJoin->contigNo2,newJoin->likelihood);
                newJoin->count=extraJoin->count;
                countsReversed[newJoin->contigNo1].set(newJoin->contigNo2,newJoin->count);

                //  set gap
                gapsReversed[newJoin->contigNo1].set(newJoin->contigNo2,extraJoin->gap);
            }
        }
//        cout<<"adding to joins "<<contigNames[newJoin->contigNo1]<<" "<<contigNames[newJoin->contigNo2]<<" "<<newJoin->joinType<<" "<<newJoin->likelihood<<" "<<newJoin->gap<<endl;
	if(newJoin->count > MIN_READ_JOIN) 
       	insertJoin(newJoin);
        
    }
    
    
}

void joinPath(long int contigNo1, long int contigNo2, int joinType, int gap, double likelihoodIn)
{
    
    Join *j;
    
//    cout<<contigNames[contigNo1]<<endl;
//    cout<<contigNames[contigNo2]<<endl;

    
    int isDummy=0;
    
    if(intraJoins1.size()>0 || intraJoins2.size()>0)
    {
        isDummy=1;
    }
    
    int joined=joinContigs(contigNo1, contigNo2, joinType,MIN_GAP-1,isDummy,likelihoodIn);

    if(joined==1)
    {
        
        addExtraJoins();
    }
    
    while(extraJoins1.size()!=0)
    {
        j=extraJoins1.front();
//        cout<<contigNames[j->contigNo1]<<" "<<contigNames[j->contigNo2]<<" "<<j->joinType<<" "<<j->gap<<" "<<j->likelihood<<endl;
        delete j;
        extraJoins1.erase(extraJoins1.begin());
    }
//    cout<<endl;
    while(extraJoins2.size()!=0)
    {
        j=extraJoins2.front();
//        cout<<contigNames[j->contigNo1]<<" "<<contigNames[j->contigNo2]<<" "<<j->joinType<<" "<<j->gap<<" "<<j->likelihood<<endl;
        delete j;
        extraJoins2.erase(extraJoins2.begin());
    }
//    cout<<endl;
    
    if(joined==0)
    {
        while(intraJoins1.size()!=0)
        {
            j=intraJoins1.front();
            delete j;
            intraJoins1.erase(intraJoins1.begin());
        }

        while(intraJoins2.size()!=0)
        {
            j=intraJoins2.front();
            delete j;
            intraJoins2.erase(intraJoins2.begin());
        }
        return;
    }
    
    
    long int contigRoot1, contigRoot2;
    long int contigLength1, contigLength2;
    
    contigRoot1=getContigRoot(contigNo1);
    contigRoot2=getContigRoot(contigNo2);

    contigLength1=contigLengths[contigRoot1];
    contigLength2=contigLengths[contigRoot2];
    
    vector<long int> intraContigs;

    long int contigIn;
    long int contigOut;
    
//    intraContigs.push_back(contigRoot1);
//    intraContigs.push_back(contigRoot2);

    long int tempContigNo, tempContigRoot, tempContigLength,tempContigStart, tempContigEnd;
    
    int orderReversed=0;
    
//    cout<<intraJoins1.size()<<" "<<intraJoins2.size()<<endl;
    
    for(int i=0;i<intraJoins1.size();i++)
    {
        j=intraJoins1[i];
        
        if(j->contigNo1==contigNo1)
        {
            tempContigNo=j->contigNo2;
        }
        else
        {
            tempContigNo=j->contigNo1;
        }
        tempContigRoot=getContigRoot(tempContigNo);
        
        if(isMember(intraContigs, tempContigRoot)==0)
        {
            intraContigs.push_back(tempContigRoot);
        }
    }
    for(int i=0;i<intraJoins2.size();i++)
    {
        j=intraJoins2[i];
        
        if(j->contigNo1==contigNo2)
        {
            tempContigNo=j->contigNo2;
        }
        else
        {
            tempContigNo=j->contigNo1;
        }
        tempContigRoot=getContigRoot(tempContigNo);
        
        if(isMember(intraContigs, tempContigRoot)==0)
        {
            intraContigs.push_back(tempContigRoot);
        }
    }
    
    long int gapIn, gapOut;
    
    for(int i=0;i<intraJoins1.size();i++)
    {
        j=intraJoins1[i];
        
        if(j->contigNo1==contigNo1)
        {
            tempContigNo=j->contigNo2;
        }
        else
        {
            tempContigNo=j->contigNo1;
        }
        tempContigRoot=getContigRoot(tempContigNo);
        tempContigLength=contigLengths[tempContigRoot];
        tempContigStart=contigStarts[tempContigRoot];
        tempContigEnd=contigEnds[tempContigRoot];
        
        
        if(j->joinType==0)
        {
            if(j->contigNo1==contigNo1)
            {
                gapIn=j->gap+contigLength1;
                gapOut=gap-j->gap-tempContigLength+contigLength2;
                
                contigIn=contigRoot1;
                contigOut=contigRoot2;
            }
            else
            {
                gapOut=j->gap+contigLength1;
                gapIn=gap-j->gap-tempContigLength+contigLength2;
                
                contigIn=contigRoot2;
                contigOut=contigRoot1;
            }
        }
        else if(j->joinType==1)
        {
            gapOut=j->gap+contigLength1;
            gapIn=gap-j->gap-tempContigLength+contigLength2;
            
            contigIn=contigRoot2;
            contigOut=contigRoot1;
        }
        else if(j->joinType==2)
        {
            gapIn=j->gap+contigLength1;
            gapOut=gap-j->gap-tempContigLength+contigLength2;
            
            contigIn=contigRoot1;
            contigOut=contigRoot2;

        }
        
//        cout<<contigNames[tempContigNo]<<" "<<tempContigLength<<" "<<gapIn<<" "<<gapOut<<endl;
        
        //in degree
        
        for(long int c=0;c<noContigs;c++)
        {
            if(likelihoods[c][tempContigStart]<=0 || checkLikelihood(j->likelihood,likelihoods[c][tempContigStart])==1)
                continue;
            if(isMember(intraContigs, getContigRoot(c))==1 || getContigRoot(c)==contigIn)
                continue;
            if(gaps[c][tempContigStart]*stretchFactor<gapIn)
            {
                j->joinType=-1;
                break;
            }
        }
        
        //out degree
        
        for(long int c=0;c<noContigs;c++)
        {
            if(j->joinType==-1)
                break;
            if(likelihoods[tempContigEnd][c]<=0 || checkLikelihood(j->likelihood,likelihoods[tempContigEnd][c])==1)
                continue;
            if(isMember(intraContigs, getContigRoot(c))==1 || getContigRoot(c)==contigOut)
                continue;
            if(gaps[tempContigEnd][c]*stretchFactor<gapOut)
            {
                j->joinType=-1;
                break;
            }
        }
        
        //out reversed
        
        for(long int c=0;c<tempContigEnd;c++)
        {
            if(j->joinType==-1)
                break;
            if(likelihoodsReversed[c][tempContigEnd]<=0 || checkLikelihood(j->likelihood,likelihoodsReversed[c][tempContigEnd])==1)
                continue;
            if(isMember(intraContigs, getContigRoot(c))==1  || getContigRoot(c)==contigOut)
                continue;
            if(gapsReversed[c][tempContigEnd]*stretchFactor<gapOut)
            {
                j->joinType=-1;
                break;
            }
        }
        for(long int c=tempContigEnd+1;c<noContigs;c++)
        {
            if(j->joinType==-1)
                break;
            if(likelihoodsReversed[tempContigEnd][c]<=0 || checkLikelihood(j->likelihood,likelihoodsReversed[tempContigEnd][c])==1)
                continue;
            if(isMember(intraContigs, getContigRoot(c))==1  || getContigRoot(c)==contigOut)
                continue;
            if(gapsReversed[tempContigEnd][c]*stretchFactor<gapOut)
            {
                j->joinType=-1;
                break;
            }
        }
        
        //in reversed
        
        for(long int c=0;c<tempContigStart;c++)
        {
            if(j->joinType==-1)
                break;
            if(likelihoodsReversed[tempContigStart][c]<=0 || checkLikelihood(j->likelihood,likelihoodsReversed[tempContigStart][c])==1)
                continue;
            if(isMember(intraContigs, getContigRoot(c))==1  || getContigRoot(c)==contigIn)
                continue;
            if(gapsReversed[tempContigStart][c]*stretchFactor<gapIn)
            {
                j->joinType=-1;
                break;
            }
        }
        for(long int c=tempContigStart+1;c<noContigs;c++)
        {
            if(j->joinType==-1)
                break;
            if(likelihoodsReversed[c][tempContigStart]<=0 || checkLikelihood(j->likelihood,likelihoodsReversed[c][tempContigStart])==1)
                continue;
            if(isMember(intraContigs, getContigRoot(c))==1 || getContigRoot(c)==contigIn)
                continue;
            if(gapsReversed[c][tempContigStart]*stretchFactor<gapIn)
            {
                j->joinType=-1;
                break;
            }
        }

        if(j->joinType!=-1)
        {
            
            Interval *newInterval= new Interval();
            newInterval->contigNo=tempContigNo;
            newInterval->gap=j->gap;
            newInterval->end=j->gap+tempContigLength;
            newInterval->selected=0;
            newInterval->likelihood=j->likelihood;
            newInterval->count=j->count;
            newInterval->realGap=j->gap;

            
            if((joinType==0 && j->joinType==0) || (joinType==1 && j->joinType==0)  || (joinType==2 && j->joinType==2))
            {
                newInterval->contigOrientation=0;
            }
            else if((joinType==0 && j->joinType==1) || (joinType==1 && j->joinType==1) || (joinType==2 && j->joinType==0))
            {
                newInterval->contigOrientation=1;
            }
            
            insertInterval(newInterval);
        }
        else
        {
            QEntry *qe=new QEntry();
            qe->contigNo=tempContigRoot;
            qe->likelihood=j->likelihood;
            quarantine.push_back(qe);
            
        }
    }
    
    for(int i=0;i<intraJoins2.size();i++)
    {
        j=intraJoins2[i];
        
        if(j->contigNo1==contigNo2)
        {
            tempContigNo=j->contigNo2;
        }
        else
        {
            tempContigNo=j->contigNo1;
        }
        tempContigRoot=getContigRoot(tempContigNo);
        tempContigLength=contigLengths[tempContigRoot];
        tempContigStart=contigStarts[tempContigRoot];
        tempContigEnd=contigEnds[tempContigRoot];
       
        if(j->joinType==0)
        {
            if(j->contigNo1==contigNo2)
            {
                gapIn=j->gap+contigLength2;
                gapOut=gap-j->gap-tempContigLength+contigLength1;
                
                contigIn=contigRoot2;
                contigOut=contigRoot1;
            }
            else
            {
                gapOut=j->gap+contigLength2;
                gapIn=gap-j->gap-tempContigLength+contigLength1;
                
                contigIn=contigRoot1;
                contigOut=contigRoot2;
                
            }
        }
        else if(j->joinType==1)
        {
            gapOut=j->gap+contigLength2;
            gapIn=gap-j->gap-tempContigLength+contigLength1;
            
            contigIn=contigRoot1;
            contigOut=contigRoot2;
        }
        else if(j->joinType==2)
        {
            gapIn=j->gap+contigLength2;
            gapOut=gap-j->gap-tempContigLength+contigLength1;
            
            contigIn=contigRoot2;
            contigOut=contigRoot1;
            
        }
        
//        cout<<contigNames[tempContigNo]<<" "<<tempContigLength<<" "<<gapIn<<" "<<gapOut<<endl;
        
        //in
        
        for(long int c=0;c<noContigs;c++)
        {
            if(likelihoods[c][tempContigStart]<=0 || checkLikelihood(j->likelihood,likelihoods[c][tempContigStart])==1)
                continue;
            if(isMember(intraContigs, getContigRoot(c))==1 || getContigRoot(c)==contigIn)
                continue;
            if(gaps[c][tempContigStart]*stretchFactor<gapIn)
            {
                j->joinType=-1;
                break;
            }
        }
        
        //out
        
        for(long int c=0;c<noContigs;c++)
        {
            if(j->joinType==-1)
                break;
            if(likelihoods[tempContigEnd][c]<=0 || checkLikelihood(j->likelihood,likelihoods[tempContigEnd][c])==1)
                continue;
            if(isMember(intraContigs, getContigRoot(c))==1 || getContigRoot(c)==contigOut)
                continue;
            if(gaps[tempContigEnd][c]*stretchFactor<gapOut)
            {
                j->joinType=-1;
                break;
            }
        }
        
        //out reversed
        
        for(long int c=0;c<tempContigEnd;c++)
        {
            if(j->joinType==-1)
                break;
            if(likelihoodsReversed[c][tempContigEnd]<=0 || checkLikelihood(j->likelihood,likelihoodsReversed[c][tempContigEnd])==1)
                continue;
            if(isMember(intraContigs, getContigRoot(c))==1 || getContigRoot(c)==contigOut)
                continue;
            if(gapsReversed[c][tempContigEnd]*stretchFactor<gapOut)
            {
                j->joinType=-1;
                break;
            }
        }
        for(long int c=tempContigEnd+1;c<noContigs;c++)
        {
            if(j->joinType==-1)
                break;
            if(likelihoodsReversed[tempContigEnd][c]<=0 || checkLikelihood(j->likelihood,likelihoodsReversed[tempContigEnd][c])==1)
                continue;
            if(isMember(intraContigs, getContigRoot(c))==1 || getContigRoot(c)==contigOut)
                continue;
            if(gapsReversed[tempContigEnd][c]*stretchFactor<gapOut)
            {
                j->joinType=-1;
                break;
            }
        }
        
        //in reversed
        
        for(long int c=0;c<tempContigStart;c++)
        {
            if(j->joinType==-1)
                break;
            if(likelihoodsReversed[tempContigStart][c]<=0 || checkLikelihood(j->likelihood,likelihoodsReversed[tempContigStart][c])==1)
                continue;
            if(isMember(intraContigs, getContigRoot(c))==1 || getContigRoot(c)==contigIn)
                continue;
            if(gapsReversed[tempContigStart][c]*stretchFactor<gapIn)
            {
                j->joinType=-1;
                break;
            }
        }
        for(long int c=tempContigStart+1;c<noContigs;c++)
        {
            if(j->joinType==-1)
                break;
            if(likelihoodsReversed[c][tempContigStart]<=0 || checkLikelihood(j->likelihood,likelihoodsReversed[c][tempContigStart])==1)
                continue;
            if(isMember(intraContigs, getContigRoot(c))==1 || getContigRoot(c)==contigIn)
                continue;
            if(gapsReversed[c][tempContigStart]*stretchFactor<gapIn)
            {
                j->joinType=-1;
                break;
            }
        }
        
        if(j->joinType!=-1)
        {
            
            Interval *newInterval= new Interval();
            newInterval->contigNo=tempContigNo;
            newInterval->gap=gap-j->gap-tempContigLength;
            newInterval->end=gap-j->gap;
            newInterval->selected=0;
            newInterval->likelihood=j->likelihood;
            newInterval->count=j->count;
            newInterval->realGap=j->gap;
            
            
            if((joinType==0 && j->joinType==0) || (joinType==1 && j->joinType==1)  || (joinType==2 && j->joinType==0))
            {
                newInterval->contigOrientation=0;
            }
            else if((joinType==0 && j->joinType==2) || (joinType==1 && j->joinType==0) || (joinType==2 && j->joinType==2))
            {
                newInterval->contigOrientation=1;
            }

            
            insertInterval(newInterval);
        }
        else
        {
            QEntry *qe=new QEntry();
            qe->contigNo=tempContigRoot;
            qe->likelihood=j->likelihood;
            quarantine.push_back(qe);
            
        }

        
    }
    
    
    selectIntraContigs();
    
    long int curContig=contigNo1;
    
    int curOrientation;
    
    if(joinType==0)
    {
        curOrientation=0;
    }
    else if(joinType==1)
    {
        curOrientation=0;
    }
    else if(joinType==2)
    {
        curOrientation=1;
    }
    
    int curPosition=0;
    int tempJoinType,tempJoined;
    double curLikelihood=0;
    int tempGap;
    
    for(int i=0;i<intervals.size();i++)
    {
//        cout<<intervals[i]->contigNo<<" "<<contigNames[intervals[i]->contigNo]<<" "<<intervals[i]->contigOrientation<<" "<<intervals[i]->gap<<" "<<intervals[i]->end<<" "<<intervals[i]->selected<<endl;

       
        
        if(intervals[i]->selected==1)
        {
            orderReversed=0;
            if(curOrientation==0 && intervals[i]->contigOrientation==0)
            {
                tempJoinType=0;
            }
            else if(curOrientation==0 && intervals[i]->contigOrientation==1)
            {
                tempJoinType=1;
            }
            else if(curOrientation==1 && intervals[i]->contigOrientation==0)
            {
                tempJoinType=2;
            }
            else if(curOrientation==1 && intervals[i]->contigOrientation==1)
            {
                tempJoinType=0;
                orderReversed=1;
            }

            if(tempJoinType!=-1)
            {
                tempContigRoot=getContigRoot(intervals[i]->contigNo);
                tempContigStart=contigStarts[tempContigRoot];
                tempContigEnd=contigEnds[tempContigRoot];
                
                if(orderReversed==0)
                {
                    if(intervals[i]->gap-curPosition>10)
                    {
                        tempJoined=joinContigs(curContig, intervals[i]->contigNo, tempJoinType, intervals[i]->gap-curPosition, 0, curLikelihood>intervals[i]->likelihood?curLikelihood:intervals[i]->likelihood);
                    }
                    else
                    {
                        tempGap=-getOverlap(curContig, intervals[i]->contigNo, tempJoinType);
                        if(tempGap<0)
                        {
                            tempJoined=joinContigs(curContig, intervals[i]->contigNo, tempJoinType, tempGap, 0, curLikelihood>intervals[i]->likelihood?curLikelihood:intervals[i]->likelihood);
                            
                        }
                        else
                        {
                            tempJoined=joinContigs(curContig, intervals[i]->contigNo, tempJoinType, 10, 0, curLikelihood>intervals[i]->likelihood?curLikelihood:intervals[i]->likelihood);

                        }
                    }
                    
                    
                }
                else
                {
                    if(intervals[i]->gap-curPosition>10)
                    {
                        tempJoined=joinContigs(intervals[i]->contigNo,curContig, tempJoinType, intervals[i]->gap-curPosition, 0,curLikelihood>intervals[i]->likelihood?curLikelihood:intervals[i]->likelihood);
                    }
                    else
                    {
                        tempGap=-getOverlap(intervals[i]->contigNo, curContig, tempJoinType);
                        if(tempGap<0)
                        {
                            tempJoined=joinContigs(intervals[i]->contigNo,curContig, tempJoinType, tempGap, 0,curLikelihood>intervals[i]->likelihood?curLikelihood:intervals[i]->likelihood);
                        }
                        else
                        {
                            tempJoined=joinContigs(intervals[i]->contigNo,curContig, tempJoinType, 10, 0,curLikelihood>intervals[i]->likelihood?curLikelihood:intervals[i]->likelihood);
                        }
                        
                    }
                }
                if(tempJoined==1)
                {
                    curOrientation=intervals[i]->contigOrientation;
                    curLikelihood=intervals[i]->likelihood;
                    if(curOrientation==0)
                    {
                //      curContig=contigEnds[getContigRoot(intervals[i]->contigNo)];
                        curContig=tempContigEnd;
                    }
                    else if(curOrientation==1)
                    {
                //       curContig=contigStarts[getContigRoot(intervals[i]->contigNo)];
                         curContig=tempContigStart;
                        
                    }
                    curPosition=intervals[i]->end;
                    intervals[i]->contigStart=tempContigStart;
                    intervals[i]->contigEnd=tempContigEnd;
                    
                    
                }
                else
                {
                    intervals[i]->selected=0;
                }
                
            }
            else
            {
                QEntry *qe=new QEntry();
                qe->contigNo=getContigRoot(intervals[i]->contigNo);
                qe->likelihood=intervals[i]->likelihood;
                quarantine.push_back(qe);
                
            }
            
        }
       
        //       delete intervals[i];
        
    }
    
    int endOrientation;
    
    if(joinType==0)
    {
        endOrientation=0;
    }
    else if(joinType==1)
    {
        endOrientation=1;
    }
    else if(joinType==2)
    {
        endOrientation=0;
    }

    orderReversed=0;
    
    if(curOrientation==0 && endOrientation==0)
    {
        tempJoinType=0;
    }
    else if(curOrientation==0 && endOrientation==1)
    {
        tempJoinType=1;
    }
    else if(curOrientation==1 && endOrientation==0)
    {
        tempJoinType=2;
    }
    else if(curOrientation==1 && endOrientation==1)
    {
        tempJoinType=0;
        orderReversed=1;
        
    }
    if(orderReversed==0)
    {
        if(gap-curPosition>10)
        {
            tempJoined=joinContigs(curContig, contigNo2, tempJoinType, gap-curPosition, 0, likelihoodIn);

        }
        else
        {
            tempGap=-getOverlap(curContig, contigNo2, tempJoinType);
            if(tempGap<0)
            {
                tempJoined=joinContigs(curContig, contigNo2, tempJoinType, tempGap, 0, likelihoodIn);

            }
            else
            {
                tempJoined=joinContigs(curContig, contigNo2, tempJoinType, 10, 0, likelihoodIn);
            }
        }
    }
    else
    {
        if(gap-curPosition>10)
        {
            tempJoined=joinContigs(contigNo2,curContig,tempJoinType, gap-curPosition, 0, likelihoodIn);

        }
        else
        {
            tempGap=-getOverlap(contigNo2, curContig, tempJoinType);
            if(tempGap<0)
            {
                tempJoined=joinContigs(contigNo2,curContig,tempJoinType, tempGap, 0, likelihoodIn);
            }
            else
            {
                tempJoined=joinContigs(contigNo2,curContig,tempJoinType, 10, 0, likelihoodIn);
            }
        }
    }
    
    int gap1, gap2;
    long int start,end;
    
    for(int i=0;i<intervals.size();i++)
    {
        if(intervals[i]->selected==1)
        {
            
            gap1=intervals[i]->gap+contigLength1;
            gap2=gap-intervals[i]->end+contigLength2;
            tempContigStart=intervals[i]->contigStart;
            tempContigEnd=intervals[i]->contigEnd;
            
            if(intervals[i]->contigOrientation==0)
            {
                for(int c=0;c<noContigs;c++)
                {
                    if(contigLengths[c]==-1)
                        continue;
                
                    if(getContigRoot(c)==getContigRoot(intervals[i]->contigNo))
                        continue;
                    
                    
                    start=contigStarts[c];
                    end=contigEnds[c];
                    
                    
                    if(likelihoods[end][tempContigStart]>0)
                    {
                        if(gaps[end][tempContigStart]*stretchFactor>gap1)
                        {
                            if(joinType==0 || joinType==1)
                            {
                                Join *tempJoin=new Join();
                                tempJoin->contigNo1=end;
                                tempJoin->contigNo2=contigNo1;
                                tempJoin->likelihood=likelihoods[end][tempContigStart];
                                tempJoin->joinType=0;
                                tempJoin->gap=gaps[end][tempContigStart]-gap1;
                                extraJoins2.push_back(tempJoin);
                                
                            }
                            else if(joinType==2)
                            {
                                Join *tempJoin=new Join();
                                tempJoin->contigNo1=end<contigNo1?end:contigNo1;
                                tempJoin->contigNo2=end>contigNo1?end:contigNo1;
                                tempJoin->likelihood=likelihoods[end][tempContigStart];
                                tempJoin->joinType=1;
                                tempJoin->gap=gaps[end][tempContigStart]-gap1;
                                extraJoins2.push_back(tempJoin);

                            }
                        }
                    }
                    if(likelihoods[tempContigEnd][start]>0)
                    {
                        if(gaps[tempContigEnd][start]*stretchFactor>gap2)
                        {
                            if(joinType==0 || joinType==2)
                            {
                                Join *tempJoin=new Join();
                                tempJoin->contigNo1=contigNo2;
                                tempJoin->contigNo2=start;
                                tempJoin->likelihood=likelihoods[tempContigEnd][start];
                                tempJoin->joinType=0;
                                tempJoin->gap=gaps[tempContigEnd][start]-gap2;
                                extraJoins1.push_back(tempJoin);
                                
                            }
                            else if(joinType==1)
                            {
                                Join *tempJoin=new Join();
                                tempJoin->contigNo1=start>contigNo2?start:contigNo2;
                                tempJoin->contigNo2=start<contigNo2?start:contigNo2;
                                tempJoin->likelihood=likelihoods[tempContigEnd][start];
                                tempJoin->joinType=2;
                                tempJoin->gap=gaps[tempContigEnd][start]-gap2;
                                extraJoins1.push_back(tempJoin);
                            }
                        }
                    }
                    if(start>tempContigStart)
                    {
                        if(likelihoodsReversed[start][tempContigStart]>0)
                        {
                            if(gapsReversed[start][tempContigStart]*stretchFactor>gap1)
                            {
                                if(joinType==0 || joinType==1)
                                {
                                    Join *tempJoin=new Join();
                                    tempJoin->contigNo1=start>contigNo1?start:contigNo1;
                                    tempJoin->contigNo2=start<contigNo1?start:contigNo1;
                                    tempJoin->likelihood=likelihoodsReversed[start][tempContigStart];
                                    tempJoin->joinType=2;
                                    tempJoin->gap=gapsReversed[start][tempContigStart]-gap1;
                                    extraJoins2.push_back(tempJoin);
                                    
                                }
                                else if(joinType==2)
                                {
                                    Join *tempJoin=new Join();
                                    tempJoin->contigNo1=contigNo1;
                                    tempJoin->contigNo2=start;
                                    tempJoin->likelihood=likelihoodsReversed[start][tempContigStart];
                                    tempJoin->joinType=0;
                                    tempJoin->gap=gapsReversed[start][tempContigStart]-gap1;
                                    extraJoins2.push_back(tempJoin);
                                }
                            }
                        }
                    }
                    else if(start<tempContigStart)
                    {
                        if(likelihoodsReversed[tempContigStart][start]>0)
                        {
                            if(gapsReversed[tempContigStart][start]*stretchFactor>gap1)
                            {
                                if(joinType==0 || joinType==1)
                                {
                                    Join *tempJoin=new Join();
                                    tempJoin->contigNo1=start>contigNo1?start:contigNo1;
                                    tempJoin->contigNo2=start<contigNo1?start:contigNo1;
                                    tempJoin->likelihood=likelihoodsReversed[tempContigStart][start];
                                    tempJoin->joinType=2;
                                    tempJoin->gap=gapsReversed[tempContigStart][start]-gap1;
                                    extraJoins2.push_back(tempJoin);
                                    
                                }
                                else if(joinType==2)
                                {
                                    Join *tempJoin=new Join();
                                    tempJoin->contigNo1=contigNo1;
                                    tempJoin->contigNo2=start;
                                    tempJoin->likelihood=likelihoodsReversed[tempContigStart][start];
                                    tempJoin->joinType=0;
                                    tempJoin->gap=gapsReversed[tempContigStart][start]-gap1;
                                    extraJoins2.push_back(tempJoin);
                                }
                            }
                        }
                        
                    }
                    if(end<tempContigEnd)
                    {
                        if(likelihoodsReversed[end][tempContigEnd]>0)
                        {
                            if(gapsReversed[end][tempContigEnd]*stretchFactor>gap2)
                            {
                                if(joinType==0 || joinType==2)
                                {
                                    Join *tempJoin=new Join();
                                    tempJoin->contigNo1=end<contigNo2?end:contigNo2;
                                    tempJoin->contigNo2=end>contigNo2?end:contigNo2;
                                    tempJoin->likelihood=likelihoodsReversed[end][tempContigEnd];
                                    tempJoin->joinType=1;
                                    tempJoin->gap=gapsReversed[end][tempContigEnd]-gap2;
                                    extraJoins1.push_back(tempJoin);
                                    
                                }
                                else if(joinType==1)
                                {
                                    Join *tempJoin=new Join();
                                    tempJoin->contigNo1=end;
                                    tempJoin->contigNo2=contigNo2;
                                    tempJoin->likelihood=likelihoodsReversed[end][tempContigEnd];
                                    tempJoin->joinType=0;
                                    tempJoin->gap=gapsReversed[end][tempContigEnd]-gap2;
                                    extraJoins1.push_back(tempJoin);
                                }
                            }
                        }
                    }
                    else if(end>tempContigEnd)
                    {
                        if(likelihoodsReversed[tempContigEnd][end]>0)
                        {
                            if(gapsReversed[tempContigEnd][end]*stretchFactor>gap2)
                            {
                                if(joinType==0 || joinType==2)
                                {
                                    Join *tempJoin=new Join();
                                    tempJoin->contigNo1=end<contigNo2?end:contigNo2;
                                    tempJoin->contigNo2=end>contigNo2?end:contigNo2;
                                    tempJoin->likelihood=likelihoodsReversed[tempContigEnd][end];
                                    tempJoin->joinType=1;
                                    tempJoin->gap=gapsReversed[tempContigEnd][end]-gap2;
                                    extraJoins1.push_back(tempJoin);
                                    
                                }
                                else if(joinType==1)
                                {
                                    Join *tempJoin=new Join();
                                    tempJoin->contigNo1=end;
                                    tempJoin->contigNo2=contigNo2;
                                    tempJoin->likelihood=likelihoodsReversed[tempContigEnd][end];
                                    tempJoin->joinType=0;
                                    tempJoin->gap=gapsReversed[tempContigEnd][end]-gap2;
                                    extraJoins1.push_back(tempJoin);
                                }
                            }
                        }
                    }
                
                }
            }
            else if(intervals[i]->contigOrientation==1)
            {
                for(int c=0;c<noContigs;c++)
                {
                    if(contigLengths[c]==-1)
                        continue;
                    
                    if(getContigRoot(c)==getContigRoot(intervals[i]->contigNo))
                        continue;
                    
                    
                    start=contigStarts[c];
                    end=contigEnds[c];
                    
                    
                    if(likelihoods[end][tempContigStart]>0)
                    {
                        if(gaps[end][tempContigStart]*stretchFactor>gap2)
                        {
                            if(joinType==0 || joinType==2)
                            {
                                Join *tempJoin=new Join();
                                tempJoin->contigNo1=end<contigNo2?end:contigNo2;
                                tempJoin->contigNo2=end>contigNo2?end:contigNo2;
                                tempJoin->likelihood=likelihoods[end][tempContigStart];
                                tempJoin->joinType=1;
                                tempJoin->gap=gaps[end][tempContigStart]-gap2;
                                extraJoins1.push_back(tempJoin);
                                
                            }
                            else if(joinType==1)
                            {
                                Join *tempJoin=new Join();
                                tempJoin->contigNo1=end;
                                tempJoin->contigNo2=contigNo2;
                                tempJoin->likelihood=likelihoods[end][tempContigStart];
                                tempJoin->joinType=0;
                                tempJoin->gap=gaps[end][tempContigStart]-gap2;
                                extraJoins1.push_back(tempJoin);
                                
                            }
                        }
                    }
                    if(likelihoods[tempContigEnd][start]>0)
                    {
                        if(gaps[tempContigEnd][start]*stretchFactor>gap1)
                        {
                            if(joinType==2)
                            {
                                Join *tempJoin=new Join();
                                tempJoin->contigNo1=contigNo1;
                                tempJoin->contigNo2=start;
                                tempJoin->likelihood=likelihoods[tempContigEnd][start];
                                tempJoin->joinType=0;
                                tempJoin->gap=gaps[tempContigEnd][start]-gap1;
                                extraJoins2.push_back(tempJoin);
                                
                            }
                            else if(joinType==0 || joinType==1)
                            {
                                Join *tempJoin=new Join();
                                tempJoin->contigNo1=start>contigNo1?start:contigNo1;
                                tempJoin->contigNo2=start<contigNo1?start:contigNo1;
                                tempJoin->likelihood=likelihoods[tempContigEnd][start];
                                tempJoin->joinType=2;
                                tempJoin->gap=gaps[tempContigEnd][start]-gap1;
                                extraJoins2.push_back(tempJoin);
                            }
                        }
                    }
                    if(start>tempContigStart)
                    {
                        if(likelihoodsReversed[start][tempContigStart]>0)
                        {
                            if(gapsReversed[start][tempContigStart]*stretchFactor>gap2)
                            {
                                if(joinType==1)
                                {
                                    Join *tempJoin=new Join();
                                    tempJoin->contigNo1=start>contigNo2?start:contigNo2;
                                    tempJoin->contigNo2=start<contigNo2?start:contigNo2;
                                    tempJoin->likelihood=likelihoodsReversed[start][tempContigStart];
                                    tempJoin->joinType=2;
                                    tempJoin->gap=gapsReversed[start][tempContigStart]-gap2;
                                    extraJoins1.push_back(tempJoin);
                                    
                                }
                                else if(joinType==0 || joinType==2)
                                {
                                    Join *tempJoin=new Join();
                                    tempJoin->contigNo1=contigNo2;
                                    tempJoin->contigNo2=start;
                                    tempJoin->likelihood=likelihoodsReversed[start][tempContigStart];
                                    tempJoin->joinType=0;
                                    tempJoin->gap=gapsReversed[start][tempContigStart]-gap2;
                                    extraJoins1.push_back(tempJoin);
                                }
                            }
                        }
                    }
                    else if(start<tempContigStart)
                    {
                        if(likelihoodsReversed[tempContigStart][start]>0)
                        {
                            if(gapsReversed[tempContigStart][start]*stretchFactor>gap2)
                            {
                                if(joinType==1)
                                {
                                    Join *tempJoin=new Join();
                                    tempJoin->contigNo1=start>contigNo2?start:contigNo2;
                                    tempJoin->contigNo2=start<contigNo2?start:contigNo2;
                                    tempJoin->likelihood=likelihoodsReversed[tempContigStart][start];
                                    tempJoin->joinType=2;
                                    tempJoin->gap=gapsReversed[tempContigStart][start]-gap2;
                                    extraJoins1.push_back(tempJoin);
                                    
                                }
                                else if(joinType==0 || joinType==2)
                                {
                                    Join *tempJoin=new Join();
                                    tempJoin->contigNo1=contigNo2;
                                    tempJoin->contigNo2=start;
                                    tempJoin->likelihood=likelihoodsReversed[tempContigStart][start];
                                    tempJoin->joinType=0;
                                    tempJoin->gap=gapsReversed[tempContigStart][start]-gap2;
                                    extraJoins1.push_back(tempJoin);
                                }
                            }
                        }
                        
                    }
                    if(end<tempContigEnd)
                    {
                        if(likelihoodsReversed[end][tempContigEnd]>0)
                        {
                            if(gapsReversed[end][tempContigEnd]*stretchFactor>gap1)
                            {
                                if(joinType==2)
                                {
                                    Join *tempJoin=new Join();
                                    tempJoin->contigNo1=end<contigNo1?end:contigNo1;
                                    tempJoin->contigNo2=end>contigNo1?end:contigNo1;
                                    tempJoin->likelihood=likelihoodsReversed[end][tempContigEnd];
                                    tempJoin->joinType=1;
                                    tempJoin->gap=gapsReversed[end][tempContigEnd]-gap1;
                                    extraJoins2.push_back(tempJoin);
                                    
                                }
                                else if(joinType==0 || joinType==1)
                                {
                                    Join *tempJoin=new Join();
                                    tempJoin->contigNo1=end;
                                    tempJoin->contigNo2=contigNo1;
                                    tempJoin->likelihood=likelihoodsReversed[end][tempContigEnd];
                                    tempJoin->joinType=0;
                                    tempJoin->gap=gapsReversed[end][tempContigEnd]-gap1;
                                    extraJoins2.push_back(tempJoin);
                                }
                            }
                        }
                    }
                    else if(end>tempContigEnd)
                    {
                        if(likelihoodsReversed[tempContigEnd][end]>0)
                        {
                            if(gapsReversed[tempContigEnd][end]*stretchFactor>gap1)
                            {
                                if(joinType==2)
                                {
                                    Join *tempJoin=new Join();
                                    tempJoin->contigNo1=end<contigNo1?end:contigNo1;
                                    tempJoin->contigNo2=end>contigNo1?end:contigNo1;
                                    tempJoin->likelihood=likelihoodsReversed[tempContigEnd][end];
                                    tempJoin->joinType=1;
                                    tempJoin->gap=gapsReversed[tempContigEnd][end]-gap1;
                                    extraJoins2.push_back(tempJoin);
                                    
                                }
                                else if(joinType==0 || joinType==1)
                                {
                                    Join *tempJoin=new Join();
                                    tempJoin->contigNo1=end;
                                    tempJoin->contigNo2=contigNo1;
                                    tempJoin->likelihood=likelihoodsReversed[tempContigEnd][end];
                                    tempJoin->joinType=0;
                                    tempJoin->gap=gapsReversed[tempContigEnd][end]-gap1;
                                    extraJoins2.push_back(tempJoin);
                                }
                            }
                        }
                    }
                    
                }
                
            }
            
            
            
            addExtraJoins();
            
            while(extraJoins1.size()!=0)
            {
                j=extraJoins1.front();
                delete j;
                extraJoins1.erase(extraJoins1.begin());
            }
            while(extraJoins2.size()!=0)
            {
                j=extraJoins2.front();
                delete j;
                extraJoins2.erase(extraJoins2.begin());
            }
        }
        
        delete intervals[i];
    }
    
    intervals.clear();
    
    
    while(intraJoins1.size()!=0)
    {
        j=intraJoins1.front();
//        cout<<contigNames[j->contigNo1]<<" "<<contigNames[j->contigNo2]<<" "<<j->joinType<<" "<<j->gap<<" "<<j->likelihood<<endl;
        delete j;
        intraJoins1.erase(intraJoins1.begin());
    }
//    cout<<endl;

    while(intraJoins2.size()!=0)
    {
        j=intraJoins2.front();
//        cout<<contigNames[j->contigNo1]<<" "<<contigNames[j->contigNo2]<<" "<<j->joinType<<" "<<j->gap<<" "<<j->likelihood<<endl;
        delete j;
        intraJoins2.erase(intraJoins2.begin());
    }
//    cout<<endl;

    
  

}



int selectJoin(Join *join)
{
 
    long int contigNo1, contigNo2;
    long long int contigLength1, contigLength2;
    long int contigRoot1, contigRoot2;
    long int contigStart1, contigStart2, contigEnd1, contigEnd2;
    int toJoin=1;
  
    
    int gap;
    
    long int start, end;
    
    contigNo1=join->contigNo1;
    contigNo2=join->contigNo2;
    
    contigRoot1=getContigRoot(contigNo1);
    contigRoot2=getContigRoot(contigNo2);
    
    contigLength1=contigLengths[contigRoot1];
    contigLength2=contigLengths[contigRoot2];

    contigStart1=contigStarts[contigRoot1];
    contigStart2=contigStarts[contigRoot2];
    
    contigEnd1=contigEnds[contigRoot1];
    contigEnd2=contigEnds[contigRoot2];

    
//    cout<<"--------------------------------------------------------"<<endl;
    
    if(join->joinType==0)
    {
        gap=gaps[contigNo1][contigNo2];
        
//        cout<<contigLength1<<" "<<contigLength2<<" "<<gap<<endl;
        
        for(int i=0;i<noContigs;i++)
        {
            if(contigLengths[i]==-1)
                continue;
            
            start=contigStarts[i];
            end=contigEnds[i];
      
            //outdegree of contig 1
      
            if(start!=contigNo2 && likelihoods[contigNo1][start]>0)
            {
                if(gaps[contigNo1][start]+contigLengths[i]<=gap*stretchFactor)
                {
                    
                    Join *tempJoin=new Join();
                    tempJoin->contigNo1=contigNo1;
                    tempJoin->contigNo2=start;
                    tempJoin->likelihood=likelihoods[contigNo1][start];
                    tempJoin->count=counts[contigNo1][start];
                    tempJoin->joinType=0;
                    tempJoin->gap=gaps[contigNo1][start];
                    intraJoins1.push_back(tempJoin);
                    
                }
                else if(gap+contigLength2<=gaps[contigNo1][start]*stretchFactor)
                {
                    Join *tempJoin=new Join();
                    tempJoin->contigNo1=contigEnds[getContigRoot(contigNo2)];
                    tempJoin->contigNo2=start;
                    tempJoin->likelihood=likelihoods[contigNo1][start];
                    tempJoin->count=counts[contigNo1][start];
                    tempJoin->joinType=0;
                    tempJoin->gap=gaps[contigNo1][start]-gap-contigLength2;
                    extraJoins1.push_back(tempJoin);

                    
                }
                else if(checkLikelihood(join->likelihood,likelihoods[contigNo1][start])==0)
                {
                    toJoin=0;
                    break;
                }
                
            }
            
            //indegree of contig 2

            if(end!=contigNo1 && likelihoods[end][contigNo2]>0)
            {
                if(gaps[end][contigNo2]+contigLengths[i]<=gap*stretchFactor)
                {
                    Join *tempJoin=new Join();
                    tempJoin->contigNo1=end;
                    tempJoin->contigNo2=contigNo2;
                    tempJoin->likelihood=likelihoods[end][contigNo2];
                    tempJoin->count=counts[end][contigNo2];
                    tempJoin->joinType=0;
                    tempJoin->gap=gaps[end][contigNo2];
                    intraJoins2.push_back(tempJoin);

                    
                }
                else if(gap+contigLength1<=gaps[end][contigNo2]*stretchFactor)
                {
                    Join *tempJoin=new Join();
                    tempJoin->contigNo1=end;
                    tempJoin->contigNo2=contigStarts[getContigRoot(contigNo1)];;
                    tempJoin->likelihood=likelihoods[end][contigNo2];
                    tempJoin->count=counts[end][contigNo2];
                    tempJoin->joinType=0;
                    tempJoin->gap=gaps[end][contigNo2]-gap-contigLength1;
                    extraJoins2.push_back(tempJoin);

                    
                }
                else if(checkLikelihood(join->likelihood,likelihoods[end][contigNo2])==0)
                {
                    toJoin=0;
                    break;
                }
                
            }
            
            //outdegree of contig 1 with other one reversed

            if(end<contigNo1)
            {
                if(likelihoodsReversed[end][contigNo1]>0)
                {
                    if(gapsReversed[end][contigNo1]+contigLengths[i]<=gap*stretchFactor)
                    {
                        Join *tempJoin=new Join();
                        tempJoin->contigNo1=end;
                        tempJoin->contigNo2=contigNo1;
                        tempJoin->likelihood=likelihoodsReversed[end][contigNo1];
                        tempJoin->count=countsReversed[end][contigNo1];
                        tempJoin->joinType=1;
                        tempJoin->gap=gapsReversed[end][contigNo1];
                        intraJoins1.push_back(tempJoin);

                    }
                    else if(gap+contigLength2<=gapsReversed[end][contigNo1]*stretchFactor)
                    {
                        Join *tempJoin=new Join();
                        tempJoin->contigNo1=end<contigEnd2?end:contigEnd2;
                        tempJoin->contigNo2=end>contigEnd2?end:contigEnd2;
                        tempJoin->likelihood=likelihoodsReversed[end][contigNo1];
                        tempJoin->count=countsReversed[end][contigNo1];
                        tempJoin->joinType=1;
                        tempJoin->gap=gapsReversed[end][contigNo1]-gap-contigLength2;
                        extraJoins1.push_back(tempJoin);
                        
                    }
                    else if(checkLikelihood(join->likelihood,likelihoodsReversed[end][contigNo1])==0)
                    {
                        toJoin=0;
                        break;
                    }
                
                }
            }
            else if(end>contigNo1)
            {
                if(likelihoodsReversed[contigNo1][end]>0)
                {
                    if(gapsReversed[contigNo1][end]+contigLengths[i]<=gap*stretchFactor)
                    {
                        Join *tempJoin=new Join();
                        tempJoin->contigNo1=contigNo1;
                        tempJoin->contigNo2=end;
                        tempJoin->likelihood=likelihoodsReversed[contigNo1][end];
                        tempJoin->count=countsReversed[contigNo1][end];
                        tempJoin->joinType=1;
                        tempJoin->gap=gapsReversed[contigNo1][end];
                        intraJoins1.push_back(tempJoin);
                    }
                    else if(gap+contigLength2<=gapsReversed[contigNo1][end]*stretchFactor)
                    {
                        Join *tempJoin=new Join();
                        tempJoin->contigNo1=end<contigEnd2?end:contigEnd2;
                        tempJoin->contigNo2=end>contigEnd2?end:contigEnd2;
                        tempJoin->likelihood=likelihoodsReversed[contigNo1][end];
                        tempJoin->count=countsReversed[contigNo1][end];
                        tempJoin->joinType=1;
                        tempJoin->gap=gapsReversed[contigNo1][end]-gap-contigLength2;
                        extraJoins1.push_back(tempJoin);
                        
                    }
                    else if(checkLikelihood(join->likelihood,likelihoodsReversed[contigNo1][end])==0)
                    {
                        toJoin=0;
                        break;
                    }
                    
                }
                
            }
            //indegree of contig 2 with other one reversed
            if(start<contigNo2)
            {
                if(likelihoodsReversed[contigNo2][start]>0)
                {
                    if(gapsReversed[contigNo2][start]+contigLengths[i]<=gap*stretchFactor)
                    {
                        Join *tempJoin=new Join();
                        tempJoin->contigNo1=contigNo2;
                        tempJoin->contigNo2=start;
                        tempJoin->likelihood=likelihoodsReversed[contigNo2][start];
                        tempJoin->count=countsReversed[contigNo2][start];
                        tempJoin->joinType=2;
                        tempJoin->gap=gapsReversed[contigNo2][start];
                        intraJoins2.push_back(tempJoin);
                    }
                    else if(gap+contigLength1<=gapsReversed[contigNo2][start]*stretchFactor)
                    {
                        Join *tempJoin=new Join();
                        tempJoin->contigNo1=start>contigStart1?start:contigStart1;
                        tempJoin->contigNo2=start<contigStart1?start:contigStart1;
                        tempJoin->likelihood=likelihoodsReversed[contigNo2][start];
                        tempJoin->count=countsReversed[contigNo2][start];
                        tempJoin->joinType=2;
                        tempJoin->gap=gapsReversed[contigNo2][start]-gap-contigLength1;
                        extraJoins2.push_back(tempJoin);
                        
                    }
                    else if(checkLikelihood(join->likelihood,likelihoodsReversed[contigNo2][start])==0)
                    {
                        toJoin=0;
                        break;
                    }
                    
                }
            }
            else if(start>contigNo2)
            {
                if(likelihoodsReversed[start][contigNo2]>0)
                {
                    if(gapsReversed[start][contigNo2]+contigLengths[i]<=gap*stretchFactor)
                    {
                        Join *tempJoin=new Join();
                        tempJoin->contigNo1=start;
                        tempJoin->contigNo2=contigNo2;
                        tempJoin->likelihood=likelihoodsReversed[start][contigNo2];
                        tempJoin->count=countsReversed[start][contigNo2];
                        tempJoin->joinType=2;
                        tempJoin->gap=gapsReversed[start][contigNo2];
                        intraJoins2.push_back(tempJoin);
                    }
                    else if(gap+contigLength1<=gapsReversed[start][contigNo2]*stretchFactor)
                    {
                        Join *tempJoin=new Join();
                        tempJoin->contigNo1=start>contigStart1?start:contigStart1;
                        tempJoin->contigNo2=start<contigStart1?start:contigStart1;
                        tempJoin->likelihood=likelihoodsReversed[start][contigNo2];
                        tempJoin->count=countsReversed[start][contigNo2];
                        tempJoin->joinType=2;
                        tempJoin->gap=gapsReversed[start][contigNo2]-gap-contigLength1;
                        extraJoins2.push_back(tempJoin);

                        
                    }
                    else if(checkLikelihood(join->likelihood,likelihoodsReversed[start][contigNo2])==0)
                    {
                        toJoin=0;
                        break;
                    }
                    
                }
                
            }
            
        }
        
        
    }
    
    if(join->joinType==1)
    {
        gap=gapsReversed[contigNo1][contigNo2];
     
        for(int i=0;i<noContigs;i++)
        {
            if(contigLengths[i]==-1)
                continue;
            
            start=contigStarts[i];
            end=contigEnds[i];
            
            //outdegree of contig 1
            
            if(likelihoods[contigNo1][start]>0)
            {
                if(gaps[contigNo1][start]+contigLengths[i]<=gap*stretchFactor)
                {
                    Join *tempJoin=new Join();
                    tempJoin->contigNo1=contigNo1;
                    tempJoin->contigNo2=start;
                    tempJoin->likelihood=likelihoods[contigNo1][start];
                    tempJoin->count=counts[contigNo1][start];
                    tempJoin->joinType=0;
                    tempJoin->gap=gaps[contigNo1][start];
                    intraJoins1.push_back(tempJoin);
                    
                }
                else if(gap+contigLength2<=gaps[contigNo1][start]*stretchFactor)
                {
                    Join *tempJoin=new Join();
                    tempJoin->contigNo1=start>contigStart2?start:contigStart2;
                    tempJoin->contigNo2=start<contigStart2?start:contigStart2;
                    tempJoin->likelihood=likelihoods[contigNo1][start];
                    tempJoin->count=counts[contigNo1][start];
		      tempJoin->joinType=2;
                    tempJoin->gap=gaps[contigNo1][start]-gap-contigLength2;
                    extraJoins1.push_back(tempJoin);
                    
                }
                else if(checkLikelihood(join->likelihood,likelihoods[contigNo1][start])==0)
                {
                    toJoin=0;
                    break;
                }
                
            }
            
            //outdegree of contig 2
            
            if(likelihoods[contigNo2][start]>0)
            {
                if(gaps[contigNo2][start]+contigLengths[i]<=gap*stretchFactor)
                {
                    Join *tempJoin=new Join();
                    tempJoin->contigNo1=contigNo2;
                    tempJoin->contigNo2=start;
                    tempJoin->likelihood=likelihoods[contigNo2][start];
                    tempJoin->count=counts[contigNo2][start];
                    tempJoin->joinType=0;
                    tempJoin->gap=gaps[contigNo2][start];
                    intraJoins2.push_back(tempJoin);
                    
                }
                else if(gap+contigLength1<=gaps[contigNo2][start]*stretchFactor)
                {
                    Join *tempJoin=new Join();
                    tempJoin->contigNo1=start>contigStart1?start:contigStart1;
                    tempJoin->contigNo2=start<contigStart1?start:contigStart1;
                    tempJoin->likelihood=likelihoods[contigNo2][start];
                    tempJoin->count=counts[contigNo2][start];
                    tempJoin->joinType=2;
                    tempJoin->gap=gaps[contigNo2][start]-gap-contigLength1;
                    extraJoins2.push_back(tempJoin);
                    
                }
                else if(checkLikelihood(join->likelihood,likelihoods[contigNo2][start])==0)
                {
                    toJoin=0;
                    break;
                }
                
            }
            
            //outdegree of contig 1 with other one reversed
            
            if(end<contigNo1)
            {
                if(end !=contigNo2 && likelihoodsReversed[end][contigNo1]>0)
                {
                    if(gapsReversed[end][contigNo1]+contigLengths[i]<=gap*stretchFactor)
                    {
                        Join *tempJoin=new Join();
                        tempJoin->contigNo1=end;
                        tempJoin->contigNo2=contigNo1;
                        tempJoin->likelihood=likelihoodsReversed[end][contigNo1];
                        tempJoin->count=countsReversed[end][contigNo1];
                        tempJoin->joinType=1;
                        tempJoin->gap=gapsReversed[end][contigNo1];
                        intraJoins1.push_back(tempJoin);
                    }
                    else if(gap+contigLength2<=gapsReversed[end][contigNo1]*stretchFactor)
                    {
                        Join *tempJoin=new Join();
                        tempJoin->contigNo1=end;
                        tempJoin->contigNo2=contigStart2;
                        tempJoin->likelihood=likelihoodsReversed[end][contigNo1];
                        tempJoin->count=countsReversed[end][contigNo1];
                        tempJoin->joinType=0;
                        tempJoin->gap=gapsReversed[end][contigNo1]-gap-contigLength2;
                        extraJoins1.push_back(tempJoin);
                        
                    }
                    else if(checkLikelihood(join->likelihood,likelihoodsReversed[end][contigNo1])==0)
                    {
                        toJoin=0;
                        break;
                    }
                    
                }
            }
            else if(end>contigNo1)
            {
                if(end !=contigNo2 && likelihoodsReversed[contigNo1][end]>0)
                {
                    if(gapsReversed[contigNo1][end]+contigLengths[i]<=gap*stretchFactor)
                    {
                        Join *tempJoin=new Join();
                        tempJoin->contigNo1=contigNo1;
                        tempJoin->contigNo2=end;
                        tempJoin->likelihood=likelihoodsReversed[contigNo1][end];
                        tempJoin->count=countsReversed[contigNo1][end];
                        tempJoin->joinType=1;
                        tempJoin->gap=gapsReversed[contigNo1][end];
                        intraJoins1.push_back(tempJoin);
                    }
                    else if(gap+contigLength2<=gapsReversed[contigNo1][end]*stretchFactor)
                    {
                        Join *tempJoin=new Join();
                        tempJoin->contigNo1=end;
                        tempJoin->contigNo2=contigStart2;
                        tempJoin->likelihood=likelihoodsReversed[contigNo1][end];
                        tempJoin->count=countsReversed[contigNo1][end];
                        tempJoin->joinType=0;
                        tempJoin->gap=gapsReversed[contigNo1][end]-gap-contigLength2;
                        extraJoins1.push_back(tempJoin);
                        
                    }
                    else if(checkLikelihood(join->likelihood,likelihoodsReversed[contigNo1][end])==0)
                    {
                        toJoin=0;
                        break;
                    }
                    
                }
                
            }
            //outdegree of contig 2 with other one reversed
            if(end<contigNo2)
            {
                if(end !=contigNo1 && likelihoodsReversed[end][contigNo2]>0)
                {
                    if(gapsReversed[end][contigNo2]+contigLengths[i]<=gap*stretchFactor)
                    {
                        Join *tempJoin=new Join();
                        tempJoin->contigNo1=end;
                        tempJoin->contigNo2=contigNo2;
                        tempJoin->likelihood=likelihoodsReversed[end][contigNo2];
                        tempJoin->count=countsReversed[end][contigNo2];
                        tempJoin->joinType=1;
                        tempJoin->gap=gapsReversed[end][contigNo2];
                        intraJoins2.push_back(tempJoin);
                    }
                    else if(gap+contigLength1<=gapsReversed[end][contigNo2]*stretchFactor)
                    {
                        Join *tempJoin=new Join();
                        tempJoin->contigNo1=end;
                        tempJoin->contigNo2=contigStart1;
                        tempJoin->likelihood=likelihoodsReversed[end][contigNo2];
                        tempJoin->count=countsReversed[end][contigNo2];			   
                        tempJoin->joinType=0;
                        tempJoin->gap=gapsReversed[end][contigNo2]-gap-contigLength1;
                        extraJoins2.push_back(tempJoin);
                        
                    }
                    else if(checkLikelihood(join->likelihood,likelihoodsReversed[end][contigNo2])==0)
                    {
                        toJoin=0;
                        break;
                    }
                    
                }
            }
            else if(end>contigNo2)
            {
                if(end !=contigNo1 && likelihoodsReversed[contigNo2][end]>0)
                {
                    if(gapsReversed[contigNo2][end]+contigLengths[i]<=gap*stretchFactor)
                    {
                        Join *tempJoin=new Join();
                        tempJoin->contigNo1=contigNo2;
                        tempJoin->contigNo2=end;
                        tempJoin->likelihood=likelihoodsReversed[contigNo2][end];
                        tempJoin->count=countsReversed[contigNo2][end];
                        tempJoin->joinType=1;
                        tempJoin->gap=gapsReversed[contigNo2][end];
                        intraJoins2.push_back(tempJoin);
                    }
                    else if(gap+contigLength1<=gapsReversed[contigNo2][end]*stretchFactor)
                    {
                        Join *tempJoin=new Join();
                        tempJoin->contigNo1=end;
                        tempJoin->contigNo2=contigStart1;
                        tempJoin->likelihood=likelihoodsReversed[contigNo2][end];
                        tempJoin->count=countsReversed[contigNo2][end];
                        tempJoin->joinType=0;
                        tempJoin->gap=gapsReversed[contigNo2][end]-gap-contigLength1;
                        extraJoins2.push_back(tempJoin);
                        
                    }
                    else if(checkLikelihood(join->likelihood,likelihoodsReversed[contigNo2][end])==0)
                    {
                        toJoin=0;
                        break;
                    }
                    
                }
                
            }
            
        }
        
    }

    if(join->joinType==2)
    {
        gap=gapsReversed[contigNo1][contigNo2];
        for(int i=0;i<noContigs;i++)
        {
            if(contigLengths[i]==-1)
                continue;
            
            start=contigStarts[i];
            end=contigEnds[i];
            
            //indegree of contig 1
            
            if(likelihoods[end][contigNo1]>0)
            {
                if(gaps[end][contigNo1]+contigLengths[i]<=gap*stretchFactor)
                {
                    Join *tempJoin=new Join();
                    tempJoin->contigNo1=end;
                    tempJoin->contigNo2=contigNo1;
                    tempJoin->likelihood=likelihoods[end][contigNo1];
                    tempJoin->count=counts[end][contigNo1];
                    tempJoin->joinType=0;
                    tempJoin->gap=gaps[end][contigNo1];
                    intraJoins1.push_back(tempJoin);
                }
                else if(gap+contigLength2<=gaps[end][contigNo1]*stretchFactor)
                {
                    Join *tempJoin=new Join();
                    tempJoin->contigNo1=end<contigEnd2?end:contigEnd2;
                    tempJoin->contigNo2=end>contigEnd2?end:contigEnd2;
                    tempJoin->likelihood=likelihoods[end][contigNo1];
                    tempJoin->count=counts[end][contigNo1];
                    tempJoin->joinType=1;
                    tempJoin->gap=gaps[end][contigNo1]-gap-contigLength2;
                    extraJoins1.push_back(tempJoin);
                    
                }
                else if(checkLikelihood(join->likelihood,likelihoods[end][contigNo1])==0)
                {
                    toJoin=0;
                    break;
                }
                
            }
            
            //indegree of contig 2
            
            if(likelihoods[end][contigNo2]>0)
            {
                if(gaps[end][contigNo2]+contigLengths[i]<=gap*stretchFactor)
                {
                    Join *tempJoin=new Join();
                    tempJoin->contigNo1=end;
                    tempJoin->contigNo2=contigNo2;
                    tempJoin->likelihood=likelihoods[end][contigNo2];
                    tempJoin->count=counts[end][contigNo2];
                    tempJoin->joinType=0;
                    tempJoin->gap=gaps[end][contigNo2];
                    intraJoins2.push_back(tempJoin);
                    
                }
                else if(gap+contigLength1<=gaps[end][contigNo2]*stretchFactor)
                {
                    Join *tempJoin=new Join();
                    tempJoin->contigNo1=end<contigEnd1?end:contigEnd1;
                    tempJoin->contigNo2=end>contigEnd1?end:contigEnd1;
                    tempJoin->likelihood=likelihoods[end][contigNo2];
                    tempJoin->count=counts[end][contigNo2];
   			tempJoin->joinType=1;
                    tempJoin->gap=gaps[end][contigNo2]-gap-contigLength1;
                    extraJoins2.push_back(tempJoin);
                    
                }
                else if(checkLikelihood(join->likelihood,likelihoods[end][contigNo2])==0)
                {
                    toJoin=0;
                    break;
                }
                
            }
            
            //indegree of contig 1 with other one reversed
            
            if(start<contigNo1)
            {
                if(start !=contigNo2 && likelihoodsReversed[contigNo1][start]>0)
                {
                    if(gapsReversed[contigNo1][start]+contigLengths[i]<=gap*stretchFactor)
                    {
                        Join *tempJoin=new Join();
                        tempJoin->contigNo1=contigNo1;
                        tempJoin->contigNo2=start;
                        tempJoin->likelihood=likelihoodsReversed[contigNo1][start];
                        tempJoin->count=countsReversed[contigNo1][start];
                        tempJoin->joinType=2;
                        tempJoin->gap=gapsReversed[contigNo1][start];
                        intraJoins1.push_back(tempJoin);

                    }
                    else if(gap+contigLength2<=gapsReversed[contigNo1][start]*stretchFactor)
                    {
                        Join *tempJoin=new Join();
                        tempJoin->contigNo1=contigEnd2;
                        tempJoin->contigNo2=start;
                        tempJoin->likelihood=likelihoodsReversed[contigNo1][start];
                        tempJoin->count=countsReversed[contigNo1][start];
                        tempJoin->joinType=0;
                        tempJoin->gap=gapsReversed[contigNo1][start]-gap-contigLength2;
                        extraJoins1.push_back(tempJoin);
                        
                    }
                    else if(checkLikelihood(join->likelihood,likelihoodsReversed[contigNo1][start])==0)
                    {
                        toJoin=0;
                        break;
                    }
                    
                }
            }
            else if(start>contigNo1)
            {
                if(start !=contigNo2 && likelihoodsReversed[start][contigNo1]>0)
                {
                    if(gapsReversed[start][contigNo1]+contigLengths[i]<=gap*stretchFactor)
                    {
                        Join *tempJoin=new Join();
                        tempJoin->contigNo1=start;
                        tempJoin->contigNo2=contigNo1;
                        tempJoin->likelihood=likelihoodsReversed[start][contigNo1];
                        tempJoin->count=countsReversed[start][contigNo1];
                        tempJoin->joinType=2;
                        tempJoin->gap=gapsReversed[start][contigNo1];
                        intraJoins1.push_back(tempJoin);
                    }
                    else if(gap+contigLength2<=gapsReversed[start][contigNo1]*stretchFactor)
                    {
                        Join *tempJoin=new Join();
                        tempJoin->contigNo1=contigEnd2;
                        tempJoin->contigNo2=start;
                        tempJoin->likelihood=likelihoodsReversed[start][contigNo1];
                        tempJoin->count=countsReversed[start][contigNo1];
                        tempJoin->joinType=0;
                        tempJoin->gap=gapsReversed[start][contigNo1]-gap-contigLength2;
                        extraJoins1.push_back(tempJoin);

                        
                    }
                    else if(checkLikelihood(join->likelihood,likelihoodsReversed[start][contigNo1])==0)
                    {
                        toJoin=0;
                        break;
                    }
                    
                }
                
            }
            
            //indegree of contig 2 with other one reversed
            if(start<contigNo2)
            {
                if(start !=contigNo1 && likelihoodsReversed[contigNo2][start]>0)
                {
                    if(gapsReversed[contigNo2][start]+contigLengths[i]<=gap*stretchFactor)
                    {
                        Join *tempJoin=new Join();
                        tempJoin->contigNo1=contigNo2;
                        tempJoin->contigNo2=start;
                        tempJoin->likelihood=likelihoodsReversed[contigNo2][start];
                        tempJoin->count=countsReversed[contigNo2][start];
                        tempJoin->joinType=2;
                        tempJoin->gap=gapsReversed[contigNo2][start];
                        intraJoins2.push_back(tempJoin);
                    }
                    else if(gap+contigLength1<=gapsReversed[contigNo2][start]*stretchFactor)
                    {
                        Join *tempJoin=new Join();
                        tempJoin->contigNo1=contigEnd1;
                        tempJoin->contigNo2=start;
                        tempJoin->likelihood=likelihoodsReversed[contigNo2][start];
                        tempJoin->count=countsReversed[contigNo2][start];
                        tempJoin->joinType=0;
                        tempJoin->gap=gapsReversed[contigNo2][start]-gap-contigLength1;
                        extraJoins2.push_back(tempJoin);
                        
                    }
                    else if(checkLikelihood(join->likelihood,likelihoodsReversed[contigNo2][start])==0)
                    {
                        toJoin=0;
                        break;
                    }
                    
                }
            }
            else if(start>contigNo2)
            {
                if(start !=contigNo1 && likelihoodsReversed[start][contigNo2]>0)
                {
                    if(gapsReversed[start][contigNo2]+contigLengths[i]<=gap*stretchFactor)
                    {
                        Join *tempJoin=new Join();
                        tempJoin->contigNo1=start;
                        tempJoin->contigNo2=contigNo2;
                        tempJoin->likelihood=likelihoodsReversed[start][contigNo2];
                        tempJoin->count=countsReversed[start][contigNo2];
                        tempJoin->joinType=2;
                        tempJoin->gap=gapsReversed[start][contigNo2];
                        intraJoins2.push_back(tempJoin);
                    }
                    else if(gap+contigLength1<=gapsReversed[start][contigNo2]*stretchFactor)
                    {
                        Join *tempJoin=new Join();
                        tempJoin->contigNo1=contigEnd1;
                        tempJoin->contigNo2=start;
                        tempJoin->likelihood=likelihoodsReversed[start][contigNo2];
                        tempJoin->count=countsReversed[start][contigNo2];
                        tempJoin->joinType=0;
                        tempJoin->gap=gapsReversed[start][contigNo2]-gap-contigLength1;
                        extraJoins2.push_back(tempJoin);
                        
                    }
                    else if(checkLikelihood(join->likelihood,likelihoodsReversed[start][contigNo2])==0)
                    {
                        toJoin=0;
                        break;
                    }
                    
                }
                
            }
            
        }
        
        
    }

    if(toJoin==1)
    {
        
        joinPath(contigNo1, contigNo2, join->joinType,gap,join->likelihood);
        return 1;
        
    }
    else
    {
        
        for(int i=0;i<intraJoins1.size();i++)
        {
            delete intraJoins1[i];
        }
        for(int i=0;i<intraJoins2.size();i++)
        {
            delete intraJoins2[i];
        }

        for(int i=0;i<extraJoins1.size();i++)
        {
            delete extraJoins1[i];
        }
        for(int i=0;i<extraJoins2.size();i++)
        {
            delete extraJoins2[i];
        }
        intraJoins1.clear();
        intraJoins2.clear();

        extraJoins1.clear();
        extraJoins2.clear();

    }
    
    return 0;
}

void joinPair(long int contigNo1, long int contigNo2, int joinType, int gap, double likelihoodIn)
{
    
    Join *j;
    
//    cout<<contigNames[contigNo1]<<endl;
//    cout<<contigNames[contigNo2]<<endl;

    
    int isDummy=0;
    
    
    int joined=joinContigs(contigNo1, contigNo2, joinType,MIN_GAP-1,isDummy,likelihoodIn);

    if(joined==0)
    {
        while(intraJoins1.size()!=0)
        {
            j=intraJoins1.front();
            delete j;
            intraJoins1.erase(intraJoins1.begin());
        }

        while(intraJoins2.size()!=0)
        {
            j=intraJoins2.front();
            delete j;
            intraJoins2.erase(intraJoins2.begin());
        }
        return;
    }

	long int tempContigNo,tempContigRoot;

    for(int i=0;i<intraJoins1.size();i++)
    {
        j=intraJoins1[i];
        
        if(j->contigNo1==contigNo1)
        {
            tempContigNo=j->contigNo2;
        }
        else
        {
            tempContigNo=j->contigNo1;
        }
        tempContigRoot=getContigRoot(tempContigNo);
                
        QEntry *qe=new QEntry();
        qe->contigNo=tempContigRoot;
        qe->likelihood=j->likelihood;
        quarantine.push_back(qe);
            
    }
	
    for(int i=0;i<intraJoins2.size();i++)
    {
        j=intraJoins2[i];
        
        if(j->contigNo1==contigNo2)
        {
            tempContigNo=j->contigNo2;
        }
        else
        {
            tempContigNo=j->contigNo1;
        }
        tempContigRoot=getContigRoot(tempContigNo);

        QEntry *qe=new QEntry();
        qe->contigNo=tempContigRoot;
        qe->likelihood=j->likelihood;
        quarantine.push_back(qe);
            
    }
    while(intraJoins1.size()!=0)
    {
        j=intraJoins1.front();
//        cout<<contigNames[j->contigNo1]<<" "<<contigNames[j->contigNo2]<<" "<<j->joinType<<" "<<j->gap<<" "<<j->likelihood<<endl;
        delete j;
        intraJoins1.erase(intraJoins1.begin());
    }
//    cout<<endl;

    while(intraJoins2.size()!=0)
    {
        j=intraJoins2.front();
//        cout<<contigNames[j->contigNo1]<<" "<<contigNames[j->contigNo2]<<" "<<j->joinType<<" "<<j->gap<<" "<<j->likelihood<<endl;
        delete j;
        intraJoins2.erase(intraJoins2.begin());
    }
//    cout<<endl;




}

int selectJoin_Conservative(Join *join)
{
 
    long int contigNo1, contigNo2;
    long long int contigLength1, contigLength2;
    long int contigRoot1, contigRoot2;
    long int contigStart1, contigStart2, contigEnd1, contigEnd2;
    int toJoin=1;
  
    
    int gap;
    
    long int start, end;
    
    contigNo1=join->contigNo1;
    contigNo2=join->contigNo2;
    
    contigRoot1=getContigRoot(contigNo1);
    contigRoot2=getContigRoot(contigNo2);
    
    contigLength1=contigLengths[contigRoot1];
    contigLength2=contigLengths[contigRoot2];

    contigStart1=contigStarts[contigRoot1];
    contigStart2=contigStarts[contigRoot2];
    
    contigEnd1=contigEnds[contigRoot1];
    contigEnd2=contigEnds[contigRoot2];

    
//    cout<<"--------------------------------------------------------"<<endl;
    
    if(join->joinType==0)
    {
        gap=gaps[contigNo1][contigNo2];
        
//        cout<<contigLength1<<" "<<contigLength2<<" "<<gap<<endl;
        
        for(int i=0;i<noContigs;i++)
        {
            if(contigLengths[i]==-1)
                continue;
            
            start=contigStarts[i];
            end=contigEnds[i];
      
            //outdegree of contig 1
      
            if(start!=contigNo2 && likelihoods[contigNo1][start]>0)
            {
               if(checkLikelihood(join->likelihood,likelihoods[contigNo1][start])==0)
                {
                    toJoin=0;
                    break;
                }
  
               if(gaps[contigNo1][start]+contigLengths[i]<=gap*stretchFactor)
                {
                    
                    Join *tempJoin=new Join();
                    tempJoin->contigNo1=contigNo1;
                    tempJoin->contigNo2=start;
                    tempJoin->likelihood=likelihoods[contigNo1][start];
                    tempJoin->joinType=0;
                    tempJoin->gap=gaps[contigNo1][start];
                    intraJoins1.push_back(tempJoin);
                    
                }
                else if(gap+contigLength2<=gaps[contigNo1][start]*stretchFactor)
                {
                     
                }
                
            }
            
            //indegree of contig 2

            if(end!=contigNo1 && likelihoods[end][contigNo2]>0)
            {
                if(checkLikelihood(join->likelihood,likelihoods[end][contigNo2])==0)
                {
                    toJoin=0;
                    break;
                }

                if(gaps[end][contigNo2]+contigLengths[i]<=gap*stretchFactor)
                {
                    Join *tempJoin=new Join();
                    tempJoin->contigNo1=end;
                    tempJoin->contigNo2=contigNo2;
                    tempJoin->likelihood=likelihoods[end][contigNo2];
                    tempJoin->joinType=0;
                    tempJoin->gap=gaps[end][contigNo2];
                    intraJoins2.push_back(tempJoin);

                    
                }
                else if(gap+contigLength1<=gaps[end][contigNo2]*stretchFactor)
                {
                    
                }
                
            }
            
            //outdegree of contig 1 with other one reversed

            if(end<contigNo1)
            {
                if(likelihoodsReversed[end][contigNo1]>0)
                {
                    if(checkLikelihood(join->likelihood,likelihoodsReversed[end][contigNo1])==0)
                    {
                        toJoin=0;
                        break;
                    }

                    if(gapsReversed[end][contigNo1]+contigLengths[i]<=gap*stretchFactor)
                    {
                        Join *tempJoin=new Join();
                        tempJoin->contigNo1=end;
                        tempJoin->contigNo2=contigNo1;
                        tempJoin->likelihood=likelihoodsReversed[end][contigNo1];
                        tempJoin->joinType=1;
                        tempJoin->gap=gapsReversed[end][contigNo1];
                        intraJoins1.push_back(tempJoin);

                    }
                    else if(gap+contigLength2<=gapsReversed[end][contigNo1]*stretchFactor)
                    {

                    }
                
                }
            }
            else if(end>contigNo1)
            {
                if(likelihoodsReversed[contigNo1][end]>0)
                {
                    if(checkLikelihood(join->likelihood,likelihoodsReversed[contigNo1][end])==0)
                    {
                        toJoin=0;
                        break;
                    }

                    if(gapsReversed[contigNo1][end]+contigLengths[i]<=gap*stretchFactor)
                    {
                        Join *tempJoin=new Join();
                        tempJoin->contigNo1=contigNo1;
                        tempJoin->contigNo2=end;
                        tempJoin->likelihood=likelihoodsReversed[contigNo1][end];
                        tempJoin->joinType=1;
                        tempJoin->gap=gapsReversed[contigNo1][end];
                        intraJoins1.push_back(tempJoin);
                    }
                    else if(gap+contigLength2<=gapsReversed[contigNo1][end]*stretchFactor)
                    {
                        
                    }
                    
                }
                
            }
            //indegree of contig 2 with other one reversed
            if(start<contigNo2)
            {
                if(likelihoodsReversed[contigNo2][start]>0)
                {
                    if(checkLikelihood(join->likelihood,likelihoodsReversed[contigNo2][start])==0)
                    {
                        toJoin=0;
                        break;
                    }

                    if(gapsReversed[contigNo2][start]+contigLengths[i]<=gap*stretchFactor)
                    {
                        Join *tempJoin=new Join();
                        tempJoin->contigNo1=contigNo2;
                        tempJoin->contigNo2=start;
                        tempJoin->likelihood=likelihoodsReversed[contigNo2][start];
                        tempJoin->joinType=2;
                        tempJoin->gap=gapsReversed[contigNo2][start];
                        intraJoins2.push_back(tempJoin);
                    }
                    else if(gap+contigLength1<=gapsReversed[contigNo2][start]*stretchFactor)
                    {
                        
                    }
                    
                }
            }
            else if(start>contigNo2)
            {
                if(likelihoodsReversed[start][contigNo2]>0)
                {
                    if(checkLikelihood(join->likelihood,likelihoodsReversed[start][contigNo2])==0)
                    {
                        toJoin=0;
                        break;
                    }

                    if(gapsReversed[start][contigNo2]+contigLengths[i]<=gap*stretchFactor)
                    {
                        Join *tempJoin=new Join();
                        tempJoin->contigNo1=start;
                        tempJoin->contigNo2=contigNo2;
                        tempJoin->likelihood=likelihoodsReversed[start][contigNo2];
                        tempJoin->joinType=2;
                        tempJoin->gap=gapsReversed[start][contigNo2];
                        intraJoins2.push_back(tempJoin);
                    }
                    else if(gap+contigLength1<=gapsReversed[start][contigNo2]*stretchFactor)
                    {
                        
                    }
                    
                }
                
            }
            
        }
        
        
    }
    
    if(join->joinType==1)
    {
        gap=gapsReversed[contigNo1][contigNo2];
     
        for(int i=0;i<noContigs;i++)
        {
            if(contigLengths[i]==-1)
                continue;
            
            start=contigStarts[i];
            end=contigEnds[i];
            
            //outdegree of contig 1
            
            if(likelihoods[contigNo1][start]>0)
            {
                if(checkLikelihood(join->likelihood,likelihoods[contigNo1][start])==0)
                {
                    toJoin=0;
                    break;
                }

                if(gaps[contigNo1][start]+contigLengths[i]<=gap*stretchFactor)
                {
                    Join *tempJoin=new Join();
                    tempJoin->contigNo1=contigNo1;
                    tempJoin->contigNo2=start;
                    tempJoin->likelihood=likelihoods[contigNo1][start];
                    tempJoin->joinType=0;
                    tempJoin->gap=gaps[contigNo1][start];
                    intraJoins1.push_back(tempJoin);
                    
                }
                else if(gap+contigLength2<=gaps[contigNo1][start]*stretchFactor)
                {
                    
                }
                
            }
            
            //outdegree of contig 2
            
            if(likelihoods[contigNo2][start]>0)
            {
                if(checkLikelihood(join->likelihood,likelihoods[contigNo2][start])==0)
                {
                    toJoin=0;
                    break;
                }

                if(gaps[contigNo2][start]+contigLengths[i]<=gap*stretchFactor)
                {
                    Join *tempJoin=new Join();
                    tempJoin->contigNo1=contigNo2;
                    tempJoin->contigNo2=start;
                    tempJoin->likelihood=likelihoods[contigNo2][start];
                    tempJoin->joinType=0;
                    tempJoin->gap=gaps[contigNo2][start];
                    intraJoins2.push_back(tempJoin);
                    
                }
                else if(gap+contigLength1<=gaps[contigNo2][start]*stretchFactor)
                {
                    
                }
                
            }
            
            //outdegree of contig 1 with other one reversed
            
            if(end<contigNo1)
            {
                if(end !=contigNo2 && likelihoodsReversed[end][contigNo1]>0)
                {
                    if(checkLikelihood(join->likelihood,likelihoodsReversed[end][contigNo1])==0)
                    {
                        toJoin=0;
                        break;
                    }

                    if(gapsReversed[end][contigNo1]+contigLengths[i]<=gap*stretchFactor)
                    {
                        Join *tempJoin=new Join();
                        tempJoin->contigNo1=end;
                        tempJoin->contigNo2=contigNo1;
                        tempJoin->likelihood=likelihoodsReversed[end][contigNo1];
                        tempJoin->joinType=1;
                        tempJoin->gap=gapsReversed[end][contigNo1];
                        intraJoins1.push_back(tempJoin);
                    }
                    else if(gap+contigLength2<=gapsReversed[end][contigNo1]*stretchFactor)
                    {
                        
                    }
                    
                }
            }
            else if(end>contigNo1)
            {

                if(end !=contigNo2 && likelihoodsReversed[contigNo1][end]>0)
                {
                    if(checkLikelihood(join->likelihood,likelihoodsReversed[contigNo1][end])==0)
                    {
                        toJoin=0;
                        break;
                    }

                    if(gapsReversed[contigNo1][end]+contigLengths[i]<=gap*stretchFactor)
                    {
                        Join *tempJoin=new Join();
                        tempJoin->contigNo1=contigNo1;
                        tempJoin->contigNo2=end;
                        tempJoin->likelihood=likelihoodsReversed[contigNo1][end];
                        tempJoin->joinType=1;
                        tempJoin->gap=gapsReversed[contigNo1][end];
                        intraJoins1.push_back(tempJoin);
                    }
                    else if(gap+contigLength2<=gapsReversed[contigNo1][end]*stretchFactor)
                    {
                        
                    }
                    
                }
                
            }
            //outdegree of contig 2 with other one reversed
            if(end<contigNo2)
            {
                if(end !=contigNo1 && likelihoodsReversed[end][contigNo2]>0)
                {
                    if(checkLikelihood(join->likelihood,likelihoodsReversed[end][contigNo2])==0)
                    {
                        toJoin=0;
                        break;
                    }

                    if(gapsReversed[end][contigNo2]+contigLengths[i]<=gap*stretchFactor)
                    {
                        Join *tempJoin=new Join();
                        tempJoin->contigNo1=end;
                        tempJoin->contigNo2=contigNo2;
                        tempJoin->likelihood=likelihoodsReversed[end][contigNo2];
                        tempJoin->joinType=1;
                        tempJoin->gap=gapsReversed[end][contigNo2];
                        intraJoins2.push_back(tempJoin);
                    }
                    else if(gap+contigLength1<=gapsReversed[end][contigNo2]*stretchFactor)
                    {
                        
                    }
                    
                }
            }
            else if(end>contigNo2)
            {
                if(end !=contigNo1 && likelihoodsReversed[contigNo2][end]>0)
                {
                    if(checkLikelihood(join->likelihood,likelihoodsReversed[contigNo2][end])==0)
                    {
                        toJoin=0;
                        break;
                    }

                    if(gapsReversed[contigNo2][end]+contigLengths[i]<=gap*stretchFactor)
                    {
                        Join *tempJoin=new Join();
                        tempJoin->contigNo1=contigNo2;
                        tempJoin->contigNo2=end;
                        tempJoin->likelihood=likelihoodsReversed[contigNo2][end];
                        tempJoin->joinType=1;
                        tempJoin->gap=gapsReversed[contigNo2][end];
                        intraJoins2.push_back(tempJoin);
                    }
                    else if(gap+contigLength1<=gapsReversed[contigNo2][end]*stretchFactor)
                    {
                        
                    }
                    
                }
                
            }
            
        }
        
    }

    if(join->joinType==2)
    {
        gap=gapsReversed[contigNo1][contigNo2];
        for(int i=0;i<noContigs;i++)
        {
            if(contigLengths[i]==-1)
                continue;
            
            start=contigStarts[i];
            end=contigEnds[i];
            
            //indegree of contig 1
            
            if(likelihoods[end][contigNo1]>0)
            {
                if(checkLikelihood(join->likelihood,likelihoods[end][contigNo1])==0)
                {
                    toJoin=0;
                    break;
                }

                if(gaps[end][contigNo1]+contigLengths[i]<=gap*stretchFactor)
                {
                    Join *tempJoin=new Join();
                    tempJoin->contigNo1=end;
                    tempJoin->contigNo2=contigNo1;
                    tempJoin->likelihood=likelihoods[end][contigNo1];
                    tempJoin->joinType=0;
                    tempJoin->gap=gaps[end][contigNo1];
                    intraJoins1.push_back(tempJoin);
                }
                else if(gap+contigLength2<=gaps[end][contigNo1]*stretchFactor)
                {
                    
                }
                
            }
            
            //indegree of contig 2
            
            if(likelihoods[end][contigNo2]>0)
            {
                if(checkLikelihood(join->likelihood,likelihoods[end][contigNo2])==0)
                {
                    toJoin=0;
                    break;
                }

                if(gaps[end][contigNo2]+contigLengths[i]<=gap*stretchFactor)
                {
                    Join *tempJoin=new Join();
                    tempJoin->contigNo1=end;
                    tempJoin->contigNo2=contigNo2;
                    tempJoin->likelihood=likelihoods[end][contigNo2];
                    tempJoin->joinType=0;
                    tempJoin->gap=gaps[end][contigNo2];
                    intraJoins2.push_back(tempJoin);
                    
                }
                else if(gap+contigLength1<=gaps[end][contigNo2]*stretchFactor)
                {
                    
                }
                
            }
            
            //indegree of contig 1 with other one reversed
            
            if(start<contigNo1)
            {
                if(start !=contigNo2 && likelihoodsReversed[contigNo1][start]>0)
                {
                    if(checkLikelihood(join->likelihood,likelihoodsReversed[contigNo1][start])==0)
                    {
                        toJoin=0;
                        break;
                    }

                    if(gapsReversed[contigNo1][start]+contigLengths[i]<=gap*stretchFactor)
                    {
                        Join *tempJoin=new Join();
                        tempJoin->contigNo1=contigNo1;
                        tempJoin->contigNo2=start;
                        tempJoin->likelihood=likelihoodsReversed[contigNo1][start];
                        tempJoin->joinType=2;
                        tempJoin->gap=gapsReversed[contigNo1][start];
                        intraJoins1.push_back(tempJoin);

                    }
                    else if(gap+contigLength2<=gapsReversed[contigNo1][start]*stretchFactor)
                    {
                        
                    }
                    
                }
            }
            else if(start>contigNo1)
            {
                if(start !=contigNo2 && likelihoodsReversed[start][contigNo1]>0)
                {
                    if(checkLikelihood(join->likelihood,likelihoodsReversed[start][contigNo1])==0)
                    {
                        toJoin=0;
                        break;
                    }

                    if(gapsReversed[start][contigNo1]+contigLengths[i]<=gap*stretchFactor)
                    {
                        Join *tempJoin=new Join();
                        tempJoin->contigNo1=start;
                        tempJoin->contigNo2=contigNo1;
                        tempJoin->likelihood=likelihoodsReversed[start][contigNo1];
                        tempJoin->joinType=2;
                        tempJoin->gap=gapsReversed[start][contigNo1];
                        intraJoins1.push_back(tempJoin);
                    }
                    else if(gap+contigLength2<=gapsReversed[start][contigNo1]*stretchFactor)
                    {
                        
                    }
                    
                }
                
            }
            
            //indegree of contig 2 with other one reversed
            if(start<contigNo2)
            {
                if(start !=contigNo1 && likelihoodsReversed[contigNo2][start]>0)
                {
                    if(checkLikelihood(join->likelihood,likelihoodsReversed[contigNo2][start])==0)
                    {
                        toJoin=0;
                        break;
                    }

                    if(gapsReversed[contigNo2][start]+contigLengths[i]<=gap*stretchFactor)
                    {
                        Join *tempJoin=new Join();
                        tempJoin->contigNo1=contigNo2;
                        tempJoin->contigNo2=start;
                        tempJoin->likelihood=likelihoodsReversed[contigNo2][start];
                        tempJoin->joinType=2;
                        tempJoin->gap=gapsReversed[contigNo2][start];
                        intraJoins2.push_back(tempJoin);
                    }
                    else if(gap+contigLength1<=gapsReversed[contigNo2][start]*stretchFactor)
                    {
                        
                    }
                    
                }
            }
            else if(start>contigNo2)
            {
                if(start !=contigNo1 && likelihoodsReversed[start][contigNo2]>0)
                {
                    if(checkLikelihood(join->likelihood,likelihoodsReversed[start][contigNo2])==0)
                    {
                        toJoin=0;
                        break;
                    }

                    if(gapsReversed[start][contigNo2]+contigLengths[i]<=gap*stretchFactor)
                    {
                        Join *tempJoin=new Join();
                        tempJoin->contigNo1=start;
                        tempJoin->contigNo2=contigNo2;
                        tempJoin->likelihood=likelihoodsReversed[start][contigNo2];
                        tempJoin->joinType=2;
                        tempJoin->gap=gapsReversed[start][contigNo2];
                        intraJoins2.push_back(tempJoin);
                    }
                    else if(gap+contigLength1<=gapsReversed[start][contigNo2]*stretchFactor)
                    {
                        
                    }
                    
                }
                
            }
            
        }
        
        
    }

    if(toJoin==1)
    {
        
        joinPair(contigNo1, contigNo2, join->joinType,gap,join->likelihood);
        return 1;
        
    }
    else
    {
        
        for(int i=0;i<intraJoins1.size();i++)
        {
            delete intraJoins1[i];
        }
        for(int i=0;i<intraJoins2.size();i++)
        {
            delete intraJoins2[i];
        }

        intraJoins1.clear();
        intraJoins2.clear();

    }
    
    return 0;
}

void findAmbiguousJoins()
{
    
    for(int i=0;i<noContigs;i++)
    {
        if(contigLengths[i]==-1)
            continue;
        
        for(int j=0;j<noContigs;j++)
        {
            if(contigLengths[j]==-1)
                continue;
            
            if(likelihoods[contigEnds[i]][contigStarts[j]]>0 && counts[contigEnds[i]][contigStarts[j]]>MIN_READ_JOIN)
            {
                Join *join=new Join();
                join->contigNo1=contigEnds[i];
                join->contigNo2=contigStarts[j];
                join->likelihood=likelihoods[contigEnds[i]][contigStarts[j]];
                join->joinType=0;
                join->isExtra=0;
			join->count=counts[contigEnds[i]][contigStarts[j]];
                insertJoin(join);
            }
            
        }
        for(int j=i+1;j<noContigs;j++)
        {
            if(contigLengths[j]==-1)
                continue;
            
            if(likelihoodsReversed[contigEnds[i]][contigEnds[j]]>0 && countsReversed[contigEnds[i]][contigEnds[j]]>MIN_READ_JOIN)
            {
                Join *join=new Join();
                join->contigNo1=contigEnds[i];
                join->contigNo2=contigEnds[j];
                join->likelihood=likelihoodsReversed[contigEnds[i]][contigEnds[j]];
                join->joinType=1;
                join->isExtra=0;
			join->count=countsReversed[contigEnds[i]][contigEnds[j]];
                insertJoin(join);
            }
        }

        for(int j=0;j<i;j++)
        {
            if(contigLengths[j]==-1)
                continue;
            
            if(likelihoodsReversed[contigStarts[i]][contigStarts[j]]>0 && countsReversed[contigStarts[i]][contigStarts[j]]>MIN_READ_JOIN)
            {
                Join *join=new Join();
                join->contigNo1=contigStarts[i];
                join->contigNo2=contigStarts[j];
                join->likelihood=likelihoodsReversed[contigStarts[i]][contigStarts[j]];
                join->joinType=2;
                join->isExtra=0;
			join->count=countsReversed[contigStarts[i]][contigStarts[j]];
                //check the order of contigs
                
                insertJoin(join);
            }
        }

        
    }


    
    for(int i=0;i<joins.size();i++)
    {
//        cout<<joins[i]->contigNo1<<" "<<joins[i]->contigNo2<<" "<<joins[i]->likelihood<<" "<<joins[i]->joinType<<endl;
        
    }
    
    Join *j;
    int joined=0;
    while(joins.size()!=0)
    {
        j=joins.front();
        joins.erase(joins.begin());
    /*
        if(j->isExtra==0 || checkFailedJoins(j)==1)
        {
            joined=selectJoin(j);
            if(joined==0)
            {
                failedJoins.push_back(j);
            }
            else
            {
                delete j;
            }
        }
        else
            delete j;
     */

	if(isConservative==0)
	{
        	selectJoin(j);
	}
	else
	{
		selectJoin_Conservative(j);
	}
        delete j;
    }
    
}

void generateScaffolds()
{
    
    scaffoldFile=fopen("scaffolds.fa","w");
    
    for(int i=0;i<noContigs;i++)
    {
        if(contigLengths[i]!=-1)
        {
            
            fprintf(scaffoldFile,">%d\n",i);
            fprintf(scaffoldFile,"%s\n",contigs[i]);
//            cout<<contigLengths[i]<<" "<<i<<" "<<contigNames[contigStarts[i]]<<" "<<contigNames[contigEnds[i]]<<endl;
        }

    }
    
    fclose(scaffoldFile);
/*    
    for(int i=0;i<noContigs;i++)
    {
        if(contigLengths[i]==-1)
            continue;
        
        for(int j=0;j<noContigs;j++)
        {
            if(contigLengths[j]==-1)
                continue;
            
            if(likelihoods[contigEnds[i]][contigStarts[j]]>0)
            {
                cout<<i<<" "<<j<<" "<<likelihoods[contigEnds[i]][contigStarts[j]]<<" "<<gaps[contigEnds[i]][contigStarts[j]]<<endl;
                cout<<contigNames[contigEnds[i]]<<" "<<contigNames[contigStarts[j]]<<endl;
            }
            if(likelihoodsReversed[contigEnds[i]][contigEnds[j]]>0)
            {
                cout<<i<<" "<<j<<" "<<likelihoodsReversed[contigEnds[i]][contigEnds[j]]<<" "<<gapsReversed[contigEnds[i]][contigStarts[j]]<<endl;
                cout<<contigNames[contigEnds[i]]<<" "<<contigNames[contigEnds[j]]<<endl;
            }
        }
        
    }
 */   
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
     
    
    srand (SEED);
//    srand (rand());
    
    

	contigFileName=argv[1];
    
	for(int i=2;i<argc;i++)
	{
		if(strcmp(argv[i],"--help")==0 || strcmp(argv[i],"-h")==0)
			printHelp();

		if(strcmp(argv[i],"--jump")==0 || strcmp(argv[i],"-j")==0)
		{
			isJump=1;
		}
		if(strcmp(argv[i],"--conservative")==0 || strcmp(argv[i],"-c")==0)
		{
			isConservative=1;
		}

	}

    
	contigFile=fopen(contigFileName, "r");
	outFile=fopen("out.txt", "w");
    	joinsFile=fopen("joins.txt", "w");
    
    
	if (contigFile == NULL)
	{
		printf("Can't open contig file\n");
		exit(1);
	}
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
				contigs.push_back(contig);
				contigLengths.push_back(contigLength);
				totalContigLength+=contigLength;
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
	contigs.push_back(contig);
	contigLengths.push_back(contigLength);
	totalContigLength+=contigLength;
    
    
	fclose(contigFile);
    
	for(long int i=0;i<noContigs;i++)
	{
		for(long int j=0;j<contigLengths[i];j++)
		{
			contigs[i][j]=toupper(contigs[i][j]);
		}
	}
    
    //	cout<<"after reading contig file"<<endl;
    
	/*
     use bfast or some other tool to map reads and save mapping
     */

	
	readFromFile();
 
    findUniqueJoins();

    computeLikelihoodDist();

    cout<<likelihoodMean<<" "<<likelihoodSD<<endl;
    
    findAmbiguousJoins();
    
    generateScaffolds();

    fclose(joinsFile);
    
	return 0;
}

inline double Likelihoods::get(int index)
{
	for(int i=0;i<size;i++)
	{
		if(indexes[i]==index)
		{
			return vals[i];
		}
	}
	return getScore(0,0,valUnmapped,contigLength,contigLengthsOriginals[index]);
 }


