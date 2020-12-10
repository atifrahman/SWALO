#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#define NUM_THREADS 32
#define SEED 94703
#define VERSION "0.9.8-beta"

void itoa(int value, char *s, int base)
{
	sprintf(s,"%d",value);
}

char * toupper(char *contig)
{
	long int i=0;
	char c;
	while (contig[i])
  	{
    		c=contig[i];
    		contig[i]=toupper(c);
    		i++;
  	}
	return contig;
}


