/***********************************************************************************
***********************************************************************************/

#include "stdio.h"
#include "math.h"
#include "string.h"
#include "stdlib.h"
#include "float.h"

#define MAX_CHROM_NUM 24       //maximum chromosome number



typedef struct
{
	char chromName[100];
	int chromSize;
}CHROM_INFO;

//Structure defining the histone modification regions

CHROM_INFO chromInfo[MAX_CHROM_NUM];


//total number of bins
int totalBinNum;

int chromNum;

int ChromToIndex(char *chrom);
char *IndexToChrom(int index, char *chrom, int len);
int GetChromInfo(char *fileName);
float **GetBinCount(char *l1FileName, char *l2FileName,int binsize);
int ReadHistoneFile(char *histoneFilename,char *cell_name,int binsize);


int main(int argc, char* argv[])
{
	char histoneFileName[1000],chromFileName[1000],projectName[1000];
        int binsize;
	int count=0;
	if (argc!=5)
	{
		printf("Usage: <Histone File Name> <Chromosome Description File Name> <bin size> <cell-type name>\n");
		return 0;
	}

	strcpy(histoneFileName,argv[1]);
	strcpy(chromFileName,argv[2]);
        binsize=atoi(argv[3]);
        strcpy(projectName,argv[4]);


	//Read chromosome description file
	if (!GetChromInfo(chromFileName))
	{
		printf("chromosome description file is not valid\n");
		return 0;
	}

	//Read tag files of L1 and L2, and define the histone modification regions
float **binrpkm;	
int i,j;
ReadHistoneFile(histoneFileName,projectName,binsize);
	


	return 0;
}

//transform chromosome name the chromosome index
//
//
//
//
//
int ChromToIndex(char *chrom)
{
	int i;

	for (i=0;i<chromNum;i++)
	{
		if (!strcmp(chrom, chromInfo[i].chromName))
		{
			break;
		}
	}

	if (i<chromNum)
	{
		return i;
	}
	else
	{
		return -1;
	}
}

//transfrom chromosome index to chromosome name
char *IndexToChrom(int index, char *chrom, int len)
{
	if ((index<0)||(index>=chromNum))
	{
		return 0;
	}

	if (strlen(chromInfo[index].chromName)>=len)
	{
		return 0;
	}

	strcpy(chrom, chromInfo[index].chromName);

	return chrom;
}

//Read the chromosome description file
int GetChromInfo(char *fileName)
{
	FILE *fh;
	char tmpStr[1000];

	fh = (FILE *)fopen(fileName, "r");

	if (!fh)
	{
		return -1;
	}

	chromNum = 0;

	fscanf(fh, "%s", tmpStr);

	while (!feof(fh))
	{
		strcpy(chromInfo[chromNum].chromName, tmpStr);
		fscanf(fh, "%s", tmpStr);
		chromInfo[chromNum].chromSize = atoi(tmpStr);
		fscanf(fh, "%s", tmpStr);
		chromNum++;
	}

	fclose(fh);
	return chromNum;
}

int ReadHistoneFile(char *histoneFileName,char *cell_name,int BIN_SIZE)
{
FILE *fh,*fh2;
float **avgBinRPKM,**binRPKM;
int i,j,k,num_reps;
char word[100],hist_mod[20],l1FileName[1000],l2FileName[1000];
binRPKM =(float **)malloc(chromNum*sizeof(float *));
avgBinRPKM=(float **)malloc(chromNum*sizeof(float *));
totalBinNum = 0;

        for (i=0;i<chromNum;i++)
        {
                totalBinNum += chromInfo[i].chromSize/BIN_SIZE+1;
		binRPKM[i] = (float *)malloc((chromInfo[i].chromSize/BIN_SIZE+1)*sizeof(float));
 		memset(binRPKM[i], 0, (chromInfo[i].chromSize/BIN_SIZE+1)*sizeof(float));
		avgBinRPKM[i] = (float *)malloc((chromInfo[i].chromSize/BIN_SIZE+1)*sizeof(float));
 		memset(avgBinRPKM[i], 0, (chromInfo[i].chromSize/BIN_SIZE+1)*sizeof(float));
 
       }
fh = (FILE *)fopen(histoneFileName, "r");
 char binFileName[1000];

        if (!fh)
        {
                return -1;
        }

        
        while (!feof(fh))

        { fscanf(fh, "%s",word);
          if((word[0]=='E')&&(word[1]=='O')&&(word[2]=='F')){break;}
          if((word[0]=='M')&&(word[1]=='a'))
          { fscanf(fh, "%s",hist_mod);fscanf(fh,"%d",&num_reps);
          printf("Mark:%s\n",hist_mod);
          sprintf(binFileName,"%s.%s.rpkm",cell_name,hist_mod);
          fh2  = (FILE *)fopen(binFileName, "w");}
 
         else{sprintf(l1FileName,"%s",word);
          fscanf(fh, "%s",l2FileName);
          printf("%s\n%s\n",l1FileName,l2FileName);
          binRPKM=GetBinCount(l1FileName,l2FileName,BIN_SIZE);
for(j=0;j<chromNum;j++){
for(k=0;k<(chromInfo[j].chromSize/BIN_SIZE);k++){
avgBinRPKM[j][k]=binRPKM[j][k];}}
printf("Done with replicate 1!\n");
          for(i=1;i<num_reps;i++){ 
          fscanf(fh, "%s",l1FileName);
          fscanf(fh, "%s",l2FileName);
binRPKM=GetBinCount(l1FileName,l2FileName,BIN_SIZE);

for(j=0;j<chromNum;j++){
for(k=0;k<(chromInfo[j].chromSize/BIN_SIZE);k++){
avgBinRPKM[j][k]=avgBinRPKM[j][k]+binRPKM[j][k];}}
printf("Done with Replicate %d\n",i+1);
}
for (i=0;i<chromNum;i++)
        {
                free(binRPKM[i]);
}
free(binRPKM);

for(j=0;j<chromNum;j++){
for(k=1;k<(chromInfo[j].chromSize/BIN_SIZE);k++){
if(avgBinRPKM[j][k]!=0.0000){fprintf(fh2,"%d\t%d\t%f\n",j+1,k*BIN_SIZE,(avgBinRPKM[j][k]/num_reps));}}}

fclose(fh2);
}
}
fclose(fh);
return 0;
}

//Read the tag files of L1 and L2, and determine the histone modification sites
float **GetBinCount(char *l1FileName, char *l2FileName,int BIN_SIZE)
{
        int l1TagNum,l2TagNum;
	//binCounts and binCounts2 store the fragment counts in each bin. mask=1 flags histone modification site
	int **binCounts1, **binCounts2;
        float **binRPKM; 
	char tmpStr[1000];
	char tmpChrom[100];
	int tmpPos,tmpNeg,tmpIndex;
	char tmpStrand;
	int tmpStart, tmpEnd;
	int i,j,k;
        FILE *fh;
	//Allocate Memory
	binCounts1 = (int **)malloc(chromNum*sizeof(int *));
	binCounts2 = (int **)malloc(chromNum*sizeof(int *));
        binRPKM =(float **)malloc(chromNum*sizeof(float *));
	//Initialize the arrays
	totalBinNum = 0;

	for (i=0;i<chromNum;i++)
	{
		totalBinNum += chromInfo[i].chromSize/BIN_SIZE+1;

		binCounts1[i] = (int *)malloc((chromInfo[i].chromSize/BIN_SIZE+1)*sizeof(int));
		binCounts2[i] = (int *)malloc((chromInfo[i].chromSize/BIN_SIZE+1)*sizeof(int));
binRPKM[i] = (float *)malloc((chromInfo[i].chromSize/BIN_SIZE+1)*sizeof(float));
		memset(binCounts1[i], 0, (chromInfo[i].chromSize/BIN_SIZE+1)*sizeof(int));
		memset(binCounts2[i], 0, (chromInfo[i].chromSize/BIN_SIZE+1)*sizeof(int));
 memset(binRPKM[i], 0, (chromInfo[i].chromSize/BIN_SIZE+1)*sizeof(float));
	}

	// Read library 1 tag file
	l1TagNum = 0;
        char a[100],b[100];
	fh = (FILE *)fopen(l1FileName, "r");


	fscanf(fh, "%s", tmpChrom);

	while (!feof(fh))
	{
		tmpIndex = ChromToIndex(tmpChrom);
		if (tmpIndex<0)
		{
			fscanf(fh, "%d", &tmpPos); fscanf(fh, "%d", &tmpNeg);
	                fscanf(fh, "%s",a);
                        //fscanf(fh, "%s",b);
                        fscanf(fh, "%s",tmpStr);    		
			fscanf(fh, "%s", tmpChrom);
			continue;
		}
               

		fscanf(fh, "%d", &tmpPos);fscanf(fh, "%d", &tmpNeg);
                fscanf(fh, "%s",a);
                //fscanf(fh, "%s",b);
                fscanf(fh, "%s",tmpStr);
                tmpStrand=tmpStr[0];

		if ((tmpPos>=chromInfo[tmpIndex].chromSize)||(tmpPos<0))
		{
			fscanf(fh, "%s", tmpChrom);
			continue;
		}

		l1TagNum++;
                if(tmpStrand=='+'){
		binCounts1[tmpIndex][tmpPos/BIN_SIZE]++;}
                if(tmpStrand=='-'){
                binCounts1[tmpIndex][tmpNeg/BIN_SIZE]++;}

		fscanf(fh,"%s", tmpChrom);
	}
printf("Done reading file1! Number of tags=%d\n",l1TagNum);          
	fclose(fh);

	//read library 2 tag file
	l2TagNum = 0;
	fh = (FILE *)fopen(l2FileName, "r");


	fscanf(fh, "%s", tmpChrom);

	while (!feof(fh))
	{
		tmpIndex = ChromToIndex(tmpChrom);

		if (tmpIndex<0)
		{
			fscanf(fh, "%d", &tmpPos);fscanf(fh, "%d", &tmpNeg);
			fscanf(fh, "%s",a);
			//fscanf(fh, "%s",b);
			fscanf(fh, "%s",tmpStr);
			fscanf(fh, "%s", tmpChrom);
			continue;
		}

		fscanf(fh, "%d", &tmpPos);fscanf(fh, "%d", &tmpNeg);
		fscanf(fh, "%s",a);
		//fscanf(fh, "%s",b);
		fscanf(fh, "%s",tmpStr);
		 tmpStrand=tmpStr[0];
		if ((tmpPos>=chromInfo[tmpIndex].chromSize)||(tmpPos<0))
		{
			fscanf(fh, "%s", tmpChrom);
			continue;
		}


		l2TagNum++;
		 if(tmpStrand=='+'){
                binCounts2[tmpIndex][tmpPos/BIN_SIZE]++;}
                if(tmpStrand=='-'){
                binCounts2[tmpIndex][tmpNeg/BIN_SIZE]++;}
		fscanf(fh,"%s", tmpChrom);
	}
	fclose(fh);
 printf("Done reading file2! Number of tags=%d\n",l2TagNum);
        float binnum1,binnum2;
        for (i=0;i<chromNum;i++)
        {
        for (j=0;j<=(chromInfo[i].chromSize/BIN_SIZE);j++)
        {    
            binnum1=10000000*binCounts1[i][j];
if((j>=2)||(j<(chromInfo[i].chromSize/BIN_SIZE)-2))          
{  binnum2=2000000*(binCounts2[i][j-2]+binCounts2[i][j-1]+binCounts2[i][j]+binCounts2[i][j+1]+binCounts2[i][j+2]);
}
else
{ binnum2=10000000*binCounts2[i][j];}
float abc1 = (binnum1/l1TagNum);
float abc2 = (binnum2/l2TagNum);
binRPKM[i][j]=abc1-abc2;
  
       
        }
        }
	//Free memory

	for (i=0;i<chromNum;i++)
	{
		free(binCounts1[i]);
		free(binCounts2[i]);
		
	}

	free(binCounts1);
	free(binCounts2);
	

	return binRPKM;
}




