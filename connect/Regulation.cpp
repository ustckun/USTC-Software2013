// Regulation.cpp
// Get regualtion matrix and 166 genes' names
//
// version 2.0
// change:
// 1.add output to txt module
// 2.add notation
//
//**********************************************************************************
//this module read TF-TF regualtion from RegulonDB database
//transform + - to +1 -1 and 0 if not relationships it will be 2
//filter the name to class TFIM::geneName
//**********************************************************************************
//
//interface:
//int getGeneAmount():all gene number
//float	originalMartix[200][200]:genes' regualtion
//A to B 's regulation is originalMatrix[B][A]


#include"TFIM.h"
#include"Calculate.h"
#include"Regulation.h"
#include"EasytoDebug.h"
#include"Sequence.h"
#include"GRN.h"


//input:objexts of TFIM, read name and also get regulation, if file is not open int openFileError will be 1
void Regulation::readName(TFIM geneFirst[])
{
	ifstream data("TF-TF");
	if(!data)
	{
		openFileError=1;
	}
	else
		openFileError=0;
	char ch;
	int i,num=0;
	int geneAM;
	char *Name;
	char STRING1[N][10];
	float originalMA[N][N];
	data.get(ch);
	for(i=0;i<N;)
	{
		int a,b;
		//data.get(ch);
		while(ch=='#')//filter RegulonDB line
		{
			string noUse;
			getline(data,noUse);
			data.get(ch);
		}
		while(ch!='	')
		{
			STRING1[i][num]=ch;
			data.get(ch);
			num++;
		}
		STRING1[i][num]='\0';
		Name=STRING1[i];
		int j;
		for(j=0;j<i;j++)
		{
			strlwr(geneFirst[j].name);
			strlwr(Name);
			if(strcmp(geneFirst[j].name,Name)==0)
			{
				a=j;
				j=i+1;
			}
		}
		if(j==i)
		{
			geneFirst[i].name=Name;
			geneFirst[i].putName();
			geneFirst[i].geneNumber=i;
			a=i;
			i++;
		}
		data.get(ch);
		num=0;
		while(ch!='	')
		{
			STRING1[i][num]=ch;
			data.get(ch);
			num++;
		}
		STRING1[i][num]='\0';
		Name=STRING1[i];
		for(j=0;j<i;j++)
		{
			strlwr(geneFirst[j].name);
			strlwr(Name);
			if(strcmp(geneFirst[j].name,Name)==0)
			{
				b=j;
				j=i+1;
			}
		}
		if(j==i)
		{
			geneFirst[i].name=Name;
			geneFirst[i].putName();
			geneFirst[i].geneNumber=i;
			b=i;
			i++;
		}
		//geneFirst[i].putName();
		data.get(ch);
		if(ch=='-')
		{
			if(originalMA[b][a]==1||originalMA[b][a]==0)
			{
				originalMA[b][a]=0;
			}
			else
				originalMA[b][a]=-1;
		}
		else if(ch=='+')
		{
			data.get(ch);
			if(ch=='-')
				originalMA[b][a]=0;
			else
				originalMA[b][a]=1;
		}
		else
			originalMA[b][a]=2;
		string noUse;
		getline(data,noUse);
		num=0;
		if(!data.get(ch))
		{
			geneAM=i;
			i=N;
		}
	}
	memcpy((char *)originalMatrix,(char *)originalMA,sizeof(float)*N*N);
	//originalMatrix=originalMA;
	geneAmount=geneAM;
}

//full fill the matrix with 2
void Regulation::fullFill()
{
	float a;
	int n;
//	a=originalMatrix[0][0];
	for(n=0;n<N;n++)
	{
		for(int m=0;m<N;m++)
		{
			a=originalMatrix[m][n];
			if(a!=1.0&&a!=0.0&&a!=-1.0&&a!=2.0)
			{
				originalMatrix[m][n]=2;
			}
		}
	}
}

//get amount of genes
int Regulation::getGeneAmount()
{
	return geneAmount;
}

//know is file opened
int Regulation::getOpenError()
{
	return openFileError;
}