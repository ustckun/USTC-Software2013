#include"TFIM.h"
#include"Regulation.h"
#define N 10
void Regulation::readName(TFIM *geneFirst)
{
	ifstream data("TF-TF.txt");
	if(!data)
	{
		openFileError=1;
	}
	char ch;
	int i,m,num=0;
	int geneAM;
	char *Name;
	char STRING1[10];
	float originalMA[N][N];
	for(i=1;i<N;)
	{
		int a,b;
		if(!data.get(ch))
		{
			geneAM=i;
			i=N;
		}
		while(ch=='#')//filter RegulonDB line
		{
			string noUse;
			getline(data,noUse);
			data.get(ch);
		}
		while(ch!='	')
		{
			STRING1[num]=ch;
			data.get(ch);
			num++;
		}
		STRING1[num]='\0';
		Name=STRING1;
		int j;
		for(j=1;j<i;j++)
		{
			strlwr(geneFirst[j].name);
			strlwr(Name);
			if(strcmp(geneFirst[j].name,Name)==0)
			{
				a=j;
				j=i+1;
			}
		}
		if(j!=i+1)
		{
			geneFirst[i].name=Name;
			geneFirst[i].geneNumber=i;
			a=i;
			i++;
		}
		data.get(ch);
		num=0;
		while(ch!='	')
		{
			STRING1[num]=ch;
			data.get(ch);
			num++;
		}
		STRING1[num]='\0';
		Name=STRING1;
		for(j=1;j<i;j++)
		{
			strlwr(geneFirst[j].name);
			strlwr(Name);
			if(strcmp(geneFirst[j].name,Name)==0)
			{
				b=j;
				j=i+1;
			}
		}
		if(j!=i+1)
		{
			geneFirst[i].name=Name;
			geneFirst[i].geneNumber=i;
			a=i;
			i++;
		}
		data.get(ch);
		if(ch=='-')
		{
			if(originalMA[b][a]==1||originalMA[b][a]==0.5)
			{
				originalMA[b][a]=0.5;
			}
			else
				originalMA[b][a]=-1;
		}
		else if(ch=='+')
		{
			data.get(ch);
			if(ch=='-')
				originalMA[b][a]=0.5;
			else
				originalMA[b][a]=1;
		}
		else
			originalMA[b][a]=0;
		string noUse;
		getline(data,noUse);
	}
	originalMatrix=&originalMA[0][0];
	geneAmount=geneAM;
}

void Regulation::fullFill()
{
	float a;
	int n;
	a=*originalMatrix;
	for(n=0;n<N*N;n++)
	{
		if(a!=1.0||a!=0.0||a!=0.5||a!=-1.0)
		{
			*originalMatrix=0;
		}
		a=*originalMatrix++;
	}
}

int Regulation::getGeneAmount()
{
	return geneAmount;
}

int Regulation::getOpenError()
{
	return openFileError;
}