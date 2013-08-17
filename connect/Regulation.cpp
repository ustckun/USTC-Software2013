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
#include"ReadDNA.h"
#include"Sequence.h"
#include"GRN.h"
#include"PSOPredict.h"
#include"SBOL.h"
#include"RandSeq.h"


void Regulation::getRegulationMatrix(TFIM geneFirst[])
{
	readTFTF(geneFirst,originalGRN);
	addTF(geneFirst);
	TFAmount=geneAmount;
	readTFGene(geneFirst,originalGRN);
}

void Regulation::addTF(TFIM geneFirst[])
{
	ifstream data("TF-Gene");
	if(!data)
	{
		openFileError=1;
	}
	else
		openFileError=0;
	char ch;
	int i=geneAmount;
	data.get(ch);
	char STRING1[TFScale][10];
	char *Name;
	for(;i<TFScale;)
	{
		int num=0;
		while(ch=='#')
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
			strlwr(Name);
			if(strcmp(geneFirst[j].geneName,Name)==0)
			{
				j=i+1;
			}
		}
		if(j==i)
		{
			geneFirst[i].name=Name;
			geneFirst[i].putName();
			geneFirst[i].geneNumber=i;
			i++;
			geneAmount++;
		}
		string noUse;
		getline(data,noUse);
		num=0;
		if(!data.get(ch))
		{
			i=TFScale;
		}
	}
}
void Regulation::readTFGene(TFIM geneFirst[],double **GRN)
{
	ifstream data("TF-Gene");
	if(!data)
	{
		openFileError=1;
	}
	else
		openFileError=0;
	char ch;
	int i,num=0;
	//int geneAM,TFAM=0;
	char *Name;
	char STRING1[GENEAM][10];
	//double originalMA[N][N];
	int unknow=0;
	data.get(ch);
	for(i=geneAmount;i<GENEAM;)
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
			//strlwr(geneFirst[j].name);
			strlwr(Name);
			if(strcmp(geneFirst[j].geneName,Name)==0)
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
			//strlwr(geneFirst[j].name);
			strlwr(Name);
			if(strcmp(geneFirst[j].geneName,Name)==0)
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
			if(GRN[b][a]==1||GRN[b][a]==2)
			{
				GRN[b][a]=2;
				uncertain1.push_back(b);
				uncertain2.push_back(a);
				unknow++;
			}
			else
				GRN[b][a]=-1;
		}
		else if(ch=='+')
		{
			data.get(ch);
			if(ch=='-')
			{
				GRN[b][a]=2;
				uncertain1.push_back(b);
				uncertain2.push_back(a);
				unknow++;
			}
			else
				GRN[b][a]=1;
		}
		else
			GRN[b][a]=0;
		string noUse;
		getline(data,noUse);
		num=0;
		if(!data.get(ch))
		{
			geneAmount=i;
			i=GENEAM;
		}
		//if(TFAM<a)
			//TFAM=a;
		//cout<<j<<endl;
	}
	//memcpy((char *)originalMatrix,(char *)originalMA,sizeof(double)*N*N);
	//originalMatrix=originalMA;
	//cout<<TFAM<<endl;
	//geneAmount=geneAM;
}

//input:objexts of TFIM, read name and also get regulation, if file is not open int openFileError will be 1
void Regulation::readTFTF(TFIM geneFirst[],double **GRN)
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
	//int geneAM,TFAM=0;
	char *Name;
	char STRING1[TFScale][10];
	//double originalMA[N][N];
	int unknow=0;
	data.get(ch);
	for(i=0;i<TFScale;)
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
			if(GRN[b][a]==1||GRN[b][a]==2)
			{
				GRN[b][a]=2;
				uncertain1.push_back(b);
				uncertain2.push_back(a);
				unknow++;
			}
			else
				GRN[b][a]=-1;
		}
		else if(ch=='+')
		{
			data.get(ch);
			if(ch=='-')
			{
				GRN[b][a]=2;
				uncertain1.push_back(b);
				uncertain2.push_back(a);
				unknow++;
			}
			else
				GRN[b][a]=1;
		}
		else
			GRN[b][a]=0;
		string noUse;
		getline(data,noUse);
		num=0;
		if(!data.get(ch))
		{
			TFAmount=i;
			i=TFScale;
		}
		//if(TFAM<a)
			//TFAM=a;
		//cout<<j<<endl;
	}
	//memcpy((char *)originalMatrix,(char *)originalMA,sizeof(double)*N*N);
	//originalMatrix=originalMA;
	//cout<<TFAM<<endl;
	geneAmount=TFAmount;
}

//full fill the matrix with 2
Regulation::Regulation()
{
	originalGRN=new double*[GENEAM];
	for(int i=0;i<GENEAM;i++)
		originalGRN[i]=new double[TFScale];
	geneAmount=0;
	openFileError=0;
//	a=originalMatrix[0][0];
	for(int n=0;n<TFScale;n++)
	{
		for(int m=0;m<GENEAM;m++)
		{
			originalGRN[m][n]=0;
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

map<string,string> Regulation::mapTFIM()
{
	map<string,string> dictTFIM;
	ifstream data("TFIM");
	//ofstream text("Sequence");
	char ch;
	int i;
	for(i=1;data.get(ch);i++)
	{
		while(ch=='#')//filter RegulonDB line
		{
			string noUse;
			getline(data,noUse);
			data.get(ch);
		}
		string line;
		char sline[30];
		//char cline[30];
		getline(data,line);
		int j;
		for(j=0;line[j]!='(';j++)
		{
			sline[j]=line[j];
		}
		sline[j]='\0';
		strlwr(sline);
		int p;
		for(p=0;sline[p]!='	';p++)
		{
			//cline[p]=line[p];
		}
		char tempName[20];
		string name;
		for(int q=0;p<j;q++)
		{
			p++;
			tempName[q]=sline[p];
		}
		name=tempName;
		dictTFIM[name]=line;
		//cout<<i<<endl;
		//tempName[0]='\0';
	}
	return dictTFIM;
	//tempName="";
}

void Regulation::readTUPosition()
{
	ifstream data("TUposition");
	char ch;
	for(;data.get(ch);)
	{
		while(ch=='#')//filter RegulonDB line
		{
			string noUse;
			getline(data,noUse);
			data.get(ch);
		}
		string line;
		getline(data,line);
		int i;
		for(i=0;line[i]!='	';i++){}
		i++;
		for(;line[i]!='	';i++){}
		i++;
		char ProName[10];
		int j;
		for(j=0;line[i]!='	';j++)
		{
			ProName[j]=line[i];
			i++;
		}
		i++;
		ProName[j]='\0';
		strlwr(ProName);
		char TUPos[10];
		for(j=0;line[i]!='	';j++)
		{
			TUPos[j]=line[i];
			i++;
		}
		TUPos[j]='\0';
		string tempTUPos=TUPos;
		promoterNameLib.push_back(ProName);
		TUPositon.push_back(atoi(tempTUPos.c_str()));
	}
}

int Regulation::getTFAmount()
{
	return TFAmount;
}

void Regulation::getGenePromoter(TFIM geneFirst[])
{
	int TUAmount;
	TUAmount=TUPositon.size();
	int target;
	for(int i=0;i<getGeneAmount();i++)
	{
		int max=TUAmount,min=0;
		for(int j=(max+min)/2;(max-min)!=1;j=(max+min)/2)
		{
			if(geneFirst[i].getLeftPosition()<=TUPositon[j])
				max=j;
			else
				min=j;
		}
		geneFirst[i].putInPromoterName(promoterNameLib[min]);
	}
}

map<string,string> Regulation::mapPromoter()
{
	map<string,string> dictPromoter;
	ifstream data("Promoters");
	//ofstream text("Sequence");
	char ch;
	for(;data.get(ch);)
	{
		while(ch=='#')//filter RegulonDB line
		{
			string noUse;
			getline(data,noUse);
			data.get(ch);
		}
		string line;
		char nameLine[20];
		//char cline[30];
		getline(data,line);
		int j,i;
		for(j=0;line[j]!='	';j++){}
		j++;
		for(i=0;line[j]!='	';i++)
		{
			nameLine[i]=line[j];
			j++;
		}
		nameLine[i]='\0';
		strlwr(nameLine);
		j++;
		for(;line[j]!='	';j++){}
		j++;
		for(;line[j]!='	';j++){}
		j++;
		for(;line[j]!='	';j++){}
		j++;
		char tempSequence[100];
		for(i=0;line[j]!='	';i++)
		{
			tempSequence[i]=line[j];
			j++;
		}
		tempSequence[i]='\0';
		strupr(tempSequence);
		string name=nameLine;
		dictPromoter[name]=tempSequence;
		//cout<<i<<endl;
		//tempName[0]='\0';
	}
	return dictPromoter;
	//tempName="";
}