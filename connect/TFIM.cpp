// TFIM.cpp
// use name to find gene information in database
//
//version 2.0
//change:
//1 output gene's sequence to txt, improve our effecient of debug
//2 add notation
//
//***********************************************************************************************************
//Regulation::readName has get all the gene's name from TF-TF regulation
//search all the Genes information in about 4400 genes information
//it contains gene numbers, gene left position, gene right positon, gene ID in regulonDB and gene sequence
//***********************************************************************************************************
//
//interface:
//char* getGeneSequence()
//char* getGeneName()
//int geneNumber
//

#include"TFIM.h"
#include"Calculate.h"
#include"Sequence.h"
#include"GRN.h"
#include"Regulation.h"
#include"ReadDNA.h"
#include"PSOPredict.h"
#include"SBOL.h"

//#include"Regulation.h"
#define aM 5000

//a function to find string b in string a
int easyFind(char *a,char *b)
{
	int p=0,q=0;
	int c=strlen(b);
	while(p<60)
	{
		if(a[p]==b[0])
		{
			for(q=1;q<c;)
			{
				if(a[p+q]==b[q])
					q++;
				else
					q=100;
			}
			if(q==c)
				if(a[p-1]=='\t')
					if(a[p+c]=='(')
					return 1;
		}
		p++;
	}
	return -1;
}

//get 166 genes information in database
void TFIM::getGeneInformation(FILE *fp,map<string,string> dict)
{
	ifstream data("TFIM");
	//ofstream text("Sequence");
	string allLine;
	if(dict.find(geneName)!=dict.end())
	{
		allLine=dict[geneName];
	char fullLine[10000];
	int h=0;
	while(allLine[h]!='\0')
	{
		fullLine[h]=allLine[h];
		h++;
	}
	fullLine[h]='\0';
	const char *delims="	";
	char *p;
	p=strtok(fullLine,delims);
	iD=p;
	p=strtok(NULL,delims);
	p=strtok(NULL,delims);
	leftPosition=p;
	p=strtok(NULL,delims);
	rightPosition=p;
	p=strtok(NULL,delims);
	p=strtok(NULL,delims);
	p=strtok(NULL,delims);
	geneDescription=p;
	p=strtok(NULL,delims);
	p=strtok(NULL,delims);
	p=strtok(NULL,delims);
	geneSequence=p;
	//i=aM+1;
	//text.seekp(text.end);
	//text<<geneNumber<<"	"<<geneSequence<<endl;
	fprintf(fp,"%s	%s\n",geneName,geneSequence);}
	//cout<<i<<endl;
	//if(i==aM)
	//	cout<<"can't find gene sequence"<<endl;
	else
		cout<<geneName<<endl;
}

//char *name has memory leaks so put name to private geneName
void TFIM::putName()
{
	for(int i=0;i<10;i++)
	{
		geneName[i]=name[i];
	}
}

//get Gene ID of regulonDB
string TFIM::getID()
{
	return iD;
}

//get gene left end position in genome
string TFIM::getLeftPosition()
{
	return leftPosition;
}

//get gene right end position in genome
string TFIM::getRightPosition()
{
	return rightPosition;
}

//get Gene sequence in database
string TFIM::getGeneSequence()
{
	return geneSequence;
}

//if file not open it will return 1
int TFIM::getFileError()
{
	return getError;
}

//because geneName is put into private this function return gene name
char *TFIM::getGeneName()
{
	return geneName;
}

/*int TFIM::getRNA()
{
	cout<<RNA<<endl;
	return RNA;
}*/

string TFIM::getGeneDescription()
{
	return geneDescription;
}