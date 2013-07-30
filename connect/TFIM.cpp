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
#include"EasytoDebug.h"

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
					return 1;
		}
		p++;
	}
	return -1;
}

//get 166 genes information in database
void TFIM::getGeneInformation(FILE *fp)
{
	ifstream data("TFIM");
	//ofstream text("Sequence");
	if(!data)
	{
		getError=1;
	}
	else
		getError=0;
	char ch;
	int i;
	for(i=1;i<aM;i++)
	{
		if(!data.get(ch))
		{
			i=aM;
		}
		while(ch=='#')//filter RegulonDB line
		{
			string noUse;
			getline(data,noUse);
			data.get(ch);
		}
		string line;
		char GN[10];
		getline(data,line);
		for(int h=0;h<10;h++)
			GN[h]=geneName[h];
		char cLine[60];
		for(int H=0;H<60;H++)
			cLine[H]=line[H];
		strlwr(cLine);
		strlwr(GN);
		int haveFound=easyFind(cLine,GN);
		if(haveFound==1)
		{
			char fullLine[10000];
			int h=0;
			while(line[h]!='\0')
			{
				fullLine[h]=line[h];
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
			if(*p!='-')RNA=1;
			else RNA=0;
			p=strtok(NULL,delims);
			p=strtok(NULL,delims);
			p=strtok(NULL,delims);
			p=strtok(NULL,delims);
			geneSequence=p;
			i=aM+1;
			//text.seekp(text.end);
			//text<<geneNumber<<"	"<<geneSequence<<endl;
			fprintf(fp,"%d	%s\n",geneNumber,geneSequence);
		}
		//cout<<i<<endl;
	}
	if(i==aM)
		cout<<"can't find gene sequence"<<endl;
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
char *TFIM::getID()
{
	return iD;
}

//get gene left end position in genome
char *TFIM::getLeftPosition()
{
	return leftPosition;
}

//get gene right end position in genome
char *TFIM::getRightPosition()
{
	return rightPosition;
}

//get Gene sequence in database
char *TFIM::getGeneSequence()
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

int TFIM::getRNA()
{
	cout<<RNA<<endl;
	return RNA;
}