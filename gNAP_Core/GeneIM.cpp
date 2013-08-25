// GeneIM.cpp
// use name to find gene information in database
//
//version 2.0
//change:
//1 output gene's sequence to txt, improve our effecient of debug
//2 add notation
//
//***********************************************************************************************************
//GetReady::readName has get all the gene's name from TF-TF regulation
//search all the Genes information in about 4400 genes information
//it contains gene numbers, gene left position, gene right positon, gene ID in regulonDB and gene sequence
//***********************************************************************************************************
//
//interface:
//char* getGeneSequence()
//char* getGeneName()
//int gene_number
//

#include"GeneIM.h"

//#include"GetReady.h"

//a function to find string b in string a

//get 166 genes information in database
void GeneIM::getGeneInformation(map<string,string> dict)
{
	//ifstream data("GeneIM");
	//ofstream text("Sequence");
	string allLine;
	string tempLeft;
	string tempRight;
	if(dict.find(gene_name)!=dict.end())
	{
		allLine=dict[gene_name];
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
	tempLeft=p;
	left_position=atoi(tempLeft.c_str());
	p=strtok(NULL,delims);
	tempRight=p;
	right_position=atoi(tempRight.c_str());
	p=strtok(NULL,delims);
	p=strtok(NULL,delims);
	p=strtok(NULL,delims);
	gene_description=p;
	p=strtok(NULL,delims);
	p=strtok(NULL,delims);
	p=strtok(NULL,delims);
	gene_sequence=p;
	//i=aM+1;
	//text.seekp(text.end);
	//text<<gene_number<<"	"<<gene_sequence<<endl;
	//fprintf(fp,"%s	%s\n",gene_name,gene_sequence);
	//cout<<i<<endl;
	//if(i==aM)
	//	cout<<"can't find gene sequence"<<endl;
	}
	else
		cout<<gene_name<<endl;
}

//char *name has memory leaks so put name to private gene_name
void GeneIM::putName()
{
	for(int i=0;i<10;i++)
	{
		gene_name[i]=name[i];
	}
}

//get Gene ID of regulonDB
string GeneIM::getID()
{
	return iD;
}

//get gene left end position in genome
int GeneIM::getLeftPosition()
{
	return left_position;
}

//get gene right end position in genome
int GeneIM::getRightPosition()
{
	return right_position;
}

//get Gene sequence in database
string GeneIM::getGeneSequence()
{
	return gene_sequence;
}

//if file not open it will return 1
int GeneIM::getFileError()
{
	return get_error;
}

//because gene_name is put into private this function return gene name
char *GeneIM::getGeneName()
{
	return gene_name;
}

int GeneIM::getRNA()
{
	cout<<RNA<<endl;
	return RNA;
}

string GeneIM::getGeneDescription()
{
	return gene_description;
}

void GeneIM::putInPromoterName(string promoter)
{
	promoter_name=promoter;
}

void GeneIM::getPromoterIF(FILE *fp,map<string,string>dict)
{
	if(dict.find(promoter_name)!=dict.end())
	{
		promoter_sequence=dict[promoter_name];
		//DATA<<gene_number<<"	"<<promoter_sequence<<"	"<<gene_sequence<<endl;
		fprintf(fp,"%d	%s	%s\n",gene_number,promoter_sequence.c_str(),gene_sequence.c_str());
	}
	else
		cout<<promoter_name<<endl;
}