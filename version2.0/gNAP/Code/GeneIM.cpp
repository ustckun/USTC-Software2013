////////////////////////////////////////////////////////////////////////////////
/// COPYRIGHT NOTICE\n
/// Distribute under BSD License\n
/// Copyright (c) 2013, iGEM Software Team of University of Science and
/// Technology of China\n
/// All rights reserved.
///
/// \file GeneIM.cpp
/// \brief Statments of funcions of the class GeneIM.
/// \version 1.0
/// \author Wang Chenkun
/// \date September 2nd, 2013
////////////////////////////////////////////////////////////////////////////////

#include"GeneIM.h"


void GeneIM::getGeneInformation(map<string,string> dict)
{
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
    h=0;
    char *p=&fullLine[h];
    for(;fullLine[h]!='\t';h++)
    {
    }
    fullLine[h]='\0';
    iD=p;
    h++;
    p=&fullLine[h];
    for(;fullLine[h]!='(';h++)
    {
    }
    fullLine[h]='\0';
    true_name=p;
    h++;
    for(;fullLine[h]!='\t';h++)
    {
    }
    h++;
    p=&fullLine[h];
    for(;fullLine[h]!='\t';h++)
    {
    }
    fullLine[h]='\0';
    tempLeft=p;
    left_position=atoi(tempLeft.c_str());
    h++;
    p=&fullLine[h];
    for(;fullLine[h]!='\t';h++)
    {
    }
    fullLine[h]='\0';
    tempRight=p;
    right_position=atoi(tempRight.c_str());
    h++;
    for(;fullLine[h]!='\t';h++)
    {
    }
    h++;
    for(;fullLine[h]!='\t';h++)
    {
    }
    h++;
    p=&fullLine[h];
    for(;fullLine[h]!='\t';h++)
    {
    }
    fullLine[h]='\0';
    gene_description=p;
    h++;
    for(;fullLine[h]!='\t';h++)
    {
    }
    h++;
    for(;fullLine[h]!='\t';h++)
    {
    }
    h++;
    p=&fullLine[h];
    for(;fullLine[h]!='\0';h++)
    {
    }
    fullLine[h]='\0';
    gene_sequence=p;
	}
	else
		cout<<gene_name<<endl;
}

void GeneIM::putName()
{
	for(int i=0;i<10;i++)
	{
		gene_name[i]=name[i];
	}
}

string GeneIM::getID()
{
	return iD;
}

int GeneIM::getLeftPosition()
{
	return left_position;
}

int GeneIM::getRightPosition()
{
	return right_position;
}

string GeneIM::getGeneSequence()
{
	return gene_sequence;
}

char *GeneIM::getGeneName()
{
	return gene_name;
}

string GeneIM::getGeneTrueName()
{
    return true_name;
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

void GeneIM::getPromoterIF(map<string,string>dict)
{
	if(dict.find(promoter_name)!=dict.end())
	{
		promoter_sequence=dict[promoter_name];
	}
}

string GeneIM::getPromoterName()
{
    return promoter_name;
}

string GeneIM::getPromoterSequence()
{
    return promoter_sequence;
}
