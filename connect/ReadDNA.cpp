//ReadDNA.cpp
//get gene number. sequence and length of sequence
//

#include"TFIM.h"
//#include"calculation.h"
#include"Calculate.h"
#include"Sequence.h"
#include"GRN.h"
#include"Regulation.h"
#include"ReadDNA.h"
#include"PSOPredict.h"
#include"SBOL.h"


void ReadDNA::fetch(FILE *temp)
{
	//ifstream temp("Sequence");
	char s[2700];
	fscanf(temp,"%d	%s",&number,s);
	sequence=s;
	length=sequence.length();
}

int ReadDNA::getLength()
{
	return length;
}

int ReadDNA::getNumber()
{
	return number;
}

string ReadDNA::getSequence()
{
	return sequence;
}