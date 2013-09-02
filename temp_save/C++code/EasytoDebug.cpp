//EasytoDebug.cpp
//get gene number. sequence and length of sequence
//

#include"TFIM.h"
#include"Regulation.h"
#include"EasytoDebug.h"

void EasytoDebug::fetch(FILE *temp)
{
	//ifstream temp("Sequence");
	char *s=new char;
	fscanf(temp,"%d	%s",&number,s);
	sequence=s;
	length=sequence.length();
}

int EasytoDebug::getLength()
{
	return length;
}

int EasytoDebug::getNumber()
{
	return number;
}

string EasytoDebug::getSequence()
{
	return sequence;
}