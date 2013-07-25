//EasytoDebug.cpp
//get gene number. sequence and length of sequence
//

#include"TFIM.h"
#include"Regulation.h"
#include"EasytoDebug.h"

void EasytoDebug::fetch()
{
	ifstream temp("sequence");
	temp>>number>>sequence;
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