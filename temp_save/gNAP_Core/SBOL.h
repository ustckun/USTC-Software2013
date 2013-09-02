#ifndef __USTC_Software__SBOL__
#define __USTC_Software__SBOL__

#include <fstream>
#include <iostream>


class SBOL
{
public:
	void CreatSBOL(string position,GeneIM IM);
private:
	string head;
	string Combine(string title,string detail);
	string FormartStart(string a);
	string FormartEnd(string b);
};


#endif