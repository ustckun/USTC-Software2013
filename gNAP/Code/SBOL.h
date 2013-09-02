#ifndef __USTC_Software__SBOL__
#define __USTC_Software__SBOL__

#include <fstream>
#include <iostream>
#include <string>

using namespace std;


class SBOL
{
public:
    void CreatSBOL(string gene_name,string ID,string left,string right,string description,string seq);
private:
    string head;
    string Combine(string title,string detail);
    string FormartStart(string a);
    string FormartEnd(string b);
};


#endif
