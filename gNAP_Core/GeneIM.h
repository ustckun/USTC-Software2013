//class GeneIM
//search imformation about the Gene.

#ifndef __USTC_Software__GeneIM__
#define __USTC_Software__GeneIM__

#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <map>

#include"define.h"

using namespace std;


class GeneIM
{
public:
	GeneIM()
	{
		get_error=0;
		left_position=0;
		right_position=0;
	}
	int gene_number;//get from GetReady
	char *name;//get from GetReady
	void getGeneInformation(map<string,string> dict);//fp is the position of output file which contains gene number and DNA sequence
	void getPromoterIF(FILE *fp , map<string,string> dict);
	string getID();
	string getGeneSequence();
	int getLeftPosition();
	int getRightPosition();
	int getFileError();
	char *getGeneName();
	void putName();
	void putInPromoterName(string promoter);
	int getRNA();
	string getGeneDescription();
	char gene_name[10];
private:
	string iD;
	string gene_sequence;
	//string left_position;
	int left_position;
	int right_position;
	//string right_position;
	string gene_description;
	string promoter_name;
	string promoter_sequence;
	int get_error;
	int RNA;
};

#endif