//class TFIM
//search imformation about the Gene.


#include<iostream>
#include<fstream>
#include<string>
#include<vector>
#include<algorithm>
using namespace std;

class TFIM
{
public:
	int geneNumber;//get from Regulation
	char *name;//get from Regulation
	void getGeneInformation();
	char *getID();
	char *getGeneSequence();
	char *getLeftPosition();
	char *getRightPosition();
	int getFileError();
	char *getGeneName();
	void putName();
private:
	char *iD;
	char *geneSequence;
	char *leftPosition;
	char *rightPosition;
	int getError;
	char geneName[10];
};