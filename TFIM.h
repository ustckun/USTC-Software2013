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
	void getGeneImformation();
	char *getID();
	char *getGeneSequence();
	char *getLeftPosition();
	char *getRightPosition();
	int getFileError();
private:
	char *iD;
	char *geneSequence;
	char *leftPosition;
	char *rightPosition;
	int getError;
};