//class TFIM
//search imformation about the Gene.

#include <ctime>
#include <cmath>
#include <iomanip>
#include <cstdlib>
#include<iostream>
#include<fstream>
#include<string>
#include<vector>
#include<algorithm>
#include<map>
#define NN 100                                   //count number
#define M 100										//PSO number
#define PETS 128                                
#define STEP (1.0/PETS)                           //step length
#define MAXTIME 100                            //maxtime
#define INITIALVALUE 2.5                          //initial value
#define DIMENS 170    
#define N 170
#define scale 170
#define gap 0
#define aM 4700
using namespace std;

class TFIM
{
public:
	TFIM()
	{
		getError=0;
	}
	int geneNumber;//get from Regulation
	char *name;//get from Regulation
	void getGeneInformation(FILE *fp,map<string,string> dict);//fp is the position of output file which contains gene number and DNA sequence
	string getID();
	string getGeneSequence();
	string getLeftPosition();
	string getRightPosition();
	int getFileError();
	char *getGeneName();
	void putName();
	//int getRNA();
	string getGeneDescription();
private:
	string iD;
	string geneSequence;
	string leftPosition;
	string rightPosition;
	string geneDescription;
	int getError;
	char geneName[10];
	//int RNA;
};