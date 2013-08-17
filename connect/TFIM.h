//class TFIM
//search imformation about the Gene.

#include <ctime>
#include <cmath>
#include <iomanip>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <map>

#define NN 100                                   //count number
#define M 30										//PSO number
#define PETS 128                                
#define STEP (1.0/PETS)                           //step length
#define MAXTIME 100                            //maxtime
#define INITIALVALUE 2.5                          //initial value
#define DIMENS 100    
#define N 100
#define scale 100
#define gap 0
#define aM 4700
#define TFScale 220
#define GENEAM 1800
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
	void getGeneInformation(map<string,string> dict);//fp is the position of output file which contains gene number and DNA sequence
	void getPromoterIF(FILE *fp,map<string,string>dict);
	string getID();
	string getGeneSequence();
	int getLeftPosition();
	int getRightPosition();
	int getFileError();
	char *getGeneName();
	void putName();
	void putInPromoterName(string promoter);
	//int getRNA();
	string getGeneDescription();
	char geneName[10];
private:
	string iD;
	string geneSequence;
	//string leftPosition;
	int leftPosition;
	int rightPosition;
	//string rightPosition;
	string geneDescription;
	string promoterName;
	string promoterSequence;
	int getError;
	
	//int RNA;
};