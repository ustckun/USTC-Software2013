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
//#include <iomanip>
//#include <fstream>
//#include <cmath>
//#include <cstdlib>
//#include <ctime>
#define NN 1                                   //count number
#define PETS 128                                
#define STEP (1.0/PETS)                           //step length
#define MAXTIME 100                            //maxtime
#define INITIALVALUE 2.5                          //initial value
#define DIMENS 170    
#define N 170
#define scale 170
using namespace std;

class TFIM
{
public:
	int geneNumber;//get from Regulation
	char *name;//get from Regulation
	void getGeneInformation(FILE *fp);
	char *getID();
	char *getGeneSequence();
	char *getLeftPosition();
	char *getRightPosition();
	int getFileError();
	char *getGeneName();
	void putName();
	int getRNA();
private:
	char *iD;
	char *geneSequence;
	char *leftPosition;
	char *rightPosition;
	int getError;
	char geneName[10];
	int RNA;
};