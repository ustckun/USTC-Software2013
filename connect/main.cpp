#include"TFIM.h"
//#include"calculation.h"

#include"Regulation.h"
#include"EasytoDebug.h"
#include"Sequence.h"
//#include"GRN.h"

int main()
{
	TFIM test[N];
	Regulation Hope;
	Hope.readName(test);
	Hope.fullFill();
	int A,C;
	string B;
	FILE *fp;
	//FILE *fp1;
	//fp=fopen("DNA","w");
	/*for(int i=160;i<Hope.getGeneAmount();i++)
	{
		test[i].getGeneInformation(fp);
		cout<<i<<endl;
		//cout<<test[i].getGeneSequence()<<endl;
	}
	fclose(fp);*/
	//EasytoDebug read[100];
	EasytoDebug luck[N];
	Sequence good[N];
	fp=fopen("DNA","r");
	//fp1=fopen("AAS","w");
	ofstream data("AAS");
	for(int i=0;i<Hope.getGeneAmount();i++)
	{
		//test[i].getGeneInformation(fp);
		luck[i].fetch(fp);
		A=luck[i].getNumber();
		B=luck[i].getSequence();
		C=luck[i].getLength();
		//B=test[i].getGeneSequence();
		//C=B.length();
		good[i].initializeGeneSequence(B,A,C);
		B=good[i].getAminoAcidSequence();
		data<<A<<"	"<<B<<endl;
	}
	//test[30].getGeneInformation();
	//test[31].getGeneInformation();
	//int a;
	//a=test[31].getRNA();
	//a=test[30].getRNA();
	getchar();
}