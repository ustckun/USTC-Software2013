#include"TFIM.h"
#include"Calculate.h"
#include"Regulation.h"
#include"ReadDNA.h"
#include"Sequence.h"
#include"GRN.h"
#include"PSOPredict.h"


int main()
{
	srand((unsigned)time(0));
	double matrix[170][170];
	ReadRegulation(matrix);
	//double reg[170][170];	
	Calculate Hope;
	//Hope.RandMatrix(matrix,reg,166);
	Hope.Network_1(matrix,166);
	//for(int i=0;i<170;i++)cout<<Hope.consistence[i]<<endl;
	/*TFIM test[N];
	Regulation Hope;
	Hope.readName(test);
	Hope.fullFill();
	Calculate xinyu;
	//xinyu.Network_1(Hope.originalMatrix,166);
	//xinyu.Network_2(Hope.originalMatrix,166);
	//cout<<xinyu.nong[90];
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
	fclose(fp);
	//ReadDNA read[100];
	ReadDNA luck[N];
	Sequence good[N];
	fp=fopen("DNA","r");
	//fp1=fopen("AAS","w");
	//ofstream data("AAS");
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
		good[i].translation();
		B=good[i].aminoAcidSequence;
		//cout<<A<<"	"<<B<<endl;
	}
	fclose(fp);
	//data.close();
	ofstream RN("newGRN");
	string insertGene="MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTFAYGVQCFSRYPDHMKRHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSVLSKDPNEKRDHMVLLEFVTAAGITHGMDELYK";
	GRN newGRN;
	newGRN.initializeGRN(Hope.originalMatrix,Hope.getGeneAmount());
	good[166].getNewASS(insertGene);
	newGRN.constructNewGRN(good);
	for(int j=0;j<N;j++)
	{
		for(int i=0;i<N;i++)
		{
			RN<<newGRN.newGRNCorrelation[j][i]<<"	";
		}
		RN<<endl;
	}
	//test[30].getGeneInformation();
	//test[31].getGeneInformation();
	//int a;
	//a=test[31].getRNA();
	//a=test[30].getRNA();*/
	getchar();
}