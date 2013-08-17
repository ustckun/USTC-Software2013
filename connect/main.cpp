#include"TFIM.h"
#include"Calculate.h"
#include"Regulation.h"
#include"ReadDNA.h"
#include"Sequence.h"
#include"GRN.h"
#include"PSOPredict.h"
#include"ReadRegulation.h"
#include"SBOL.h"

int main()
{
	//srand((unsigned)time(0));
	//double matrix[170][170];
	//ReadRegulation(matrix);
	//double reg[170][170];	
	//Calculate Hope;
	//Hope.RandMatrix(matrix,reg,166);
	//Hope.Network_1(matrix,166);
	//for(int i=0;i<170;i++)cout<<Hope.consistence[i]<<endl;
	TFIM test[GENEAM];
	Regulation luck;
	//SBOL happy;
	luck.readTUPosition();
	/*for(int i=0;i<luck.TUPositon.size();i++)
	{
		cout<<luck.promoterNameLib[i]<<"	"<<luck.TUPositon[i]<<endl;
	}*/
	luck.getRegulationMatrix(test);
	//luck.fullFill();
	//ofstream dna("Old_GRN");
	//ofstream DATA("NO&PRO&GENE");
	FILE *fp;
	fp=fopen("NO&PRE&DNA","w");
	map<string,string> dict;
	dict=luck.mapTFIM();
	/*for(int i=0;i<GENEAM;i++)
	{
		for(int j=0;j<TFScale;j++)
		{
			dna<<luck.originalGRN[i][j]<<"	";
		}
		dna<<endl;
	}*/
	for(int i=0;i<luck.getGeneAmount();i++)
	{
		test[i].getGeneInformation(dict);
		cout<<i<<endl;
	}
	//dict.erase(dict.begin(),dict.end());
	luck.readTUPosition();
	luck.getGenePromoter(test);
	dict=luck.mapPromoter();
	for(int i=0;i<luck.getGeneAmount();i++)
		test[i].getPromoterIF(fp,dict);
	getchar();
	/*for(int i=0;i<luck.uncertain1.size();i++)
	{
		cout<<luck.uncertain1[i]<<"	";
		cout<<luck.uncertain2[i]<<endl;
	}
	cout<<luck.getGeneAmount()<<endl;
	getchar();*/
	/*map<string,string> good;
	good=luck.mapTFIM();
	FILE *fp;
	fp=fopen("DNA","w");
	for(int i=0;i<luck.getGeneAmount();i++)
	{
		test[i].getGeneInformation(fp,good);
		//happy.CreatSBOL("c:\\Users\\kun\\Desktop\\igem\\Code\\USTC-Software\\All\\All\\All\\test\\",test[i]);
	//	cout<<test[i].getLeftPosition();
		//cout<<i<<endl;
	}
	cout << (double)clock()/CLOCKS_PER_SEC << endl;
	getchar();*/
	//PSOPredict happy;
	//Calculate cal;
	/*for(int i=0;i<N;i++)
	{
		happy.randomMax[i]=0;
		happy.randomMin[i]=100;
	}
	for(int i=0;i<100;i++)
		happy.getRange(luck.originalMatrix,luck.getGeneAmount(),cal);
	ofstream file("Range");
	for(int i=0;i<luck.getGeneAmount();i++)
	{
		file<<happy.randomMin[i]<<"	"<<happy.randomMax[i]<<endl;
	}*/
	/*happy.targetGene[156]=5.6;
	happy.algorithm_PSO(luck.originalMatrix,happy.targetGene,luck.getGeneAmount(),cal);
	for(int i=0;i<luck.getGeneAmount();i++)
	{
		cout<<happy.edPick[i]<<endl;
		cout<<happy.toPick[i]<<endl;
	}*/
	/*ofstream file("Old_GRN");
	for(int i=0;i<N;i++)
	{
		for(int j=0;j<N;j++)
		{
			file<<luck.originalMatrix[i][j]<<"	";
		}
		file<<endl;
	}*/
	//Hope.fullFill();
	//Calculate xinyu;
	//xinyu.Network_1(Hope.originalMatrix,166);
	//xinyu.Network_2(Hope.originalMatrix,166);
	//cout<<xinyu.nong[90];
	//int A,C;
	//string B;
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
	//getchar();
}