
#include"TFIM.h"
#include"Regulation.h"
#include"EasytoDebug.h"

int main()
{
	TFIM test[N];
	Regulation Hope;
	Hope.readName(test);
	Hope.fullFill();
	cout<<Hope.getGeneAmount()<<endl;
	cout<<Hope.originalMatrix[11][11]<<endl;
	FILE *fp;
	fp=fopen("sequence","w");
	for(int i=0;i<Hope.getGeneAmount();i++)
	{
		test[i].getGeneInformation(fp);
		cout<<i<<endl;
		//cout<<test[i].getGeneSequence()<<endl;
	}
	fclose(fp);
	//rewind(fp);
	//cout<<i<<endl;
	//cout<<test[164].getGeneName()<<endl;
	//cout<<test[164].getLeftPosition()<<endl;
	//cout<<test[164].getRightPosition()<<endl;
	//cout<<test[164].getGeneSequence()<<endl;
	/*EasytoDebug test1[200];
	test1[100].fetch();
	cout<<test1[100].getLength()<<endl;
	cout<<test1[100].getNumber()<<endl;
	cout<<test1[100].getSequence()<<endl;*/
	//EasytoDebug read[100];
	//fp=fopen("Sequence","r");
	/*for(int j=0;j<1;j++)
	{
		read[j].fetch(fp);
		//cout<<read[j].getNumber()<<endl;
		//cout<<read[j].getSequence()<<endl;
	}*/
	getchar();
}
