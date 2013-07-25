
#include"TFIM.h"
#include"Regulation.h"
#include"EasytoDebug.h"

int main()
{
	TFIM test[200];
	Regulation Hope;
	Hope.readName(test);
	Hope.fullFill();
	cout<<Hope.getGeneAmount()<<endl;
	cout<<Hope.originalMatrix[11][11]<<endl;
	test[100].getGeneInformation();
	//cout<<i<<endl;
	cout<<test[100].getGeneName()<<endl;
	cout<<test[100].getLeftPosition()<<endl;
	cout<<test[100].getRightPosition()<<endl;
	cout<<test[100].getGeneSequence()<<endl;
	EasytoDebug test1[200];
	test1[100].fetch();
	cout<<test1[100].getLength()<<endl;
	cout<<test1[100].getNumber()<<endl;
	cout<<test1[100].getSequence()<<endl;
	getchar();
}
