//#include"stdio.h"
#include"TFIM.h"
#include"Regulation.h"

void main()
{
	TFIM test1[5];
	Regulation test;
	test.readName(test1);
	int a=test.getGeneAmount();
	cout<<"gene amount:"<<a<<endl;
	int i;
	for(i=0;i<5*5;i++)
	{
		cout<<"regulation:"<<test.originalMatrix<<endl;
		test.originalMatrix=test.originalMatrix+1;
	}
	getchar();
}