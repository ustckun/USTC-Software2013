#include"TFIM.h"
#include"Regulation.h"
#define aM 100

int easyFind(char *a,char *b)
{
	int p=0,q;
	int c=strlen(b);
	while(p<20)
	{
		if(a[p]==b[0])
		{
			for(q=1;q<c;)
			{
				if(a[p+q]==b[q])
					q++;
				else
					q=100;
			}
			if(q==c)
				return 1;
		}
		p++;
	}
	return -1;
}

void TFIM::getGeneImformation()
{
	ifstream data("TFIM.txt");
	if(!data)
	{
		getError=1;
	}
	char ch;
	int i;
	for(i=1;i<aM;)
	{
		if(!data.get(ch))
		{
			i=aM;
		}
		while(ch=='#')//filter RegulonDB line
		{
			string noUse;
			getline(data,noUse);
			data.get(ch);
		}
		string line;
		char GN[10];
		getline(data,line);
		for(int h=0;h<10;h++)
			GN[h]=name[h];
		char cLine[60];
		for(int H=0;H<60;H++)
			cLine[H]=line[H];
		strlwr(cLine);
		strlwr(GN);
		int haveFound=easyFind(cLine,GN);
		if(haveFound==1)
		{
			char fullLine[10000];
			int h=0;
			while(line[h]!='\n')
			{
				fullLine[h]=line[h];
				h++;
			}
			fullLine[h]='\0';
			const char *delims="	";
			char *p;
			p=strtok(fullLine,delims);
			iD=p;
			p=strtok(NULL,delims);
			p=strtok(NULL,delims);
			leftPosition=p;
			p=strtok(NULL,delims);
			rightPosition=p;
			p=strtok(NULL,delims);
			p=strtok(NULL,delims);
			p=strtok(NULL,delims);
			p=strtok(NULL,delims);
			p=strtok(NULL,delims);
			p=strtok(NULL,delims);
			geneSequence=p;
			i=aM;
		}
	}
}

char *TFIM::getID()
{
	return iD;
}

char *TFIM::getLeftPosition()
{
	return leftPosition;
}

char *TFIM::getRightPosition()
{
	return rightPosition;
}

char *TFIM::getGeneSequence()
{
	return geneSequence;
}

int TFIM::getFileError()
{
	return getError;
}