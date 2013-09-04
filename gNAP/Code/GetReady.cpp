
#include"GeneIM.h"
#include"GetReady.h"


void GetReady::getRegulationMatrix(GeneIM temp_gene_IM[],string TF_TF_address,string TF_Gene_address)
{
    readTFTF(temp_gene_IM,originalGRN,TF_TF_address);
    addTF(temp_gene_IM,TF_Gene_address);
	TF_amount=gene_amount;
    readTFGene(temp_gene_IM,originalGRN,TF_Gene_address);
}

void GetReady::addTF(GeneIM temp_gene_IM[],string TF_Gene_address)
{
    ifstream data(TF_Gene_address.c_str());
	if(!data)
	{
		open_file_error=1;
	}
	else
		open_file_error=0;
	char ch;
	int i=gene_amount;
	data.get(ch);
	char STRING1[TFScale][10];
	char *Name;
	for(;i<TFScale;)
	{
		int num=0;
		while(ch=='#')
		{
			string noUse;
			getline(data,noUse);
			data.get(ch);
		}
		while(ch!='	')
		{
			STRING1[i][num]=ch;
			data.get(ch);
			num++;
		}
		STRING1[i][num]='\0';
		Name=STRING1[i];
		int j;
		for(j=0;j<i;j++)
		{
            strlwr(Name);
			if(strcmp(temp_gene_IM[j].gene_name,Name)==0)
			{
				j=i+1;
			}
		}
		if(j==i)
		{
			temp_gene_IM[i].name=Name;
			temp_gene_IM[i].putName();
			temp_gene_IM[i].gene_number=i;
			i++;
			gene_amount++;
		}
		string noUse;
		getline(data,noUse);
		num=0;
		if(!data.get(ch))
		{
			i=TFScale;
		}
	}
}
void GetReady::readTFGene(GeneIM temp_gene_IM[],double **old_GRN,string TF_Gene_address)
{
    ifstream data(TF_Gene_address.c_str());
	if(!data)
	{
		open_file_error=1;
	}
	else
		open_file_error=0;
	char ch;
	int i,num=0;
	char *Name;
	char STRING1[GENEAM][10];
	data.get(ch);
	for(i=gene_amount;i<GENEAM;)
	{
		int a,b;
		while(ch=='#')
		{
			string noUse;
			getline(data,noUse);
			data.get(ch);
		}
		while(ch!='	')
		{
			STRING1[i][num]=ch;
			data.get(ch);
			num++;
		}
		STRING1[i][num]='\0';
		Name=STRING1[i];
		int j;
		for(j=0;j<i;j++)
		{
            strlwr(Name);
			if(strcmp(temp_gene_IM[j].gene_name,Name)==0)
			{
				a=j;
				j=i+1;
			}
		}
		if(j==i)
		{
			temp_gene_IM[i].name=Name;
			temp_gene_IM[i].putName();
			temp_gene_IM[i].gene_number=i;
			a=i;
			i++;
		}
		data.get(ch);
		num=0;
		while(ch!='	')
		{
			STRING1[i][num]=ch;
			data.get(ch);
			num++;
		}
		STRING1[i][num]='\0';
		Name=STRING1[i];
		for(j=0;j<i;j++)
		{
            strlwr(Name);
			if(strcmp(temp_gene_IM[j].gene_name,Name)==0)
			{
				b=j;
				j=i+1;
			}
		}
		if(j==i)
		{
			temp_gene_IM[i].name=Name;
			temp_gene_IM[i].putName();
			temp_gene_IM[i].gene_number=i;
			b=i;
			i++;
		}
		data.get(ch);
		if(ch=='-')
		{
			if(old_GRN[b][a]==1||old_GRN[b][a]==2)
			{
                if(old_GRN[b][a]!=2)
                    uncertain<<"?\t"<<temp_gene_IM[b].getGeneName()<<"->"<<temp_gene_IM[a].getGeneName()<<endl;
                old_GRN[b][a]=2;
				uncertain_row.push_back(b);
				uncertain_column.push_back(a);
				unknow++;
			}
			else
				old_GRN[b][a]=-1;
		}
		else if(ch=='+')
		{
			data.get(ch);
			if(ch=='-')
			{
				old_GRN[b][a]=2;
				uncertain_row.push_back(b);
				uncertain_column.push_back(a);
				uncertain<<"?\t"<<temp_gene_IM[b].getGeneName()<<"->"<<temp_gene_IM[a].getGeneName()<<endl;
				unknow++;
			}
			else
				old_GRN[b][a]=1;
		}
		else
			old_GRN[b][a]=0;
		string noUse;
		getline(data,noUse);
		num=0;
		if(!data.get(ch))
		{
			gene_amount=i;
			i=GENEAM;
		}
	}
}

void GetReady::readTFTF(GeneIM temp_gene_IM[],double **old_GRN,string TF_TF_address)
{
    ifstream data(TF_TF_address.c_str());
	if(!data)
	{
		open_file_error=1;
	}
	else
		open_file_error=0;
	char ch;
	int i,num=0;
	char *Name;
	char STRING1[TFScale][10];
	data.get(ch);
	for(i=0;i<TFScale;)
	{
		int a,b;
		while(ch=='#')
		{
			string noUse;
			getline(data,noUse);
			data.get(ch);
		}
		while(ch!='	')
		{
			STRING1[i][num]=ch;
			data.get(ch);
			num++;
		}
		STRING1[i][num]='\0';
		Name=STRING1[i];
		int j;
		for(j=0;j<i;j++)
		{
            strlwr(temp_gene_IM[j].name);
            strlwr(Name);
			if(strcmp(temp_gene_IM[j].name,Name)==0)
			{
				a=j;
				j=i+1;
			}
		}
		if(j==i)
		{
			temp_gene_IM[i].name=Name;
			temp_gene_IM[i].putName();
			temp_gene_IM[i].gene_number=i;
			a=i;
			i++;
		}
		data.get(ch);
		num=0;
		while(ch!='	')
		{
			STRING1[i][num]=ch;
			data.get(ch);
			num++;
		}
		STRING1[i][num]='\0';
		Name=STRING1[i];
		for(j=0;j<i;j++)
		{
            strlwr(temp_gene_IM[j].name);
            strlwr(Name);
			if(strcmp(temp_gene_IM[j].name,Name)==0)
			{
				b=j;
				j=i+1;
			}
		}
		if(j==i)
		{
			temp_gene_IM[i].name=Name;
			temp_gene_IM[i].putName();
			temp_gene_IM[i].gene_number=i;
			b=i;
			i++;
		}
		data.get(ch);
		if(ch=='-')
		{
			if(old_GRN[b][a]==1||old_GRN[b][a]==2)
			{
                if(old_GRN[b][a]!=2)
                    uncertain<<"?\t"<<temp_gene_IM[b].getGeneName()<<"->"<<temp_gene_IM[a].getGeneName()<<endl;
                old_GRN[b][a]=2;
				uncertain_row.push_back(b);
				uncertain_column.push_back(a);
				unknow++;
			}
			else
				old_GRN[b][a]=-1;
		}
		else if(ch=='+')
		{
			data.get(ch);
			if(ch=='-')
			{
				old_GRN[b][a]=2;
				uncertain_row.push_back(b);
				uncertain_column.push_back(a);
				uncertain<<"?\t"<<temp_gene_IM[b].getGeneName()<<"->"<<temp_gene_IM[a].getGeneName()<<endl;
				unknow++;
			}
			else
				old_GRN[b][a]=1;
		}
		else
			old_GRN[b][a]=0;
		string noUse;
		getline(data,noUse);
		num=0;
		if(!data.get(ch))
		{
			TF_amount=i;
			i=TFScale;
		}
	}
	gene_amount=TF_amount;
}

GetReady::GetReady()
{
	originalGRN=new double*[GENEAM];
	for(int i=0;i<GENEAM;i++)
		originalGRN[i]=new double[TFScale];
	gene_amount=0;
	open_file_error=0;
    unknow=0;
    uncertain.open("uncertain_database");
	for(int n=0;n<TFScale;n++)
	{
		for(int m=0;m<GENEAM;m++)
		{
			originalGRN[m][n]=0;
		}
	}
}

int GetReady::getGeneAmount()
{
	return gene_amount;
}

int GetReady::getOpenError()
{
	return open_file_error;
}

map<string,string> GetReady::mapTFIM(string Gene_IM_address)
{
	map<string,string> dictTFIM;
    ifstream data(Gene_IM_address.c_str());
	char ch;
	int i;
	for(i=1;data.get(ch);i++)
	{
		while(ch=='#')
		{
			string noUse;
			getline(data,noUse);
			data.get(ch);
		}
		string line;
		char sline[30];
		getline(data,line);
		int j;
		for(j=0;line[j]!='(';j++)
		{
			sline[j]=line[j];
		}
		sline[j]='\0';
        strlwr(sline);
		int p;
		for(p=0;sline[p]!='	';p++)
		{
		}
		char tempName[20];
		string name;
		for(int q=0;p<j;q++)
		{
			p++;
			tempName[q]=sline[p];
		}
		name=tempName;
		dictTFIM[name]=line;
	}
	return dictTFIM;
}

void GetReady::readTUPosition(string TU_position_address)
{
    ifstream data(TU_position_address.c_str());
	char ch;
	for(;data.get(ch);)
	{
		while(ch=='#')
		{
			string noUse;
			getline(data,noUse);
			data.get(ch);
		}
		string line;
		getline(data,line);
		int i;
		for(i=0;line[i]!='	';i++){}
		i++;
		for(;line[i]!='	';i++){}
		i++;
		char ProName[10];
		int j;
		for(j=0;line[i]!='	';j++)
		{
			ProName[j]=line[i];
			i++;
		}
		i++;
		ProName[j]='\0';
        strlwr(ProName);
		char TUPos[10];
		for(j=0;line[i]!='	';j++)
		{
			TUPos[j]=line[i];
			i++;
		}
		TUPos[j]='\0';
		string tempTUPos=TUPos;
		promoter_name_dict.push_back(ProName);
		TU_position.push_back(atoi(tempTUPos.c_str()));
	}
}

int GetReady::getTFAmount()
{
	return TF_amount;
}

void GetReady::getGenePromoter(GeneIM temp_gene_IM[])
{
	int TUAmount;
    TUAmount=(int)TU_position.size();
	for(int i=0;i<getGeneAmount();i++)
	{
		int max=TUAmount,min=0;
		for(int j=(max+min)/2;(max-min)!=1;j=(max+min)/2)
		{
			if(temp_gene_IM[i].getLeftPosition()<=TU_position[j])
				max=j;
			else
				min=j;
		}
		temp_gene_IM[i].putInPromoterName(promoter_name_dict[min]);
	}
}

map<string,string> GetReady::mapPromoter(string promoters_address)
{
	map<string,string> dictPromoter;
    ifstream data(promoters_address.c_str());
	char ch;
	for(;data.get(ch);)
	{
		while(ch=='#')
		{
			string noUse;
			getline(data,noUse);
			data.get(ch);
		}
		string line;
		char nameLine[20];
		getline(data,line);
		int j,i;
		for(j=0;line[j]!='	';j++){}
		j++;
		for(i=0;line[j]!='	';i++)
		{
			nameLine[i]=line[j];
			j++;
		}
		nameLine[i]='\0';
		strlwr(nameLine);
		j++;
		for(;line[j]!='	';j++){}
		j++;
		for(;line[j]!='	';j++){}
		j++;
		for(;line[j]!='	';j++){}
		j++;
		char tempSequence[100];
		for(i=0;line[j]!='	';i++)
		{
			tempSequence[i]=line[j];
			j++;
		}
		tempSequence[i]='\0';
		strupr(tempSequence);
		string name=nameLine;
		dictPromoter[name]=tempSequence;
	}
	return dictPromoter;
}

void GetReady::inputUncertainGene()
{
	ifstream uncertain("uncertain_database");
	char ch;
	string del;
	for(int i=0;uncertain.get(ch);i++)
	{
		if(ch=='+')
			originalGRN[uncertain_row[i]][uncertain_column[i]]=1;
		else if(ch=='-')
			originalGRN[uncertain_row[i]][uncertain_column[i]]=-1;
		else
			originalGRN[uncertain_row[i]][uncertain_column[i]]=0;
		getline(uncertain,del);
	}
}
