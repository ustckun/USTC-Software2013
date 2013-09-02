#include "PSO.h"

PSO::PSO(ModleNetwork New, int row,int column)
{
    ifstream matrix("old_GRN");
    temp_GRN=new double*[GENEAM];
	for(int i=0;i<GENEAM;i++)
		temp_GRN[i]=new double[TFScale];
    for(int i=0;i<GENEAM;i++)
	{
        for(int j=0;j<TFScale;j++)
		{
            matrix>>temp_GRN[i][j];
		}
	}
	New.Network_2(temp_GRN,column,row);
	for(int i=0;i<GENEAM;i++)
	{
		random_max[i]=0;
		random_min[i]=100;
		target[i]=New.value[i];
	}
}

double PSO::getFitness(vector<double> row_column_matrix,ModleNetwork New,int row,int column)
{
	double row_regulation[TFScale];
	double column_regulation[GENEAM];
    for(int i=0;i<TFScale;i++)
		row_regulation[i]=row_column_matrix[i];
    for(int i=0;i<GENEAM;i++)
		column_regulation[i]=row_column_matrix[TFScale+GENEAM-i-1];
	for(int i=0;i<TFScale;i++)
	{
		temp_GRN[row][i]=row_regulation[i];
	}
	for(int i=0;i<GENEAM;i++)
	{
		temp_GRN[i][column]=column_regulation[i];
	}
	New.Network_2(temp_GRN,column+1,row+1);
	return getVariance(New.value);
}

double PSO::getVariance(double A[GENEAM])
{
	double Variance=0;
	for(int i=0;i<GENEAM;i++)
	{
		Variance=Variance+(A[i]-target[i])*(A[i]-target[i]);
	}
	return Variance;
}

double PSO::random(double min,double max)
{
	double ran;
    ran = ((double)rand() / RAND_MAX) * (max - min) + min;
    return ran;
}

void PSO::getRange(int row,int column,ModleNetwork cal)
{
	for(int i=0;i<row;i++)
	{
		if(random(0,1)>0.5)
		{
			temp_GRN[column][i]=1;
		}
		else
			temp_GRN[column][i]=-1;
		//cout<<random(0,1)<<endl;
	}
	for(int j=0;j<column;j++)
	{
		if(random(0,1)>0.5)
		{
			temp_GRN[j][row]=1;
		}
		else
			temp_GRN[j][row]=-1;
	}
	cal.Network_2(temp_GRN,column,row);
	for(int m=0;m<GENEAM;m++)
	{
		if(random_max[m]<cal.value[m])
			random_max[m]=cal.value[m];
		if(random_min[m]>cal.value[m])
			random_min[m]=cal.value[m];
	}
}

int PSO::getMinLine(double A[GENEAM],int column)
{
	double Min=A[0];
	int B=0;
	for(int i=1;i<column;i++)
	{
		if(A[i]<Min)
		{
			B=i;
			Min=A[i];
		}
	}
	return B;
}

void PSO::Filter(int n)
{
	/*for(int i=0;i<n;i++)
	{
		if(toBest[i]>0&&toBest[i]<0.2)
			toPick[i]=1;
		else if(toBest[i]>0.2&&toBest[i]<0.4)
			toPick[i]=2;
		else if(toBest[i]>0.4&&toBest[i]<0.6)
			toPick[i]=3;
		else if(toBest[i]>0.6&&toBest[i]<0.8)
			toPick[i]=4;
		else if(toBest[i]>0.8&&toBest[i]<1)
			toPick[i]=5;
	}
	for(int i=0;i<n;i++)
	{
		if(edBest[i]>0&&edBest[i]<0.2)
			edPick[i]=1;
		else if(edBest[i]>0.2&&edBest[i]<0.4)
			edPick[i]=2;
		else if(edBest[i]>0.4&&edBest[i]<0.6)
			edPick[i]=3;
		else if(edBest[i]>0.6&&edBest[i]<0.8)
			edPick[i]=4;
		else if(edBest[i]>0.8&&edBest[i]<1)
			edPick[i]=5;
	}*/
}

void PSO::getPrediction(ModleNetwork New,int row,int column)
{
	vector<double> temp_row_column(TFScale+GENEAM);
	vector<double> temp_row_column_V(TFScale+GENEAM);
    vector<vector<double> > row_column_matrix(PARTICLENUM,vector<double>(TFScale+GENEAM,0));
    vector<vector<double> > PBest(PARTICLENUM,vector<double>(TFScale+GENEAM,0));
    vector<vector<double> > row_column_V_matrix(PARTICLENUM,vector<double>(TFScale+GENEAM,0));
	double fitness[PARTICLENUM];
	double tempFitness[PARTICLENUM];
	double bestFitness;
	for(int j=0;j<PARTICLENUM;j++)
	{
		for(int i=0;i<TFScale+GENEAM;i++)
		{
			temp_row_column[i]=random(Pmin,Pmax);
			temp_row_column_V[i]=random(Vmin,Vmax);
		}
		/*for(int i=0;i<TFScale+GENEAM;i++)
		{
			if(fabs(temp_row_column[i])<0.8)
			{
				temp_row_column[i]=0;
			}
		}*/
		row_column_matrix[j]=temp_row_column;
		row_column_V_matrix[j]=temp_row_column_V;
		PBest[j]=row_column_matrix[j];
		fitness[j]=getFitness(row_column_matrix[j],New,row,column);
	}
	int bestLine=getMinLine(fitness,PARTICLENUM);
	bestFitness=tempFitness[bestLine];
	int runTime=0;
	int w=0;
	while(bestFitness>minAccu&&runTime<NN)
	{
		for(int i=0;i<PARTICLENUM;++i)
		{
			for(int j=0;j<TFScale+GENEAM;++j)
			{
				w=1-(runTime/NN*0.5);
				row_column_V_matrix[i][j]=w*row_column_V_matrix[i][j]+random(0,2)*(PBest[i][j]-row_column_matrix[i][j])+random(0,2)*(PBest[bestLine][j]-row_column_matrix[i][j]);
				if(row_column_V_matrix[i][j]>Vmax)
					row_column_V_matrix[i][i]=Vmax;
				if(row_column_V_matrix[i][j]<Vmin)
					row_column_V_matrix[i][j]=Vmin;
				row_column_matrix[i][j]+=row_column_V_matrix[i][j];
				if(row_column_matrix[i][j]>Pmax)
					row_column_matrix[i][j]=Pmax;
				if(row_column_matrix[i][j]<Pmin)
					row_column_matrix[i][j]=Pmin;
			}
			tempFitness[i]=getFitness(row_column_matrix[i],New,row,column);
			if(tempFitness[i]>fitness[i])
			{
				fitness[i]=tempFitness[i];
				PBest[i]=row_column_matrix[i];
			}
		}
		bestLine=getMinLine(fitness,PARTICLENUM);
		runTime++;
	}
	for(int i=0;i<column;i++)
		edPick.push_back(PBest[bestLine][i]);
	for(int i=0;i<row;i++)
		toPick.push_back(PBest[bestLine][TFScale+GENEAM-i-1]);
}
