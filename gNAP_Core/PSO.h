#ifndef __USTC_Software__PSO__
#define __USTC_Software__PSO__


#include <vector>
#include <cstdlib>

#include "define.h"
#include "ModleNetwork.h"

class PSO
{
public:
	PSO(ModleNetwork New,double **matrix, int row,int column);
	double target[GENEAM];
	void getPrediction(ModleNetwork New,int row,int column);
	vector<double> toPick;
	vector<double> edPick;
	void getRange(int row,int column,ModleNetwork cal);
	double random_max[GENEAM];
	double random_min[GENEAM];
	void Filter(int n);
private:
	double **temp_GRN;
	int getMinLine(double A[GENEAM],int column);
	double getFitness(vector<double> row_column_matrix,ModleNetwork New,int row,int column);
	double getVariance(double A[GENEAM]);
	double random(double min,double max);
};

#endif