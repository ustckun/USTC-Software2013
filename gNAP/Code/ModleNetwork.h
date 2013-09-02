#ifndef __USTC_Software__ModleNetwork__
#define __USTC_Software__ModleNetwork__

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <ctime>

#include"define.h"

using namespace std;

class ModleNetwork
{
public:
    ModleNetwork()
    {
        int ik;
        p=new double[GENEAM];
        q=new double[GENEAM];
        r=new double[GENEAM];
        nn=new double[GENEAM];
        for(ik=0;ik<GENEAM;++ik)
        {
            p[ik]=18;q[ik]=9;r[ik]=2;nn[ik]=6;
        }
        value=new double[GENEAM];
        MaxMa=new double*[GENEAM];
        for(ik=0;ik<GENEAM;++ik)MaxMa[ik]=new double[TFScale];
    }
    double **MaxMa;
    double *value;
    void Network_1(double **ReguMatrix,int nx,int ny);
    void Network_2(double **Matr,int nx,int ny);
private:
    double *p,*q,*r,*nn;
    void RandMatrix(double **a,double **b,const int nx,const int ny);
    double FaNexVal(double **Matr,double a[],const int nx,const int i,double p[],double q[],double nn[],double r[]);
};

#endif
