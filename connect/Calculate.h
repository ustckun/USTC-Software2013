
class Calculate
{
public:
	Calculate()
    {
        int ik;
        for(ik=0;ik<DIMENS;++ik)
        {
            p[ik]=18;q[ik]=9;r[ik]=2;nn[ik]=6;
        }
    }
    double MaxMa[DIMENS][DIMENS];
    double consistence[DIMENS];
    //New regulation Matrix and number of genes into this function. It will putout the score and gene's changing into documents
	void Network_1(double ReguMatrix[][DIMENS],int n);
	//New Matrix with random regulation and gene number into this function. It will calculate the stable consistence to consistence array. 
    void Network_2(double Matr[][DIMENS],int n);
	void RandMatrix(double a[][DIMENS],double b[][DIMENS],const int n);
private:
    double p[DIMENS],q[DIMENS],r[DIMENS],nn[DIMENS];
	
	double FaNexVal(double Matr[][DIMENS],double a[],const int n,const int i,double p[],double q[],double nn[],double r[]);
	double NextValue(double Matr[][DIMENS],double a[],const int n,const int i,double step,double p[],double q[],double nn[],double r[]);
};
