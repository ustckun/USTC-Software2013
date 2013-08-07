class PSOPredict
{
public:
	PSOPredict();
	double targetGene[N];
	void algorithm_PSO(double target_gene[],int GA);
	vector<double> best;
private:
	int min_Line(double fitness[], int M);
	double fitness_Function (vector<double> x, double GRN[scale][scale], double target_gene[N],int GA);
	double compareConcentration (double c_new[], double c_target[]);
	double random(double min, double max);
};

PSOPredict::PSOPredict()
{
	for(int i=0;i<N;i++)
		targetGene[i]=-1;
}
