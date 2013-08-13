class PSOPredict
{
public:
	PSOPredict();
	double targetGene[N];
	void algorithm_PSO(double GRN_Origin[scale][scale],double target_gene[N],int n, Calculate NEW);
	vector<double> toPick;
	vector<double> edPick;
	void getRange(double GRN_Origin[scale][scale],int n,Calculate cal);
	double randomMax[N];
	double randomMin[N];
	void Filter(int n);
private:
	int min_Line(double fitness[]);
	double fitness_Function (vector<double> x, double GRN[scale][scale], double target_gene[scale], int n, Calculate NEW);
	double compareConcentration (double c_new[], double c_target[],int n);
	double random(double min, double max);
	vector<double> toBest;
	vector<double> edBest;
};


