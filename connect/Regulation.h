class Regulation
{
public:
	Regulation();
	void getRegulationMatrix(TFIM geneFirst[]);
	int getGeneAmount();
	int getTFAmount();
	void fullFill();
	int getOpenError();
	//double originalMatrix[N][N];//0:no regulation 1:+active -1:-negative 2:NULL 0.5:+ -
	double **originalGRN;
	map<string,string> mapTFIM();
	map<string,string> mapPromoter();
	vector<int>uncertain1;
	vector<int>uncertain2;
	void readTUPosition();
	vector<int> TUPositon;
	vector<string> promoterNameLib;
	void getGenePromoter(TFIM geneFirst[]);
private:
	void readTFTF(TFIM geneFirst[],double **GRN);//read Name in TF-TF.txt to TFIM, intput is TFIM projects array
	void readTFGene(TFIM geneFirst[],double **GRN);
	void addTF(TFIM geneFirst[]);
	int geneAmount;
	int TFAmount;
	int openFileError;
};