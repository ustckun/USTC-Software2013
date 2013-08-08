//class Regulation
//get regulation into a matrix

class Regulation
{
public:
	Regulation()
	{
		geneAmount=0;
		openFileError=0;
	}
	void readName(TFIM geneFirst[]);//read Name in TF-TF.txt to TFIM, intput is TFIM projects array
	int getGeneAmount();
	void fullFill();
	int getOpenError();
	double originalMatrix[N][N];;//0:no regulation 1:+active -1:-negative 2:NULL 0.5:+ -
private:
	int geneAmount;
	int openFileError;
};