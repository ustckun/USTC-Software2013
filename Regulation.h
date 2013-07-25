//class Regulation
//get regulation into a matrix

class Regulation
{
public:
	void readName(TFIM geneFirst[]);//read Name in TF-TF.txt to TFIM
	int getGeneAmount();
	void fullFill();
	int getOpenError();
	float originalMatrix[250][250];;//0:no regulation 1:+active -1:-negative 2:NULL 0.5:+ -
private:
	int geneAmount;
	int openFileError;
};