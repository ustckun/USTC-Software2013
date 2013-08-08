//EasytoDebug.h
//Class EasytoDebug
//read number, sequence, length of sequence in "Sequecne"
//

class ReadDNA
{
public:
	ReadDNA()
	{
		number=0;
		length=0;
	}
	int getNumber();
	int getLength();
	string getSequence();
	void fetch(FILE *temp);
private:
	int number;
	int length;
	string sequence;
};