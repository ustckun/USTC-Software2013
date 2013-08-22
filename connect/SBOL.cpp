#include"GeneIM.h"
#include"SBOL.h"

string SBOL::Combine(string title,string detail)
{
	string all;
	all=FormartStart(title)+detail+FormartEnd(title);
	return all;
}

string SBOL::FormartStart(string a)
{
	string all;
	string start="<";
	string end=">";
	all=start+a+end;
	return all;
}

string SBOL::FormartEnd(string b)
{
	string all;
	string start="</";
	string end=">";
	all=start+b+end;
	return all;
}

void SBOL::CreatSBOL(string position,GeneIM IM)
{
	string FP=".xml";
	string temp=IM.getGeneName();
	string Seq=IM.getGeneSequence();
	//vector <char>seq;
	//for(int i=0;i<Seq.length();i++)
	//seq[i]=Seq[i];
	FP=position+temp+FP;
	ofstream SBOLfile(FP);
	//head='<?xml version="1.0"?>\n<rdf:RDF\n	xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#"\n	xmlns:rdfs="http://www.w3.org/2000/01/rdf-schema#" \n	xmlns:so="http://purl.obolibrary.org/obo/" \n	xmlns="http://sbols.org/v1#"\n	\n	>\n\n';
	SBOLfile<<"<?xml version=\"1.0\"?>"<<endl;
	SBOLfile<<"<rdf:RDF"<<endl;
	SBOLfile<<"	xmlns:rdf=\"http://www.w3.org/1999/02/22-rdf-syntax-ns#\""<<endl;
	SBOLfile<<"	xmlns:rdfs=\"http://www.w3.org/2000/01/rdf-schema#\" "<<endl;
	SBOLfile<<"	xmlns:so=\"http://purl.obolibrary.org/obo/\" "<<endl;
	SBOLfile<<"	xmlns=\"http://sbols.org/v1#\""<<endl<<endl<<"	>"<<endl<<endl;
	SBOLfile<<"	<DnaComponent rdf:about=\"http://example.com/dc1\">"<<endl;
	SBOLfile<<"		"<<Combine("displayId",IM.getID())<<endl;
	SBOLfile<<"		"<<Combine("name",IM.getGeneName())<<endl;
	SBOLfile<<"		"<<Combine("description",IM.getGeneDescription())<<endl;
	SBOLfile<<"		"<<"<dnaSequence>"<<endl;
	SBOLfile<<"			"<<"<DnaSequence rdf:about=\"http://example.com/ds1\">"<<endl;
	//SBOLfile<<"				"<<Combine("nucleotides",Seq)<<endl;
	SBOLfile<<"				"<<"<nucleotides>"<<Seq<<"</nucleotides>"<<endl;
	SBOLfile<<"			"<<"</DnaSequence>"<<endl;
	SBOLfile<<"		"<<"</dnaSequence>"<<endl;
	SBOLfile<<"		"<<"<annotation>"<<endl;
	SBOLfile<<"			<SequenceAnnotation rdf:about=\"http://example.com/sa1\">"<<endl;
	SBOLfile<<"				"<<FormartStart("biostart")<<IM.getLeftPosition()<<FormartEnd("biostart")<<endl;
	SBOLfile<<"				"<<FormartStart("biostart")<<IM.getRightPosition()<<FormartEnd("biostart")<<endl;
	SBOLfile<<"			</SequenceAnnotation>"<<endl;
	SBOLfile<<"		"<<"</annotation>"<<endl;
	SBOLfile<<"	</DnaComponent>"<<endl;
	SBOLfile<<"</rdf:RDF>";
	SBOLfile.close();
}