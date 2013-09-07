////////////////////////////////////////////////////////////////////////////////
/// COPYRIGHT NOTICE\n
/// Distribute under BSD License\n
/// Copyright (c) 2013, iGEM Software Team of University of Science and
/// Technology of China\n
/// All rights reserved.
///
/// \file SBOL.cpp
/// \brief Statments of funcions of the class SBOL.
/// \version 1.0
/// \author Wang Chenkun
/// \date September 2nd, 2013
////////////////////////////////////////////////////////////////////////////////

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

void SBOL::CreatSBOL(string gene_name,string ID,string left,string right,string description,string seq)
{
	string FP=".xml";
    string temp=gene_name;
	//vector <char>seq;
	//for(int i=0;i<Seq.length();i++)
	//seq[i]=Seq[i];
    FP=temp+FP;
    ofstream SBOLfile(FP.c_str());
	//head='<?xml version="1.0"?>\n<rdf:RDF\n	xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#"\n	xmlns:rdfs="http://www.w3.org/2000/01/rdf-schema#" \n	xmlns:so="http://purl.obolibrary.org/obo/" \n	xmlns="http://sbols.org/v1#"\n	\n	>\n\n';
	SBOLfile<<"<?xml version=\"1.0\"?>"<<endl;
	SBOLfile<<"<rdf:RDF"<<endl;
	SBOLfile<<"	xmlns:rdf=\"http://www.w3.org/1999/02/22-rdf-syntax-ns#\""<<endl;
	SBOLfile<<"	xmlns:rdfs=\"http://www.w3.org/2000/01/rdf-schema#\" "<<endl;
	SBOLfile<<"	xmlns:so=\"http://purl.obolibrary.org/obo/\" "<<endl;
	SBOLfile<<"	xmlns=\"http://sbols.org/v1#\""<<endl<<endl<<"	>"<<endl<<endl;
    SBOLfile<<"	<DnaComponent rdf:about=\"http://2013.igem.org/Team:USTC-Software/SBOL/E"+ID+"\">"<<endl;
    SBOLfile<<"		"<<Combine("displayId",ID)<<endl;
    SBOLfile<<"		"<<Combine("name",gene_name)<<endl;
    SBOLfile<<"		"<<Combine("description",description)<<endl;
	SBOLfile<<"		"<<"<dnaSequence>"<<endl;
    SBOLfile<<"			"<<"<DnaSequence rdf:about=\"http://2013.igem.org/Team:USTC-Software/SBOL/E"+ID+"\">"<<endl;
	//SBOLfile<<"				"<<Combine("nucleotides",Seq)<<endl;
    SBOLfile<<"				"<<"<nucleotides>"<<seq<<"</nucleotides>"<<endl;
	SBOLfile<<"			"<<"</DnaSequence>"<<endl;
	SBOLfile<<"		"<<"</dnaSequence>"<<endl;
	SBOLfile<<"		"<<"<annotation>"<<endl;
    SBOLfile<<"			<SequenceAnnotation rdf:about=\"http://2013.igem.org/Team:USTC-Software/SBOL/E"+ID+"\">"<<endl;
    SBOLfile<<"				"<<FormartStart("biostart")<<left<<FormartEnd("biostart")<<endl;
    SBOLfile<<"				"<<FormartStart("biostart")<<right<<FormartEnd("biostart")<<endl;
	SBOLfile<<"			</SequenceAnnotation>"<<endl;
	SBOLfile<<"		"<<"</annotation>"<<endl;
	SBOLfile<<"	</DnaComponent>"<<endl;
	SBOLfile<<"</rdf:RDF>";
	SBOLfile.close();
}
