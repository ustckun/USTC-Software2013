//
//  Sequence.cpp
//  GRN
//
//  Created by jinyang on 13-7-26.
//  Copyright (c) 2013Äê Li Jinyang. All rights reserved.
//


#include"TFIM.h"
#include"Calculate.h"
#include"Sequence.h"
#include"GRN.h"
#include"Regulation.h"
#include"EasytoDebug.h"

void Sequence::initializeGeneSequence( std::string sequence,int number, int size ){
    geneSequence += sequence;
    geneNumber = number;
    DNASize = size;
    aminoASSize = size / 3 - 1;
    //RNASequence +=sequence;
}
//get gene amino acid sequence;
std::string Sequence::getAminoAcidSequence(){
        //transcription();
    translation();
    return aminoAcidSequence;
}
    
//****************************************************************************************
//  Some explain of the transcription and translation process:
//    1.Actually, the protein expression process is:
//      DNA --> mRNA (i.e. transcription);
//      mRNA --> protein (i.e. translation).
//    2.DNA has double strands, but only one strand takes part in transcription.
//    3.Codons are the sequence messages carried by mRNA;
//      Take initiation codon "AUG" for example:
//      ---ATG--- : DNA strand which doesn't take part in transcription process;
//      ---TAC--- : DNA strand which exactlly takes part in transcription proess;
//      ---AUG--- : mRNA strand which carries codons. In this case, it carries initiation
//                 codon, i.e. "AUG";
//      ---TAC--- : tRNA which also carries amino acid Methionine(M);
//    4.Owing to the DNA sequences that our database provides are the UNEXPRESSION
//      strands, the translation process of the program can just use DNA sequence without
//      transcription. It should be noted that users have to enable the transcription
//      module and modify appropriate lines if their databases use different strands of
//      DNA sequences.
//****************************************************************************************
    
//transcript gene sequence(DNA) into RNA sequence;
/*void transcription(){
    for(int i = 1; i != geneSequence.size(); i++){
    RNASequence[i] = transcript(geneSequence[i]);
    }
     //return RNASequence;
}
char transcript(char s){
    if(s == 'G') return 'C';
    else if(s == 'A') return 'U';
    else if(s == 'T') return 'A';
    else if(s == 'C') return 'G';
}*/
    
//translate RNA sequence into amino acid sequence;
void Sequence::translation(){
    //amino acid codon;
    //U or T = 0; C = 1; A = 2; G = 3;
    char codon[4][4][4] = {
        {
            {'F', 'F', 'L', 'L'},
            {'S', 'S', 'S', 'S'},
            {'Y', 'Y', '#', '#'},
            {'C', 'C', '#', 'W'}
        },
        {
            {'L', 'L', 'L', 'L'},
            {'P', 'P', 'P', 'P'},
            {'H', 'H', 'Q', 'Q'},
            {'R', 'R', 'R', 'R'}
        },
        {
            {'I', 'I', 'I', 'M'},
            {'T', 'T', 'T', 'T'},
            {'N', 'N', 'K' ,'K'},
            {'S', 'S', 'R', 'R'}
        },
        {
            {'V', 'V', 'V', 'V'},
            {'A', 'A', 'A', 'A'},
            {'D', 'D', 'E', 'E'},
            {'G', 'G', 'G', 'G'}
        }
    };
    char mRNA[3];
    int index[3] = { 0 };
    aminoAcidSequence = " ";
    //If need to use transcription module, replace geneSequence with RNASquence;
    for (int i = 1; i != geneSequence.size(); i += 3) {
        mRNA[0] = geneSequence[i];
        mRNA[1] = geneSequence[i + 1];
        mRNA[2] = geneSequence[i + 2];
        //change A, U, C, G into number;
        for (int j = 0; j != 3; ++j) {
            index[j] = translate(mRNA[j]);
        }
        if (codon[index[0]][index[1]][index[2]] != '#') {
            aminoAcidSequence += codon[index[0]][index[1]][index[2]];
        }
    }
}
int Sequence::translate(char s){
    if (s == 'T') {
        return 0;
    }
    if (s == 'U') {
        return 0;
    }
    if (s == 'C') {
        return 1;
    }
    if (s == 'A') {
        return 2;
    }
    if (s == 'G') {
        return 3;
    }
}

void Sequence::getNewASS(string insertGene)
{
	geneNumber=165+1;
	DNASize=insertGene.length()*3+3;
	aminoASSize=insertGene.length();
	geneSequence=insertGene;
}

