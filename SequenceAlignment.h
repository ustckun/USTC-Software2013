//
//  SequenceAlignment.h
//  SequenceAlignment
//
//  Created by jinyang on 13-7-25.
//  Copyright (c) 2013年 jinyang. All rights reserved.
//

#ifndef SequenceAlignment_SequenceAlignment_h
#define SequenceAlignment_SequenceAlignment_h


#include <iostream>
#include <vector>
class Sequence{
public:
    Sequence()
    {
        initializeGeneSequence(" "," ",0);
    }
    void initializeGeneSequence( std::string sequence,std::string RNAsequence,int number ){
        geneNumber = number;
        geneSequence += sequence;
        //RNASequence +=RNAsequence;
    }
    void initializeGeneSequence( std::string sequence,int number ){
        geneNumber = number;
        geneSequence += sequence;
        //RNASequence +=sequence;
    }
    //get gene amino acid sequence;
    std::string getAminoAcidSequence(){
        //transcription();
        translation();
        return aminoAcidSequence;
    }
private:
    std::string geneSequence;
    //std::string RNASequence;
    std::string aminoAcidSequence;
    int geneNumber;
    
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
    void translation(){
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
    int translate(char s){
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
    
};

class Alignment {
public:
    Alignment(){
        correlation = 0;
        gene_number = 0;
    }
    double gene_number;
    double getCorrelation(std::string refAAS, int refAASSize, std::string inqAAS, int inqAASSize){
        aminoASAlignment(refAAS, refAASSize, inqAAS, inqAASSize);
        return correlation;
    }
    int getGeneNumber(){
        return gene_number;
    }
private:
    double correlation;
    //AAS is the short for Amino Acid Sequence;
    //ref: reference;
    //inq: inquire;
    void aminoASAlignment(std::string refAAS, int refAASSize, std::string inqAAS, int inqAASSize ){
        int G_ = 0;
        int _G = 0;
        std::vector<std::vector<int>> alignMatirx(refAASSize, std::vector<int> (inqAASSize));
        //initialize alignMatrix;
        alignMatirx[0][0] = 0;
        for (int i = 1; i != refAASSize; ++i) {
            alignMatirx[i][0] = alignMatirx[i - 1][0] + G_;
        }
        for (int j = 1; j != inqAASSize; ++j) {
            alignMatirx[0][j] = alignMatirx[0][j - 1] + _G;
        }
        for (int i = 1; i != refAASSize; ++i) {
            for (int j = 1; j != inqAASSize; ++j) {
                alignMatirx[i][j] =
                maxValue(
                         alignMatirx[i-1][j-1] + alignScore(inqAAS[j] ,refAAS[i]),
                         alignMatirx[i-1][j] + alignScore(inqAAS[j], ' '),
                         alignMatirx[i][j-1] + alignScore(' ', refAAS[i])
                         );
            }
        }
        correlation = (double)alignMatirx[refAASSize - 1][inqAASSize - 1] * 2 / (refAASSize + inqAASSize);
        
    }
    int alignScore (char s, char t)
    {
        int score = 0;
        int GG = 1;
        int _G = 0;
        int GT = 0;
        if (s == t)
            score = GG;
        else if (s == ' ')
            score = _G;
        else if (t == ' ')
            score = _G;
        else if (s != t)
            score = GT;
        return score;
    }
    int maxValue (int a, int b, int c)
    {
        if(a<b)
            a=b;
        if(a<c)
            a=c;
        return a;
    }
    
};

#endif /* defined(__SequenceAlignment__transcription__) */
