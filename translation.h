//
//  transcription.h
//  SequenceAlignment
//
//  Created by jinyang on 13-7-24.
//  Copyright (c) 2013 jinyang. All rights reserved.
//

#ifndef __SequenceAlignment__transcription__
#define __SequenceAlignment__transcription__

#include <iostream>
class Sequence{
public:
    Sequence()
    {
       initializeGeneSequence(" "," ",0);
    }
    void initializeGeneSequence( std::string sequence,std::string RNAsequence,int number ){
        geneNumber = number;
        geneSequence += sequence;
        RNASequence +=RNAsequence;
    }
    void initializeGeneSequence( std::string sequence,int number ){
        geneNumber = number;
        geneSequence += sequence;
        RNASequence +=sequence;
    }
    //get gene amino acid sequence;
    std::string getAminoAcidSequence(){
        transcription();
        translation();
        return aminoAcidSequence;
    }
private:
    std::string geneSequence;
    std::string RNASequence;
    std::string aminoAcidSequence;
    int geneNumber;
    //transcript gene sequence(DNA) into RNA sequence;
    void transcription(){
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
    }
    //translate RNA sequence into amino acid sequence;
    void translation(){
        //amino acid codon;
        //U = 0; C = 1; A = 2; G = 3;
        char codon[4][4][4] = {{{'F', 'F', 'L', 'L'}, {'S', 'S', 'S', 'S'}, {'Y', 'Y', '#', '#'}, {'C', 'C', '#', 'W'}}, {{'L', 'L', 'L', 'L'}, {'P', 'P', 'P', 'P'}, {'H', 'H', 'Q', 'Q'}, {'R', 'R', 'R', 'R'}}, {{'I', 'I', 'I', 'M'}, {'T', 'T', 'T', 'T'}, {'N', 'N', 'K' ,'K'}, {'S', 'S', 'R', 'R'}}, {{'V', 'V', 'V', 'V'}, {'A', 'A', 'A', 'A'}, {'D', 'D', 'E', 'E'}, {'G', 'G', 'G', 'G'}}};
        char tRNA[3];
        int index[3] = { 0 };
        aminoAcidSequence = " ";
        for (int i = 1; i != RNASequence.size(); i += 3) {
            tRNA[0] = RNASequence[i];
            tRNA[1] = RNASequence[i + 1];
            tRNA[2] = RNASequence[i + 2];
            //change A, U, C, G into number;
            for (int j = 0; j != 3; ++j) {
                index[j] = translate(tRNA[j]);
            }
        if (codon[index[0]][index[1]][index[2]] != '#') {
            aminoAcidSequence += codon[index[0]][index[1]][index[2]];
        }
        }
    }
    int translate(char s){
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

#endif /* defined(__SequenceAlignment__transcription__) */
