//
//  Sequence.h
//  GRN
//
//  Created by jinyang on 13-7-26.
//  Copyright (c) 2013年 Li Jinyang. All rights reserved.
//

#ifndef __GRN__Sequence__
#define __GRN__Sequence__

#include <iostream>
class Sequence{
public:
    Sequence()
    {
        geneNumber = 0;
        DNASize = 0;
        aminoASSize = 0;
        geneSequence = " ";
        aminoAcidSequence = " ";
        //RNASequence = " ";
    }
    void initializeGeneSequence( std::string sequence,int number, int size );
    std::string aminoAcidSequence;
    int geneNumber;
    int DNASize;
    int aminoASSize;
    //translate RNA sequence into amino acid sequence;
    void translation();
private:
    std::string geneSequence;
    //std::string RNASequence;
    //transcript gene sequence(DNA) into RNA sequence;
    //void transcription();
    
    int translate(char s);
};



#endif /* defined(__GRN__Sequence__) */
