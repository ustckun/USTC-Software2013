//
//  Sequence.h
//  GRN
//
//  Created by jinyang on 13-7-26.
//  Copyright (c) 2013Äê Li Jinyang. All rights reserved.
//

#ifndef __GRN__Sequence__
#define __GRN__Sequence__

//#include <iostream>
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
    std::string geneSequence;
    std::string aminoAcidSequence;
    int geneNumber;
    int DNASize;
    int aminoASSize;
    //translate RNA sequence into amino acid sequence;
    void translation();
	void getNewASS(string insertGene,int n);
private:
    //std::string RNASequence;
    //transcript gene sequence(DNA) into RNA sequence;
    //void transcription();
    
    int translate(char s);
};



#endif /* defined(__GRN__Sequence__) */
