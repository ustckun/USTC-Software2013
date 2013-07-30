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
    //get gene amino acid sequence;
    std::string getAminoAcidSequence();
    int geneNumber;
    int DNASize;
    int aminoASSize;
	void getNewASS(string insertGene);
private:
    std::string geneSequence;
    //std::string RNASequence;
    std::string aminoAcidSequence;
    //transcript gene sequence(DNA) into RNA sequence;
    //void transcription();
    
    //translate RNA sequence into amino acid sequence;
    void translation();
    int translate(char s);
};



#endif /* defined(__GRN__Sequence__) */
