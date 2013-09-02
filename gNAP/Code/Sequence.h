//
//  Sequence.h
//  GRN
//
//  Created by jinyang on 13-7-26.
//  Copyright (c) 2013 Li Jinyang. All rights reserved.
//

#ifndef __GRN__Sequence__
#define __GRN__Sequence__

#include <iostream>
class Sequence{
public:
    Sequence()
    {
        regulation_unit_number = 0;
        gene_size = 0;
        amino_acid_sequence_size = 0;
        gene_sequence = " ";
        promoter_sequence = " ";
        amino_acid_sequence = " ";
        //RNASequence = " ";
    }
    std::string gene_sequence;
    std::string promoter_sequence;
    std::string amino_acid_sequence;
    int regulation_unit_number;
    int gene_size;
    int promoter_size;
    int amino_acid_sequence_size;
    void initialize_Sequence(int RU_number, std::string promoter, int p_size,
                            std::string gene, int g_size);
    //translate RNA sequence into amino acid sequence;
    void Translation();
    void GetAAS(std::string aas, int ass_size,
                std::string promoter, int p_size){
        amino_acid_sequence = aas;
        amino_acid_sequence_size = ass_size;
        promoter_sequence += promoter;
        promoter_size = p_size;
    }
private:
    //std::string RNASequence;
    //transcript gene sequence(DNA) into RNA sequence;
    //void transcription();
    //char transcript(char s)'
    int Translate(char s);
};

#endif /* defined(__GRN__Sequence__) */
