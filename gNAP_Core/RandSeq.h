//
//  RandSeq.h
//  Random Sequence
//
//  Created by jinyang on 13-8-9.
//  Copyright (c) 2013Äê Li Jinyang. All rights reserved.
//

#ifndef __Random_Sequence__RandSeq__
#define __Random_Sequence__RandSeq__

#include <iostream>
class RandomSequence{
public:
    RandomSequence(){
        random_amino_acid_sequence = " ";
    }
    void generate_random_amino_acid_sequence(int length);
    std::string random_amino_acid_sequence;
private:
    char GenerateRandomAminoAcid();
};
#endif /* defined(__Random_Sequence__RandSeq__) */
