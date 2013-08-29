////////////////////////////////////////////////////////////////////////////////
/// COPYRIGHT NOTICE\n
/// Distribute under BSD License\n
/// Copyright (c) 2013, iGEM Software Team of University of Science and
/// Technology of China\n
/// All rights reserved.
///
/// \file RandSeq.h
/// \brief Define the class RandSeq.
///
///    Generate a random amino acid sequence at a specific length.
/// \version 1.0
/// \author Li Jinyang
/// \date Aug. 9, 2013
////////////////////////////////////////////////////////////////////////////////
#ifndef __Random_Sequence__RandSeq__
#define __Random_Sequence__RandSeq__

#include <iostream>
///     Generate a random amino aicd sequence at a specific length.
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
