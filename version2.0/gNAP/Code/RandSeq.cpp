////////////////////////////////////////////////////////////////////////////////
/// COPYRIGHT NOTICE\n
/// Distribute under BSD License\n
/// Copyright (c) 2013, iGEM Software Team of University of Science and
/// Technology of China\n
/// All rights reserved.
///
/// \file RandSeq.cpp
/// \brief Statements of functions of class RandSeq
/// \version 1.0
/// \author Li Jinyang
/// \date Aug. 9, 2013
////////////////////////////////////////////////////////////////////////////////
#include "RandSeq.h"
#include <iostream>
#include <ctime>
#include "stdlib.h"
void RandomSequence::generate_random_amino_acid_sequence(int length){
    for (int i = 0; i != length; ++i) {
        random_amino_acid_sequence += GenerateRandomAminoAcid();
    }
}
char RandomSequence::GenerateRandomAminoAcid(){
    switch (rand() % 20) {
        case 0:
            return 'A';
            break;
        case 1:
            return 'R';
            break;
        case 2:
            return 'N';
            break;
        case 3:
            return 'D';
            break;
        case 4:
            return 'C';
            break;
        case 5:
            return 'Q';
            break;
        case 6:
            return 'E';
            break;
        case 7:
            return 'G';
            break;
        case 8:
            return 'H';
            break;
        case 9:
            return 'I';
            break;
        case 10:
            return 'L';
            break;
        case 11:
            return 'K';
            break;
        case 12:
            return 'M';
            break;
        case 13:
            return 'F';
            break;
        case 14:
            return 'P';
            break;
        case 15:
            return 'S';
            break;
        case 16:
            return 'T';
            break;
        case 17:
            return 'W';
            break;
        case 18:
            return 'Y';
            break;
        case 19:
            return 'V';
            break;
    }
}

