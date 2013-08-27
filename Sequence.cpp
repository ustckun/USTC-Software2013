////////////////////////////////////////////////////////////////////////////////
/// COPYRIGHT NOTICE\n
/// Distribute under BSD License\n
/// Copyright (c) 2013, iGEM Software Team of University of Science and
/// Technology of China\n
/// All rights reserved.
///
/// \file Sequence.cpp
/// \brief Statments of funcions of the class Sequence.
/// \version 1.0
/// \author Li Jinyang
/// \date July 26, 2013
////////////////////////////////////////////////////////////////////////////////
#include "Sequence.h"
void Sequence::initialize_Sequence(int RU_number, std::string promoter, int p_size,
                                  std::string gene, int g_size) {
    regulation_unit_number = RU_number;
    gene_sequence += gene;
    gene_size = g_size;
    promoter_sequence += promoter;
    promoter_size = p_size;
    Translation();
}
void Sequence::Translation(){
    // Amino acid codon;
    // U or T = 0; C = 1; A = 2; G = 3;
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
    amino_acid_sequence = " ";
    for (int i = 1; i < gene_size; i += 3) {
        mRNA[0] = gene_sequence[i];
        mRNA[1] = gene_sequence[i + 1];
        mRNA[2] = gene_sequence[i + 2];
        // Change A, T, C, G into number;
        for (int j = 0; j != 3; ++j) {
            index[j] = Translate(mRNA[j]);
        }
        if (codon[index[0]][index[1]][index[2]] != '#') {
            amino_acid_sequence += codon[index[0]][index[1]][index[2]];
        }
    }
    amino_acid_sequence_size = (int)amino_acid_sequence.size() - 1;
}
int Sequence::Translate(char s){
    switch (s) {
        case 'T':
            return 0;
            break;
        case 'C':
            return 1;
            break;
        case 'A':
            return 2;
            break;
        case 'G':
            return 3;
            break;
        default:
            break;
    }
}
    

