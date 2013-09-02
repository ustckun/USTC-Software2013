//
//  GRN.h
//  GRN
//
//  Created by jinyang on 13-7-26.
//  Copyright (c) 2013Äê Li Jinyang. All rights reserved.
//

#ifndef __GRN__GRN__
#define __GRN__GRN__
#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>
#include "Sequence.h"
class GRN{
public:
    GRN() {
        new_GRN = new double*[1800];
        for (int i = 0; i != 1800; i++) {
            new_GRN[i] = new double[220];
        }
        for (int i = 0; i != 20; ++i) {
            for (int j = 0; j != 20; ++j) {
                BLOSUM[i][j] = 0;
            }
        }
        number_row = 0;
        number_column = 0;
    }
    double **new_GRN;
    void initialize_GRN(double **old_GRN, int num_row, int num_column);
    void construct_new_GRN(Sequence reg_unit[]);
    double AminoAcidSeqAlignment(std::string query, int query_size,
                                 std::string subject, int subject_size);
    double DNASeqAlignment(std::string  query, int query_size,
                                std::string subject, int subject_size);
    void load_matrix_BLOSUM();
private:
    // The number of rows.
    int number_row;
    // The number of columns.
    int number_column;
    int BLOSUM[20][20];
    int AminoAcidSequenceAlignScore (char t, char s);
    int DNASequenceAlignScore(char t, char s);
    double get_max_value (double a, double b, double c);
    int get_index_of_BLOSUM50(char s);
};

#endif /* defined(__GRN__GRN__) */
