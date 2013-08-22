//
//  SeqAlign.h
//  alignment2.0
//
//  Created by jinyang on 13-8-9.
//  Copyright (c) 2013年 Li Jinyang. All rights reserved.
//

#ifndef __alignment2_0__SequenceAlignment__
#define __alignment2_0__SequenceAlignment__

#include <iostream>
using namespace std;
class SequenceAlignment{
public:
    SequenceAlignment() {
        for (int i = 0; i != 20; ++i) {
            for (int j = 0; j != 20; ++j) {
                BLOSUM[i][j] = 0;
            }
        }
        AAS_similarity_ = 0;
        DNA_similarity_ = 0;
        load_matrix_BLOSUM();
    }
    void AASAlignment(string query, int query_size,
                                 string subject, int subject_size);
    void DNASeqAlignment(string  query, int query_size,
                           string subject, int subject_size);
    double get_AAS_similarity_(){
        return AAS_similarity_;
    }
    double get_DNA_similarity_(){
        return DNA_similarity_;
    }
private:
    double AAS_similarity_;
    double DNA_similarity_;
    int BLOSUM[20][20];
    int AASAlignScore (char t, char s);
    int DNASequenceAlignScore(char t, char s);
    double get_max_value (double a, double b, double c);
    int get_index_of_BLOSUM50(char s);
    void load_matrix_BLOSUM();
};

#endif /* defined(__alignment2_0__SeqAlign__) */
