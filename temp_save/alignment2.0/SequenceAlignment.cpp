//
//  SeqAlign.cpp
//  alignment2.0
//
//  Created by jinyang on 13-8-9.
//  Copyright (c) 2013年 Li Jinyang. All rights reserved.
//

#include "SeqAlign.h"
#include <vector>
#include <fstream>
#define GAP -8
#define GAP_2 -1
using namespace std;
void SequenceAlignment::AASAlignment(string query, int query_size,
                                     string subject, int subject_size) {
    double G_ = GAP;
    double _G = GAP;
    int query_length = 0;
    int subject_length = 0;
    query_length = query_size;
    subject_length = subject_size;
    string common_subsequence;
    string query_sub;
    string subject_sub;
    vector< vector<double> > align_matrix
        (subject_size + 1, vector<double>(query_size + 1));
    //initialize alignMatrix;
    align_matrix[0][0] = 0;
    for (int i = 1; i != subject_size + 1; ++i) {
        align_matrix[i][0] = align_matrix[i - 1][0] + G_;
    }
    for (int j = 1; j != query_size + 1; ++j) {
        align_matrix[0][j] = align_matrix[0][j - 1] + _G;
    }
    for (int i = 1; i != subject_size + 1; ++i) {
        for (int j = 1; j != query_size + 1; ++j) {
            double st = 0;
            double s_ = 0;
            double _t = 0;
            st = align_matrix[i - 1][j - 1] +
            AASAlignScore(subject[i], query[j]);
            _t = align_matrix[i - 1][j] +
            AASAlignScore(subject[i], ' ');
            s_ = align_matrix[i][j - 1] +
            AASAlignScore(' ', query[j]);
            align_matrix[i][j] = get_max_value(st, s_, _t);
        }
    }
    for (int i = subject_size, j = query_size; (i > 0 || j > 0);) {
        if (
            (i * j > 0) && (
                            align_matrix[i][j] ==
                            align_matrix[i - 1][j - 1] +
                            AASAlignScore(subject[i], query[j])
                            )
            ) {
            if (query[j] == subject[i]) {
                common_subsequence += query[j];
                query_sub += query[j];
                subject_sub += subject[i];
            } else {
                int similarity_test = 0;
                similarity_test = BLOSUM[get_index_of_BLOSUM50(query[j])]
                [get_index_of_BLOSUM50(subject[i])];
                if (similarity_test >= 0) {
                    common_subsequence += '+';
                }
                query_sub += query[j];
                subject_sub += subject[i];
            }
            i--;
            j--;
        } else if ((i > 0) &&
                   (align_matrix[i][j] == align_matrix[i - 1][j] +
                    AASAlignScore(subject[i], ' '))) {
                       //insert a space into string_s;
                       query_length += 1;
                       query_sub += '-';
                       subject_sub += subject[i];
                       i--;
                   } else if ((j > 0) &&
                              (align_matrix[i][j] == align_matrix [i][j - 1] +
                               AASAlignScore(' ', query[j] ))){
                                  //insert a space into string_s;
                                  subject_length += 1;
                                  subject_sub += '-';
                                  query_sub += query[j];
                                  j--;
                              }
    }
    AAS_similarity_ = (double)common_subsequence.size() / query_length;
}
void SequenceAlignment::DNASeqAlignment(string  query, int query_size,
                                        string subject, int subject_size) {
    double G_ = GAP_2;
    double _G = GAP_2;
    int query_length = 0;
    int subject_length = 0;
    query_length = query_size;
    subject_length = subject_size;
    string common_subsequence;
    // Records result of alignment.
    string query_sub;
    string subject_sub;
    vector< vector<double> > align_matrix
    (subject_size + 1, vector<double>(query_size + 1));
    //initialize alignMatrix;
    align_matrix[0][0] = 0;
    for (int i = 1; i != subject_size + 1; ++i) {
        align_matrix[i][0] = align_matrix[i - 1][0] + G_;
    }
    for (int j = 1; j != query_size + 1; ++j) {
        align_matrix[0][j] = align_matrix[0][j - 1] + _G;
    }
    for (int i = 1; i != subject_size + 1; ++i) {
        for (int j = 1; j != query_size + 1; ++j) {
            double st = 0;
            double s_ = 0;
            double _t = 0;
            st = align_matrix[i - 1][j - 1] +
            DNASequenceAlignScore(subject[i], query[j]);
            _t = align_matrix[i - 1][j] +
            DNASequenceAlignScore(subject[i], ' ');
            s_ = align_matrix[i][j - 1] +
            DNASequenceAlignScore(' ', query[j]);
            align_matrix[i][j] = get_max_value(st, s_, _t);
        }
    }
    for (int i = subject_size, j = query_size; (i > 0 || j > 0);) {
        if ((i * j > 0) && (align_matrix[i][j] == align_matrix[i - 1][j - 1] +
                            DNASequenceAlignScore(subject[i], query[j]))) {
            if (query[j] == subject[i]) {
                common_subsequence += query[j];
                query_sub += query[j];
                subject_sub += subject[i];
            } else {
                query_sub += query[j];
                subject_sub += subject[i];
            }
            i--;
            j--;
        } else if ((i > 0) && (align_matrix[i][j] == align_matrix[i - 1][j] +
                               DNASequenceAlignScore(subject[i], ' '))) {
            //insert a space into string_s;
            query_length += 1;
            query_sub += '-';
            subject_sub += subject[i];
            i--;
        } else if ((j > 0) && (align_matrix[i][j] == align_matrix[i][j - 1] +
                               DNASequenceAlignScore(' ', query[j]))) {
            //insert a space into string_s;
            subject_length += 1;
            subject_sub += '-';
            query_sub += query[j];
            j--;
        }
    }
    DNA_similarity_ = (double)common_subsequence.size() / query_length;
}
int SequenceAlignment::DNASequenceAlignScore(char t, char s) {
    int score = 0;
    int GG = 9;
    int _G = GAP_2;
    int GT = GAP_2;
    if (s == t)
        score = GG;
    else if (s == ' ')
        score = _G;
    else if (t == ' ')
        score = _G;
    else if (s != t)
        score = GT;
    return score;
}
int SequenceAlignment::AASAlignScore (char t, char s) {
    int score = 0;
    int G_ = GAP;
    int index_s = 0;
    int index_t = 0;
    if (s == ' ') {
        score = G_;
    } else if ( t == ' ') {
        score = G_;
    } else if (s == t || s != t) {
        index_s = get_index_of_BLOSUM50(s);
        index_t = get_index_of_BLOSUM50(t);
        score = BLOSUM[index_s][index_t];
    }
    return score;
}
int SequenceAlignment::get_index_of_BLOSUM50 (char s) {
    switch (s) {
        case 'A':
            return 0;
            break;
        case 'R':
            return 1;
            break;
        case 'N':
            return 2;
            break;
        case 'D':
            return 3;
            break;
        case 'C':
            return 4;
            break;
        case 'Q':
            return 5;
            break;
        case 'E':
            return 6;
            break;
        case 'G':
            return 7;
            break;
        case 'H':
            return 8;
            break;
        case 'I':
            return 9;
            break;
        case 'L':
            return 10;
            break;
        case 'K':
            return 11;
            break;
        case 'M':
            return 12;
            break;
        case 'F':
            return 13;
            break;
        case 'P':
            return 14;
            break;
        case 'S':
            return 15;
            break;
        case 'T':
            return 16;
            break;
        case 'W':
            return 17;
            break;
        case 'Y':
            return 18;
            break;
        case 'V':
            return 19;
            break;
        default:
            break;
    }
}
void SequenceAlignment::load_matrix_BLOSUM() {
    ifstream infile;
    //make sure to fill in the correct path;
    infile.open("/Users/jinyang/Documents/iGEM_Programmes/GRN/BLOSUM_50");
    for (int i = 0; i != 20; ++i) {
        for (int j = 0; j != 20; ++j) {
            infile >> BLOSUM[i][j];
        }
    }
    infile.close();
}
double SequenceAlignment::get_max_value (double a, double b, double c) {
    if(a<b)
        a=b;
    if(a<c)
        a=c;
    return a;
}