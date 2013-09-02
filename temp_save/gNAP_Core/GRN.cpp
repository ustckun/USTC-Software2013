//
//  GRN.cpp
//  GRN
//
//  Created by jinyang on 13-7-26.
//  Copyright (c) 2013Äê Li Jinyang. All rights reserved.
//

#include "GRN.h"
#include "RandSeq.h"
#include <vector>
#include <string>
#include <fstream>
#include <ctime> 
#include <cmath>
#include <stdlib.h>

#define GAP -8 //AAS Gap penalty;
#define GAP_2 -1 //Promoter sequence gap penalty;
#define RAND_SCALE 100 //Filtering module: the number of random sequence;
#define SIGMA_NUM 0.2//Filtering module: numbers of sigma;

void GRN::initialize_GRN(double **old_GRN, int num_row, int num_column) {
    srand((unsigned)time(0));
    for (int i = 0; i != 1800; ++i) {
        for (int j = 0; j != 220; ++j) {
            //Chooses direction of regulation if it is +/-(2);
            if (old_GRN[i][j] == 2) {
                if (((rand() % 100) / 100.0) < 0.5) {
                    new_GRN[i][j] = 0;
                } else {
                    new_GRN[i][j] = 0;
                }
            }else{
                new_GRN[i][j] = old_GRN[i][j];
            }
        }
    }
    number_row = num_row;
    number_column = num_column;
}
void GRN::construct_new_GRN(Sequence reg_unit[]) {
    double promoter_similarity[1800] = { 0 };
    double gene_similarity[220] = { 0 };
    load_matrix_BLOSUM();
    // Generates promoter similarity matrix.
    for (int RU_number = 0; RU_number != number_row; ++RU_number) {
        promoter_similarity[RU_number] =
        DNASeqAlignment(reg_unit[number_row].promoter_sequence,
                        reg_unit[number_row].promoter_size,
                        reg_unit[RU_number].promoter_sequence,
                        reg_unit[RU_number].promoter_size);
        std::cout << "promoter simi: " << RU_number << '\t' <<
        promoter_similarity[RU_number] << std::endl;
    }
    // Generates gene similarity matrix.
    for (int RU_number = 0; RU_number != number_column; ++RU_number) {
        gene_similarity[RU_number] =
        AminoAcidSeqAlignment(
                              reg_unit[number_row].amino_acid_sequence,
                              reg_unit[number_row].amino_acid_sequence_size,
                              reg_unit[RU_number].amino_acid_sequence,
                              reg_unit[RU_number].amino_acid_sequence_size);
        std::cout << "gene simi: " << RU_number  << '\t' <<
        gene_similarity[RU_number] << std::endl;
    }
    // Filters promoter similarity matrix.
    /*for (int RU_number = 0; RU_number != number_row; ++RU_number) {
        if (promoter_similarity[RU_number] < 0.53) {
            promoter_similarity[RU_number] = 0;
        }
    }*/
    // Filters gene similarity matrix.
    for (int RU_number = 0; RU_number != number_column; ++RU_number) {
        RandomSequence random_sequence[RAND_SCALE];
        double sizer = 0;
        double sigma = 0;
        double random_similarity[RAND_SCALE] = { 0 };
        // Initialize random sizer at the same size as query sequence.
        for (int i = 0; i != RAND_SCALE; ++i) {
            random_sequence[i].generate_random_amino_acid_sequence(
                reg_unit[number_column].amino_acid_sequence_size);
            random_similarity[i] =
            AminoAcidSeqAlignment(random_sequence[i].random_amino_acid_sequence,
                             reg_unit[number_column].amino_acid_sequence_size,
                             reg_unit[RU_number].amino_acid_sequence,
                             reg_unit[RU_number].amino_acid_sequence_size);
            sizer += random_similarity[i];
        }
        sizer = sizer / RAND_SCALE; 
        for (int i = 0; i != RAND_SCALE; ++i) {
            sigma += pow(random_similarity[i]- sizer, 2);
        }
        // Standard deviation.
        sigma = pow(sigma, 0.5);
        // Set (SIGMA_NUM * sigma) as threshold.
        sizer = sizer + sigma * SIGMA_NUM;
        if (gene_similarity[RU_number] < sizer) {
            gene_similarity[RU_number] = 0;
        }
        std::cout << "filtering: " << RU_number << std::endl;
    }
    // Insert new correlations to (number_row + 1) row.
    for (int j_RU_number = 0; j_RU_number != number_column; ++j_RU_number) {
        int counter_positive = 0;
        int counter_negative = 0;
        double average_similarity_positive = 0;
        double regulation_positive = 0;
        double average_similarity_negative = 0;
        double regulation_negative = 0;
        for (int i_RU_number = 0; i_RU_number != number_row; ++i_RU_number) {
            if (reg_unit[i_RU_number].promoter_sequence ==
                reg_unit[number_row].promoter_sequence) {
                if (new_GRN[i_RU_number][j_RU_number] == 1) {
                    average_similarity_positive +=
                        promoter_similarity[i_RU_number];
                    regulation_positive += new_GRN[i_RU_number][j_RU_number] *
                        promoter_similarity[i_RU_number];
                    if (promoter_similarity[i_RU_number] != 0) {
                        counter_positive += 1;
                    }
                } else if (new_GRN[i_RU_number][j_RU_number] == -1) {
                    average_similarity_negative +=
                        promoter_similarity[i_RU_number];
                    regulation_negative += new_GRN[i_RU_number][j_RU_number] *
                        promoter_similarity[i_RU_number];
                    if (promoter_similarity[i_RU_number] != 0) {
                        counter_negative += 1;
                    }
                }
            }
        }
        if (counter_positive != 0) {
            average_similarity_positive = average_similarity_positive /
                                          counter_positive;
            regulation_positive = regulation_positive / counter_positive;
        } else {
            average_similarity_positive = 0;
            regulation_positive = 0;
        }
        if (counter_negative != 0) {
            average_similarity_negative = average_similarity_negative /
                                          counter_negative;
            regulation_negative = regulation_negative / counter_negative;
        } else {
            average_similarity_negative = 0;
            regulation_negative = 0;
        }
        if (average_similarity_positive > average_similarity_negative) {
            new_GRN[number_row][j_RU_number] = regulation_positive;
        } else if (average_similarity_positive < average_similarity_negative) {
            new_GRN[number_row][j_RU_number] = regulation_negative;
        } else {// If equal, chooses random;
            switch (rand() % 2) {
                case 0:
                    new_GRN[number_row][j_RU_number] = regulation_positive;
                    break;
                case 1:
                    new_GRN[number_row][j_RU_number] = regulation_negative;
                    break;
            }
        }
        std::cout << j_RU_number << std::endl;
    }
    // Insert new correlations to (number_cloumn + 1) column;
    for (int i_RU_number = 0; i_RU_number != number_row; ++i_RU_number) {
        int counter_positive = 0;
        int counter_negative = 0;
        double average_similarity_positive = 0;
        double regulation_positive = 0;
        double average_similarity_negative = 0;
        double regulation_negative = 0;
        for (int j_RU_number = 0; j_RU_number != number_column; ++j_RU_number) {
            if (new_GRN[i_RU_number][j_RU_number] == 1) {
                average_similarity_positive += gene_similarity[j_RU_number];
                regulation_positive += new_GRN[i_RU_number][j_RU_number] *
                gene_similarity[j_RU_number];
                if (gene_similarity[j_RU_number] != 0) {
                    counter_positive += 1;
                }
            }else if (new_GRN[i_RU_number][j_RU_number] == -1){
                average_similarity_negative += gene_similarity[j_RU_number];
                regulation_negative += new_GRN[i_RU_number][j_RU_number] *
                gene_similarity[j_RU_number];
                if (gene_similarity[j_RU_number] != 0) {
                    counter_negative += 1;
                }
            }
        }
        if (counter_positive != 0) {
            average_similarity_positive = average_similarity_positive /
            counter_positive;
            regulation_positive = regulation_positive / counter_positive;
        } else {
            average_similarity_positive = 0;
            regulation_positive = 0;
        }
        if (counter_negative != 0) {
            average_similarity_negative = average_similarity_negative /
            counter_negative;
            regulation_negative = regulation_negative / counter_negative;
        } else {
            average_similarity_negative = 0;
            regulation_negative = 0;
        }
        if (average_similarity_positive > average_similarity_negative) {
            new_GRN[i_RU_number][number_column] = regulation_positive;
        } else if (average_similarity_positive < average_similarity_negative) {
            new_GRN[i_RU_number][number_column] = regulation_negative;
        } else {//If equal, chooses random;
            switch (rand() % 2) {
                case 0:
                    new_GRN[i_RU_number][number_column] = regulation_positive;
                    break;
                case 1:
                    new_GRN[i_RU_number][number_column] = regulation_negative;
                    break;
            }
        }
        std::cout << i_RU_number << std::endl;
    }
    new_GRN[number_row][number_column] = 0;
    // TODO(jinyangustc@gmail.com): If can't get ideal result, wirte a module to
    //                              filter low correlations in the column.
    //use system time as file name;
    std::ofstream outfile;
    time_t nowtime = time(NULL);
    struct tm *p;
    p = gmtime(&nowtime);
    char filename_1[256] ={ 0 };
    std::string filename = "/Users/jinyang/Desktop/Parameter Data Test/";
    sprintf(filename_1, "%d-%d %d%02d", 1 + p -> tm_mon, p -> tm_mday,
                                        8 + p -> tm_hour, p -> tm_min);
    filename += filename_1;
    filename += " p_simi & g_simi & para.txt";
    outfile.open(filename);
    for (int i = 0; i != 1800; ++i) {
        outfile << i + 1 << '\t' << promoter_similarity[i] << std::endl;
    }
    outfile << std::endl;
    for (int i = 0; i != 220; ++i) {
        outfile << i + 1 << '\t' << gene_similarity[i] << std::endl;
    }
    outfile << std::endl;
    outfile << "Gap: " << GAP << std::endl;
    outfile << "Numbers of random sequence: " << RAND_SCALE << std::endl;
    outfile << "Sigma control: " << SIGMA_NUM << std::endl;
    outfile.close();
}
double GRN::AminoAcidSeqAlignment(std::string query, int query_size,
                                  std::string subject, int subject_size) {
    double G_ = GAP;
    double _G = GAP;
    int query_length = 0;
    int subject_length = 0;
    double positives = 0;
    query_length = query_size;
    subject_length = subject_size;
    std::string common_subsequence;
    std::string query_sub;
    std::string subject_sub;
    std::vector< std::vector<double> > align_matrix
        (subject_size + 1, std::vector<double>(query_size + 1));
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
                 AminoAcidSequenceAlignScore(subject[i], query[j]);
            _t = align_matrix[i - 1][j] +
                 AminoAcidSequenceAlignScore(subject[i], ' ');
            s_ = align_matrix[i][j - 1] +
                 AminoAcidSequenceAlignScore(' ', query[j]);
            align_matrix[i][j] = get_max_value(st, s_, _t);
        }
    }
    for (int i = subject_size, j = query_size; (i > 0 || j > 0);) {
            if (
                (i * j > 0) && (
                                align_matrix[i][j] ==
                                align_matrix[i - 1][j - 1] +
                                AminoAcidSequenceAlignScore(subject[i], query[j]
                                                            )
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
                      AminoAcidSequenceAlignScore(subject[i], ' '))) {
                //insert a space into string_s;
                query_length += 1;
                query_sub += '-';
                subject_sub += subject[i];
                i--;
            } else if ((j > 0) &&
                      (align_matrix[i][j] == align_matrix [i][j - 1] +
                       AminoAcidSequenceAlignScore(' ', query[j] ))){
                //insert a space into string_s;
                subject_length += 1;
                subject_sub += '-';
                query_sub += query[j];
                j--;
            }
    }
    positives = (double)common_subsequence.size() / query_length;
    return positives;
}
double GRN::DNASeqAlignment(std::string  query, int query_size,
                            std::string subject, int subject_size) {
    double G_ = GAP_2;
    double _G = GAP_2;
    int query_length = 0;
    int subject_length = 0;
    double positives = 0;
    query_length = query_size;
    subject_length = subject_size;
    std::string common_subsequence;
    // Records result of alignment.
    std::string query_sub;
    std::string subject_sub;
    std::vector< std::vector<double> > align_matrix
        (subject_size + 1, std::vector<double>(query_size + 1));
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
    positives = (double)common_subsequence.size() / query_length;
    return positives;
}
int GRN::DNASequenceAlignScore(char t, char s) {
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
int GRN::AminoAcidSequenceAlignScore (char t, char s) {
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
int GRN::get_index_of_BLOSUM50 (char s) {
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
void GRN::load_matrix_BLOSUM() {
    std::ifstream infile;
    //make sure to fill in the correct path;
    infile.open("BLOSUM_50");
    for (int i = 0; i != 20; ++i) {
        for (int j = 0; j != 20; ++j) {
            infile >> BLOSUM[i][j];
        }
    }
    infile.close();
}
double GRN::get_max_value (double a, double b, double c) {
    if(a<b)
        a=b;
    if(a<c)
        a=c;
    return a;
}