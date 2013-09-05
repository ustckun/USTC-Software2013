////////////////////////////////////////////////////////////////////////////////
/// COPYRIGHT NOTICE\n
/// Distribute under BSD License\n
/// Copyright (c) 2013, iGEM Software Team of University of Science and
/// Technology of China\n
/// All rights reserved.
///
/// \file GRN.cpp
/// \brief Statements of functions of the class GRN.
/// \version 1.0
/// \author Li Jinyang
/// \date July 26, 2013
////////////////////////////////////////////////////////////////////////////////

#include "GRN.h"
#include "RandSeq.h"
#include <vector>
#include <string>
#include <fstream>
#include <ctime>
#include <cmath>
#include <stdlib.h>
/// Linear gap penalty of amino acid sequence alignment.
#define GAP -8
/// Linear gap penalty of DNA sequence alignment.
#define GAP_2 -1
/// The number of generated random sequences.
#define RAND_SCALE 100
/// Filter control determins the range of similary to be filtered.
#define SIGMA_NUM 0.2

void GRN::initialize_GRN(double **old_GRN, int num_row, int num_column) {
    srand((unsigned)time(0));
// FIXME(jinyangustc@gmail.com): Delete random selection lines when finish
//                              environment parameter input function.
    for (int i = 0; i != 1800; ++i) {
        for (int j = 0; j != 220; ++j) {
            // Chooses direction of regulation if it is +/-(2);
            if (old_GRN[i][j] == 2) {
                if (((rand() % 100) / 100.0) < 0.5) {
                    new_GRN[i][j] = 1;
                } else {
                    new_GRN[i][j] = -1;
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
    std::vector<int> same_promoter_position;
    int counter_same_promoter = 0;
    // Align query RU DNA sequence with all original RU DNA sequences and get
    // promoter similarities.
    for (int RU_number = 0; RU_number != number_row; ++RU_number) {
        promoter_similarity[RU_number] =
        DNASeqAlignment(reg_unit[number_row].promoter_sequence,
                        reg_unit[number_row].promoter_size,
                        reg_unit[RU_number].promoter_sequence,
                        reg_unit[RU_number].promoter_size);
        // Record the position of the RU which has the same promoter as the
        // query one
        if (promoter_similarity[RU_number] == 1) {
            same_promoter_position.push_back(RU_number);
            //std::cout << same_promoter_position[0] << std::endl;
            counter_same_promoter += 1;
        }
        //std::cout << "promoter simi: " << RU_number << '\t' <<
        //promoter_similarity[RU_number] << std::endl;
    }
    // Align query RU amino acid sequence with all orginal RU amino acid
    // sequences and get gene similarities.
    for (int RU_number = 0; RU_number != number_column; ++RU_number) {
        gene_similarity[RU_number] =
        AminoAcidSeqAlignment(
                              reg_unit[number_row].amino_acid_sequence,
                              reg_unit[number_row].amino_acid_sequence_size,
                              reg_unit[RU_number].amino_acid_sequence,
                              reg_unit[RU_number].amino_acid_sequence_size);
        //std::cout << "gene simi: " << RU_number  << '\t' <<
        //gene_similarity[RU_number] << std::endl;
    }
    // Filter gene similarities.
    //
    //    Generate 100 random sequences at the same
    //    length of the query gene sequence. Get the random simialrity
    //    distribution by aligning the subject sequence and the random
    //    sequences. Analyse the distribution and filter the gene similarities
    //    with standard deviation.
    for (int RU_number = 0; RU_number != number_column; ++RU_number) {
        RandomSequence random_sequence[RAND_SCALE];
        double sizer = 0;
        double sigma = 0;
        double random_similarity[RAND_SCALE] = { 0 };
        // Generated random sequences at the same length of query sequence.
        // Align random sequences with the subject sequence which is under
        // alignment with the query one.
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
        // Calculate standard deviation.
        sigma = pow(sigma, 0.5);
        // Set (SIGMA_NUM * sigma) as threshold.
        sizer = sizer + sigma * SIGMA_NUM;
        if (gene_similarity[RU_number] < sizer) {
            gene_similarity[RU_number] = 0;
        }
        //std::cout << "filtering: " << RU_number << std::endl;
    }
    // Predict new gene's regulated behavior and construct a new row in GRN.
    //
    //    There are two conditions:
    //      (1) There are some RUs in original GRN having the same promoter as
    //         the query one.
    //      (2) There is no RU in original GRN having the same promoter as the
    //         query one.
    if (counter_same_promoter != 0) {  // Condition (1).
        // Traverses all regulator RUs.
        for (int j_RU_number = 0; j_RU_number != number_column; ++j_RU_number) {
            int counter_positive = 0;
            int counter_negative = 0;
            double average_similarity_positive = 0;
            // Sum of all positive regulations.
            double regulation_positive = 0;
            double average_similarity_negative = 0;
            // Sum of all negative regulations.
            double regulation_negative = 0;
            // Traverses target RUs which have the same promoter as query one.
            for (int i_same_promoter = 0;
                 i_same_promoter != counter_same_promoter; ++i_same_promoter) {
                double promoter_followed_gene_similarity = 0;
                promoter_followed_gene_similarity =
                DNASeqAlignment(reg_unit[number_row].gene_sequence,
                                reg_unit[number_row].gene_size,
                                reg_unit[
                                         same_promoter_position[i_same_promoter]
                                         ].gene_sequence,
                                reg_unit[
                                         same_promoter_position[i_same_promoter]
                                         ].gene_size);
                // Counts and calculates new regulated relationship
                // in different types.
                if (new_GRN[same_promoter_position[i_same_promoter]]
                           [j_RU_number] == 1) {
                    average_similarity_positive +=
                        promoter_followed_gene_similarity;
                    regulation_positive +=
                        new_GRN[same_promoter_position[i_same_promoter]]
                               [j_RU_number] *
                        promoter_followed_gene_similarity;
                    counter_positive += 1;
                } else if (new_GRN[same_promoter_position[i_same_promoter]]
                           [j_RU_number] == -1) {
                    average_similarity_negative +=
                        promoter_followed_gene_similarity;
                    regulation_negative +=
                        new_GRN[same_promoter_position[i_same_promoter]]
                               [j_RU_number] *
                        promoter_followed_gene_similarity;
                    counter_negative += 1;
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
            // Chooses the regulated direction.
            if (average_similarity_positive > average_similarity_negative) {
                new_GRN[number_row][j_RU_number] = regulation_positive;
            } else if (average_similarity_positive <
                       average_similarity_negative) {
                new_GRN[number_row][j_RU_number] = regulation_negative;
            } else {  // If average similarities are equal, chooses randomly.
                switch (rand() % 2) {
                    case 0:
                        new_GRN[number_row][j_RU_number] =regulation_positive;
                        break;
                    case 1:
                        new_GRN[number_row][j_RU_number] = regulation_negative;
                        break;
                }
            }
        }
    } else {  // Condition (2).
        // Traverses all regulator RUs.
        for (int j_RU_number = 0; j_RU_number != number_column; ++j_RU_number) {
            int counter_positive = 0;
            int counter_negative = 0;
            double average_similarity_positive = 0;
            // Sum of positive regulations.
            double regulation_positive = 0;
            double average_similarity_negative = 0;
            // Sum of negtive regulations.
            double regulation_negative = 0;
            // Traverses all target RUs.
            for (int i_RU_number = 0; i_RU_number != number_row;
                 ++i_RU_number) {
                // Counts and calculates new regulated relationships in
                // different types.
                if (new_GRN[i_RU_number][j_RU_number] == 1) {
                    average_similarity_positive += gene_similarity[i_RU_number];
                    regulation_positive += new_GRN[i_RU_number][j_RU_number] *
                    gene_similarity[i_RU_number];
                    if (gene_similarity[i_RU_number] != 0) {
                        counter_positive += 1;
                    }
                } else if (new_GRN[i_RU_number][j_RU_number] == -1){
                    average_similarity_negative += gene_similarity[i_RU_number];
                    regulation_negative += new_GRN[i_RU_number][j_RU_number] *
                    gene_similarity[i_RU_number];
                    if (gene_similarity[i_RU_number] != 0) {
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
            // Chooses the regulated direction.
            if (average_similarity_positive > average_similarity_negative) {
                new_GRN[number_row][j_RU_number] = regulation_positive;
            } else if (average_similarity_positive <
                       average_similarity_negative) {
                new_GRN[number_row][j_RU_number] = regulation_negative;
            } else {//If average similarities are equal, chooses randomly;
                switch (rand() % 2) {
                    case 0:
                        new_GRN[number_row][j_RU_number] = regulation_positive;
                        break;
                    case 1:
                        new_GRN[number_row][j_RU_number] = regulation_negative;
                        break;
                }
            }
        }
    }
    // Predicts new RU's regulating behavior and constructs a new column in GRN.
    //
    //    Uses the same method as the one used to prediction of new RU's
    //    regulated behavior in condition (2).
    //
    // Traverses all target RUs.
    for (int i_RU_number = 0; i_RU_number != number_row; ++i_RU_number) {
        int counter_positive = 0;
        int counter_negative = 0;
        double average_similarity_positive = 0;
        double regulation_positive = 0;
        double average_similarity_negative = 0;
        double regulation_negative = 0;
        // Traverses all regulator RUs.
        for (int j_RU_number = 0; j_RU_number != number_column; ++j_RU_number) {
            // Counts and calculates new regualting relationships in different
            // types.
            if (new_GRN[i_RU_number][j_RU_number] == 1) {
                average_similarity_positive += gene_similarity[j_RU_number];
                regulation_positive += new_GRN[i_RU_number][j_RU_number] *
                gene_similarity[j_RU_number];
                if (gene_similarity[j_RU_number] != 0) {
                    counter_positive += 1;
                }
            } else if (new_GRN[i_RU_number][j_RU_number] == -1){
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
        // Chooses regualting direction.
        if (average_similarity_positive > average_similarity_negative) {
            new_GRN[i_RU_number][number_column] = regulation_positive;
        } else if (average_similarity_positive < average_similarity_negative) {
            new_GRN[i_RU_number][number_column] = regulation_negative;
        } else {//If average similarities are equal, chooses randomly;
            switch (rand() % 2) {
                case 0:
                    new_GRN[i_RU_number][number_column] = regulation_positive;
                    break;
                case 1:
                    new_GRN[i_RU_number][number_column] = regulation_negative;
                    break;
            }
        }
        //std::cout << i_RU_number << std::endl;
    }
    new_GRN[number_row][number_column] = 0;
    // TODO(jinyangustc@gmail.com): If can't get ideal result, wirte a module to
    //                              filter low correlations in the column.
    //use system time as file name;
    /*std::ofstream outfile;
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
    outfile.close();*/
}
double GRN::AminoAcidSeqAlignment(std::string query, int query_size,
                                  std::string subject, int subject_size) {
    // Linear gap penalty when inserts a space in query sequence.
    double G_ = GAP;
    // Linear gep penalty when inserts a space in subject sequence.
    double _G = GAP;
    // Records the length of query sequence after alignment.
    int query_length = 0;
    // Records the length of subjuect sequence after alignment.
    int subject_length = 0;
    // The number of matches and similar matches in alignment. When it's divided
    // by query_length or subject_length, it means the percentage similariy.
    double positives = 0;
    query_length = query_size;
    subject_length = subject_size;
    // Records longest common subsequence.
    std::string common_subsequence;
    // Records query sequence after alignment.
    std::string query_sub;
    // Records subject sequence afer alignemt.
    std::string subject_sub;
    // Score matrix of dynamic planning.
    //
    //    Sets query sequence to be X axis.
    //    Sets subject sequence to be Y axis.
    std::vector< std::vector<double> > align_matrix
        (subject_size + 1, std::vector<double>(query_size + 1));
    // Initialize score matrix.
    //
    //    Sets (0, 0) to be 0.
    //    (i, 0) = (i - 1, 0) + gap penalty.
    //    (0, j) = (0, j - 1) + gap penalty.
    align_matrix[0][0] = 0;
    for (int i = 1; i != subject_size + 1; ++i) {
        align_matrix[i][0] = align_matrix[i - 1][0] + G_;
    }
    for (int j = 1; j != query_size + 1; ++j) {
        align_matrix[0][j] = align_matrix[0][j - 1] + _G;
    }
    // Dynamic planning.
    //
    //    There are 3 possible paths moving to (i, j) element in score matrix:
    //    (1) Diagonal:
    //      (i, j) = (i - 1, j - 1) + (subject[i] vs. query[j]).
    //      It means a match or a similar match according to BLOSUM_50.
    //    (2) Vertical:
    //      (i, j) = (i - 1, j) + (subject[i] vs. a space in query[j]).
    //      It means inserting a space in query sequence.
    //    (3) Horizontal:
    //      (i, j) = (i, j - 1) + (a space in subject[i] vs. query[j]).
    //      It means inserting a space in subject sequence.
    //    Chooses the maximum of the three as the actual path.
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
    // Trace back.
    //
    //    Finds the longest common sequence and records the length of query
    //    sequence and the length of subject sequence after alignment.
    //    The trace begins at the last element in the score matrix and ends at
    //    the first element in the score matrix.
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
                // Insert a space in the query sequence.
                query_length += 1;
                query_sub += '-';
                subject_sub += subject[i];
                i--;
            } else if ((j > 0) &&
                      (align_matrix[i][j] == align_matrix [i][j - 1] +
                       AminoAcidSequenceAlignScore(' ', query[j] ))){
                // Insert a space in the subject sequence.
                subject_length += 1;
                subject_sub += '-';
                query_sub += query[j];
                j--;
            }
    }
    // Calculates the percentage simialrity.
    positives = (double)common_subsequence.size() / query_length;
    return positives;
}
double GRN::DNASeqAlignment(std::string  query, int query_size,
                            std::string subject, int subject_size) {
    // Linear gap penalty when inserts a space in query sequence.
    double G_ = GAP_2;
    // Linear gap penalty when inserts a space in subject sequence.
    double _G = GAP_2;
    // Records the length of query sequence after alignment.
    int query_length = 0;
    // Records the length of subject sequence after aligment.
    int subject_length = 0;
    // The number of matches and similar matches in alignment. When it's divided
    // by query_length or subject_length, it means the percentage similariy.
    double positives = 0;
    query_length = query_size;
    subject_length = subject_size;
    // Records longest common subsequence.
    std::string common_subsequence;
    // Records query sequence after alignment.
    std::string query_sub;
    // Records subject sequence after alignment.
    std::string subject_sub;
    // Score matrix of dynamic planning.
    //
    //    Sets query sequence to be X axis.
    //    Sets subject sequence to be Y axis.
    std::vector< std::vector<double> > align_matrix
        (subject_size + 1, std::vector<double>(query_size + 1));
    // Initialize score matrix.
    //
    //    Sets (0, 0) to be 0.
    //    (i, 0) = (i - 1, 0) + gap penalty.
    //    (0, j) = (0, j - 1) + gap penalty.
    align_matrix[0][0] = 0;
    for (int i = 1; i != subject_size + 1; ++i) {
        align_matrix[i][0] = align_matrix[i - 1][0] + G_;
    }
    for (int j = 1; j != query_size + 1; ++j) {
        align_matrix[0][j] = align_matrix[0][j - 1] + _G;
    }
    // Dynamic planning.
    //
    //    There are 3 possible paths moving to (i, j) element in score matrix:
    //    (1) Diagonal:
    //      (i, j) = (i - 1, j - 1) + (subject[i] vs. query[j]).
    //      It means a match or a similar match according to BLOSUM_50.
    //    (2) Vertical:
    //      (i, j) = (i - 1, j) + (subject[i] vs. a space in query[j]).
    //      It means inserting a space in query sequence.
    //    (3) Horizontal:
    //      (i, j) = (i, j - 1) + (a space in subject[i] vs. query[j]).
    //      It means inserting a space in subject sequence.
    //    Chooses the maximum of the three as the actual path.
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
    // Trace back.
    //
    //    Finds the longest common sequence and records the length of query
    //    sequence and the length of subject sequence after alignment.
    //    The trace begins at the last element in the score matrix and ends at
    //    the first element in the score matrix.
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
            // Insert a space in query sequence.
            query_length += 1;
            query_sub += '-';
            subject_sub += subject[i];
            i--;
        } else if ((j > 0) && (align_matrix[i][j] == align_matrix[i][j - 1] +
                               DNASequenceAlignScore(' ', query[j]))) {
            // Insert a space in subject sequence.
            subject_length += 1;
            subject_sub += '-';
            query_sub += query[j];
            j--;
        }
    }
    // Calculates the percentage simialrity.
    positives = (double)common_subsequence.size() / query_length;
    return positives;
}
int GRN::DNASequenceAlignScore(char t, char s) {
    int score = 0;
    // Match score.
    int GG = 9;
    // Linear gap penalty.
    int _G = GAP_2;
    // Unmatch score.
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
    // Linear gap penalty.
    int G_ = GAP;
    // The index in BLOSUM_50 of character s.
    int index_s = 0;
    // The index in BOLSUM_50 of character t.
    int index_t = 0;
    if (s == ' ') {
        score = G_;
    } else if ( t == ' ') {
        score = G_;
    } else if (s == t || s != t) {
        index_s = get_index_of_BLOSUM50(s);
        index_t = get_index_of_BLOSUM50(t);
        // Looks up socre in BLOSUM_50.
        score = BLOSUM[index_s][index_t];
    }
    return score;
}
int GRN::get_index_of_BLOSUM50 (char s) {
    switch (s) {  // Every case represents a kind of amino acid.
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
    // Make sure to fill in the correct path;
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
