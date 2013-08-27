////////////////////////////////////////////////////////////////////////////////
/// COPYRIGHT NOTICE\n
/// Distribute under BSD License\n
/// Copyright (c) 2013, iGEM Software Team of University of Science and
/// Technology of China\n
/// All rights reserved.
///
/// \file GRN.h
/// \brief Define the class GRN.
/// \version 1.0
/// \author Li Jinyang
/// \date July 26, 2013
////////////////////////////////////////////////////////////////////////////////

#ifndef __GRN__GRN__
#define __GRN__GRN__
#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>
#include "Sequence.h"
///    Calculate sequence simialrity and construct new GRN.
///
///    1. Get original Gene Regulatory Network matrix.\n
///       Get original GRN matrix from the object of class [FIXME] and add a
///       new row and column in the end to be filled in the new relationship.\n
///    2. Get sequence simialrity.\n
///       Get sequence similarity by sequence alignment using dynamic planning
///       with the substituion matrix BLOSUM_50.\n
///    3. Predict exogenous gene regulatory behavior.\n
///       Using simialrity vector and regulatory vectors predict the behavior of
///       exogenous gene. And fill the correlations in GRN.\n
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
    /// New GRN matrix.
    double **new_GRN;
    /// Initialize the object.
    ///
    ///    \param old_GRN Original GRN.
    ///    \param num_row The numbers of rows of original GRN.
    ///    \param num_column The numbers of column of orginal GRN.
    ///    \see [FIXME]
    void initialize_GRN(double **old_GRN, int num_row, int num_column);
    /// Construct the new GRN with exogenous gene's row and column filled.
    ///
    ///    \param reg_unit The object array of class Sequence. Contains original
    ///                    RU sequences and the query sequences.
    ///    \see Sequence
    void construct_new_GRN(Sequence reg_unit[]);
    /// Align amino aicd sequence.
    ///
    ///    \param query The query amino acid sequence.
    ///    \param query_size The length of query amino acid sequence.
    ///    \param subject The subject amino acid sequence.
    ///    \param subject_size The length of subject amino acid sequence.
    ///    \return Precentage similarity of the two amino acid sequences.
    ///    \see DNASeqAlignment
    double AminoAcidSeqAlignment(std::string query, int query_size,
                                 std::string subject, int subject_size);
    /// Align DNA sequence.
    ///
    ///    \param query The query DNA sequence.
    ///    \param query_size The length of query DNA sequence.
    ///    \param subject The subject DNA sequence.
    ///    \param subject_size The length of subject DNA sequence.
    ///    \return Percentage simialrity of the two DNA sequence.
    ///    \see AminoAcidSeqAlignment
    double DNASeqAlignment(std::string  query, int query_size,
                                std::string subject, int subject_size);
    /// Read the substiturion matrix BLOSUM_50.
    ///
    void load_matrix_BLOSUM();
private:
    /// The number of rows of original GRN.
    int number_row;
    /// The number of columns of original GRN.
    int number_column;
    /// The substiturion matrix.
    int BLOSUM[20][20];
    /// Score a aligment of two amino acids.
    ///
    ///    One amino acid comes from the query sequence. Another comes from the
    ///    subject seqence. The socre will be filled in the socre matrix of
    ///    dynamic planning.
    ///    \param t An amino acid comes from the subject sequence.
    ///    \param s An amino acid comes from the query sequence.
    ///    \return The score of the alignment.
    ///    \see DNASequenceAlignScore
    ///    \see AminoAcidSeqAlignment
    ///    \note The alignment score is dependent on the substitution matrix.
    int AminoAcidSequenceAlignScore(char t, char s);
    /// Score a algnment of two DNAs.
    ///    One DNA comes from the query sequence. Another comes from the subject
    ///    seqence. The socre will be filled in the socre matrix of dynamic
    ///    planning.
    ///    \param t A DNA comes from the subject sequence.
    ///    \param s A DNA comes from the query sequence.
    ///    \return The score of the alignment.
    ///    \see AminoAcidSequenceAlignScore
    ///    \see DNASeqAlignment
    int DNASequenceAlignScore(char t, char s);
    /// Find the biggest value.
    ///
    ///    \return The biggest value of the input.
    double get_max_value (double a, double b, double c);
    /// Get the index of BLOSUM_50.
    ///
    ///    \param s An amino acid.
    ///    \return The index of the amino acid in BLOSUM_50.
    int get_index_of_BLOSUM50(char s);
};

#endif /* defined(__GRN__GRN__) */
