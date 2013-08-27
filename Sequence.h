////////////////////////////////////////////////////////////////////////////////
/// COPYRIGHT NOTICE\n
/// Distribute under BSD License\n
/// Copyright (c) 2013, iGEM Software Team of University of Science and
/// Technology of China\n
/// All rights reserved.
///
/// \file Sequence.h
/// \brief Define the class Sequence.
///
///    Achieve the construction of a regulation unit. The object contains a
///    a promoter sequence, the length of the promoter sequence, a protein
///    coding sequence, the length of the promoter sequence, the amino acid
///    sequence translated from the protein coding sequence, and the length of
///    the amino acid sequence. The regulation unit is identified by a number
///    which is alse stored in the object.
/// \version 1.0
/// \author Li Jinyang
/// \date July 26, 2013
////////////////////////////////////////////////////////////////////////////////
#ifndef __GRN__Sequence__
#define __GRN__Sequence__

#include <iostream>
///    Store promoter and protein coding sequence and construct regulation unit.
///
///    An object of class Sequence is a "regualtion unit". It contains a
///    promoter sequence, a protein coding sequence, the corresponding amino
///    acid sequence,and thier lengths. An RU is identified by a number which is
///    also stored in the object.
class Sequence{
public:
    Sequence()
    {
        regulation_unit_number = 0;
        gene_size = 0;
        amino_acid_sequence_size = 0;
        gene_sequence = " ";
        promoter_sequence = " ";
        amino_acid_sequence = " ";
        //RNASequence = " ";
    }
    /// The protein coding sequence(DNA) of an regulation unit(RU).
    std::string gene_sequence;
    /// The promoter sequence of the regulation unit(RU).
    std::string promoter_sequence;
    /// The translation product(amino acid sequence) of the RU.
    std::string amino_acid_sequence;
    /// Number of the RU.
    int regulation_unit_number;
    /// The length of protein coding DNA sequence.
    int gene_size;
    /// The length of promoter sequence.
    int promoter_size;
    /// The length of amino acid sequence.
    int amino_acid_sequence_size;
    /// Initializes an object.
    ///
    ///    Initialize an object and translates the gene sequence into amino acid
    ///    sequence. Get the length of the amino acid sequence.
    ///    \param RU_number The number of the RU.
    ///    \param promoter The promoter sequence of the RU.
    ///    \param p_size The length of the promoter sequence.
    ///    \param gene The protein coding sequence.
    ///    \param g_size The length of the protein coding sequence.
    ///    \see GRN
    void initialize_Sequence(int RU_number, std::string promoter, int p_size,
                            std::string gene, int g_size);
    /// Translates gene sequence into amino acid sequence.
    ///
    ///    Some explain of the transcription and translation process:\n
    ///    1. Actually, the protein expression process is:\n
    ///       DNA --> mRNA (i.e. transcription);\
    ///       mRNA --> protein (i.e. translation).\n
    ///    2.DNA has double strands, but only one strand takes part in
    ///      transcription.\n
    ///    3.Codons are the sequence messages carried by mRNA;\n
    ///     Take initiation codon "AUG" for example:\n
    ///     ---ATG--- : DNA strand which doesn't take part in transcription
    ///                 process;\n
    ///     ---TAC--- : DNA strand which exactlly takes part in transcription
    ///                 proess;\n
    ///     ---AUG--- : mRNA strand which carries codons. In this case, it
    ///                 carries initiation codon, i.e. "AUG";\n
    ///     ---TAC--- : tRNA which also carries amino acid Methionine(M);\n
    ///   4.Owing to the DNA sequences that our database provided are the
    ///     UNEXPRESSION strands, the translation process of the program can
    ///     just use DNA sequence without the simulation of transcription
    ///     process.
    void Translation();
private:
    int Translate(char s);
};

#endif /* defined(__GRN__Sequence__) */
