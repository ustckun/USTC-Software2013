////////////////////////////////////////////////////////////////////////////////
/// COPYRIGHT NOTICE\n
/// Distribute under BSD License\n
/// Copyright (c) 2013, iGEM Software Team of University of Science and
/// Technology of China\n
/// All rights reserved.
///
/// \file GeneIM.h
/// \brief Define the class GeneIM.
/// \version 1.0
/// \author Wang Chenkun
/// \date September 2nd, 2013
////////////////////////////////////////////////////////////////////////////////

#ifndef __USTC_Software__GeneIM__
#define __USTC_Software__GeneIM__

#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <map>

#include"define.h"

using namespace std;

///    A class which contain one gene information such as gene name, gene position,
///    gene ID, gene sequence, gene description, promoter name and promoter sequence.
///
///    1. Get gene information.\n
///       Get gene information in map of gene info and input them into corresponding
///       variable.\n
///    2. Get promoter information.\n
///       Get promoter information in map of promoter which is constructed by gene
///       position in TU.
///    3. Find RNA gene.\n
///       Some genes are not expressed to amino acid but RNA or tRNA. Avoid of aligning
///       the AAS of RNA sequence, we find out the RNA gene.\n

class GeneIM
{
public:
    GeneIM()
    {
        left_position=0;
        right_position=0;
    }
/// The number of gene in GRN
    int gene_number;
/// Name used to store name in TF-TF file temporary
    char *name;
/// Get gene information from map of gene info constructed in class GetReady
///    \param map of gene info
///    \see GetReady
    void getGeneInformation(map<string,string> dict);
/// Get promoter information from map of promoter sequence also constructed in class
/// GetReady
///    \param Map of promoter sequence
///    \see GetReady
    void getPromoterIF(map<string,string> dict);
/// Get ID of gene in RegulonDB
///    \return ID of gene
    string getID();
/// Get gene sequence
///    \return gene sequence
    string getGeneSequence();
/// Get promoter sequence
///    \return promoter sequence
    string getPromoterSequence();
/// Get promoter name
///    \return promoter name
    string getPromoterName();
/// Get gene name which distinct capital and small letter
///    \return Gene name
    string getGeneTrueName();
/// Get gene left position which mean the position of gene
///    \return left position of gene
    int getLeftPosition();
/// Get gene right position
///    \return right position of gene
    int getRightPosition();
/// Put gene name into gene_name[10]
    void putName();
/// Put promoter name into promoter_name
    void putInPromoterName(string promoter);
/// Get those genes which are not expressed into amino acid but RNA
    int getRNA();
/// Get gene name for finding
///    \return Gene name
    char *getGeneName();
/// Get gene description
///    \return gene description
    string getGeneDescription();
/// contain gene name which not distinct capital and small letter
///    It is used to find the right gene in map
    char gene_name[10];
private:
/// ID in RegulonDB
    string iD;
/// Gene sequence
    string gene_sequence;
/// Gene left position
    int left_position;
/// Gene right position
    int right_position;
/// Gene description which contains the gene expression products
    string gene_description;
/// Promoter name
    string promoter_name;
/// Promoter sequence
    string promoter_sequence;
/// Gene name which distinct capital and small letter
    string true_name;
/// Represent to RNA or not by yes(1), no(0)
    int RNA;
};

#endif
