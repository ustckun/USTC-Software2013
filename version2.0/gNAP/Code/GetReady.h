////////////////////////////////////////////////////////////////////////////////
/// COPYRIGHT NOTICE\n
/// Distribute under BSD License\n
/// Copyright (c) 2013, iGEM Software Team of University of Science and
/// Technology of China\n
/// All rights reserved.
///
/// \file GetReady.h
/// \brief Define the class GetReady.
/// \version 1.0
/// \author Wang Chenkun
/// \date September 2nd, 2013
////////////////////////////////////////////////////////////////////////////////

#ifndef __USTC_Software__Regulation__
#define __USTC_Software__Regulaiton__


#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <map>

#include <cstring>

#include"define.h"
#include"strlwr.h"

using namespace std;

///    Input the database files downloaded and get ready to fullfill all 
///    information needed to calculate.
///
///    1. Get gene regulatory network.\n
///       Get gene regulatory network from gene to gene interaction files like
///       TF-TF and TF-Gene database on RegulonDB. Build a matrix which contains
///       active(1),repressive(-1),uncertain(2),unknow or no interaction(0).\n
///    2. Map genes' and promoters' information.\n
///       Use map function to build a map of genes' information and promoters' 
///       detail preparing for getting sequence, position and so on.\n
///    3. Ensure the uncertain genes.\n
///       When build regulatory matrix, there will be some uncertain interactions 
///       such as "ada->ada" which has both active and repressive interaction based
///       on the enviroument outside. An "uncertain" is output for users to make 
///       sure those uncertain interactions as needed.\n

class GetReady
{
public:
    GetReady();
    /// Build GRN matrix
    ///
    ///    \param array of GeneIM's objects
    ///    \param file address of TF-TF file
    ///    \param file address of TF-Gene file
    ///    \see GeneIM
    ///    \see readTFTF
    ///    \see readTFGene
    ///    \see addTF
    void getRegulationMatrix(GeneIM temp_gene_IM[],string TF_TF_address,string TF_Gene_address);
    /// Get the amount of all genes in GRN
    ///
    ///    \return number of genes
    int getGeneAmount();
    /// Get the amount of TFs in GRN
    ///
    ///    \return number of transcription factors
    int getTFAmount();
    /// Original GRN matrix
    double **originalGRN;
    /// Construct genes' information map
    ///
    ///    \param file address of gene info file
    ///    \return map of gene info whose flag is gene name
    map<string,string> mapTFIM(string Gene_IM_address);
    /// Get transcription unit position
    ///
    ///    This position is used to ensure the promoter to each gene.
    ///    \param file address of TU info file
    map<string,string> mapPromoter(string promoters_address);
    /// a vector contains the position of each TU
    void readTUPosition(string TU_position_address);
    /// a vector contains the promoter name of each promoter
    vector<int> TU_position;
    /// a vector contains the promoter name of each promoter
    vector<string> promoter_name_dict;
    /// Get promoter name and sequence
    ///
    ///    Use gene position to confirm the TU which contains it.
    ///    Search promoter name in promoter info map and get its sequence.
    ///    \param array of GeneIM's object
    ///    \see GeneIM
    void getGenePromoter(GeneIM temp_gene_IM[]);
    /// Ensure uncertain genes interaction
    ///
    ///    File named "uncertain" having been output in getRegulationMatrix function is read to change the original matrix.
    ///    Users make sure the interaction and change those uncertain genes in that file.
    void inputUncertainGene();

private:
    /// Build TF-TF GRN and get TF name
    ///
    ///    \param array of GeneIM's object
    ///    \param original GRN matrix
    ///    \param file address of TF-TF regulation file 
    void readTFTF(GeneIM temp_gene_IM[],double **old_GRN,string TF_TF_address);
    /// Build TF-Gene GRN
    ///
    ///    \param array of GeneIM's object
    ///    \param original GRN matrix
    ///    \param file address of TF-Gene regulation file
    void readTFGene(GeneIM temp_gene_IM[],double **old_GRN,string TF_Gene_address);
    /// Add TF not included in TF-TF regualtion
    ///
    ///    Some trandcription factors are not included in TF-TF regualtion but included in TF-Gene regulation.
    ///    \param array of GeneIM's object
    ///    \param file address of TF-Gene regualtion file
    void addTF(GeneIM temp_gene_IM[],string TF_Gene_address);
    /// Get row number of uncertain genes
    vector<int>uncertain_row;
    /// Get column number of uncertain genes
    vector<int>uncertain_column;
    /// Gene number of GRN
    int gene_amount;
    /// Transcription Factor number of GRN
    int TF_amount;
    /// Uncertain gene number
    int unknow;
    /// Output stream of uncertain genes
    ofstream uncertain;
};


#endif
