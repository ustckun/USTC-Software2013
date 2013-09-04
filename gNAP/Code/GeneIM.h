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

class GeneIM
{
public:
    GeneIM()
    {
        get_error=0;
        left_position=0;
        right_position=0;
    }

    int gene_number;
    char *name;
    void getGeneInformation(map<string,string> dict);
    void getPromoterIF(map<string,string> dict);
    string getID();
    string getGeneSequence();
    string getPromoterSequence();
    string getPromoterName();
    string getGeneTrueName();
    int getLeftPosition();
    int getRightPosition();
    int getFileError();
    char *getGeneName();
    void putName();
    void putInPromoterName(string promoter);
    int getRNA();
    string getGeneDescription();
    char gene_name[10];
private:
    string iD;
    string gene_sequence;
    int left_position;
    int right_position;
    string gene_description;
    string promoter_name;
    string promoter_sequence;
    string true_name;
    int get_error;
    int RNA;
};

#endif
