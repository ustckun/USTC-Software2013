#ifndef __USTC_Software__Regulation__
#define __USTC_Software__Regulaiton__


#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <map>

#include <string.h>

#include <cstring>

#include"define.h"

using namespace std;

class GetReady
{
public:
    GetReady();
    void getRegulationMatrix(GeneIM temp_gene_IM[],string TF_TF_address,string TF_Gene_address);
    int getGeneAmount();
    int getTFAmount();
    //void fullFill();
    int getOpenError();
    //double originalMatrix[N][N];//0:no regulation 1:+active -1:-negative 2:NULL 0.5:+ -
    double **originalGRN;
    map<string,string> mapTFIM(string Gene_IM_address);
    map<string,string> mapPromoter(string promoters_address);
    vector<int>uncertain_row;
    vector<int>uncertain_column;
    void readTUPosition(string TU_position_address);
    vector<int> TU_position;
    vector<string> promoter_name_dict;
    void getGenePromoter(GeneIM temp_gene_IM[]);

    void debug();

private:
    void readTFTF(GeneIM temp_gene_IM[],double **old_GRN,string TF_TF_address);//read Name in TF-TF.txt to GeneIM, intput is GeneIM projects array
    void readTFGene(GeneIM temp_gene_IM[],double **old_GRN,string TF_Gene_address);
    void addTF(GeneIM temp_gene_IM[],string TF_Gene_address);
    int gene_amount;
    int TF_amount;
    int open_file_error;
};


#endif
