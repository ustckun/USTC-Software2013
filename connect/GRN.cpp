//
//  GRN.cpp
//  GRN
//
//  Created by jinyang on 13-7-26.
//  Copyright (c) 2013Äê Li Jinyang. All rights reserved.
//

#include"TFIM.h"
//#include"calculation.h"
#include"Sequence.h"
#include"GRN.h"
#include"Regulation.h"
#include"EasytoDebug.h"

void GRN::initializeGRN(double oldGRN[][scale], int mSize){
    for (int i = 0; i != scale; ++i) {
        for (int j = 0; j != scale; ++j) {
            newGRNCorrelation[i][j] = oldGRN[i][j];
        }
    }
    matrixSize = mSize;
}
void GRN::constructNewGRN(Sequence seqArry[]){
    //insert new correlations to (matrixSize + 1) row;
    for (int j_geneNum = 0; j_geneNum != matrixSize; ++j_geneNum) {
        newGRNCorrelation[matrixSize + 1][j_geneNum] = aminoASAlignment(seqArry[matrixSize + 1].getAminoAcidSequence(), seqArry[matrixSize + 1].aminoASSize, seqArry[j_geneNum].getAminoAcidSequence(), seqArry[j_geneNum].aminoASSize);
    }
    //insert new correlations to (matrixSize + 1) column;
    for (int i_geneNum = 0; i_geneNum; ++i_geneNum) {
        newGRNCorrelation[i_geneNum][matrixSize + 1] = aminoASAlignment(seqArry[matrixSize + 1].getAminoAcidSequence(), seqArry[matrixSize + 1].aminoASSize, seqArry[i_geneNum].getAminoAcidSequence(), seqArry[i_geneNum].aminoASSize);
    }
    newGRNCorrelation[matrixSize + 1][matrixSize + 1] = 0;
}

double GRN::aminoASAlignment(std::string s, int s_size, std::string t, int t_size){
    int G_ = 0;
    int _G = 0;
    std::vector< std::vector<double> > alignMatrix (t_size, std::vector<double>(s_size));
    //initialize alignMatrix;
    alignMatrix[0][0] = 0;
    for (int i = 1; i != t_size; ++i) {
        alignMatrix[i][0] = alignMatrix[i - 1][0] + G_;
    }
    for (int j = 1; j != s_size; ++j) {
        alignMatrix[0][j] = alignMatrix[0][j - 1] + _G;
    }
    for (int i = 1; i != t_size; ++i) {
        for (int j = 1; j != s_size; ++j) {
            alignMatrix[i][j] = maxValue(alignMatrix[i-1][j-1] + alignScore(s[j],t[i]), alignMatrix[i-1][j] + alignScore(s[j], ' '), alignMatrix[i][j-1] + alignScore(' ', t[i]));
        }
    }
    return alignMatrix[t_size - 1][s_size - 1];
}

int GRN::alignScore (char s, char t)
{
    int score = 0;
    int GG = 1;
    int _G = 0;
    int GT = 0;
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

int GRN::maxValue (int a, int b, int c)
{
    if(a<b)
        a=b;
    if(a<c)
        a=c;
    return a;
}