//
//  GRN.cpp
//  GRN
//
//  Created by jinyang on 13-7-26.
//  Copyright (c) 2013年 Li Jinyang. All rights reserved.
//

#include "GRN.h"
#include <vector>
#include <string>
void GRN::initializeGRN(double oldGRN[][scale], int mSize){
    for (int i = 0; i != scale; ++i) {
        for (int j = 0; j != scale; ++j) {
            newGRNCorrelation[i][j] = oldGRN[i][j];
        }
    }
    matrixSize = mSize;
}
void GRN::constructNewGRN(Sequence seqArry[]){
    loadMatrixBLOSUM_62();
    //insert new correlations to (matrixSize + 1) row;
    for (int j_geneNum = 0; j_geneNum != matrixSize; ++j_geneNum) {
        int counter = 0;
        newGRNCorrelation[matrixSize][j_geneNum] = 0;
        for (int i_geneNum = 0; i_geneNum != matrixSize; ++i_geneNum) {
            newGRNCorrelation[matrixSize][j_geneNum] += newGRNCorrelation[i_geneNum][j_geneNum] * aminoASAlignment(seqArry[matrixSize ].aminoAcidSequence, seqArry[matrixSize].aminoASSize, seqArry[j_geneNum].aminoAcidSequence, seqArry[j_geneNum].aminoASSize);
            if (newGRNCorrelation[i_geneNum][j_geneNum] != 2) {
                counter += 1;
            }
        }
        newGRNCorrelation[matrixSize][j_geneNum] = newGRNCorrelation[matrixSize][j_geneNum] / counter;
    }
    //insert new correlations to (matrixSize + 1) column;
    for (int i_geneNum = 0; i_geneNum != matrixSize; ++i_geneNum) {
        int counter = 0;
        newGRNCorrelation[i_geneNum][matrixSize] = 0;
        for (int j_geneNum = 0; j_geneNum != matrixSize; ++j_geneNum) {
            newGRNCorrelation[i_geneNum][matrixSize] += newGRNCorrelation[matrixSize][j_geneNum] * aminoASAlignment(seqArry[matrixSize].aminoAcidSequence, seqArry[matrixSize].aminoASSize, seqArry[i_geneNum].aminoAcidSequence, seqArry[i_geneNum].aminoASSize);
            if (newGRNCorrelation[i_geneNum][j_geneNum] != 2) {
                counter += 1;
            }
        }
        newGRNCorrelation[i_geneNum][matrixSize] = newGRNCorrelation[i_geneNum][matrixSize] / counter;
    }
    newGRNCorrelation[matrixSize][matrixSize] = 0;
}

double GRN::aminoASAlignment(std::string s, int s_size, std::string t, int t_size){
    double G_ = -5;
    double _G = -5;
    double normalization_s = 0;
    double normalization_t = 0;
    double similarity = 0;
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
            //alignMatrix[i][j] = maxValue(alignMatrix[i-1][j-1] + alignScore(s[j],t[i]), alignMatrix[i-1][j] + alignScore(s[j], ' '), alignMatrix[i][j-1] + alignScore(' ', t[i]));
            double st = 0;
            double s_ = 0;
            double _t = 0;
            st = alignMatrix[i - 1][j - 1] + alignScore(s[j], t[i]);
            s_ = alignMatrix[i - 1][j] + alignScore(s[j], ' ');
            _t = alignMatrix[i][j - 1] + alignScore(' ', t[i]);
            alignMatrix[i][j] = maxValue(st, s_, _t);
        }
    }
    //normalize alignment score;
    for (int i = 1; i != s_size; ++i) {
        normalization_s += alignScore(s[i], s[i]);
    }
    for (int i = 0; i != t_size; ++i) {
        normalization_t += alignScore(t[i], t[i]);
    }
    similarity = fabs(2 * alignMatrix[t_size - 1][s_size - 1]) / (normalization_s + normalization_t);
    return similarity;
}

int GRN::alignScore (char s, char t)
{
    /*int score = 0;
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
        score = GT;*/
    int score = 0;
    int G_ = -5;
    int index_s = 0;
    int index_t = 0;

    if (s == ' ') {
        score = G_;
    }
    else if ( t == ' '){
        score = G_;
    }
    else if (s == t || s != t){
        index_s = alignIndex_BLOSUM62(s);
        index_t = alignIndex_BLOSUM62(t);
        score = matrixBLOSUM_62[index_s][index_t];
    }
    return score;
}

int GRN::alignIndex_BLOSUM62 (char s){
    switch (s) {
        case 'C':
            return 0;
            break;
        case 'S':
            return 1;
            break;
        case 'T':
            return 2;
            break;
        case 'P':
            return 3;
            break;
        case 'A':
            return 4;
            break;
        case 'G':
            return 5;
            break;
        case 'N':
            return 6;
            break;
        case 'D':
            return 7;
            break;
        case 'E':
            return 8;
            break;
        case 'Q':
            return 9;
            break;
        case 'H':
            return 10;
            break;
        case 'R':
            return 11;
            break;
        case 'K':
            return 12;
            break;
        case 'M':
            return 13;
            break;
        case 'I':
            return 14;
            break;
        case 'L':
            return 15;
            break;
        case 'V':
            return 16;
            break;
        case 'F':
            return 17;
            break;
        case 'Y':
            return 18;
            break;
        case 'W':
            return 19;
            break;
    }
}

void GRN::loadMatrixBLOSUM_62(){
    std::ifstream infile;
    //make sure to fill in the correct path;
    infile.open("/Users/jinyang/Documents/iGEM_Programmes/GRN/BLOSUM_62");
    for (int i = 0; i != 20; ++i) {
        for (int j = 0; j != 20; ++j) {
            infile >> matrixBLOSUM_62[i][j];
        }
    }
    infile.close();
}

double GRN::maxValue (double a, double b, double c)
{
    if(a<b)
        a=b;
    if(a<c)
        a=c;
    return a;
}