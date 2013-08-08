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
#include <fstream>
#include <ctime>
#include <cmath>
#include <stdlib.h>
#define gap -8
void GRN::initializeGRN(double oldGRN[][scale], int mSize){
    srand((unsigned)time(0));
    for (int i = 0; i != scale; ++i) {
        for (int j = 0; j != scale; ++j) {
            //newGRNCorrelation[i][j] = oldGRN[i][j];
            if (oldGRN[i][j] == 2) {
                if (((rand() % 100) / 100.0) < 0.5) {
                    newGRNCorrelation[i][j] = 0;
                }else{
                    newGRNCorrelation[i][j] = 0;
                }
            }else{
                newGRNCorrelation[i][j] = oldGRN[i][j];
            }
        }
    }
    matrixSize = mSize;
}
void GRN::constructNewGRN(Sequence seqArry[]){
    double simiMatrix[170] = { 0 };
    loadMatrixBLOSUM_XX();
//*************************************************************************************；
// method_1:
//          Default new gene has regulatory relation with all genes.
//*************************************************************************************；
    //generate similarity matrix;
    for (int geneNum = 0; geneNum != matrixSize; ++geneNum) {
        simiMatrix[geneNum] = aminoASAlignment(seqArry[geneNum].aminoAcidSequence, seqArry[geneNum].aminoASSize, seqArry[matrixSize].aminoAcidSequence, seqArry[matrixSize].aminoASSize);
    }
    //artificial filtering;
    for (int geneNum = 0; geneNum != matrixSize; ++geneNum) {
        /*if (simiMatrix[geneNum] < 0.4) {
            simiMatrix[geneNum] = 0;
        }else if (simiMatrix[geneNum] > 0.8){
            simiMatrix[geneNum] = (1 - 0.8) * simiMatrix[geneNum] + 0.8;
        }*/
        simiMatrix[geneNum] = pow(simiMatrix[geneNum] + 0.4, 20.0)/(pow(simiMatrix[geneNum] + 0.4, 20.0) + 1);
    }
    //insert new correlations to (matrixSize + 1) row;
    for (int j_geneNum = 0; j_geneNum != matrixSize; ++j_geneNum) {
        int counter = 0;
        newGRNCorrelation[matrixSize][j_geneNum] = 0;
        for (int i_geneNum = 0; i_geneNum != matrixSize; ++i_geneNum) {
            if (newGRNCorrelation[i_geneNum][j_geneNum] != 0) {
                newGRNCorrelation[matrixSize][j_geneNum] += newGRNCorrelation[i_geneNum][j_geneNum] * simiMatrix[i_geneNum];
                counter += 1;
            }
        }
        if (counter != 0) {
            newGRNCorrelation[matrixSize][j_geneNum] = newGRNCorrelation[matrixSize][j_geneNum] / 166;
        }
        else
            newGRNCorrelation[matrixSize][j_geneNum] = 0;
    }
    //insert new correlations to (matrixSize + 1) column;
    for (int i_geneNum = 0; i_geneNum != matrixSize; ++i_geneNum) {
        int counter = 0;
        newGRNCorrelation[i_geneNum][matrixSize] = 0;
        for (int j_geneNum = 0; j_geneNum != matrixSize; ++j_geneNum) {
            if (newGRNCorrelation[i_geneNum][j_geneNum] != 0) {
                newGRNCorrelation[i_geneNum][matrixSize] += newGRNCorrelation[i_geneNum][j_geneNum] * simiMatrix[j_geneNum];
                counter += 1;
            }
        }
        if (counter != 0) {
            newGRNCorrelation[i_geneNum][matrixSize] = newGRNCorrelation[i_geneNum][matrixSize] / 166;
        }
        else
            newGRNCorrelation[i_geneNum][matrixSize] = 0;
    }
    newGRNCorrelation[matrixSize][matrixSize] = 0;
    
    std::ofstream outfile;
    outfile.open("/Users/jinyang/Documents/iGEM_Programmes/GRN/similarity");
    for (int i = 0; i != 170; ++i) {
        outfile << i + 1 << '\t' << simiMatrix[i] << std::endl;
    }
    outfile.close();

//*************************************************************************************；
//method_2:
//         Default new gene has RR with top 3 similar gene.
//*************************************************************************************；
/*
    double simiMatrix[170] = { 0 };
    int top_3_index[3] = { 0 };
    //insert new correlations '0' to (matrixSize + 1) row;
    for (int j_geneNum = 0; j_geneNum != matrixSize; ++j_geneNum){
        newGRNCorrelation[matrixSize][j_geneNum] = 0;
    }
    //insert new correlations '0' to (matrixSize + 1) column;
    for (int i_geneNum = 0; i_geneNum != matrixSize; ++i_geneNum){
        newGRNCorrelation[i_geneNum][matrixSize] = 0;
    }
    newGRNCorrelation[matrixSize][matrixSize] = 0;
    //generate similarity matrix;
    for (int geneNum = 0; geneNum != matrixSize; ++geneNum) {
        simiMatrix[geneNum] = aminoASAlignment(seqArry[geneNum].aminoAcidSequence, seqArry[geneNum].aminoASSize, seqArry[matrixSize].aminoAcidSequence, seqArry[matrixSize].aminoASSize);
    }
    //find out top 3 similar gene;
    for (int geneNum = 0; geneNum != matrixSize; ++geneNum) {
        if (simiMatrix[top_3_index[0]] < simiMatrix[geneNum]) {
            top_3_index[0] = geneNum;//find out 1st similar gene;
        }
    }
    for (int geneNum = 0; geneNum != matrixSize; ++geneNum) {
        if (simiMatrix[top_3_index[1]] < simiMatrix[geneNum]) {
            if (geneNum != top_3_index[0]) {
                top_3_index[1] = geneNum;//find out 2nd similar gene;
            }
        }
    }
    for (int geneNum = 0; geneNum != matrixSize; ++geneNum) {
        if (simiMatrix[top_3_index[2]] < simiMatrix[geneNum]) {
            if (geneNum != top_3_index[0] && geneNum != top_3_index[1]) {
                top_3_index[2] = geneNum;//find out 2nd similar gene;
            }
        }
    }
    //insert new correlations to (matrixSize + 1) row;
    for (int i = 0; i != 3; ++i) {
        int counter = 0;
        for (int i_geneNum = 0; i_geneNum != matrixSize; ++i_geneNum) {
            if (newGRNCorrelation[i_geneNum][top_3_index[i]] != 0) {
                newGRNCorrelation[matrixSize][top_3_index[i]] += newGRNCorrelation[i_geneNum][top_3_index[i]] * simiMatrix[i_geneNum];
                counter += 1;
            }
        }
        if (counter != 0) {
            newGRNCorrelation[matrixSize][top_3_index[i]]= newGRNCorrelation[matrixSize][top_3_index[i]] / counter;
        }
        else
            newGRNCorrelation[matrixSize][top_3_index[i]] = 0;
    }
    //insert new correlations to (matrixSize + 1) column;
    for (int i = 0; i != 3; ++i) {
        int counter = 0;
        for (int j_geneNum = 0; j_geneNum != matrixSize; ++j_geneNum) {
            if (newGRNCorrelation[top_3_index[i]][j_geneNum] != 0) {
                newGRNCorrelation[top_3_index[i]][matrixSize] += newGRNCorrelation[top_3_index[i]][j_geneNum] * simiMatrix[j_geneNum];
                counter += 1;
            }
        }
        if (counter != 0) {
            newGRNCorrelation[top_3_index[i]][matrixSize]= newGRNCorrelation[top_3_index[i]][matrixSize] / counter;
        }
        else
            newGRNCorrelation[top_3_index[i]][matrixSize] = 0;
    }
    //test module;
    //std::cout << top_3_index[0] << std::endl;
    //std::cout << top_3_index[1] << std::endl;
    //std::cout << top_3_index[2] << std::endl;
    std::ofstream outfile;
    outfile.open("/Users/jinyang/Documents/iGEM_Programmes/GRN/similarity");
    for (int i = 0; i != 170; ++i) {
        outfile << i + 1 << '\t' << simiMatrix[i] << std::endl;
    }
    outfile.close();
*/
}

double GRN::aminoASAlignment(std::string s, int s_size, std::string t, int t_size){
    double G_ = gap;
    double _G = gap;
    //double normalization_s = 0;
    //double normalization_t = 0;
    //double similarity = 0;
    int s_length = 0;
    int t_length = 0;
    double positives = 0;
    s_length = s_size;
    t_length = t_size;
    std::string commonSubsequence;
    std::string s_sub;
    std::string t_sub;
    std::vector< std::vector<double> > alignMatrix (t_size + 1, std::vector<double>(s_size + 1));
    //initialize alignMatrix;
    alignMatrix[0][0] = 0;
    for (int i = 1; i != t_size + 1; ++i) {
        alignMatrix[i][0] = alignMatrix[i - 1][0] + G_;
    }
    for (int j = 1; j != s_size + 1; ++j) {
        alignMatrix[0][j] = alignMatrix[0][j - 1] + _G;
    }
    for (int i = 1; i != t_size + 1; ++i) {
        for (int j = 1; j != s_size + 1; ++j) {
            //alignMatrix[i][j] = maxValue(alignMatrix[i-1][j-1] + alignScore(s[j],t[i]), alignMatrix[i-1][j] + alignScore(s[j], ' '), alignMatrix[i][j-1] + alignScore(' ', t[i]));
            double st = 0;
            double s_ = 0;
            double _t = 0;
            st = alignMatrix[i - 1][j - 1] + alignScore(t[i], s[j]);
            s_ = alignMatrix[i - 1][j] + alignScore(t[i], ' ');
            _t = alignMatrix[i][j - 1] + alignScore(' ', s[j]);
            alignMatrix[i][j] = maxValue(st, s_, _t);
        }
    }

    for (int i = t_size, j = s_size; (i > 0 || j > 0); ) {
        //for (int j = s_size; j != 0;) {
            if ((i * j > 0) && (alignMatrix[i][j] == alignMatrix[i - 1][j - 1] + alignScore(t[i], s[j]))) {
                if (s[j] == t[i]) {
                    commonSubsequence += s[j];
                    s_sub += s[j];
                    t_sub += t[i];
                }else{
                    int simi_if = 0;
                    simi_if = matrixBLOSUM_XX[alignIndex_BLOSUM50(s[j])][alignIndex_BLOSUM50(t[i])];
                    if (simi_if >= 0) {
                        commonSubsequence += '+';
                    }
                    /*else{
                        commonSubsequence += '#';
                    }*/
                    s_sub += s[j];
                    t_sub += t[i];
                }
            
            i--;
            j--;
            }else if((i > 0) && (alignMatrix[i][j] == alignMatrix[i - 1][j] +  alignScore(t[i], ' '))){
                //insert a space into string_s;
                s_length += 1;
                s_sub += '-';
                t_sub += t[i];
                i--;
                //commonSubsequence += '#';
            }else if ((j > 0) && (alignMatrix[i][j] == alignMatrix [i][j - 1] + alignScore(' ', s[j] ))){
                //insert a space into string_s;
                t_length += 1;
                t_sub += '-';
                s_sub += s[j];
                j--;
                //commonSubsequence += '#';
            }
        //}
    }
    //Give the positives;
    //std::cout << "test_s: " << s_length << std::endl;
    //std::cout << "test_t: " << t_length << std::endl;
    //std::cout << s_sub << std::endl;
    //std::cout << commonSubsequence << std::endl;
    //std::cout << t_sub << std::endl;
    positives = (double)commonSubsequence.size() / s_length;
    return positives;
}

int GRN::alignScore (char t, char s)
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
    int G_ = gap;
    int index_s = 0;
    int index_t = 0;

    if (s == ' ') {
        score = G_;
    }
    else if ( t == ' '){
        score = G_;
    }
    else if (s == t || s != t){
        index_s = alignIndex_BLOSUM50(s);
        index_t = alignIndex_BLOSUM50(t);
        score = matrixBLOSUM_XX[index_s][index_t];
    }
    return score;
}

int GRN::alignIndex_BLOSUM50 (char s){
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
    }
}

/*
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
*/
void GRN::loadMatrixBLOSUM_XX(){
    std::ifstream infile;
    //make sure to fill in the correct path;
    infile.open("/Users/jinyang/Documents/iGEM_Programmes/GRN/BLOSUM_50");
    for (int i = 0; i != 20; ++i) {
        for (int j = 0; j != 20; ++j) {
            infile >> matrixBLOSUM_XX[i][j];
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