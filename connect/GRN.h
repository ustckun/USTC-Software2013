//
//  GRN.h
//  GRN
//
//  Created by jinyang on 13-7-26.
//  Copyright (c) 2013Äê Li Jinyang. All rights reserved.
//

#ifndef __GRN__GRN__
#define __GRN__GRN__
#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>
#include "Sequence.h"

#define scale 100

class GRN{
public:
    GRN(){
        for (int i = 0; i != scale; ++i) {
            for (int j = 0; j != scale; ++j) {
                newGRNCorrelation[i][j] = 2;
            }
        }
        for (int i = 0; i != 20; ++i) {
            for (int j = 0; j != 20; ++j) {
                matrixBLOSUM_XX[i][j] = 0;
            }
        }
        matrixSize = 0;
    }
    void initializeGRN(double oldGRN[][scale], int mSize);
    void constructNewGRN(Sequence seqArray[]);
    double newGRNCorrelation[scale][scale];
    
    double aminoASAlignment(std::string s, int s_size, std::string t, int t_size);
    void loadMatrixBLOSUM_XX();
    
private:
    int matrixSize;
    int matrixBLOSUM_XX[20][20];
    //Results of alignments are correlations;
    
    int alignScore (char t, char s);
    double maxValue (double a, double b, double c);
    int alignIndex_BLOSUM50(char s);
    /*
    int alignIndex_BLOSUM62(char s);
    */
    
};

#endif /* defined(__GRN__GRN__) */
