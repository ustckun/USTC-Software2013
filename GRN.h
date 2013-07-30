//
//  GRN.h
//  GRN
//
//  Created by jinyang on 13-7-26.
//  Copyright (c) 2013年 Li Jinyang. All rights reserved.
//

#ifndef __GRN__GRN__
#define __GRN__GRN__
#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>
#include "Sequence.h"

#define scale 170

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
                matrixBLOSUM_62[i][j] = 0;
            }
        }
        matrixSize = 0;
    }
    void initializeGRN(double oldGRN[][scale], int mSize);
    void constructNewGRN(Sequence seqArray[]);
    double newGRNCorrelation[scale][scale];
private:
    int matrixSize;
    int matrixBLOSUM_62[20][20];
    //Results of alignments are correlations;
    double aminoASAlignment(std::string s, int s_size, std::string t, int t_size);
    int alignScore (char s, char t);
    double maxValue (double a, double b, double c);
    int alignIndex_BLOSUM62(char s);
    void loadMatrixBLOSUM_62();
};

#endif /* defined(__GRN__GRN__) */
