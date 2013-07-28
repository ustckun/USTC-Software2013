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
        matrixSize = 0;
    }
    void initializeGRN(double oldGRN[][scale], int mSize);
    double newGRNCorrelation[scale][scale];
private:
    //Results of alignments are correlations;
    double aminoASAlignment(std::string s, int s_size, std::string t, int t_size);
    void constructNewGRN(Sequence seqArray[]);
    int alignScore (char s, char t);
    int maxValue (int a, int b, int c);
    int matrixSize;
};

#endif /* defined(__GRN__GRN__) */
