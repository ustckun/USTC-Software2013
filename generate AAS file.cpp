//
//  main.cpp
//  generate AAS file
//
//  Created by jinyang on 13-7-29.
//  Copyright (c) 2013年 Li Jinyang. All rights reserved.
//

#include <iostream>
#include <fstream>
#include "Sequence.h"
using namespace std;

int main(){
    Sequence seqArray[166];
    int gene_number = 0;
    int size;
    string gene;
    ifstream infile;
    infile.open("/Users/jinyang/Documents/iGEM_Programmes/generate AAS file/Sequence");
    for (int i = 0; i != 166; ++i) {
        infile >> gene_number;
        infile >> gene;
        size = (int)gene.size();
        seqArray[i].initializeGeneSequence(gene, gene_number, size);
    }
    infile.close();
    ofstream outfile;
    outfile.open("/Users/jinyang/Documents/iGEM_Programmes/generate AAS file/AAS");
    for (int i = 0; i != 166; ++i) {
        outfile << seqArray[i].geneNumber;
        outfile << seqArray[i].getAminoAcidSequence() << endl;
    }
    outfile.close();
}

