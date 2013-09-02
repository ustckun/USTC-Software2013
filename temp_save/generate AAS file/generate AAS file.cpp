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
    int max_size_dna = 0;
    int line = 0;
    string gene;
    ifstream infile;
    infile.open("/Users/jinyang/Documents/iGEM_Programmes/generate AAS file/Sequence");
    for (int i = 0; i != 166; ++i) {
        infile >> gene_number;
        infile >> gene;
        size = (int)gene.size();
        //test longest sequence;
        if (max_size_dna < size) {
            max_size_dna = size;
            line = i;
        }
        seqArray[i].initializeGeneSequence(gene, gene_number, size);
        seqArray[i].translation();
    }
    infile.close();
    ofstream outfile;
    outfile.open("/Users/jinyang/Documents/iGEM_Programmes/generate AAS file/AAS");
    for (int i = 0; i != 166; ++i) {
        outfile << seqArray[i].geneNumber;
        outfile << seqArray[i].aminoAcidSequence << endl;
    }
    outfile.close();
    cout << max_size_dna << endl;
    cout << line << endl;
}

