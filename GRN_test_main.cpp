//
//  main.cpp
//  GRN
//
//  Created by jinyang on 13-7-26.
//  Copyright (c) 2013年 Li Jinyang. All rights reserved.
//

#include <iostream>
#include <fstream>
#include "Sequence.h"
#include "GRN.h"

using namespace std;

int main(){
    //the index is gengeNumber;
    /*string dna[3];
    dna[0] = "ATGAAACGCCCGGACTACAGAACATTACAGGCACTGGATGCGGTGATACGTGAACGAGGATTTGAGCGCGCGGCACAAAAGCTGTGCATTACACAATCAGCCGTCTCACAGCGCATTAAGCAACTGGAAAATATGTTCGGGCAGCCGCTGTTGGTGCGTACCGTACCGCCGCGCCCGACGGAACAAGGGCAAAAACTGCTGGCACTGCTGCGCCAGGTGGAGTTGCTGGAAGAAGAGTGGCTGGGCGATGAACAAACCGGTTCGACTCCGCTGCTGCTTTCACTGGCGGTCAACGCCGACAGTCTGGCGACGTGGTTGCTTCCTGCACTGGCTCCTGTGTTGGCTGATTCGCCTATCCGCCTCAACTTGCAGGTAGAAGATGAAACCCGCACTCAGGAACGTCTGCGCCGCGGCGAAGTGGTCGGCGCGGTGAGTATTCAACATCAGGCGCTGCCGAGTTGTCTTGTCGATAAACTTGGTGCGCTCGACTATCTGTTCGTCAGCTCAAAACCCTTTGCCGAAAAATATTTCCCTAACGGCGTAACGCGTTCGGCATTACTGAAAGCGCCAGTGGTCGCGTTTGACCATCTTGACGATATGCACCAGGCCTTTTTGCAGCAAAACTTCGATCTGCCTCCAGGCAGCGTGCCCTGCCATATCGTTAATTCTTCAGAAGCGTTCGTACAACTTGCTCGCCAGGGCACCACCTGCTGTATGATCCCGCACCTGCAAATCGAGAAAGAGCTGGCCAGCGGTGAACTGATTGACTTAACGCCTGGGCTATTTCAACGACGGATGCTCTACTGGCACCGCTTTGCTCCTGAAAGCCGCATGATGCGTAAAGTCACTGATGCGTTACTCGATTATGGTCACAAAGTCCTTCGTCAGGATTAA";
    dna[1] = "ATGACGCAGGATGAATTGAAAAAAGCAGTAGGATGGGCGGCACTTCAGTATGTTCAGCCCGGCACCATTGTTGGTGTAGGTACAGGTTCCACCGCCGCACACTTTATTGACGCGCTCGGTACAATGAAAGGCCAGATTGAAGGGGCCGTTTCCAGTTCAGATGCTTCCACTGAAAAACTGAAAAGCCTCGGCATTCACGTTTTTGATCTCAACGAAGTCGACAGCCTTGGCATCTACGTTGATGGCGCAGATGAAATCAACGGCCACATGCAAATGATCAAAGGCGGCGGCGCGGCGCTGACCCGTGAAAAAATCATTGCTTCGGTTGCAGAAAAATTTATCTGTATTGCAGACGCTTCCAAGCAGGTTGATATTCTGGGTAAATTCCCGCTGCCAGTAGAAGTTATCCCGATGGCACGTAGTGCAGTGGCGCGTCAGCTGGTGAAACTGGGCGGTCGTCCGGAATACCGTCAGGGCGTGGTGACCGATAATGGCAACGTGATCCTCGACGTCCACGGCATGGAAATCCTTGACCCGATAGCGATGGAAAACGCCATAAATGCGATTCCTGGCGTGGTGACTGTTGGCTTGTTTGCTAACCGTGGCGCGGACGTTGCGCTGATTGGCACACCTGACGGTGTCAAAACCATTGTGAAATGA";
    dna[2] = "ATGTCTTATCAGTATGTTAACGTTGTCACTATCAACAAAGTGGCGGTCATTGAGTTTAACTATGGCCGAAAACTTAATGCCTTAAGTAAAGTCTTTATTGATGATCTTATGCAGGCGTTAAGCGATCTCAACCGGCCGGAAATTCGCTGTATCATTTTGCGCGCACCGAGTGGATCCAAAGTCTTCTCCGCAGGTCACGATATTCACGAACTGCCGTCTGGCGGTCGCGATCCGCTCTCCTATGATGATCCATTGCGTCAAATCACCCGCATGATCCAAAAATTCCCGAAACCGATCATTTCGATGGTGGAAGGTAGTGTTTGGGGTGGCGCATTTGAAATGATCATGAGTTCCGATCTGATCATCGCCGCCAGTACCTCAACCTTCTCAATGACGCCTGTAAACCTCGGCGTCCCGTATAACCTGGTCGGCATTCACAACCTGACCCGCGACGCGGGCTTCCACATTGTCAAAGAGCTGATTTTTACCGCTTCGCCAATCACCGCCCAGCGCGCGCTGGCTGTCGGCATCCTCAACCATGTTGTGGAAGTGGAAGAACTGGAAGATTTCACCTTACAAATGGCGCACCACATCTCTGAGAAAGCGCCGTTAGCCATTGCCGTTATCAAAGAAGAGCTGCGTGTACTGGGCGAAGCACACACCATGAACTCCGATGAATTTGAACGTATTCAGGGGATGCGCCGCGCGGTGTATGACAGCGAAGATTACCAGGAAGGGATGAACGCTTTCCTCGAAAAACGTAAACCTAATTTCGTTGGTCATTAA";
    int dna_size[3];
    dna_size[0] = 894;
    dna_size[1] = 660;
    dna_size[2] = 786;
    
    //cout << dna[0].size() << endl;
    //cout << dna[1].size() << endl;
    //cout << dna[2].size() << endl;
    
    Sequence DNA_Array[3];
    
    //initialize DNA_Array;
    for (int geneNumber = 0; geneNumber != 3; ++geneNumber) {
        DNA_Array[geneNumber].initializeGeneSequence(dna[geneNumber], geneNumber, dna_size[geneNumber]);
        DNA_Array[geneNumber].translation();
    }
    
    int m_size = 2;
    double old_grn[170][170];
    for (int i = 0; i != 170; ++i) {
        for (int j = 0; j != 170; ++j) {
            old_grn[i][j] = 2;
        }
    }
    old_grn[0][0] = 0;
    old_grn[0][1] = -1;
    old_grn[1][0] = 1;
    old_grn[1][1] = 2;
    
    GRN test_GRN;
    test_GRN.initializeGRN(old_grn, m_size);
    test_GRN.constructNewGRN(DNA_Array);
    for (int i = 0; i != 5; ++i) {
        for (int j = 0; j != 5; ++j) {
            cout << test_GRN.newGRNCorrelation[i][j] << '\t';
        }
        cout << endl;
    }*/
    
//======================================================================================；
    string dna[166];
    int geneNum[166];
    int dna_size[166];
    double grn[170][170];
    GRN test_GRN;
    Sequence seq_array[166];
    ifstream infile;
    infile.open("/Users/jinyang/Documents/iGEM_Programmes/GRN/Sequence");
    for (int i = 0; i != 166; ++i) {
        infile >> geneNum[i];
        infile >> dna[i];
        dna_size[i] = (int)dna[i].size();
        seq_array[i].initializeGeneSequence(dna[i], geneNum[i], dna_size[i]);
        seq_array[i].translation();
    }
    infile.close();
    int m_size = 165;
    //intialize GRN;
    infile.open("/Users/jinyang/Documents/iGEM_Programmes/GRN/regulation");
    for (int i = 0; i != 170; ++i) {
        for (int j = 0; j != 170; ++j) {
            grn[i][j] = 2;
        }
    }
    for (int i = 0; i != 170; ++i) {
        for (int j = 0; j != 170 ; ++j) {
            infile >> grn[i][j];
        }
    }
    infile.close();
    test_GRN.initializeGRN(grn, m_size);
    test_GRN.constructNewGRN(seq_array);
    ofstream outfile;
    outfile.open("/Users/jinyang/Documents/iGEM_Programmes/GRN/test");
    for (int i = 0; i != 170; ++i) {
        for (int j = 0; j != 170; ++j) {
            outfile << test_GRN.newGRNCorrelation[i][j] << '\t';
        }
        outfile << endl;
    }
 
//=======================================================================================;
/*    Sequence debug;
    string s;
    int geneNum = 0;
    int size;
    
    ifstream infile;
    infile.open("/Users/jinyang/Documents/iGEM_Programmes/GRN/Sequence");
    infile >> geneNum;
    infile >> s;
    size = (int)s.size();
    debug.initializeGeneSequence(s, geneNum, size);
    debug.translation();
    cout << "Sequence_DNA size is: " << debug.geneSequence.size() << endl;
    cout << "Sequence_dna_size is:" << debug.DNASize << endl;
    cout << endl;
    cout << "Sequence_AAS size is:" << debug.aminoAcidSequence.size() << endl;
    cout << "Sequence_AAS_size is:" << debug.aminoASSize <<endl;
*/
}

