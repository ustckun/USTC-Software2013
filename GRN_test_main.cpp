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
/*
    //the index is gengeNumber;
    string dna[3];
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
        DNA_Array[geneNumber].initializeGeneSequence(dna[geneNumber], geneNumber, dna_size[geneNumber], 0);
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
    }
*/
//=====================================================================================；

    string dna[166];
    int geneNum[166];
    int dna_size[166];
    int prod_type[166] = { 0 };
    double grn[170][170];
    GRN test_GRN;
    Sequence seq_array[167];
    ifstream infile;
    infile.open("/Users/jinyang/Documents/iGEM_Programmes/GRN/Sequence");
    for (int i = 0; i != 166; ++i) {
        infile >> geneNum[i];
        infile >> dna[i];
        dna_size[i] = (int)dna[i].size();
        seq_array[i].initializeGeneSequence(dna[i], geneNum[i], dna_size[i], prod_type[i]);
        seq_array[i].translation();
    }
    infile.close();
    //seq_array[166].test(" MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTFAYGVQCFSRYPDHMKRHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSVLSKDPNEKRDHMVLLEFVTAAGITHGMDELYK", 238);
    seq_array[166].test(" VHWQTHTVFNQPIPLNNSNLYLSDGALCEAVTREGAGWDSDFLASIGQQLGTAESLELGRLANVNPPELLRYDAQGRRLDDVRFHPAWHLLMQALCTNRVHNLAWEEDARSGAFVARAARFMLHAQVEAGSLCPITMTFAATPLLLQMLPAPFQDWTTPLLSDRYDSHLLPGGQKRGLLIGMGMTEKQGGSDVMSNTTRAERLEDGSYRLVGHKWFFSVPQSDAHLVLAQTAGGLSCFFVPRFLPDGQRNAIRLERLKDKLGNRSNASCEVEFQDAIGWLLGLEGEGIRLILKMGGMTRFDCALGSHAMMRRAFSLAIYHAHQRHVFGNPLIQQPLMRHVLSRMALQLEGQTALLFRLARAWDRRADAKEALWARLFTPAAKFVICKRGMPFVAEAMEVLGGIGYCEESELPRLYREMPVNSIWEGSGNIMCLDVLRVLNKQAGVYDLLSEAFVEVKGQDRYFDRAVRRLQQQLRKPAEELGREITHQLFLLGCGAQMLKYASPPMAQAWCQVMLDTRGGVRLSEQIQNDLLLRATGGVCV", 541);
    int m_size = 166;
    //intialize GRN;
    infile.open("/Users/jinyang/Documents/iGEM_Programmes/GRN/regulation_2.0");
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
    outfile.open("/Users/jinyang/Documents/iGEM_Programmes/GRN/test_method_2");
    for (int i = 0; i != 170; ++i) {
        for (int j = 0; j != 170; ++j) {
            outfile << test_GRN.newGRNCorrelation[i][j] << '\t';
        }
        outfile << endl;
    }
    outfile << endl;
    outfile.close();
    outfile.open("/Users/jinyang/Documents/iGEM_Programmes/GRN/test.txt");
    
    for (int i = 0; i != 170; ++i) {
        outfile << i + 1 << '\t' << test_GRN.newGRNCorrelation[166][i] << endl;
    }
    outfile << endl;
    
    for (int i = 0; i != 170; ++i) {
        outfile << i + 1 << '\t' << test_GRN.newGRNCorrelation[i][166] << endl;
    }
    outfile.close();
    
    //use system time as file name;
    time_t nowtime = time(NULL);
    struct tm *p;
    p = gmtime(&nowtime);
    char filename_1[256] ={ 0 };
    string filename = "/Users/jinyang/Desktop/Parameter Data Test/row+colume ";
    sprintf(filename_1, "%d-%d %d%02d", 1 + p -> tm_mon, p -> tm_mday, 8 + p -> tm_hour, p -> tm_min);
    filename += filename_1;
    filename += ".txt";
    outfile.open(filename);
    for (int i = 0; i != 170; ++i) {
        outfile << i + 1 << '\t' << test_GRN.newGRNCorrelation[166][i] << '\t' << test_GRN.newGRNCorrelation[i][166]<< endl;
    }
    outfile << endl;
    
 
//=====================================================================================;

// rewrite regulation file, exchange 0 and 2;
/*    ifstream infile;
    int s[170][170];
    infile.open("/Users/jinyang/Documents/iGEM_Programmes/GRN/regulation");
    for (int i = 0; i != 170; ++i) {
        for (int j = 0; j != 170; ++j) {
            infile >> s[i][j];
        }
    }
    infile.close();
    for (int i = 0; i != 170; ++i) {
        for (int j = 0; j != 170; ++j) {
            switch (s[i][j]) {
                case 2:
                    s[i][j] = 0;
                    break;
                case 0:
                    s[i][j] = 2;
                    break;
            }
        }
    }
    ofstream outfile;
    outfile.open("/Users/jinyang/Documents/iGEM_Programmes/GRN/regulation_2.0");
    for (int i = 0; i != 170; ++i) {
        for (int j = 0; j != 170; ++j) {
            outfile << s[i][j] << '\t';
        }
        outfile << endl;
    }*/
   
/*
    double x = 0;
    Sequence s;
    Sequence t;
    s.aminoAcidSequence = " HEAGAWGHEE";
    t.aminoAcidSequence = " PAWHEAE";
    s.aminoASSize = 10;
    t.aminoASSize = 7;
    GRN test;
    test.loadMatrixBLOSUM_XX();
    x = test.aminoASAlignment(s.aminoAcidSequence, s.aminoASSize, t.aminoAcidSequence, t.aminoASSize);
    cout << x << endl;
*/

}


