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
    //seq_array[166].test(" MCEGYVEKPLYLLIAEWMMAENRWVIAREISIHFDIEHSKAVNTLTYILSEVTEISCEVKMIPNKLEGRGCQCQRLVKVVDIDEQIYARLRNNSREKLVGVRKTPRIPAVPLTELNREQKWQMMLSKSMRR", 131);
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
    string filename = "/Users/jinyang/Desktop/Parameter Data Test/";
    sprintf(filename_1, "%d-%d %d%02d", 1 + p -> tm_mon, p -> tm_mday, 8 + p -> tm_hour, p -> tm_min);
    filename += filename_1;
    filename += " row+colume.txt";
    outfile.open(filename);
    for (int i = 0; i != 170; ++i) {
        outfile << i + 1 << '\t' << test_GRN.newGRNCorrelation[166][i] << '\t' << test_GRN.newGRNCorrelation[i][166]<< endl;
    }
    outfile << endl;
    cout << (double)clock()/CLOCKS_PER_SEC << endl;
*/
    double **old_GRN;
    old_GRN = new double*[1800];
    for (int i = 0; i != 1800; ++i) {
        old_GRN[i] = new double[220];
    }
    ifstream infile;
    // Old_GRN is a matrix at the size 1800 * 220.
    infile.open("/Users/jinyang/Documents/iGEM_Programmes/GRN/Old_GRN");
    for (int i = 0; i != 1800; ++i) {
        for (int j = 0; j != 220; ++j) {
            infile >> old_GRN[i][j];
        }
    }
    infile.close();
    GRN test_grn;
    test_grn.initialize_GRN(old_GRN, 1749, 199);
    Sequence sequence_array[1750];
    int ru_num[1749] = { 0 };
    string pro_array[1749];
    int pro_size_array[1749] = { 0 };
    string gene_array[1749];
    int gene_size_array[1749] = { 0 };
    infile.open("/Users/jinyang/Documents/iGEM_Programmes/GRN/NO&PRE&DNA");
    for (int i = 0; i != 1749; ++i) {
        infile >> ru_num[i];
        infile >> pro_array[i];
        infile >> gene_array[i];
        pro_size_array[i] = (int)pro_array[i].size();
        gene_size_array[i] = (int)gene_array[i].size();
        sequence_array[i].initialize_Sequence(ru_num[i],
                                              pro_array[i],
                                              pro_size_array[i],
                                              gene_array[i],
                                              gene_size_array[i]);
    }
    infile.close();
    int new_num = 1749;
    string new_promoter = "CTTCACCGTCACTTCACATAGCTTCAAATTCTTCCCACATAGTCTTCGTATCCTGCTGCCATTGCAAAGGAGAAGACTATG";
    string new_gene = "ATGCAGACCCCGCACATTCTTATCGTTGAAGACGAGTTGGTAACACGCAACACGTTGAAAAGTATTTTCGAAGCGGAAGGCTATGATGTTTTCGAAGCGACAGATGGCGCGGAAATGCATCAGATCCTCTCTGAATATGACATCAACCTGGTGATCATGGATATCAATCTGCCGGGTAAGAACGGTCTTCTGTTAGCGCGTGAACTGCGCGAGCAGGCGAATGTTGCGTTGATGTTCCTGACTGGCCGTGACAACGAAGTCGATAAAATTCTCGGCCTCGAAATCGGTGCAGATGACTACATCACCAAACCGTTCAACCCGCGTGAACTGACGATTCGTGCACGCAACCTACTGTCCCGTACCATGAATCTGGGTACTGTCAGCGAAGAACGTCGTAGCGTTGAAAGCTACAAGTTCAATGGTTGGGAACTGGACATCAACAGCCGTTCGTTGATCGGCCCTGATGGCGAGCAGTACAAGCTGCCGCGCAGCGAGTTCCGCGCCATGCTTCACTTCTGTGAAAACCCAGGCAAAATTCAGTCCCGTGCTGAACTGCTGAAGAAAATGACCGGCCGTGAGCTGAAACCGCACGACCGTACTGTAGACGTGACGATCCGCCGTATTCGTAAACATTTCGAATCTACGCCGGATACGCCGGAAATCATCGCCACCATTCACGGTGAAGGTTATCGCTTCTGCGGTGATCTGGAAGATTAA";
    sequence_array[1749].initialize_Sequence(new_num,
                                             new_promoter,
                                             (int)new_promoter.size(),
                                             new_gene,
                                             (int)new_gene.size());
    test_grn.construct_new_GRN(sequence_array);
    ofstream outfile;
    //use system time as file name;
    time_t nowtime = time(NULL);
    struct tm *p;
    p = gmtime(&nowtime);
    char filename_1[256] ={ 0 };
    string filename = "/Users/jinyang/Desktop/Parameter Data Test/";
    sprintf(filename_1, "%d-%d %d%02d", 1 + p -> tm_mon, p -> tm_mday, 8 + p -> tm_hour, p -> tm_min);
    filename += filename_1;
    filename += " row+colume.txt";
    outfile.open(filename);
    for (int i = 0; i != 220; ++i) {
        outfile << i + 1 << '\t' << test_grn.new_GRN[1749][i] << endl;
    }
    outfile << endl;
    for (int i = 0; i != 1800; ++i) {
        outfile << i + 1 << '\t' << test_grn.new_GRN[i][199] << endl;
    }
    outfile << endl;
    outfile << (double)clock()/CLOCKS_PER_SEC << endl;
    outfile.close();
    cout << (double)clock()/CLOCKS_PER_SEC << endl;
}


