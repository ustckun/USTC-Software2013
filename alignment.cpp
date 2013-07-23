//
//  main.cpp
//  SequenceAlignment
//
//  Created by jinyang on 13-7-19.
//  Copyright (c) 2013年 jinyang. All rights reserved.
//
//********************************************************************************
// Some explain for GRN:
//  1.For sequence-alignment module, GRN is obtained from data-capture module;
//  2.The GRN matrix: line for REGULATOR, row for TARGET;
//  3.The elements are supposed to be 0, -1 and +1;
//  4.i.e.
//      If gene_A enhances gene_B, then the matrix records it as:
//      GRN[index_B][index_A] = 1;
//********************************************************************************
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

using namespace std;

int getSequenceSize (string s)
{
    int i;
    for (i = 0; s[i] != '\0'; ++i);
    return i;
}

int MaxValue (int a, int b, int c)
{
    if(a<b)
        a=b;
    if(a<c)
        a=c;
    return a;
}

//alignment score;
int AlignScore (char s, char t)
{
    int score = 0;
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
        score = GT;
    return score;
}

//Needleman_Wunch sequence alignment;
double alignSequence (string s, string t)
{
    int G_ = 0;
    int _G = 0;
    int m = getSequenceSize (t);
    int n = getSequenceSize (s);
    //cout << m << endl << n << endl;
    
    ofstream ofile;
    ofile.open("result_Haosen.txt");
    //int D[8][7];
    vector< vector<int> > D(m, vector<int>(n));
    //initialize the matrix D, (0, 0), 0-line, 0-row;
    D[0][0] = 0;
    for (int i = 1; i != m; ++i)
        D[i][0] = D[i-1][0] + G_;
    for (int i = 1; i != n; ++i)
        D[0][i] = D[0][i-1] + _G;
    for (int i = 1; i != m; ++i)
        for (int j = 1; j!= n; ++j)
            D[i][j] = MaxValue(D[i-1][j-1] + AlignScore(s[j],t[i]), D[i-1][j] + AlignScore(s[j], ' '), D[i][j-1] + AlignScore(' ', t[i]));
    return (double) (D[m-1][n-1] / (double) (MaxValue(m, n, 0) - 1));
}

//locate the end of matrix, and record ARRAY INDEX;
void EndIndex (double a[][250], int end[2])
{
    int linecounter = 0;
    int rowcounter = 0;
    for (int i = 0; i != 250; ++i)
    {
        for (int j = 0; j != 250; ++j)
        {
            if (a[i][j] != 2)
            {
                linecounter = i;
                rowcounter = j;
            }
            if (end[0] < linecounter)
                end[0] = linecounter;
            if (end[1] < rowcounter)
                end[1] = rowcounter;
        }
    }
}

void Prediction (double GRN[][250])
{
    ofstream testlog;
    testlog.open("/Users/jinyang/Documents/Program C++/SequenceAlignment/SequenceAlignment/test_log.txt");
    string sequence_inquire, sequence_align;
    int endindex[2] = {0, 0};//endindex[0] records end line, endindex[1] records end row;
    //double GRN[250][250];
    double similarity[250] = { 0 };
    //initialize sequence, espetially, add a ' ' in the beginning of each sequence;
    sequence_inquire = " ";
//******************************* test module ********************************************
    ifstream inquire;
    string s;
    inquire.open("/Users/jinyang/Documents/Program C++/SequenceAlignment/SequenceAlignment/test_inquire.txt");
    inquire >> s;
    //cout << s << endl;
    sequence_inquire += s;
    inquire.close();
    //testlog << "The inquire sequence: " << sequence_inquire << endl;
//****************************************************************************************
    //calculate similarities;
    ifstream sequence;
    sequence.open("/Users/jinyang/Documents/Program C++/SequenceAlignment/SequenceAlignment/test_sequence.txt");
    string t;
    for (int gene_index = 0; gene_index != 3; ++gene_index)
    {
        sequence_align = " ";
        sequence >> t;
        sequence_align += t;
        //testlog << "The " << gene_index + 1 << " sequence is: " << sequence_align << endl;
        //sequence_align += getGeneSequence(gene_index);
        similarity[gene_index] = alignSequence(sequence_inquire, sequence_align);
        testlog << similarity[gene_index] << '\t';
    }
    testlog << endl;
    sequence.close();
    //construct GRN with new gene;
//********************************* test module ******************************************
    for (int i = 0; i != 250; ++i)
    {
        for (int j = 0; j != 250; ++j)
        {
            GRN[i][j] = 2;
        }
    }
    GRN[0][2] = 1;
    GRN[1][0] = 1;
    GRN[2][0] = 0;
    GRN[2][1] = -1;
    GRN[2][2] = 0;
//****************************************************************************************
    EndIndex(GRN, endindex);
    cout << endindex[0] << '\t' << endindex[1] << endl;
    //prepare to insert new gene, make (endindex + 1) line and row 0;
    for (int i = 0; i != endindex[0] + 1; ++i)
    {
        GRN[i][endindex[1] + 1] = 0;
    }
    for (int j = 0; j != endindex[1] + 1; ++j)
    {
        GRN[endindex[0] + 1][j] = 0;
    }
    //insert new regulary value to (endindex[1] + 1) column;
    for (int i = 0; i != endindex[0] + 1; ++i)
    {
        int counter = 0;
        for (int j = 0; j != endindex[1] + 1; ++j)
        {
            if (GRN[i][j] != 2)
            {
                counter += 1;
                GRN[i][endindex[1] + 1] += similarity[j] * GRN[i][j];
            }
        }
        GRN[i][endindex[1] + 1] = GRN[i][endindex[1] + 1] / counter;
    }
    //insert new regulary value to (endindex[0] + 1) line;
    for (int j = 0; j != endindex[1] + 1; ++j)
    {
        int counter = 0;
        for (int i = 0; i != endindex[0] + 1; ++i)
        {
            if (GRN[i][j] != 2)
            {
                counter += 1;
                GRN[endindex[0] + 1][j] += similarity[i] * GRN[i][j];
            }
        }
        GRN[endindex[0] + 1][j] = GRN[endindex[0] + 1][j] / counter;
    }
    GRN[endindex[0] + 1][endindex[1] + 1] = 0;
//*********************** test module ****************************************************
    for (int i = 0; i != endindex[0] + 2; ++i)
    {
        for (int j = 0; j != endindex[1] + 2; ++j)
        {
            testlog << GRN[i][j] << '\t';
        }
        testlog << endl;
    }
//****************************************************************************************
    testlog.close();
}
int main()
{
    double GRN[250][250];
    Prediction(GRN);
    return 0;
}
