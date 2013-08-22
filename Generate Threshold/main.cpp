//
//  main.cpp
//  Generate Threshold
//
//  Created by jinyang on 13-8-9.
//  Copyright (c) 2013年 Li Jinyang. All rights reserved.
//

#include <iostream>
#include "RandSeq.h"
#include "SeqAlign.h"
#include <ctime>
#include <stdlib.h>
#include <fstream>
using namespace std;
int main()
{
    string subject_1 = " VHWQTHTVFNQPIPLNNSNLYLSDGALCEAVTREGAGWDSDFLASIGQQLGTAESLELGRLANVNPPELLRYDAQGRRLDDVRFHPAWHLLMQALCTNRVHNLAWEEDARSGAFVARAARFMLHAQVEAGSLCPITMTFAATPLLLQMLPAPFQDWTTPLLSDRYDSHLLPGGQKRGLLIGMGMTEKQGGSDVMSNTTRAERLEDGSYRLVGHKWFFSVPQSDAHLVLAQTAGGLSCFFVPRFLPDGQRNAIRLERLKDKLGNRSNASCEVEFQDAIGWLLGLEGEGIRLILKMGGMTRFDCALGSHAMMRRAFSLAIYHAHQRHVFGNPLIQQPLMRHVLSRMALQLEGQTALLFRLARAWDRRADAKEALWARLFTPAAKFVICKRGMPFVAEAMEVLGGIGYCEESELPRLYREMPVNSIWEGSGNIMCLDVLRVLNKQAGVYDLLSEAFVEVKGQDRYFDRAVRRLQQQLRKPAEELGREITHQLFLLGCGAQMLKYASPPMAQAWCQVMLDTRGGVRLSEQIQNDLLLRATGGVCV";
    int sub_1_size = 541;
    string subject_2 = " MTEVRRRGRPGQAEPVAQKGAQALERGIAILQYLEKSGGSSSVSDISLNLDLPLSTTFRLLKVLQAADFVYQDSQLGWWHIGLGVFNVGAAYIHNRDVLSVAGPFMRRLMLLSGETVNVAIRNGNEAVLIGQLECKSMVRMCAPLGSRLPLHASGAGKALLYPLAEEELMSIILQTGLQQFTPTTLVDMPTLLKDLEQARELGYTVDKEEHVVGLNCIASAIYDDVGSVVAAISISGPSSRLTEDRFVSQGELVRDTARDISTALGLKAHP";
    int sub_2_size = 271;
    srand((unsigned)time(0));
    RandSeq query[100];
    RandSeq subject;
    subject.GenRandSeq(541);
    for (int i = 0; i != 100; ++i) {
        query[i].GenRandSeq(541);
    }
    SeqAlign test[100];
    for (int i = 0; i != 100; ++i) {
        test[i].GetSimi(subject.randAAS, sub_1_size, query[i].randAAS, sub_1_size);
    }
    
    //statistic analysis;
    int stat[1000] = { 0 };
    int sizer = 0;
    for (int i = 0; i != 100; ++i) {
        sizer = (int)(test[i].simi * 1000);
        /*for (int j = 0; j != 1000; ++j) {
            if (sizer == (j + 1)) {
                stat[j] += 1;
            }
        }*/
        stat[sizer] += 1;
    }
    double sum = 0;
    for (int i = 0; i != 100; ++i) {
        sum += test[i].simi;
    }
    int stat_2[1000] = { 0 };
    for (int i = 0; i != 1000; ++i) {
        stat_2[i] = 10;
    }
    stat_2[(int)(sum * 10)] = 0;
    cout << sum * 10 <<'\t' << (int)(sum * 10) << endl;
    
   
    //use system time as file name;
    ofstream outfile;
    time_t nowtime = time(NULL);
    struct tm *p;
    p = gmtime(&nowtime);
    char filename_1[256] ={ 0 };
    string filename = "/Users/jinyang/Desktop/Parameter Data Test/random ";
    sprintf(filename_1, "%d-%d %d%02d", 1 + p -> tm_mon, p -> tm_mday, 8 + p -> tm_hour, p -> tm_min);
    filename += filename_1;
    filename += ".txt";
    outfile.open(filename);
    for (int i = 0; i != 1000; ++i) {
        outfile << i + 1 << '\t' << stat[i] << '\t' << stat_2[i] << endl;
    }
    outfile <<"Running time: " << (double)clock()/CLOCKS_PER_SEC << endl;
    outfile << endl;
}

