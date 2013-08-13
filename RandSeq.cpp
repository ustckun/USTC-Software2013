//
//  RandSeq.cpp
//  Random Sequence
//
//  Created by jinyang on 13-8-9.
//  Copyright (c) 2013年 Li Jinyang. All rights reserved.
//

#include "RandSeq.h"
#include <iostream>
#include <ctime>
#include "stdlib.h"
void RandSeq::GenRandSeq(int length){
    //srand((unsigned)time(0));
    for (int i = 0; i != length; ++i) {
        randAAS += GenRandAA();
    }
}
char RandSeq::GenRandAA(){
    switch (rand() % 20) {
        case 0:
            return 'A';
            break;
        case 1:
            return 'R';
            break;
        case 2:
            return 'N';
            break;
        case 3:
            return 'D';
            break;
        case 4:
            return 'C';
            break;
        case 5:
            return 'Q';
            break;
        case 6:
            return 'E';
            break;
        case 7:
            return 'G';
            break;
        case 8:
            return 'H';
            break;
        case 9:
            return 'I';
            break;
        case 10:
            return 'L';
            break;
        case 11:
            return 'K';
            break;
        case 12:
            return 'M';
            break;
        case 13:
            return 'F';
            break;
        case 14:
            return 'P';
            break;
        case 15:
            return 'S';
            break;
        case 16:
            return 'T';
            break;
        case 17:
            return 'W';
            break;
        case 18:
            return 'Y';
            break;
        case 19:
            return 'V';
            break;
    }
    
}

