//
//  RandSeq.h
//  Random Sequence
//
//  Created by jinyang on 13-8-9.
//  Copyright (c) 2013年 Li Jinyang. All rights reserved.
//

#ifndef __Random_Sequence__RandSeq__
#define __Random_Sequence__RandSeq__

#include <iostream>
class RandSeq{
public:
    RandSeq(){
        randAAS = " ";
    }
    void GenRandSeq(int length);
    std::string randAAS;
private:
    char GenRandAA();
};
#endif /* defined(__Random_Sequence__RandSeq__) */
