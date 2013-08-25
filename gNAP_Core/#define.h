#ifndef __USTC_Software__GeneIM__
#define __USTC_Software__GeneIM__

#include <ctime>
#include <cmath>
#include <iomanip>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <map>

using namespace std;

#define CALNO 100                                   //count number
#define PARTICLENUM 30										//PSO number
#define PETS 128                                
#define STEP (1.0/PETS)                           //step length
#define CALMAXTIME 100                            //maxtime
#define INITIALVALUE 2.5                          //initial value
#define DIMENS 100    
#define gap 0
#define aM 4700
#define TFScale 220
#define GENEAM 1800