////////////////////////////////////////////////////////////////////////////////
/// COPYRIGHT NOTICE\n
/// Distribute under BSD License\n
/// Copyright (c) 2013, iGEM Software Team of University of Science and
/// Technology of China\n
/// All rights reserved.
///
/// \file SBOL.h
/// \brief Define the class SBOL.
/// \version 1.0
/// \author Wang Chenkun
/// \date September 2nd, 2013
////////////////////////////////////////////////////////////////////////////////


#ifndef __USTC_Software__SBOL__
#define __USTC_Software__SBOL__

#include <fstream>
#include <iostream>
#include <string>

using namespace std;
///    Creat SBOL files outside based on gene information.


class SBOL
{
public:
/// Create SBOL files named by gene name
///    \param string of gene name
///    \param string of RegulonDB ID
///    \param string of left position
///    \param string of right position
///    \param string of gene description
///    \param string of gene sequence
    void CreatSBOL(string gene_name,string ID,string left,string right,string description,string seq);
private:
/// head of SBOL files
    string head;
/// Combine SBOL detail and its lable
///    \param lable of info
///    \param lable of detail about lable
///    \return string in the formart of lable and detail
    string Combine(string title,string detail);
/// Formart of Start lable
///    \param string of lable in each line
///    \return string contain both lable and start form
    string FormartStart(string a);
/// Formart of End lable
///    \param string of lable in each line
///    \return string contain both lable and end form
    string FormartEnd(string b);
};


#endif
