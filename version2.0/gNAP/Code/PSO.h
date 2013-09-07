////////////////////////////////////////////////////////////////////////////////
/// COPYRIGHT NOTICE\n
/// Distribute under BSD License\n
/// Copyright (c) 2013, iGEM Software Team of University of Science and
/// Technology of China\n
/// All rights reserved.
///
/// \file PSO.h
/// \brief Define the class PSO.
/// \version 1.0
/// \author Wang Chenkun
/// \date September 2nd, 2013
////////////////////////////////////////////////////////////////////////////////

#ifndef __USTC_Software__PSO__
#define __USTC_Software__PSO__


#include <vector>
#include <cstdlib>

#include "define.h"
#include "ModleNetwork.h"
///    Use PSO to predict interactions between gene needed to put into GRN and original network.\n
///
///    PSO is Particle Swarm Optimization which is be used to find the best regulation
///    fitting to users' goal.\n

class PSO
{
public:
/// Initialize the PSO object
///    Using class ModleNetwork to figure out the starting value of each genes.
///    \param an object of class ModleNetwork
///    \param row number of GRN
///    \param column number of GRN
///    \see ModleNetwork
    PSO(ModleNetwork New, int row,int column);
/// Target gene which needed to change
    double target[GENEAM];
/// Main function which use PSO method to predict interactions
///    This function using getMinLine(), getFitness(), getVariance() and random().
///    \param an object of class ModleNetwork
///    \param row number of GRN
///    \param column number of GRN
///    \see ModleNetwork
///    \see getMinLine
///    \see getFitness
///    \see getVariance
///    \see random
    void getPrediction(ModleNetwork New,int row,int column);
/// New gene interact to genes in original GRN
///    This vector contain the strength of interaction
    vector<double> toPick;
/// New gene is interacted by genes in original GRN
///    This vector contain the strength of interaction
    vector<double> edPick;
/// Get range of each gene's strength of expression
///    Use random regulation to figure out the Maximum and Minimum expression strength.
///    Those range have been put into random_max and random_min.
///    \param row number of GRN
///    \param column number of GRN
///    \param an object of class ModleNetwork
///    \see ModleNetwork
    void getRange(int row,int column,ModleNetwork cal);
/// Max expression value of genes in original GRN
///    These value is used to set the users' target genes which need high expression.
    double random_max[GENEAM];
/// Min expression value of genes in original GRN
///    These value is used to set the users' target genes which need low expression.
    double random_min[GENEAM];
/// Filt predicted regualtion
///    Classify the interactions to different degrees.
///    \param row number of GRN
///    \param column number of GRN
    void Filter(int row, int column);
private:
/// Store GRN in this vector and easy using
    double **temp_GRN;
/// Find out the minimum number in an array
///    \param variance for different particles in PSO method
///    \param particle number
///    \return this minimum line number
    int getMinLine(double A[GENEAM],int column);
/// Get fitness for each new regulation
///    \param a vector which contains the interactions between new gene and original genes
///    \param an object of class ModleNetwork
///    \param row number of GRN
///    \param column number of GRN
///    \return the variance of prediction
///    \see getVariance
    double getFitness(vector<double> row_column_matrix,ModleNetwork New,int row,int column);
/// Figure out the variance between target and prediction
///    \param gene expression strength array
///    \return the variance between prediction and users' goal
    double getVariance(double A[GENEAM]);
/// Produce a random "double" figure from "min" to "max"
///    \param Random's lower limit
///    \param Random's higher limit
///    \return random figure
    double random(double min,double max);
};

#endif
