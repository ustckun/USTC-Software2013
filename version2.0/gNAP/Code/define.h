////////////////////////////////////////////////////////////////////////////////
/// COPYRIGHT NOTICE\n
/// Distribute under BSD License\n
/// Copyright (c) 2013, iGEM Software Team of University of Science and
/// Technology of China\n
/// All rights reserved.
///
/// \file define.h
/// \brief Define the class define.
/// \version 1.0
/// \author Wang Chenkun
/// \date September 2nd, 2013
////////////////////////////////////////////////////////////////////////////////
///
///    This .h file is used to define some statistic value of factors in most command line.\n
///
/// The maximum TF amount which could contain in database
#define TFScale 220
/// The maximum gene amount which could contain in database
#define GENEAM 1800
/// Interactions of ModleNetwork's score
///    \see ModleNetwork
#define NN 100
/// Pets of solving differential equations
///    \see ModleNetwork
#define PETS 128
/// Step of solving differential equations
///    \see ModleNetwork
#define STEP (1.0/PETS)
/// Interactions of PSO
///    \see PSO
#define MAXTIME 100
/// Initial value of each gene in modle
///    \see ModleNetwork
#define INITIALVALUE 2.5
/// Partical number of PSO method
///    \see PSO
#define PARTICLENUM 30
/// Minimum accuracy of PSO
///    \see PSO
#define minAccu 0.01
/// Minimum position value of each partical in PSO
///    \see PSO
#define Pmin -1
/// Maximum position value of each partical in PSO
///    \see PSO
#define Pmax 1
/// Minimun velocity value of each partical in PSO
///    \see PSO
#define Vmin -0.01
/// Maximum velocity value of each partical in PSO
///    \see PSO
#define Vmax 0.01
