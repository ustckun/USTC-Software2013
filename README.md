USTC-Software 2013
=================

We are USTC-Software, a team from University of Science and Technology of China. We will be competing in iGem 2013!

###Introduction###
Our application aims to simulate genetic networks. The application analyzes the stability of genetic networks after introduction of exogenous genes. Meanwhile, given the original network and specific purposes, the application traces the regulative process back and gives possible regulative patterns.


## gNAP: Genetic Network Analyse and Predict ##

This software contains four parts, dealing with separate functions in forward and backward modeling of GRN(Genetic Regulatory Network) analyse.
#####1. Start#####
#####2. Monitor#####
#####3. Result#####
#####4. Display#####

###Start###

**Start** is used to prepare for the later analysis and prediction. In this part, users could input their database downloaded on Internet and sequences of exogenous gene which is needed to analyse. Also, if not input sequence in **Start**, users could also use the "Predict" function in next part.

###Monitor###

**Monitor** undertakes several functions of our software as the core methods of **gNAP**. First of them is **Analyse** function which figure out the network change when input an exogenous gene. In the same time a score presenting stablility of new GRN by statist stable time and value variation for lots of times. **Analyse** result could be saw intuitively in **Result** part next. Secondly, **Predict** function use target gene exprssion to figure out possible interaction whose result could also receive in **Result**.

###Result###

**Result** is a output part which contains all results of operations used. It is easy to read each gene's information and changing consequence in this part. What's more, all gene information could be output in [SBOL](http://www.sbolstandard.org/).

###Display###

**Display** 

This software can be built on Windows, Linux and MacOS operating platform.

For more information, please refer to our [wiki page](http://2013.igem.org/Team:USTC-Software).

####Source Files####
**gNAP** floder contains the command line source files in **Code** floder and GUI source files.The command line source files are written in C++ language and visualization parts are written in Java language. Both of them can be complied across platforms.

The GUI source files are written in C++ language with Qt Creator, it can also be compiled across platforms using Qt 5.1.0, which can be found [here](http://qt-project.org/downloads).

####Database####
The example database has been put into **data** floder and it can also be downloaded from RegulonDB, which can be found [here](http://regulondb.ccg.unam.mx/menu/download/datasets/index.jsp).

The data which used in **gNAP** is flexible. All database in those form could be read in our software.

###Contacts###

For any questions, feel free to contact:

Chenkun Wang(ustckun@gmail.com)

Jinyang Li(jinyangustc@gmail.com)