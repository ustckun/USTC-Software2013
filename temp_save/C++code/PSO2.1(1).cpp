#include <iostream>
#include <vector>
#include <math.h>
#include <string>
#include <algorithm>
#include <cstdlib>
#include <time.h>

using namespace std;

//16807 random numbers creator
/*int schrage(int seed)
{
    int x,a = 16807, m = 2147483647, q = 127773, r = 2836;
    x = a * (seed % q) - r * floor(seed / q);
    if(x < 0)
        x = x + m;
    return x;
}*/
/************************************* randoom ******************************************************/
double random(double min, double max)
{
    double ran;
    ran = ((double)rand() / RAND_MAX) * (max - min) + min;
    return ran;
}
/************************************** compareconcentration ***************************************/
//double compareConcentration (double c_new[], double c_target[], int n)
//{
//    double score = 0;
//    for (int i = 0; i != n; ++i)
//    {
//        if(c_target[i] <(/*number choosen for the gene concentration user do not require*/) )
//            score += pow((1 - c_new[i] / c_target[i]), 2);
//   )
//   return score;
//}


/****************************************** test PSO **********************************************/
//v contains x
double test_Function (vector<double> v, int n)
{
    double y = 122;
    double fit;
    double a = 13;
    double b = 17;
    fit = fabs(y - (a * v[0] + b));
    return fit;
}
/***************************************** fitness_Function ***************************************/

//double fitness_Function (vector<double> x, double GRN[250][250], int n)
//{
//   int fitness = 0;
//   for(int j = 0;j <= n; ++j)
//        GRN[n][j] = x[j];
//    for(int j = 0;j < n; ++j)
//        GRN[j][n] = x[2*n-j];
//    function_Liao(, /* an array carries gene concentration message */); // expect array fulfiled;
    //locate target gene, get concentration of target gene
//    fitness = compareConcentration(c_new, c_target, n);
//    return fitness;
//}

/************************************* min_Line ***************************************************/
int min_Line(double fitness[], int M)
{
    int comp = fitness[0];
    int minLine;
    for(int i = 0; i < M; ++i)
        if(comp > fitness[i])
        {
            comp = fitness[i];
            minLine = i;
        }
    return minLine;
}
/************************************ algorithm_PSO **********************************************/
void *algorithm_PSO(/*double GRN[250][250], *//*int const n*//* numbers of gene */)
{
    /*const int N = 2*n+1;*/
    const int N = 1;
    const int M = 100;
    const double c1 = 2;
    const double c2 = 2;
    const double Xmax = 1000;
    const double Xmin = 1;
    const double Vmax = 0.01;
    const double Vmin = -0.01;
    const int max_literation = 1000;
    const double min_error = 0.001;
    vector< vector<double> > x_Matrix(M,vector<double>(N,0));
    vector< vector<double> > v_Matrix(M,vector<double>(N,0));
    srand(time(NULL));
    //initialize postion and velocity matrix of paticles;

    for(int i = 0; i < M; ++i)
    {
        for(int j = 0; j < N; ++j)
        {
            x_Matrix[i][j] = random(Xmin,Xmax);
            v_Matrix[i][j] = random(Vmin,Vmax);
        }
    }

    vector< vector <double> > pbest(M,vector<double>(N,0));
    //initialize pbest;
    pbest = x_Matrix;
    vector<double> gbest(N);
    double fitness_x[M];
    double fitness_pbest[M];
    double fitness_gbest;
    //initialize 2 array: fitness_x and fitness_pbest;
    int best_Line = 0;
    for(int i = 0; i < M ; ++i)
    {
        fitness_pbest[i] = /*fitness_Function(x_Matrix[i], GRN, n);*/test_Function(x_Matrix[i],N);
        fitness_x[i] = fitness_pbest[i];
    }
    best_Line = min_Line(fitness_pbest, M);
    fitness_gbest = fitness_pbest[best_Line];
    gbest = pbest[best_Line];
    /********************************************** all initializaton done. ************************************/
    // Now:
    //     fulfill x_Matrix (i.e. paticle-position matrix) with random numbers;
    //     fulfill v_Matrix (i.e. paticle-velocity matrix) with random numbers;
    //     fulfill pbest (i.e. partial best particle-position matrix) with x_Matrix;
    //     fulfill fitness_x (i.e. the fitness of each particle-position);
    //     fulfill fitness_pbest(i.e. the fitness of each particle-position in history) with fitness_x;
    //     fulfill gbest (i.e. global best particle-position) with the best line of pbest;
    //
    /**********************************************************************************************************/

    int iter = 0;  //iterator;
    double w = 0;
    while((fitness_gbest > min_error )/*&&(iter < max_literation )*/)
    {
        //move paticles;
        for(int i = 0; i < M; ++i)
        {
            for(int j = 0; j < N; ++j)
            {
                w = 1/*0.9 - (iter / max_literation * 0.5)*/;
                v_Matrix[i][j] = w * v_Matrix[i][j] + c1 * random(0, 1) * (pbest[i][j] - x_Matrix[i][j]) + c2 * random(0, 1) * (gbest[j] - x_Matrix[i][j]);
                if(v_Matrix[i][j] > Vmax)
                    v_Matrix[i][j] = Vmax;
                else if(v_Matrix[i][j] < Vmin)
                    v_Matrix[i][j] = Vmin;
                x_Matrix[i][j] = x_Matrix[i][j] + v_Matrix[i][j];
                if(x_Matrix[i][j] > Xmax)
                    x_Matrix[i][j] = Xmax;
                else if(x_Matrix[i][j] < Xmin)
                    x_Matrix[i][j] = Xmin;
            }
        }
        //calculate fitness of x_Matrix
        for(int i = 0; i < M; ++i)
        {
            fitness_x[i] = /*fitness_Function(x_Matrix[i], GRN, n);*/test_Function(x_Matrix[i], N);
        }
        //update pbest array;
        for(int i = 0; i < M ; ++i)
        {
            if(fitness_pbest[i] > fitness_x[i])
            {
                pbest[i] = x_Matrix[i];
                fitness_pbest[i] = fitness_x[i];
            }
        }
        //update gbest array;
        best_Line = min_Line(fitness_pbest, M);
        gbest = pbest[best_Line];
        fitness_gbest = fitness_pbest[best_Line];
        iter = iter + 1;
    }
      cout << gbest[0] << "\n";
      cout << iter;
}


int main()
{
    double best_answer;
    algorithm_PSO();

   /* best_answer = *algorithm_PSO();
    cout << best_answer; */
}
