#include"TFIM.h"
#include"Calculate.h"
#include"Regulation.h"
#include"ReadDNA.h"
#include"Sequence.h"
#include"GRN.h"
#include"PSOPredict.h"






/************************************* randoom ******************************************************/
double PSOPredict::random(double min, double max)
{
    double ran;
    ran = ((double)rand() / RAND_MAX) * (max - min) + min;
    return ran;
}
/************************************** compareconcentration ****************************************/
double PSOPredict::compareConcentration (double c_new[], double c_target[])
{
    double score = 0;
    for (int i = 0; i != N; ++i)
    {
        if(c_target[i] >= 0)
            score += pow((1 - c_new[i] / c_target[i]), 2);
    }
   return score;
}


/****************************************** test PSO **********************************************/
//v contains x
//double test_Function (vector<double> v, int GA)
//{
//   double y = 122;
//  double fit;
// double a = 13;
//    double b = 17;
//    fit = fabs(y - (a * v[0] + b));
//   return fit;
//}
/***************************************** fitness_Function ***************************************/

double PSOPredict::fitness_Function (vector<double> x, double GRN[scale][scale], double target_gene[N],int GA)
{
   double fitness ;
   fitness = 0.0;
   vector<double> x1(2*GA+1);//a copy of x because x should not be changed

    //transformation of GRN for function_Liao

    for(int i = 0; i < (2*GA+1); ++i)
    {
        if((x[i] >= -0.9) && (x[i] <= 0.9))
            x1[i] = 0;
        else if((x[i] < -0.9) && (x[i] >= -1))
            x1[i] = (x[i] * 40.0 + 34.0);
        else if((x[i] > 0.9) && (x[i] <= 1))
            x1[i] = (x[i] * 40.0 - 34.0);
    }

    //add x1 to GRN

    for(int j = 0;j <= GA; ++j)
        GRN[GA][j] = x1[j];
    for(int j = 0;j < GA; ++j)
        GRN[j][GA] = x1[2*GA-j];
    /*for(int i=0;i<170;++i)
    {
        for(int j=0;j<170;++j)
        {
            cout<<GRN[i][j]<<"  ";
        }
        cout<<endl;
    }*/
    Calculate NEW;
    NEW.Network_2(GRN,GA);
   // for(int i=0;i<170;i++)cout<<NEW.nong[i]<<endl;    // expect array fulfiled;
	fitness = compareConcentration(NEW.consistence, target_gene);
    return fitness;
}

/************************************* min_Line ***************************************************/
int PSOPredict::min_Line(double fitness[], int M)
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
void PSOPredict::algorithm_PSO(double target_gene[],int GA,float GRN_Origin[N][N])
{
    const int n = 2*GA+1;//the number of the unknows
    const int M = 100;
    const double c1 = 2;
    const double c2 = 2;
    const double Xmax = 1;
    const double Xmin = -1;
    const double Vmax = 0.01;
    const double Vmin = -0.01;
    const int max_literation =100;
    const double min_error = 0.1;

    srand(time(NULL));
    double GRN_Origin[scale][scale];
    double GRN[scale][scale];
    //read GRN from regulation1.0

    /*ifstream infile;
    infile.open("regulation1.txt");
    for(int i = 0; i < scale; ++i)
        for(int j = 0; j < scale; ++j)
            {
                infile >> GRN_Origin[i][j];
                //cout<<GRN_Origin[i][j];
            }*/

    //transformation of GRN for function_Liao

    for(int i = 0; i < scale; ++i)
        for(int j = 0; j < scale; ++j)
        {
            if(GRN_Origin[i][j] == 2.0)
                GRN[i][j] = 0.0;
            else if(GRN_Origin[i][j] == 1.0)
                GRN[i][j] = 4.0;
            else if(GRN_Origin[i][j] == -1.0)
                GRN[i][j] = -4.0;
            else if(GRN_Origin[i][j] == 0)
                GRN[i][j] = (rand()%2)*8.0 - 4.0;
        }

    vector< vector<double> > x_Matrix(M,vector<double>(n,0));
    vector< vector<double> > v_Matrix(M,vector<double>(n,0));

    //initialize postion and velocity matrix of paticles;

    for(int i = 0; i < M; ++i)
    {
        for(int j = 0; j < n; ++j)
        {
            x_Matrix[i][j] = random(Xmin,Xmax);
            v_Matrix[i][j] = random(Vmin,Vmax);
        }
    }

    vector< vector <double> > pbest(M,vector<double>(n,0));

    //initialize pbest;

    pbest = x_Matrix;
    vector<double> gbest(n);
    double fitness_x[M];
    double fitness_pbest[M];
    double fitness_gbest;

    //initialize 2 array: fitness_x and fitness_pbest;

    int best_Line = 0;
    for(int i = 0; i < M ; ++i)
    {
        fitness_pbest[i] = fitness_Function(x_Matrix[i], GRN, target_gene,GA);
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
    while((fitness_gbest > min_error )&&(iter < max_literation ))
    {
        //move paticles;

        for(int i = 0; i < M; ++i)
        {
            for(int j = 0; j < n; ++j)
            {
                w = 1;//0.9 - (iter / max_literation * 0.5);
                v_Matrix[i][j] = w * v_Matrix[i][j] + c1 * random(0, 1) * (pbest[i][j] - x_Matrix[i][j]) + c2 * random(0, 1) * (gbest[j] - x_Matrix[i][j]);

                //limit v

                if(v_Matrix[i][j] > Vmax)
                    v_Matrix[i][j] = Vmax;
                else if(v_Matrix[i][j] < Vmin)
                    v_Matrix[i][j] = Vmin;
                x_Matrix[i][j] = x_Matrix[i][j] + v_Matrix[i][j];

                //limit x

                if(x_Matrix[i][j] > Xmax)
                    x_Matrix[i][j] = Xmax;
                else if(x_Matrix[i][j] < Xmin)
                    x_Matrix[i][j] = Xmin;
            }
        }

        //calculate fitness of x_Matrix

        for(int i = 0; i < M; ++i)
        {
            fitness_x[i] = fitness_Function(x_Matrix[i], GRN, target_gene,GA);
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
        //cout << (double)clock()/CLOCKS_PER_SEC << endl;
    }
	best=gbest;
    /*for(int i = 0; i < GA;++i)
    {
          cout << gbest[i];
          cout << "\n";
    }
    cout << fitness_gbest;*/
}