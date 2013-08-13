
#include"TFIM.h"
#include"Calculate.h"
#include"Regulation.h"
#include"ReadDNA.h"
#include"Sequence.h"
#include"GRN.h"
#include"PSOPredict.h"
//#include"ReadRegulation.h"
#include"SBOL.h"

void PSOPredict::getRange(double GRN_Origin[scale][scale],int n,Calculate cal)
{
	for(int i=0;i<n;i++)
	{
		if(random(0,1)>0.5)
		{
			GRN_Origin[n][i]=1;
		}
		else
			GRN_Origin[n][i]=-1;
		cout<<random(0,1)<<endl;
	}
	for(int j=0;j<n;j++)
	{
		if(random(0,1)>0.5)
		{
			GRN_Origin[j][n]=1;
		}
		else
			GRN_Origin[j][n]=-1;
	}
	/*for(int i=0;i<N;i++)
	{
		//randomMax[i]=0;
		//randomMin[i]=100;
		for(int j=0;j<N;j++)
		{
			if(GRN_Origin[j][i]==1)
				GRN_Origin[j][i]=4;
			else if(GRN_Origin[j][i]==-1)
				GRN_Origin[j][i]=-4;
			else if(GRN_Origin[j][i]==2)
			{
				if(((double)rand()/RAND_MAX)>0.5)
					GRN_Origin[j][i]=4;
				else
					GRN_Origin[j][i]=-4;
			}
			//cout<<GRN_Origin[j][i]<<"	";
		}
		//cout<<endl;
	}*/
		cal.Network_2(GRN_Origin,n+1);
		//for(int i=0;i<N;i++)
			//cout<<cal.consistence[i]<<endl;
		//getchar();
		for(int m=0;m<N;m++)
		{
			if(randomMax[m]<cal.consistence[m])
				randomMax[m]=cal.consistence[m];
			if(randomMin[m]>cal.consistence[m])
				randomMin[m]=cal.consistence[m];
		}
}

/************************************* randoom ******************************************************/
//creat random double from min to max
double PSOPredict::random(double min, double max)
{
    double ran;
    ran = ((double)rand() / RAND_MAX) * (max - min) + min;
    return ran;
}
/************************************** compareconcentration ****************************************/
//calculate variance with c_new and c_target
//n is gene Amount
double PSOPredict::compareConcentration (double c_new[], double c_target[],int n)
{
    double score = 0;
    for (int i = 0; i != n; ++i)
    {
        if(c_target[i] >= 0)
            score += pow((1 - c_new[i] / c_target[i]), 2);
    }
   return score;
}


/****************************************** test PSO **********************************************/
//v contains x
//double test_Function (vector<double> v, int n)
//{
//   double y = 122;
//  double fit;
// double a = 13;
//    double b = 17;
//    fit = fabs(y - (a * v[0] + b));
//   return fit;
//}
/***************************************** fitness_Function ***************************************/
//figure out stable gene consistence with regulation of x
//GRN is old regulation factor matrix
//target_gene[n] is use's needs
//n is gene Amount
double PSOPredict::fitness_Function (vector<double> x, double GRN[scale][scale], double target_gene[scale], int n, Calculate NEW)
{
   double fitness=0.0;
   vector<double> x1(2*n+1);//a copy of x because x should not be changed

    //transformation of GRN for function_Liao

    for(int i = 0; i < (2*n+1); ++i)
    {
        if((x[i] >= -0.9) && (x[i] <= 0.9))
            x1[i] = 0;
        else if((x[i] < -0.9) && (x[i] >= -1))
            x1[i] = (x[i] * 40.0 + 34.0);
        else if((x[i] > 0.9) && (x[i] <= 1))
            x1[i] = (x[i] * 40.0 - 34.0);
    }

    //add x1 to GRN

    for(int j = 0;j <= n; ++j)
        GRN[n][j] = x1[j];
    for(int j = 0;j < n; ++j)
        GRN[j][n] = x1[2*n-j];
    /*for(int i=0;i<170;++i)
    {
        for(int j=0;j<170;++j)
        {
            cout<<GRN[i][j]<<"  ";
        }
        cout<<endl;
    }*/
    NEW.Network_2(GRN,n);
   // for(int i=0;i<170;i++)cout<<NEW.nong[i]<<endl;    // expect array fulfiled;
    fitness = compareConcentration(NEW.consistence, target_gene,n);
    return fitness;
}

/************************************* min_Line ***************************************************/
//find the best line in fitness
int PSOPredict::min_Line(double fitness[])
{
    double comp = fitness[0];
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
//using PSO method to figure out the best gene regulation
//GRN_Origin is old regulation matrix, target gene is user's needs, n is gene amount
void PSOPredict::algorithm_PSO(double GRN_Origin[scale][scale],double target_gene[N],int n, Calculate NEW)
{
    const int P = 2*n+1;//the number of the unknows
    //const int M = 100;
    const double c1 = 2;
    const double c2 = 2;
    const double Xmax = 1;
    const double Xmin = -1;
    const double Vmax = 0.01;
    const double Vmin = -0.01;
    const int max_literation =100;
    const double min_error = 0.1;
	for(int i=0;i<N;i++)
	{
		for(int j=0;j<N;j++)
			cout<<GRN_Origin[i][j]<<"	";
		cout<<endl;
	}

    srand(time(NULL));
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
            if(GRN_Origin[i][j]!=1.0||GRN_Origin[i][j]!=-1.0)
                GRN[i][j] = 0.0;
            else if(GRN_Origin[i][j] == 1.0)
                GRN[i][j] = 1.0;
            else if(GRN_Origin[i][j] == -1.0)
                GRN[i][j] = -1.0;
            //else if(GRN_Origin[i][j] == 0)
              //  GRN[i][j] = (rand()%2)*8.0 - 4.0;
        }

    vector< vector<double> > x_Matrix(M,vector<double>(P,0));
    vector< vector<double> > v_Matrix(M,vector<double>(P,0));

    //initialize postion and velocity matrix of paticles;

    for(int i = 0; i < M; ++i)
    {
        for(int j = 0; j < P; ++j)
        {
            x_Matrix[i][j] = random(Xmin,Xmax);
            v_Matrix[i][j] = random(Vmin,Vmax);
        }
    }

    vector< vector <double> > pbest(M,vector<double>(P,0));

    //initialize pbest;

    pbest = x_Matrix;
    //vector<double> gbest(P);
	vector<double> best(P);
    double fitness_x[M];
    double fitness_pbest[M];
    double fitness_gbest;

    //initialize 2 array: fitness_x and fitness_pbest;

    int best_Line = 0;
    /*for(int i = 0; i < M ; ++i)
    {
		fitness_pbest[i] = fitness_Function(x_Matrix[i],GRN_Origin,target_gene,n,NEW);
        fitness_x[i] = fitness_pbest[i];
    }*/
    best_Line = min_Line(fitness_pbest);
    fitness_gbest = fitness_pbest[best_Line];
    best = pbest[best_Line];
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
            for(int j = 0; j < P; ++j)
            {
                w = 1;//0.9 - (iter / max_literation * 0.5);
                v_Matrix[i][j] = w * v_Matrix[i][j] + c1 * random(0, 1) * (pbest[i][j] - x_Matrix[i][j]) + c2 * random(0, 1) * (best[j] - x_Matrix[i][j]);

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
            fitness_x[i] = fitness_Function(x_Matrix[i], GRN, target_gene,n, NEW);
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

        best_Line = min_Line(fitness_pbest);
        best = pbest[best_Line];
        fitness_gbest = fitness_pbest[best_Line];
        iter = iter + 1;
        //cout << (double)clock()/CLOCKS_PER_SEC << endl;
    }
	for(int i=0;i<n;i++)
		edBest.push_back(best[i]);
	for(int i=0;i<n;i++)
		toBest.push_back(best[2*n-i]);
    /*for(int i = 0; i < P;++i)
    {
          cout << gbest[i];
          cout << "\n";
    }
    cout << fitness_gbest;
    int aha;
    cin >> aha;*/
}

PSOPredict::PSOPredict()
{
	for(int i=0;i<N;i++)
		targetGene[i]=-1;
}


void PSOPredict::Filter(int n)
{
	/*for(int i=0;i<n;i++)
	{
		if(toBest[i]>0&&toBest[i]<0.2)
			toPick[i]=1;
		else if(toBest[i]>0.2&&toBest[i]<0.4)
			toPick[i]=2;
		else if(toBest[i]>0.4&&toBest[i]<0.6)
			toPick[i]=3;
		else if(toBest[i]>0.6&&toBest[i]<0.8)
			toPick[i]=4;
		else if(toBest[i]>0.8&&toBest[i]<1)
			toPick[i]=5;
	}
	for(int i=0;i<n;i++)
	{
		if(edBest[i]>0&&edBest[i]<0.2)
			edPick[i]=1;
		else if(edBest[i]>0.2&&edBest[i]<0.4)
			edPick[i]=2;
		else if(edBest[i]>0.4&&edBest[i]<0.6)
			edPick[i]=3;
		else if(edBest[i]>0.6&&edBest[i]<0.8)
			edPick[i]=4;
		else if(edBest[i]>0.8&&edBest[i]<1)
			edPick[i]=5;
	}*/
	toPick=toBest;
	edPick=edBest;
}