//
//  main.cpp
//  RankingEDAsCEC
//
//  Created by Josu Ceberio Uribe on 11/19/13.
//  Copyright (c) 2013 Josu Ceberio Uribe. All rights reserved.
//
#include <sys/time.h>
#include <iostream>
#include <iomanip>
#include "RankingEDA.h"
#include "InfiniteRankingEDA.h"
#include "PBP.h"
#include "PFSP.h"
#include "FSSP.h"
#include "LOP.h"
#include "QAP.h"
#include "TSP.h"
#include "API.h"
#include "SDF.h"
#include "GAP.h"
#include "SOP.h"

#include "VRP.h"
// #include "VRPTW.h"

//It is the instance of the problem to optimize.
PBP * PROBLEM;

//The type of the problem to solve.
char PROBLEM_TYPE[10];

//The type of the model to use in the EDA.
char MODEL_TYPE[10];

//The type of the metric to use in the model of the EDA.
char METRIC_TYPE[10];

// The individual size.
int PROBLEM_SIZE;
int POP_SIZE;

double B_RATIO;
double CUT_POINT_COUNT;
int NUMBER_OF_EDGE;

int PREVIOUS_SAMPLED_REFERENCE_COUNT;

// Number of evaluations performed.
long int EVALUATIONS = 0;

// Convergence evaluation of the best fitness.
long int CONVERGENCE_EVALUATIONS = 0;

// Maximum number of evaluations allowed performed.
long int MAX_EVALUATIONS = 0;

// Name of the file where the result will be stored.
char RESULTS_FILENAME[50];

// Name of the file where the instances is stored.
char INSTANCE_FILENAME[50];

// The seed asigned to the process
int SEED;

//Determines if a solution needs to be inverted before evaluating.
int INVERSE;


/*
 * Get next command line option and parameter
 */
int GetOption (int argc, char** argv, char* pszValidOpts, char** ppszParam)
{
	
    static int iArg = 1;
    char chOpt;
    char* psz = NULL;
    char* pszParam = NULL;
	
    if (iArg < argc)
    {
        psz = &(argv[iArg][0]);
		
        if (*psz == '-' || *psz == '/')
        {
            // we have an option specifier
            chOpt = argv[iArg][1];
			
            if (isalnum(chOpt) || ispunct(chOpt))
            {
                // we have an option character
                psz = strchr(pszValidOpts, chOpt);
				
                if (psz != NULL)
                {
                    // option is valid, we want to return chOpt
                    if (psz[1] == ':')
                    {
                        // option can have a parameter
                        psz = &(argv[iArg][2]);
                        if (*psz == '\0')
                        {
                            // must look at next argv for param
                            if (iArg+1 < argc)
                            {
                                psz = &(argv[iArg+1][0]);
                                if (*psz == '-' || *psz == '/')
                                {
                                    // next argv is a new option, so param
                                    // not given for current option
                                }
                                else
                                {
                                    // next argv is the param
                                    iArg++;
                                    pszParam = psz;
                                }
                            }
                            else
                            {
                                // reached end of args looking for param
                            }
							
                        }
                        else
                        {
                            // param is attached to option
                            pszParam = psz;
                        }
                    }
                    else
                    {
                        // option is alone, has no parameter
                    }
                }
                else
                {
                    // option specified is not in list of valid options
                    chOpt = -1;
                    pszParam = &(argv[iArg][0]);
                }
            }
            else
            {
                // though option specifier was given, option character
                // is not alpha or was was not specified
                chOpt = -1;
                pszParam = &(argv[iArg][0]);
            }
        }
        else
        {
            // standalone arg given with no option specifier
            chOpt = 1;
            pszParam = &(argv[iArg][0]);
        }
    }
    else
    {
        // end of argument list
        chOpt = 0;
    }
	
    iArg++;
	
    *ppszParam = pszParam;
    return (chOpt);
}

/*
 * Help command output.
 */
void usage(char *progname)
{
    cout << "Algorithm -i <instance_name> -o <results_name> -s <seed> -t <problem_type> -m <model_type> -d <metric> -v <inverse> -p <population_size> -e <max_evaluations> -b <b_ratio> -c <cut_point_count> -r <previous sampled reference count>" <<endl;
    cout <<"   -i File name of the instance.\n"<<endl;
    cout <<"   -o Name of the file to store the results.\n"<<endl;
    cout <<"   -s Seed to be used for pseudo-random numbers generator.\n"<<endl;
    cout <<"   -t problem_type (TSP, QAP, LOP, PFSP, FSSP, SDF, GAP or API).\n"<<endl;
    cout <<"   -m model_type (M, GM, BW, MR, EO).\n"<<endl;
    cout <<"   -d metric (K (Kendall), C (Cayley), U (Ulam)).\n"<<endl;
    cout <<"   -v inverse (0<- no inverse, 1<- inverse)."<<endl;
    cout <<"   -p population size."<<endl;
    cout <<"   -e max_evaluation."<<endl;
    cout <<"   -b b_ratio. (virtual sampling)"<<endl;
    cout <<"   -c cut_point_count. (MR model only)"<<endl;
    cout <<"   -r previous sampled reference count. (AMR model only)"<<endl;
    cout <<"   -n number of edge. (LLMST, RMST, SRMST model only)"<<endl;

}

/*
 * Obtaint the execution parameters from the command line.
 */
bool GetParameters(int argc,char * argv[])
{
	char c;
    if(argc==1)
    {
    	usage(argv[0]);
        return false;
    }
	char** optarg;
	optarg = new char*[argc];
    while ((c = GetOption (argc, argv, "s:h:t:o:i:m:d:v:p:e:b:c:r:n:",optarg)) != '\0')
    {
    	switch (c)
    	{
                
            case 'h' :
                usage(argv[0]);
                return false;
                break;
                
            case 's' :
                SEED = atoi(*optarg);
                break;
            
            case 'p' :
                POP_SIZE = atoi(*optarg);
                break;

            case 'e' :
                MAX_EVALUATIONS = atoi(*optarg);
                break;	

            case 'b' :
                B_RATIO = atof(*optarg);
                break;	

            case 'c' :
                CUT_POINT_COUNT = atof(*optarg);
                break;	
            
            case 'r' :
                PREVIOUS_SAMPLED_REFERENCE_COUNT = atoi(*optarg);
                break;

           	case 't':
                strcpy(PROBLEM_TYPE, *optarg);
                break;
                
            case 'o' :
                strcpy(RESULTS_FILENAME, *optarg);
                break;
                
			case 'i':
                strcpy(INSTANCE_FILENAME, *optarg);
                break;

            case 'm':
                strcpy(MODEL_TYPE, *optarg);
                break;
                
           	case 'd':
                strcpy(METRIC_TYPE, *optarg);
                break;
                
           	case 'v':
                INVERSE=atoi(*optarg);
                break;
            
            case 'n':
                NUMBER_OF_EDGE=atoi(*optarg);
                break;
		
        }
    }

    delete [] optarg;
    
	return true;
}

/*
 * Writes the results of the execution.
 */
void WriteResults1(long int best_fitness, int * best, long int convergeNFE, long int evaluations, long double time_interval){

    ofstream output_file;
    output_file.open(RESULTS_FILENAME);
    output_file<<"Best fitness: "<<setprecision(15)<<best_fitness<<endl;
    output_file<<"Best solution: ";
    for (int i=0;i<PROBLEM_SIZE;i++)
        output_file<<best[i]<<" ";
    output_file<<endl;
    output_file<<"Evaluations performed: "<<setprecision(15)<<evaluations<<endl;
    output_file<<"Converge NFE: "<<setprecision(15)<<convergeNFE<<endl;
    output_file<<"Time consumed: "<<time_interval<<endl;
    output_file.close();
}

/*
 * Writes the results of the execution.
 */
void WriteResults2(long int best_fitness, int * best, long int evaluations, long double time_interval, int distance_to_optimum){
    
    ofstream output_file;
    output_file.open(RESULTS_FILENAME);
    output_file<<"Best fitness: "<<setprecision(15)<<best_fitness<<endl;
    output_file<<"Best solution: ";
    for (int i=0;i<PROBLEM_SIZE;i++)
        output_file<<best[i]<<" ";
    output_file<<endl;
    output_file<<"Distance optimum: "<<distance_to_optimum<<endl;
    output_file<<"Evaluations performed: "<<setprecision(15)<<evaluations<<endl;
    output_file<<"Time consumed: "<<time_interval<<endl;
    output_file.close();
}

/*
 * Reads the problem info of the instance set.
 */
PBP * GetProblemInfo(string problemType, string filename)
{
    PBP * problem;
    if (problemType=="PFSP")
        problem= new PFSP();
    else if (problemType=="FSSP")
        problem= new FSSP();
    else if (problemType=="TSP")
        problem= new TSP();
    else if (problemType=="QAP")
        problem= new QAP();
    else if (problemType=="LOP")
        problem= new LOP();
    else if (problemType=="API")
        problem= new API();
    else if (problemType=="SDF") // Self Define Function
        problem= new SDF();
    else if (problemType=="SOP") // Self Define Function
        problem= new SOP();
    else if (problemType=="GAP") // GA project fitness Function
        problem= new GAP();
    else if (problemType=="VRP") // Vehicle Routing Problem
        problem= new VRP();
    // else if (problemType=="VRPTW") // Vehicle Routing Problem with Time Window 
    //     problem= new VRPTW();
    else{
        cout<<"Wrong problem type was specified."<<endl;
        exit(1);
    }
    
    //Read the instance.
    int problem_size= problem->Read(filename);
    PROBLEM_SIZE=problem_size;
    // cout << "problem_size: " << problem_size << endl;
    
	return problem;
}

/*
 * Determines if the problem type of the instance to optimize is an Artifical Permutation Instance (API).
 */
bool IsAPI(string problemType){
    return (problemType=="API");
}

/*
 * Main function.
 */
int main(int argc, char * argv[])
{
    //Initialize time variables.

    struct timeval tim;
    gettimeofday(&tim, NULL);
    double t1=tim.tv_sec+(tim.tv_usec/1000000.0);

    //Get parameters
	if(!GetParameters(argc,argv)) return -1;

    //Set seed
    srand(SEED*1000);
    // mt19937 mt_rand(time(0));
    
    cout.precision(10);
    //Read the problem instance to optimize.
	PROBLEM = GetProblemInfo(PROBLEM_TYPE,INSTANCE_FILENAME);
    // MAX_EVALUATIONS=PROBLEM_SIZE*PROBLEM_SIZE*1000; //1600 000

    
    //Create the file to store the results.
    ofstream output_file;
    output_file.open(RESULTS_FILENAME);
    output_file.close();

    //Initialize the algorithm
    //cout<<"Ranking EDA..."<<endl;
    int infinite=0;
    if (infinite){
        InfiniteRankingEDA * inf=new InfiniteRankingEDA(PROBLEM,PROBLEM_SIZE, MAX_EVALUATIONS, MODEL_TYPE, METRIC_TYPE,INVERSE);
        inf->Run();
        delete inf;
        exit(1);
    }
    else{
        RankingEDA * alg= new RankingEDA(PROBLEM,PROBLEM_SIZE, POP_SIZE, MAX_EVALUATIONS, B_RATIO, CUT_POINT_COUNT, PREVIOUS_SAMPLED_REFERENCE_COUNT, NUMBER_OF_EDGE, MODEL_TYPE, METRIC_TYPE,INVERSE, SEED, RESULTS_FILENAME);
        //alg->MetricSuitability_Experiment(RESULTS_FILENAME);exit(1);
        alg->Run(PROBLEM_TYPE);
        
    
    

    //Get best solution fitness and the number of evaluations performed.
    
    // evaluation times here is evaluation times when converge
    int evaluations = alg->GetPerformedEvaluations();
    int convergeNFE = alg->GetRealConvergeNFE();
    
    CIndividual * best = alg->GetBestSolution();
    int dist=0;
    if (IsAPI(PROBLEM_TYPE)==1)
        dist=((API*)PROBLEM)->CalculateDistanceToOptimum(best->Genes());
    
    gettimeofday(&tim, NULL);
    double t2=tim.tv_sec+(tim.tv_usec/1000000.0);

    //Print results
    if (IsAPI(PROBLEM_TYPE)==1)
            WriteResults2(best->Value(),best->Genes(),evaluations,t2-t1,dist);
        else
            WriteResults1(best->Value(),best->Genes(),convergeNFE,evaluations,t2-t1);
    
        delete alg;
    }

    return 0;
}

