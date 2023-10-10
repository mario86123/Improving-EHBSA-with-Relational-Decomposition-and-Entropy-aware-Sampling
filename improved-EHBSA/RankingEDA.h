/*
 *  RankingEDA.h
 *  RankingEDAsCEC
 *
 *  Created by Josu Ceberio Uribe on 7/11/13.
 *  Copyright 2013 University of the Basque Country. All rights reserved.
 *
 */

#ifndef _RANKINGEDA_H__
#define _RANKINGEDA_H__

#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string.h>
#include <stdio.h>
#include "PBP.h"
#include "Tools.h"
#include "Population.h"
#include "RankingModel.h"

using std::istream;
using std::ostream;
using namespace std;
using std::cerr;
using std::cout;
using std::endl;
using std::ifstream;
using std::stringstream;
using std::string;

class RankingEDA
{
	
public:
	
    /*
     * Problem
     */
    PBP * m_problem;
    
    /*
     * Problem size
     */
    int m_problem_size;
    
    /*
     * Maximum number of evaluations permitted in the execution.
     */
    long int m_max_evaluations;
    
    /*
     * The number of evaluations performed at the moment.
     */
    int m_evaluations;

    /*
     * The number of the evaluations performed to achieve the best solution so far.
     */
    int m_convergence_evaluations;
    int m_real_converge_evaluations;
    
    /*
     * The best solution found so far.
     */
    CIndividual * m_best;
    
    /*
     * The size of the population.
     */
    int m_pop_size;
    
    /*
     * The size of the selection pool.
     */
    int m_sel_size;
    
    /*
     * The size of the offspring population.
     */
    int m_offspring_size;

    //The type of the metric to use in the model of the EDA.
    char m_metric_type[10];

    //The type of the model used in the EDA.
    char m_model_type[10];
    
    /*
     * Evaluates the inverse of the samples solutions.
     */
    int m_inverse;
    int m_seed;
    int m_number_of_edge;
    int m_cut_point_count;
    double m_b_ratio;

    string m_result_file_name;
    
    /*
     * The population
     */
    CPopulation * m_population;
    
    /*
     * The probalistic model.
     */
    CRankingModel * m_model;
    CRankingModel * m_model_1;
    CRankingModel * m_model_2;
    
    /*
     * The constructor.
     */
	RankingEDA(PBP * problem, int problem_size, int population_size, long int max_evaluations, double b_ratio, double cut_point_count, int previous_sampled_reference_count, int number_of_edge, char * model_type, char * metric_type, int inverse, int seed, char *mresult_file_name);
	
    /*
     * The destructor.
     */
    virtual ~RankingEDA();
	
    /*
     * Running function
     */
	int Run(char * problem_type);
    
    /*
     * Returns the number of performed evaluations.
     */
	int GetPerformedEvaluations();
    int GetRealConvergeNFE();
    
    /*
     * Returns the fitness of the best solution obtained.
     */
    long int GetBestSolutionFitness();
    
    /*
     * Returns the best solution obtained.
     */
    CIndividual * GetBestSolution();
    
    /*
     * Experiment for the Metric suitability for the different problems.
     */
    void MetricSuitability_Experiment(char * results_file);
private:
    
    /*
     * This method applies a swap of the given i,j positions in the array.
     */
    void Swap(int * array, int i, int j);
	
};
#endif
