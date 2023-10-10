//
//  RankingModel.cpp
//  RankingEDAsCEC
//
//  Created by Josu Ceberio Uribe on 11/19/13.
//  Copyright (c) 2013 Josu Ceberio Uribe. All rights reserved.
//

#include "RankingModel.h"
#include "Population.h"


/*
 * The constructor.
 */
CRankingModel::CRankingModel()
{
	
}

/*
 * The destructor.
 */
CRankingModel::~CRankingModel()
{

}

/*
 * Virtual learning function.
 */
bool CRankingModel::Learn(CPopulation * population, int size)
{
    return true;
}

bool CRankingModel::Learn(CPopulation * population, int start_idx, int end_idx)
{
    return true;
}

bool CRankingModel::Learn(CPopulation * population, int mode, int * used_distance_count_arr)
{
    return true;
}

/*
 * Virtual learning function.
 */
bool CRankingModel::Learn(CInfinitePopulation * population, float * selection_probs)
{
    return true;
}


/*
 * Virtual learning function.
 */
bool CRankingModel::Learn(int * consensus_ranking, double theta_parameter)
{
    return true;
}

bool CRankingModel::Learn(int * new_chromosome, int * templateee)
{
    return true;
}


/*
 * Virtual sampling function.
 */
void CRankingModel::Sample(int * genes)
{
    
}

void CRankingModel::Sample(int * genes, int *templateee)
{
    
}

void CRankingModel::Sample(int * genes, int *templateee, int **used_binary_relation)
{
    
}

void CRankingModel::Sample(int * genes, int *templateee, int need_to_sample_genes_count)
{
    
}

void CRankingModel::Sample(int * genes, int *templateee, int need_to_sample_genes_count, double chi_std)
{
    
}

void CRankingModel::Sample(int * genes, int number_of_edge)
{
    
}
/*
 * Calculates the probability of the solution given the probabilistic model.
 */
float CRankingModel::Probability(int * solution)
{
    return 0;
}