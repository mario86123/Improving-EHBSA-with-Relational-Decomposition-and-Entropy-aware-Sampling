//
//  RankingModel.h
//  RankingEDAsCEC
//
//  Created by Josu Ceberio Uribe on 11/19/13.
//  Copyright (c) 2013 Josu Ceberio Uribe. All rights reserved.
//

#ifndef __RankingEDA__RankingModel__
#define __RankingEDA__RankingModel__

#include <iostream>
#include "Population.h"
#include "InfinitePopulation.h"
#include "Individual.h"


class CRankingModel
{
public:
    
	/*
	 * The constructor.
	 */
	CRankingModel();
	
	/*
	 * The destructor.
	 */
	virtual ~CRankingModel();
	
	/*
	 * Virtual learning function.
	 */
	virtual bool Learn(CPopulation * population, int size);

	virtual bool Learn(CPopulation * population, int start_idx, int end_idx);

	virtual bool Learn(CPopulation * population, int mode, int * used_distance_count_arr);


    /*
     * Virtual learning function.
     */
    virtual bool Learn(CInfinitePopulation * population, float * selection_probs);
    
    /*
	 * Virtual learning function.
	 */
	virtual bool Learn(int * consensus_ranking, double theta_parameter);


	virtual bool Learn(int * new_chromosome, int * templateee);

    
    /*
	 * Virtual sampling function.
	 */
	virtual void Sample(int * genes);

	virtual void Sample(int * genes, int *templateee);

	virtual void Sample(int * genes, int *templateee, int **used_binary_relation);
	
	virtual void Sample(int * genes, int *templateee, int need_to_sample_genes_count);

	virtual void Sample(int * genes, int *templateee, int need_to_sample_genes_count, double chi_std);

	virtual void Sample(int * genes, int number_of_edge);

	/*
	 * Calculates the probability of the solution given the probabilistic model.
	 */
	virtual float Probability(int * solution);
    
private:
};

#endif /* defined(__RankingEDA__RankingModel__) */
