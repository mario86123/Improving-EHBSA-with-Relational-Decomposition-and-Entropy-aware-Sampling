
#ifndef __RankingEDA__BWModel__
#define __RankingEDA__BWModel__

#include <iostream>

#include "Individual.h"
#include "RankingModel.h"
#include "Population.h"
#include <list>
#include <vector>


class CBWModel : public CRankingModel
{
public:

	/*
	 * The constructor.
	 */
	CBWModel(int problem_size, int sel_size, double b_ratio);
	
	/*
	 * The destructor. It frees the memory allocated at the construction of the Cayley model.
	 */
	virtual ~CBWModel();

    /*
     * Given a population of samples, it learns a Mallows model from the data.
     */
    bool Learn(CPopulation * population, int size);

    /*
     * Given the consensus ranking, it samples a new individual given the s
     */
    void Sample(int * permutation);
	
private:
    
    /*
     * Problem size.
     */
    int m_problem_size;
    int m_half_relation_distance;

    double m_epsilon;

    /*
     * Sample size.
     */
    int m_sample_size;
    
    /*
     * Sample of solutions.
     */
    int ** m_samples;

    double ** m_half_relation_matrix;  
	
};

#endif /* defined(__RankingEDA__BWModel__) */