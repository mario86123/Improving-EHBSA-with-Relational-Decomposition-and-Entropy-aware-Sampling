
#ifndef __RankingEDA__MRModel__
#define __RankingEDA__MRModel__

#include <iostream>

#include "Individual.h"
#include "RankingModel.h"
#include "Population.h"
#include <list>
#include <vector>


class CMRModel : public CRankingModel
{
public:

	/*
	 * The constructor.
	 */
	CMRModel(int problem_size, int sel_size, double b_ratio);
	
	/*
	 * The destructor. It frees the memory allocated at the construction of the Cayley model.
	 */
	virtual ~CMRModel();

    /*
     * Given a population of samples, it learns a Mallows model from the data.
     */
    bool Learn(CPopulation * population, int size);

    /*
     * Given the consensus ranking, it samples a new individual given the s
     */
    void Sample(int * permutation);

    // int sample_from_array(int *arr);
	
private:
    
    /*
     * Problem size.
     */
    int m_problem_size;

    /*
     * Sample size.
     */
    int m_sample_size;
    
    /*
     * Sample of solutions.
     */
    int ** m_samples;

    double *** m_multiple_relation_matrix;  
    
    int first_node_idx;
    double * first_node_idx_arr;
    
    double m_epsilon;
};

#endif /* defined(__RankingEDA__MRModel__) */