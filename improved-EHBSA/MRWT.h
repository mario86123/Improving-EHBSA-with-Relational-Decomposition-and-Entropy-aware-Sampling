
#ifndef __RankingEDA__MRWT__
#define __RankingEDA__MRWT__

#include <iostream>

#include "Individual.h"
#include "RankingModel.h"
#include "Population.h"
#include <list>
#include <vector>


class CMRWT : public CRankingModel
{
public:

	/*
	 * The constructor.
	 */
	CMRWT(int problem_size, int sel_size, double b_ratio, double cut_point_count);
	
	/*
	 * The destructor. It frees the memory allocated at the construction of the Cayley model.
	 */
	virtual ~CMRWT();

    /*
     * Given a population of samples, it learns a Mallows model from the data.
     */
    bool Learn(CPopulation * population, int size);

    /*
     * Given the consensus ranking, it samples a new individual given the s
     */
    void Sample(int * permutation, int *templateee);
	
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

    // number of node that need to sampled
    int m_sampled_node_count;

    double *** m_multiple_relation_matrix;
    double m_epsilon;

    int first_node_idx;
    int * first_node_idx_arr;
};

#endif /* defined(__RankingEDA__EHBSAWO__) */