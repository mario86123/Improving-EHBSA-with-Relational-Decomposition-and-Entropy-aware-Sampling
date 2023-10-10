
#ifndef __RankingEDA__NNHBSAWT__
#define __RankingEDA__NNHBSAWT__

#include <iostream>

#include "Individual.h"
#include "RankingModel.h"
#include "Population.h"
#include <list>
#include <vector>


class CNNHBSAWT : public CRankingModel
{
public:

	/*
	 * The constructor.
	 */
	CNNHBSAWT(int problem_size, int sel_size, double b_ratio, int cut_point_count, CPopulation * population);
	
	/*
	 * The destructor. It frees the memory allocated at the construction of the Cayley model.
	 */
	virtual ~CNNHBSAWT();

    /*
     * Given a population of samples, it learns a Mallows model from the data.
     */
    bool Learn(int * new_chromosome, int * templateee);

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

    double ** m_NHM;
    double m_epsilon;

    int m_cut_point_count;
};

#endif /* defined(__RankingEDA__NNHBSAWT__) */