
#ifndef __RankingEDA__NHBSAWO__
#define __RankingEDA__NHBSAWO__

#include <iostream>

#include "Individual.h"
#include "RankingModel.h"
#include "Population.h"
#include <list>
#include <vector>


class CNHBSAWO : public CRankingModel
{
public:

	/*
	 * The constructor.
	 */
	CNHBSAWO(int problem_size, int sel_size, double b_ratio);
	
	/*
	 * The destructor. It frees the memory allocated at the construction of the Cayley model.
	 */
	virtual ~CNHBSAWO();

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
};

#endif /* defined(__RankingEDA__NHBSAWO__) */