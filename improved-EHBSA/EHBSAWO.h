
#ifndef __RankingEDA__EHBSAWO__
#define __RankingEDA__EHBSAWO__

#include <iostream>

#include "Individual.h"
#include "RankingModel.h"
#include "Population.h"
#include <list>
#include <vector>


class CEHBSAWO : public CRankingModel
{
public:

	/*
	 * The constructor.
	 */
	CEHBSAWO(int problem_size, int sel_size, double b_ratio);
	
	/*
	 * The destructor. It frees the memory allocated at the construction of the Cayley model.
	 */
	virtual ~CEHBSAWO();

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

    double ** m_EHM;
    double m_epsilon;
};

#endif /* defined(__RankingEDA__EHBSAWO__) */