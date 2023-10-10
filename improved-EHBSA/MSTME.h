
#ifndef __RankingEDA__MSTME__
#define __RankingEDA__MSTME__

#include <iostream>

#include "Individual.h"
#include "RankingModel.h"
#include "Population.h"
#include <list>
#include <vector>


class CMSTME : public CRankingModel
{
public:

	/*
	 * The constructor.
	 */
	CMSTME(int problem_size, int sel_size, double b_ratio);
	
	/*
	 * The destructor. It frees the memory allocated at the construction of the Cayley model.
	 */
	virtual ~CMSTME();

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

    double ** m_node_matrix;  
    double *** m_edge_matrix;
    
    double m_epsilon;

    struct edge {
        int node_a;
        int node_b;
        double length;
    };

    edge * m_mst_edge_arr;

    struct element {
        int parent;
        int rank;
    };

    int find_set(int x, element * arr);
    void link(int x, int y, element * arr);
    void union_set(int x, int y, element * arr);

    
};

#endif /* defined(__RankingEDA__MSTModel__) */