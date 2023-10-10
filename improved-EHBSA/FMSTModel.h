
#ifndef __RankingEDA__FMSTModel__
#define __RankingEDA__FMSTModel__

#include <iostream>

#include "Individual.h"
#include "RankingModel.h"
#include "Population.h"
#include <list>
#include <vector>


class CFMSTModel : public CRankingModel
{
public:

	/*
	 * The constructor.
	 */
	CFMSTModel(int problem_size, int sel_size, double b_ratio);
	
	/*
	 * The destructor. It frees the memory allocated at the construction of the Cayley model.
	 */
	virtual ~CFMSTModel();

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

    double ** m_node_matrix;  
    double *** m_edge_matrix;
    
    double m_epsilon;

    // struct edge {
    //     int node_a;
    //     int node_b;
    //     double length;
    //     int already_sampled_node_count;
    // };

    // // edge * m_mst_edge_arr;

    // struct node {
    //     int node_num;
    //     double length;
    // };

    // node * m_mst_node_arr;


    double * node_entropy_arr;
    double ** edge_entropy_arr;

    int * m_sample_order;

};

#endif /* defined(__RankingEDA__FMSTModel__) */