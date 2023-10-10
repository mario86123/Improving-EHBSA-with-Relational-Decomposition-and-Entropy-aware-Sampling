// MST2ME

#ifndef __RankingEDA__MST2ME__
#define __RankingEDA__MST2ME__

#include <iostream>

#include "Individual.h"
#include "RankingModel.h"
#include "Population.h"
#include <list>
#include <vector>


class CMST2ME : public CRankingModel
{
public:

	/*
	 * The constructor.
	 */
	CMST2ME(int problem_size, int sel_size, double b_ratio, string result_file_name);
	
	/*
	 * The destructor. It frees the memory allocated at the construction of the Cayley model.
	 */
	virtual ~CMST2ME();

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

    struct node {
        int node_num;
        double entropy;
        double sum;
    };
    node * m_mst_node_arr;


    struct edge {
        int node_a;
        int node_b;
        double length;
        double entropy;
        double sum;
        int distance;
    };

    edge * m_mst_edge_arr;

    // double * m_max_entropy_lst;
    int m_max_num_of_edge;
    // string m_result_file_name;
    ofstream m_entropy_file;

    int * useless_counter;
    
};


#endif /* defined(__RankingEDA__MST2ME__) */