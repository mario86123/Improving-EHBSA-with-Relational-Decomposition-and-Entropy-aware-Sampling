
#ifndef __RankingEDA__FERBOP__
#define __RankingEDA__FERBOP__

#include <iostream>

#include "Individual.h"
#include "RankingModel.h"
#include "Population.h"
#include <list>
#include <vector>
#include "Tools.h"




class CFERBOP : public CRankingModel
{
public:

	/*
	 * The constructor.
	 */
	CFERBOP(int problem_size, int sel_size, double b_ratio, int cut_point_count);
	
	/*
	 * The destructor. It frees the memory allocated at the construction of the Cayley model.
	 */
	virtual ~CFERBOP();

    /*
     * Given a population of samples, it learns a Mallows model from the data.
     */
    bool Learn(CPopulation * population, int size);

    /*
     * Given the consensus ranking, it samples a new individual given the s
     */
    void Sample(int * permutation, int * templateee);

    int Build_network(int * templateee, int already_sampled_idx_count, edgeee * edge_lst, int * sampled_or_not_idx, bool * is_sampled_node);

    int update_and_find_lowest_entropy_edge(int sample_idx, int sample_node_num, int total_edge_count, edgeee * edge_lst);

    void construct_probability_array(double * tmp_probability_array, bool * is_sampled_node, int lowest_entropy_edge_num, edgeee * edge_lst);

    double fast_calculate_entropy(double * arr);



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
    double ** m_multiple_relation_entropys;
    
    
    double m_epsilon;
    int m_cut_point_count;

    double * m_plogp; 

};

#endif /* defined(__RankingEDA__FERBOP__) */