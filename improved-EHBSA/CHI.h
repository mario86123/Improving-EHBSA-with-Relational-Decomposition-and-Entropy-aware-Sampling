// chi-square ACNRBOP

#ifndef __RankingEDA__CHI__
#define __RankingEDA__CHI__

#include <iostream>

#include "Individual.h"
#include "RankingModel.h"
#include "Population.h"
#include <list>
#include <vector>
#include "Tools.h"




class CCHI : public CRankingModel
{
public:

	/*
	 * The constructor.
	 */
	CCHI(int problem_size, int sel_size, double b_ratio);
	
	/*
	 * The destructor. It frees the memory allocated at the construction of the Cayley model.
	 */
	virtual ~CCHI();

    /*
     * Given a population of samples, it learns a Mallows model from the data.
     */
    bool Learn(CPopulation * population, int mode);

    void Sample(int * permutation, int * templateee, int cut_point_count);

    edgeee add_edgesss(int sampled_node_num, int sampled_node_idx, list<edgeee>&  edge_lst, bool * is_sampled_node, bool * is_sampled_idx);


    // edgeee Build_network(int * templateee, int already_sampled_idx_count, edgeee * edge_lst, int * sampled_or_not_idx, bool * is_sampled_node);
    edgeee Build_network(int * templateee, int already_sampled_idx_count, list<edgeee>&  edge_lst, int * sampled_or_not_idx, bool * is_sampled_node);

    double fast_calculate_entropy(int * arr, int summation);

    // edgeee update_and_find_lowest_entropy_edge(int sample_idx, int sample_node_num, int total_edge_count, edgeee * edge_lst);
    edgeee update_and_find_lowest_entropy_edge(int sample_idx, int sample_node_num, list<edgeee>&  edge_lst);


    // void construct_probability_array(double * tmp_probability_array, bool * is_sampled_node, int lowest_entropy_edge_num, edgeee * edge_lst);
	void construct_probability_array(double * tmp_probability_array, bool * is_sampled_node, edgeee lowest_entropy_edge);
    void construct_node_probability_array(double * tmp_probability_array, bool * is_sampled_node, int idx);

    

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

    int *** m_multiple_relation_matrix;  
    int ** m_node_matrix;  
    
    double ** m_plogp;  
    
    int * m_used_distance_count;
    int m_used_node_count;
    int m_need_to_sample_genes_count;

    double m_adj_ratio;
    
    double m_epsilon;


    // double ** m_entropy_matrix;
};

#endif /* defined(__RankingEDA__CHI__) */