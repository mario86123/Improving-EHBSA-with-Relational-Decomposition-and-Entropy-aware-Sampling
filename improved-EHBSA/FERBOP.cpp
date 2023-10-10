// Fast entropy RBOP
#include "FERBOP.h"
#include "Tools.h"
#include <errno.h>
#include <fcntl.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <list>
#include <vector>
using std::cerr;
using std::cout;
using std::endl;

/*
 * Class constructor.
 */
CFERBOP::CFERBOP(int problem_size, int sel_size, double b_ratio, int cut_point_count)
{
	m_problem_size=problem_size;
    m_sample_size=sel_size;
    m_samples= new int*[m_sample_size];

    m_multiple_relation_matrix = new double**[m_problem_size];
    for (int i = 0; i < m_problem_size; ++i) {
        m_multiple_relation_matrix[i] = new double*[m_problem_size];
        // m_multiple_relation_matrix[i][1][j] = node i is before node j for 1 position
        for (int j = 0; j < m_problem_size; ++j) {
            m_multiple_relation_matrix[i][j] = new double[m_problem_size]();
        }
    }

    // m_multiple_relation_entropys
    m_multiple_relation_entropys = new double*[m_problem_size];
    for (int j = 0; j < m_problem_size; ++j) {
        m_multiple_relation_entropys[j] = new double[m_problem_size]();
    }



    // NHBSA like
    m_epsilon = b_ratio * m_sample_size / (m_problem_size);
    cout << "m_epsilon: " << m_epsilon << endl;

    m_cut_point_count = cut_point_count;



    // (m_plogp[0] = 0) to (m_plogp[m_sample_size] = 0)
    m_plogp = new double[m_sample_size + 1]();
    // m_plogp[0] = 0;

    for (int i = 0; i < m_sample_size + 1; i ++) {
        cout << "total: " << ((m_problem_size - 1) * m_epsilon + m_sample_size) << endl;
        double p = ( (i + m_epsilon) / ((m_problem_size - 1) * m_epsilon + m_sample_size) );
        cout << "p[" << i << "]: " << p << endl;
        m_plogp[i] = p * log(p); // nature log
    }

}

CFERBOP::~CFERBOP()
{
    // PrintMRMatrix(m_multiple_relation_matrix, m_problem_size);
    // cout << "m_epsilon: " << m_epsilon << endl;

    delete [] m_samples;

    // Release matrix Memory
    for(int i = 0; i < m_problem_size; ++i){
        for(int j = 0; j < m_problem_size; ++j){
            delete[] m_multiple_relation_matrix[i][j];
        }
        delete[] m_multiple_relation_matrix[i];
    }
    delete[] m_multiple_relation_matrix;



    for(int i = 0; i < m_problem_size; ++i) {
        delete[] m_multiple_relation_entropys[i];
    }
    delete[] m_multiple_relation_entropys;
    
}

/*
 * Virtual learning function.
 */
bool CFERBOP::Learn(CPopulation * population, int size)
{   
    // population->Print();


    // clean matrix
    for(int ref_node = 0; ref_node < m_problem_size; ++ref_node){
        for(int dist = 0; dist < m_problem_size; ++dist){
            for(int sample_node = 0; sample_node < m_problem_size; ++sample_node){
                
                if (ref_node != sample_node && dist != 0) {
                    m_multiple_relation_matrix[ref_node][dist][sample_node] = m_epsilon;
                }
                else {
                    m_multiple_relation_matrix[ref_node][dist][sample_node] = 0.0;
                }
            }
        }
    }
    
    // selection ?????
    // selection ?????
    // selection ?????
    // printf("sample:\n");
    for(int k=0;k<size;k++) {
		m_samples[k] = population->m_individuals[k]->Genes();
        

        // build matrix
        // assign value to matrix
        for (int reference_node_idx = 0; reference_node_idx < m_problem_size; ++reference_node_idx) {
            for (int distance = 0; distance < m_problem_size; ++distance) { // next node
                int next_node_idx = reference_node_idx + distance;
                if (next_node_idx >= m_problem_size) {
                    next_node_idx -= m_problem_size;
                }
                m_multiple_relation_matrix[ m_samples[k][reference_node_idx] ][distance][ m_samples[k][next_node_idx] ] ++;
            }
        }
    }

    // PrintMRMatrix(m_multiple_relation_matrix, m_problem_size);


    // calculate m_multiple_relation_matrix entropys
    for(int ref_node = 0; ref_node < m_problem_size; ++ref_node){
        for(int dist = 0; dist < m_problem_size; ++dist){
            m_multiple_relation_entropys[ref_node][dist] = fast_calculate_entropy(m_multiple_relation_matrix[ref_node][dist]);
            // cout << "m_multiple_relation_entropys["<< ref_node << "][" << dist << "]: " << m_multiple_relation_entropys[ref_node][dist] << endl;
        }
    }


    return true;
}

/*
 * Virtual sampling function.
 */
void CFERBOP::Sample(int * genes, int * templateee)
{

    int already_sampled_idx_count = m_problem_size - (m_problem_size / m_cut_point_count);
    int total_edge_count = already_sampled_idx_count * (m_problem_size - already_sampled_idx_count);



    // random choose (1/5 * ell) indexs as not sampled indexs and (4/5 * ell) indexs as already sampled indexs
    int sampled_or_not_idx[m_problem_size];
	for (int i = 0; i < m_problem_size; ++i) {
		sampled_or_not_idx[i] = i;
	}
	random_shuffle(sampled_or_not_idx, sampled_or_not_idx + m_problem_size);
    
    
    
    bool is_sampled_idx[m_problem_size];
    bool is_sampled_node[m_problem_size];


    // the first already_sampled_idx_count sampled_or_not_idxs is the indexes that already sampled

    for (int i = 0; i < m_problem_size; ++i) {

        if (i < already_sampled_idx_count) {             // already sampled
		    is_sampled_idx[sampled_or_not_idx[i]] = true;
		    is_sampled_node[templateee[sampled_or_not_idx[i]]] = true;
            genes[sampled_or_not_idx[i]] = templateee[sampled_or_not_idx[i]];
        }
        else {                                           // not sampled yet
		    is_sampled_idx[sampled_or_not_idx[i]] = false;
            is_sampled_node[templateee[sampled_or_not_idx[i]]] = false;
        }
    }




    // record and calculate [ (1/5 * ell) nodes not sampled node x (4/5 * ell) nodes ] edges entropy 
        // record already sampled nodes position
        // calculate entropy
        // find a way to record which edges can be used and which cannot

    edgeee edge_lst[total_edge_count];



    //             create network edges and calculate entropy
    // time complexity:     O(n^2)       *       O(n)         =   O(n^3) 
    int all_edges_lowest_entropy_num = Build_network(templateee, already_sampled_idx_count, edge_lst, sampled_or_not_idx, is_sampled_node);




    // === sample === /
    // ---------- while still have genes not sampled yet ---------- //
    for (int kkk = 0; kkk < m_problem_size - already_sampled_idx_count - 1; kkk++) {


        int lowest_entropy_edge_num = all_edges_lowest_entropy_num;


        // === sample
        double tmp_probability_array[m_problem_size];

        for (int i = 0; i < m_problem_size; ++i) {
            
            if (is_sampled_node[i] == false) {
                tmp_probability_array[i] = m_multiple_relation_matrix[edge_lst[lowest_entropy_edge_num].node_a][edge_lst[lowest_entropy_edge_num].distance][i];
            }

            // else if that node num is already sampled
            else {
                tmp_probability_array[i] = 0;
            }
        }

        int sample_idx = edge_lst[lowest_entropy_edge_num].node_a_index + edge_lst[lowest_entropy_edge_num].distance;
        if (sample_idx >= m_problem_size) {
            sample_idx -= m_problem_size;
        }

        genes[sample_idx] = sample_from_array(tmp_probability_array, m_problem_size);

        is_sampled_idx[sample_idx] = true;
        is_sampled_node[genes[sample_idx]] = true;

        // === end of sample

        all_edges_lowest_entropy_num = update_and_find_lowest_entropy_edge(sample_idx, genes[sample_idx], total_edge_count, edge_lst);
    }
    // ----------  end of while still have genes not sampled yet ---------- //
    
    int sample_idx = 0;
    for (; sample_idx < m_problem_size; sample_idx ++) {
        if (is_sampled_idx[sample_idx] == false) break;
    }
    int sample_node = 0;
    for (; sample_node < m_problem_size; sample_node ++) {
        if (is_sampled_node[sample_node] == false) break;
    }
    genes[sample_idx] = sample_node;

    cout  << endl;


}







int CFERBOP::Build_network(int * templateee, int already_sampled_idx_count, edgeee * edge_lst, int * sampled_or_not_idx, bool * is_sampled_node) {
    int num = 0; // current edge count

    double lowest_entropy = 100000000.0;
    int lowest_entropy_edge_num = -1;


    // the first already_sampled_idx_count sampled_or_not_idxs is the indexes that already sampled

    for (int i = 0; i < already_sampled_idx_count; ++i) {                             // already sampled node
        for (int j = 0; j < (m_problem_size - already_sampled_idx_count); ++j) {      // not sampled yet node

            edge_lst[num].node_a = templateee[sampled_or_not_idx[i]];
            edge_lst[num].node_a_index = sampled_or_not_idx[i];


            edge_lst[num].distance = sampled_or_not_idx[already_sampled_idx_count + j] - sampled_or_not_idx[i];
            if (edge_lst[num].distance < 0) {
                edge_lst[num].distance += m_problem_size;
            }

            
            edge_lst[num].element_sum = 0;
            double deleted_element_sum = 0;

            
            for (int k = 0; k < m_problem_size; ++k) {
                // if that node num is not sampled yet
                if (is_sampled_node[k] == false) {
                    edge_lst[num].element_sum += m_multiple_relation_matrix[edge_lst[num].node_a][edge_lst[num].distance][k];
                }
                // else if that node num is already sampled
                else {
                    deleted_element_sum += m_multiple_relation_matrix[edge_lst[num].node_a][edge_lst[num].distance][k];
                }
            }

            cout << "edge_lst[" << num << "].element_sum: " << edge_lst[num].element_sum << endl;
            cout << "deleted_element_sum: " << deleted_element_sum << endl;
            // cout << "m_multiple_relation_entropys[edge_lst[" << num << "].node_a][edge_lst[" << num << "].distance]: " << m_multiple_relation_entropys[edge_lst[num].node_a][edge_lst[num].distance] << endl;
            // update entropy here !!
            // --- old_sum = popsize + m_epsilon * ... ;

            PrintArray(m_multiple_relation_matrix[edge_lst[num].node_a][edge_lst[num].distance], m_problem_size);
            // cout << "sum: " << sum_arr(m_multiple_relation_matrix[edge_lst[num].node_a][edge_lst[num].distance], m_problem_size) << endl;

            edge_lst[num].entropy = update_entropy(m_multiple_relation_entropys[edge_lst[num].node_a][edge_lst[num].distance], deleted_element_sum, ((m_problem_size - 1) * m_epsilon + m_sample_size));
            cout << "edge_lst[" << num << "].entropy: " << edge_lst[num].entropy << endl << endl;

            if (edge_lst[num].entropy < lowest_entropy) {
                lowest_entropy = edge_lst[num].entropy;
                lowest_entropy_edge_num = num;
            }
            


            edge_lst[num].can_be_used = true;

            num ++;
        }
    }

    cout << "lowest_entropy: " << lowest_entropy << endl;
    cout << "lowest_entropy_edge_num: " << lowest_entropy_edge_num << endl;
    
    
    return lowest_entropy_edge_num;
}




int CFERBOP::update_and_find_lowest_entropy_edge(int sample_idx, int sample_node_num, int total_edge_count, edgeee * edge_lst) {

    double lowest_entropy = 100000000.0;
    int all_edges_lowest_entropy_num = -1;


    // === update all edges entropy and edges information ===//
    for (int edge_num = 0; edge_num < total_edge_count; edge_num++) {

        if (edge_lst[edge_num].can_be_used == true) {

            int tmp_sample_idx = edge_lst[edge_num].node_a_index + edge_lst[edge_num].distance;
            if (tmp_sample_idx >= m_problem_size) {
                tmp_sample_idx -= m_problem_size;
            }
            
            // --- delete same distance edges --- //
            if (tmp_sample_idx == sample_idx) {
                edge_lst[edge_num].can_be_used = false;
            }
            // ---  end of delete same distance edges --- //

            // --- update different distance edges information --- //
            else {

                double deleted_element = m_multiple_relation_matrix[edge_lst[edge_num].node_a][edge_lst[edge_num].distance][sample_node_num];

                edge_lst[edge_num].entropy = update_entropy(edge_lst[edge_num].entropy, deleted_element, edge_lst[edge_num].element_sum);
                edge_lst[edge_num].element_sum -= deleted_element;



                if (edge_lst[edge_num].entropy < lowest_entropy) {
                    lowest_entropy = edge_lst[edge_num].entropy;
                    all_edges_lowest_entropy_num = edge_num;

                }
            }
            // --- end of update different distance edges information --- //

        }
        // --- end of if this edge can be used --- //
    }
    // === end of update all edges entropy and edges information === //


    return all_edges_lowest_entropy_num;
}



void CFERBOP::construct_probability_array(double * tmp_probability_array, bool * is_sampled_node, int lowest_entropy_edge_num, edgeee * edge_lst) {

    for (int i = 0; i < m_problem_size; ++i) {

        // if that node num is not sampled yet
        if (is_sampled_node[i] == false) {
            tmp_probability_array[i] = m_multiple_relation_matrix[edge_lst[lowest_entropy_edge_num].node_a][edge_lst[lowest_entropy_edge_num].distance][i];
        }

        // else if that node num is already sampled
        else {
            tmp_probability_array[i] = 0;
        }
    }
}



double CFERBOP::fast_calculate_entropy(double * arr) {
    
    double ans = 0;
    
    for (int i = 0; i < m_problem_size; ++i) {
    
        if (arr[i] != 0.0) {
            ans += m_plogp[ int(arr[i]) ];
            // cout << "m_plogp[ " << int(arr[i]) << " ]: " << m_plogp[ int(arr[i]) ] << endl;
        }
    }
    // cout << endl;

    return -ans;
}
