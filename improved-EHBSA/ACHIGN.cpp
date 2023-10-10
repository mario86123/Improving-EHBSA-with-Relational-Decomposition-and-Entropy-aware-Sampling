// Adaptive chi-square and gene number NRBOP

#include "ACHIGN.h"
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
CACHIGN::CACHIGN(int problem_size, int sel_size, double b_ratio)
{
	m_problem_size=problem_size;
    m_sample_size=sel_size;
    m_samples= new int*[m_sample_size];

    m_multiple_relation_matrix = new int**[m_problem_size];
    for (int i = 0; i < m_problem_size; ++i) {
        m_multiple_relation_matrix[i] = new int*[m_problem_size];
        // m_multiple_relation_matrix[i][1][j] = node i is before node j for 1 position
        for (int j = 0; j < m_problem_size; ++j) {
            m_multiple_relation_matrix[i][j] = new int[m_problem_size]();
        }
    }


    // NHM
    m_node_matrix = new int*[m_problem_size];
    for (int i = 0; i < m_problem_size; ++i) {
        m_node_matrix[i] = new int[m_problem_size]();
    }



    // NHBSA like
    m_epsilon = b_ratio * m_sample_size / (m_problem_size);
    // cout << "m_epsilon: " << m_epsilon << endl;



    // === m_plogp === //
    // initaialize the memory of m_plog_p
    m_plogp = new double*[m_sample_size + 1];

    for (int mother = 0; mother < m_sample_size + 1; mother ++) {
        m_plogp[mother] = new double[mother + 1]();

        for (int son = 0; son < mother + 1; son++) {
             double p = double (son) / (mother);
            m_plogp[mother][son] = p * log(p);
        }
    }
    // === end of m_plogp === //

    m_used_distance_count = new int[m_problem_size]();
    
    m_adj_ratio = 0.0;

}

CACHIGN::~CACHIGN()
{
    // PrintMRMatrix(m_multiple_relation_matrix, m_problem_size);

    delete [] m_samples;

    // Release matrix Memory
    for(int i = 0; i < m_problem_size; ++i){
        for(int j = 0; j < m_problem_size; ++j){
            delete[] m_multiple_relation_matrix[i][j];
        }
        delete[] m_multiple_relation_matrix[i];
    }
    delete[] m_multiple_relation_matrix;

    delete [] m_used_distance_count;


    // delete the memory of m_plog_p
    for (int mother = 0; mother < m_sample_size + 1; ++mother) {
        delete[] m_plogp[mother];
    }
    delete[] m_plogp;
}

/*
 * Virtual learning function.
 */
// bool CACHIGN::Learn(CPopulation * population, int size)
//
// mode == 1 -> upper half population
// mode == 2 -> lower half population
// mode == 3 -> whole population
bool CACHIGN::Learn(CPopulation * population, int mode)
{   
    // population->Print();
    // cout << endl;

    int sum_of_used_distance_count = sum_arr(m_used_distance_count, m_problem_size);
    // printf("\nsum_of_used_distance_count: %d\n", sum_of_used_distance_count);

    if (sum_of_used_distance_count == 0) {
        m_adj_ratio = 0.0;
    }
    else {
        // m_adj_ratio = double(m_used_distance_count[1] + m_used_distance_count[(m_problem_size-1)]) / sum_of_used_distance_count;
        m_adj_ratio = double(m_used_distance_count[1]) / sum_of_used_distance_count;
    }
    // cout << "m_adj_ratio: " << m_adj_ratio << endl;


    // PrintArray(m_used_distance_count, m_problem_size);
    // cout << "m_used_node_ratio: " << double(m_used_node_count) / m_need_to_sample_genes_count << endl<< endl;

    for(int dis = 0; dis < m_problem_size; ++dis) {
        m_used_distance_count[dis] = 0;
    }
    m_used_node_count = 0;
    m_need_to_sample_genes_count = 0;

    // clean matrix
    for(int ref_node = 0; ref_node < m_problem_size; ++ref_node){
        for(int dist = 0; dist < m_problem_size; ++dist){
            for(int sample_node = 0; sample_node < m_problem_size; ++sample_node){
                m_multiple_relation_matrix[ref_node][dist][sample_node] = 0;
            }
        }
    }


    // build NHM
    for(int i = 0; i < m_problem_size; ++i) {
        for(int j = 0; j < m_problem_size; ++j) {
            m_node_matrix[i][j] = 0;
        }    
    }    
    // end of build NHM



    // selection ?????
    // selection ?????
    // selection ?????
    // printf("sample:\n");



    if (mode == 1) {
        for (int k = 0;k < (m_sample_size/2); k++) {
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

            // build node matrix
            // m_node_matrix[i][j] = node i position is idx j
            for (int pos = 0; pos < m_problem_size; ++pos) {
                m_node_matrix[m_samples[k][pos]][pos] ++;
            }
            // end of build node matrix
        }
    }

    else if (mode == 2) {
        for (int k = (m_sample_size / 2); k < m_sample_size; k++) {
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

            // build node matrix
            // m_node_matrix[i][j] = node i position is idx j
            for (int pos = 0; pos < m_problem_size; ++pos) {
                m_node_matrix[m_samples[k][pos]][pos] ++;
            }
            // end of build node matrix
        }
    }

    else if (mode == 3) {
        for (int k = 0; k < m_sample_size; k++) {
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

            // build node matrix
            // m_node_matrix[i][j] = node i position is idx j
            for (int pos = 0; pos < m_problem_size; ++pos) {
                m_node_matrix[m_samples[k][pos]][pos] ++;
            }
            // end of build node matrix
        }
    }
    
    // PrintMRMatrix(m_multiple_relation_matrix, m_problem_size);

    return true;
}

/*
 * Virtual sampling function.
 */
void CACHIGN::Sample(int * genes, int * templateee, int need_to_sample_gene_count, double chi_std)
{
    // chi-std: the number of std use in chi-square
    int already_sampled_idx_count = m_problem_size - need_to_sample_gene_count;
    if (need_to_sample_gene_count == m_problem_size) {
        already_sampled_idx_count = 1;
    }
    int total_edge_count = already_sampled_idx_count * (m_problem_size - already_sampled_idx_count);


    // random choose (1/5 * ell) indexs as not sampled indexs and (4/5 * ell) indexs as already sampled indexs
    int sampled_or_not_idx[m_problem_size];




    // ===== random cut point start here ===== //
    if (m_adj_ratio < 0.25) {
        for (int i = 0; i < m_problem_size; ++i) {
            sampled_or_not_idx[i] = i;
        }
        random_shuffle(sampled_or_not_idx, sampled_or_not_idx + m_problem_size);
    }
    // PrintArray(sampled_or_not_idx, m_problem_size);
    // ===== end of random cut point===== //

    
    
    // ===== sequential cut point start here ===== //
    else { //adj_ratio >= 0.5

        int seeed = rand() % m_problem_size;

        for (int i = already_sampled_idx_count; i < m_problem_size; ++i) {
            sampled_or_not_idx[i] = seeed;
            if (sampled_or_not_idx[i] >= m_problem_size) {
                sampled_or_not_idx[i] -= m_problem_size;
            }
            // seeed += 5;
            seeed += 1;

        }
        // random_shuffle(sampled_or_not_idx, sampled_or_not_idx + m_problem_size);
        int tmpppp = 0;
        for (int i = 0; i < already_sampled_idx_count;) {
            bool already_used = false;
            

            for (int j = already_sampled_idx_count; j < m_problem_size; ++j) {

                if (sampled_or_not_idx[j] == tmpppp) {
                    already_used = true;
                    break;
                }
            }
            if (already_used == false) {
                sampled_or_not_idx[i] = tmpppp;
                i++;
            }
            tmpppp++;
        }
    }
    // ===== end of sequential cut point end ===== //
    



    // cout << "sampled_or_not_idx: ";
    // PrintArray(sampled_or_not_idx, m_problem_size);

    
    bool is_sampled_idx[m_problem_size];
    bool is_sampled_node[m_problem_size];

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
    // PrintArray(is_sampled_idx, m_problem_size);


    // record and calculate [ (1/5 * ell) nodes not sampled node x (4/5 * ell) nodes ] edges entropy 
        // record already sampled nodes position
        // calculate entropy
        // find a way to record which edges can be used and which cannot

    list<edgeee> edgeeeee_lst;

    //             create network edges and calculate entropy
    // time complexity:     O(n^2)       *       O(n)         =   O(n^3) 
    edgeee lowest_entropy_edge = Build_network(templateee, already_sampled_idx_count, edgeeeee_lst, sampled_or_not_idx, is_sampled_node);
    edgeee another_lowest_entropy_edge;


    // === sample === /
    // ---------- while still have genes not sampled yet ---------- //
    m_need_to_sample_genes_count += m_problem_size- already_sampled_idx_count;
    for (int kkk = already_sampled_idx_count; kkk < m_problem_size; kkk++) {


        // calculate the mean and std of entroy in edge list
        double entropy_sum = 0.0;
        for (auto it = edgeeeee_lst.begin(); it != edgeeeee_lst.end();++it) {
            entropy_sum += it->entropy;
        }
        double entropy_mean = entropy_sum / edgeeeee_lst.size();
        double tmp_sum = 0.0;
        for (auto it = edgeeeee_lst.begin(); it != edgeeeee_lst.end();++it) {
            tmp_sum += (it->entropy - entropy_mean) * (it->entropy - entropy_mean);
        }
        double entropy_std = sqrt(tmp_sum / edgeeeee_lst.size());
        // end of calculate the mean and std of entroy in edge list



        double lowest_entropy_probability_array[m_problem_size];
        int sample_idx;
        // if lowest entropy edge's entropy bigger than (mean - std), use node to sampled
        // === use node to sampled  === //
        if (lowest_entropy_edge.entropy > entropy_mean - chi_std * entropy_std) {
            
            // if use node to sample, construct another probability array
            double lowest_node_entropy = 100.0;
            int lowest_node_entropy_idx = -1;
            double tmp_probability_array[m_problem_size];

            // find lowest entropy node
            for (int idx = 0; idx < m_problem_size; ++idx) {
                if (is_sampled_idx[idx] == false) {
                    construct_node_probability_array(tmp_probability_array, is_sampled_node, idx);
                    double tmp_entropy = calculate_entropy(tmp_probability_array, m_problem_size);
                    // PrintArray(tmp_probability_array, m_problem_size);
                    if (tmp_entropy < lowest_node_entropy) {
                        tmp_entropy = lowest_node_entropy;
                        lowest_node_entropy_idx = idx;
                        copy_array(lowest_entropy_probability_array, tmp_probability_array, m_problem_size);
                    }
                }
            }
            // end of find lowest entropy node
            sample_idx = lowest_node_entropy_idx;
            // cout << "sample idx: " << sample_idx << endl;
            m_used_node_count++;
        }
        // === end of use node to sampled  === //
        // === use edge to sampled === //
        else {
            construct_probability_array(lowest_entropy_probability_array, is_sampled_node, lowest_entropy_edge);
            sample_idx = lowest_entropy_edge.node_a_index + lowest_entropy_edge.distance;
            m_used_distance_count[lowest_entropy_edge.distance] ++;
            // int sample_idx = edge_lst[lowest_entropy_edge_num].node_a_index + edge_lst[lowest_entropy_edge_num].distance;
            if (sample_idx >= m_problem_size) {
                sample_idx -= m_problem_size;
            }
        }
        // === end of use edge to sampled === //
        // end of if lowest entropy edge's entropy bigger than (mean - std), use node to sampled


        


        genes[sample_idx] = sample_from_array(lowest_entropy_probability_array, m_problem_size);
        // cout << "genes[" << sample_idx << "]: " << genes[sample_idx] << endl;

        is_sampled_idx[sample_idx] = true;
        is_sampled_node[genes[sample_idx]] = true;
        // === end of sample

        // update edges information and delete unsable edgesss
        lowest_entropy_edge = update_and_find_lowest_entropy_edge(sample_idx, genes[sample_idx], edgeeeee_lst);



        // add new edges
        another_lowest_entropy_edge = add_edgesss(genes[sample_idx], sample_idx, edgeeeee_lst, is_sampled_node, is_sampled_idx);
        // cout << "edgeeeee_lst size after adding new edges: " << edgeeeee_lst.size() << endl << endl;


        if (lowest_entropy_edge.entropy > another_lowest_entropy_edge.entropy) {
            lowest_entropy_edge = another_lowest_entropy_edge;
        }
        else if (lowest_entropy_edge.entropy == another_lowest_entropy_edge.entropy && rand() % 2 == 0) {
            lowest_entropy_edge = another_lowest_entropy_edge;
        }
    }
    // ----------  end of while still have genes not sampled yet ---------- //
    // cout << "genes: ";
    // PrintArray(genes, m_problem_size);
    // cout  << endl;
}


// similar to Build_network
edgeee CACHIGN::add_edgesss(int sampled_node_num, int sampled_node_idx, list<edgeee>& edge_lst, bool * is_sampled_node, bool * is_sampled_idx) {
    

    edgeee tmp_edge;
    edgeee lowest_entropy_edge;
    bool is_this_edge_first_updated = true;



    // for connect "the sampled node num" with all "not sampled yet index"
    for (int i = 0; i < m_problem_size; i++) {
        if (is_sampled_idx[i] == false) {
            tmp_edge.node_a = sampled_node_num;
            tmp_edge.node_a_index = sampled_node_idx;

            tmp_edge.distance = i - sampled_node_idx;
            if (tmp_edge.distance < 0) {
                tmp_edge.distance += m_problem_size;
            }

            tmp_edge.element_sum = 0;
            int tmp_probability_array[m_problem_size];

            for (int k = 0; k < m_problem_size; ++k) {
                
                // if that node num is not sampled yet
                if (is_sampled_node[k] == false) {
                    tmp_probability_array[k] = m_multiple_relation_matrix[tmp_edge.node_a][tmp_edge.distance][k];
                    tmp_edge.element_sum += m_multiple_relation_matrix[tmp_edge.node_a][tmp_edge.distance][k];
                }
                // else if that node num is already sampled
                else {
                    tmp_probability_array[k] = 0;
                }
            }
            tmp_edge.entropy = fast_calculate_entropy(tmp_probability_array, tmp_edge.element_sum);
            if (tmp_edge.entropy < lowest_entropy_edge.entropy || is_this_edge_first_updated) {
                lowest_entropy_edge = tmp_edge;
                is_this_edge_first_updated = false;
            }
            else if (tmp_edge.entropy == lowest_entropy_edge.entropy && rand() % 2 == 0) {
                lowest_entropy_edge = tmp_edge;
            }

            tmp_edge.can_be_used = true;
            edge_lst.push_back(tmp_edge);
            // num ++;
        }
    }
    return lowest_entropy_edge;
}



edgeee CACHIGN::Build_network(int * templateee, int already_sampled_idx_count, list<edgeee>&  edge_lst, int * sampled_or_not_idx, bool * is_sampled_node) {

    double lowest_entropy = 100.0;

    edgeee tmp_edge;
    edgeee lowest_entropy_edge;
    lowest_entropy_edge.entropy = 100;

    // === edge for loop === //
    for (int i = 0; i < already_sampled_idx_count; ++i) {                             // already sampled node
        for (int j = 0; j < (m_problem_size - already_sampled_idx_count); ++j) {      // not sampled yet node

            tmp_edge.node_a = templateee[sampled_or_not_idx[i]];
            tmp_edge.node_a_index = sampled_or_not_idx[i];

            // not_sampled_yet_node_idx - already_smpled_node_idx
            tmp_edge.distance = sampled_or_not_idx[already_sampled_idx_count + j] - sampled_or_not_idx[i];
            if (tmp_edge.distance < 0) {
                tmp_edge.distance += m_problem_size;
            }

            tmp_edge.element_sum = 0;
            int tmp_probability_array[m_problem_size];

            for (int k = 0; k < m_problem_size; ++k) {
                
                // if that node num is not sampled yet
                if (is_sampled_node[k] == false) {
                    tmp_probability_array[k] = m_multiple_relation_matrix[tmp_edge.node_a][tmp_edge.distance][k];
                    tmp_edge.element_sum += m_multiple_relation_matrix[tmp_edge.node_a][tmp_edge.distance][k];
                }
                // else if that node num is already sampled
                else {
                    tmp_probability_array[k] = 0;
                }
            }
            tmp_edge.entropy = fast_calculate_entropy(tmp_probability_array, tmp_edge.element_sum);
            tmp_edge.can_be_used = true;

            if (tmp_edge.entropy < lowest_entropy_edge.entropy || edge_lst.empty()) {
                lowest_entropy_edge = tmp_edge;
            }
            else if (tmp_edge.entropy == lowest_entropy_edge.entropy && rand() % 2 == 0) {
                lowest_entropy_edge = tmp_edge;
            }

            edge_lst.push_back(tmp_edge);
        }
    }
    // === end of edge for loop === //
    return lowest_entropy_edge;
}


double CACHIGN::fast_calculate_entropy(int * arr, int summation) {

    double ans = 0;
    for (int i = 0; i < m_problem_size; i++) {

        if (arr[i] != 0) {
            double p = double( arr[i] ) / summation;
            ans += double(m_plogp[summation][arr[i]]);
        }
    }
    return -ans;
}



// update edges information and delete unsable edgesss
edgeee CACHIGN::update_and_find_lowest_entropy_edge(int sample_idx, int sample_node_num, list<edgeee>&  edge_lst) {

    edgeee lowest_entropy_edge;
    lowest_entropy_edge.entropy = 100.0;
    bool is_this_edge_first_updated = true;


    // === update all edges entropy and edges information ===//
    for (auto it = edge_lst.begin(); it != edge_lst.end();) {
            
        int tmp_sample_idx = it->node_a_index + it->distance;
        if (tmp_sample_idx >= m_problem_size) {
            tmp_sample_idx -= m_problem_size;
        }
        // --- delete same distance edges --- //
        if (tmp_sample_idx == sample_idx) {
            it = edge_lst.erase(it);
        }
        // ---  end of delete same distance edges --- //

        // --- update different distance edges information --- //
        else {
            double deleted_element = m_multiple_relation_matrix[it->node_a][it->distance][sample_node_num];
            it->entropy = update_entropy(it->entropy, deleted_element, it->element_sum);


            it->element_sum -= deleted_element;
            if (it->entropy < lowest_entropy_edge.entropy || is_this_edge_first_updated) {
                lowest_entropy_edge = *it;
                is_this_edge_first_updated = false;
            }
            else if (it->entropy == lowest_entropy_edge.entropy && rand() % 2 == 0) {
                lowest_entropy_edge = *it;
            }
            ++it;
        }
        // --- end of update different distance edges information --- //
    }
    // === end of update all edges entropy and edges information === //
    return lowest_entropy_edge;
}



void CACHIGN::construct_probability_array(double * tmp_probability_array, bool * is_sampled_node, edgeee lowest_entropy_edge) {

    for (int i = 0; i < m_problem_size; ++i) {

        // if that node num is not sampled yet
        if (is_sampled_node[i] == false) {
            tmp_probability_array[i] = m_epsilon + m_multiple_relation_matrix[lowest_entropy_edge.node_a][lowest_entropy_edge.distance][i];
        }

        // else if that node num is already sampled
        else {
            tmp_probability_array[i] = 0;
        }
    }
}


void CACHIGN::construct_node_probability_array(double * tmp_probability_array, bool * is_sampled_node, int idx) {

    for (int i = 0; i < m_problem_size; ++i) {

        // if that node num is not sampled yet
        if (is_sampled_node[i] == false) {
            tmp_probability_array[i] = m_epsilon + m_node_matrix[idx][i];
        }

        // else if that node num is already sampled
        else {
            tmp_probability_array[i] = 0;
        }
    }
}