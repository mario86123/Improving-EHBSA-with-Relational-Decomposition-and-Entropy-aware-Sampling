
#include "MRModel.h"
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
CMRModel::CMRModel(int problem_size, int sel_size, double b_ratio)
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

    first_node_idx_arr = new double[m_problem_size]();

    // b_ratio = 1 ( if don't know how to set )
    // epsilon = b * ell / pop_siez
    m_epsilon = double(m_problem_size * b_ratio) / (m_sample_size); // pop_siez = sample_size

    // NHBSA like
    m_epsilon = b_ratio * m_sample_size / (m_problem_size);

    
    // printf("m_problem_size: %d \n", m_problem_size);
    // printf("m_sample_size: %d \n", m_sample_size);
}

CMRModel::~CMRModel()
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
    
    delete[] first_node_idx_arr;
}

/*
 * Virtual learning function.
 */
bool CMRModel::Learn(CPopulation * population, int size)
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
    
    for(int i = 0; i < m_problem_size; ++i) {
        first_node_idx_arr[i] = m_epsilon;
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


    // ========= count node's num distribution at first node idx ============ //



    // find lowest entropy column to decide first node idx
    first_node_idx = -1;
    double min_entropy = 100000;


    // scan every column to find out which has lowest entropy
    for (int i = 0; i < m_problem_size; ++i) {

        // count every node each column
        int *column_arr = new int[m_problem_size]();
        for (int j = 0; j < m_sample_size; ++j) {
            column_arr[m_samples[j][i]]++;
        }

        double current_entropy = calculate_entropy(column_arr, m_problem_size);
        // cout << "entropy: " << current_entropy << endl;

        if (current_entropy < min_entropy) {
            min_entropy = current_entropy;
            first_node_idx = i;
        }
        delete [] column_arr;
    }
    // ===== end of scan every column ===== // 
    // printf("first_node_idx: %d\n", first_node_idx);


    // for every chromosome in samples, count node number at first_node_idx
    for (int i = 0; i < m_sample_size; ++i){
        first_node_idx_arr[ m_samples[i][first_node_idx] ] ++;
        // first_node_idx_arr[ m_samples[i][0] ] ++;
    }
    // first_node_idx = 0;
    return true;
}

/*
 * Virtual sampling function.
 */
void CMRModel::Sample(int * genes)
{

    // record whether an idx is sampled
    bool *is_sampled_idx= new bool[m_problem_size];

    // record whether a node is sampled
    bool *is_sampled_node= new bool[m_problem_size];

    for(int i = 0; i < m_problem_size; ++i){
        is_sampled_idx[i] = false;
        is_sampled_node[i] = false;
    }


    // === sample first node === //
    int first_node_num = sample_from_array(first_node_idx_arr, m_problem_size);
    genes[first_node_idx] = first_node_num;
    is_sampled_idx[first_node_idx] = true;
    is_sampled_node[first_node_num] = true;
    // printf("first node -> genes[%d] = %d\n", first_node_idx, first_node_num);
    // === sample first node end === //

    

    // === sample remaining nodes === //
    int current_node_num = first_node_num;
    int current_node_idx = first_node_idx;
    int next_node_num;
    int next_node_idx;


    // ----- sample all nodes except for first node ----- //
    // "sample_count = 1" because first node is already sampled
    for (int sample_count = 1; sample_count < m_problem_size; ++sample_count) {

        // sample next node using distance entropy
        // find not sampled yet and lowest entropy column
        // int sampled_dist = -1;
        double min_entropy = 100000;
        double *tmp_arr = new double[m_problem_size]();
        double *best_arr = new double[m_problem_size]();

        // printf("current_node_idx: %d\n", current_node_idx);


        // scan every column and find out which column has lowest entropy
        for (int dist = 1; dist < m_problem_size; ++dist) {

            int tmp_idx = current_node_idx + dist;
            if (tmp_idx >= m_problem_size) {
                tmp_idx -= m_problem_size;
            }


            // if tmp_idx's column is not sampled yet
            if (is_sampled_idx[tmp_idx] == false) {

                // printf("idx: %d, tmp_arr: ", tmp_idx);
                for (int i = 0; i < m_problem_size; ++i){
                    if (is_sampled_node[i] == false) {
                        
                        // use tmp_arr to store the frequency of the node number 
                        // appear at (current_node_num + dist) column
                        tmp_arr[i] = m_multiple_relation_matrix[current_node_num][dist][i];
                    }
                    // printf("%d ", tmp_arr[i]);
                }


                // printf("\n");
                double current_entropy = calculate_entropy(tmp_arr, m_problem_size);
                // cout << "\n" << dist<<" entropy: " << current_entropy << endl;
                // cout << "min_entropy: " << min_entropy << endl;

                if (current_entropy < min_entropy) {
                    min_entropy = current_entropy;
                    // sampled_dist = dist;
                    next_node_idx = tmp_idx;
                    for (int i = 0; i < m_problem_size; ++i){
                        best_arr[i] = tmp_arr[i];
                    }
                }
            }
        }
        // printf("best_arr: ");
        // for (int i = 0; i < m_problem_size; ++i){
            // printf("%d ", best_arr[i]);
        // }
        // printf("\n");


        // sample next node
        // if sum of frequencies of all elements in best_arr are not zero
        next_node_num = sample_from_array(best_arr, m_problem_size);


        // printf("next_node_num = %d\n", next_node_num);
        genes[next_node_idx] = next_node_num;
        is_sampled_idx[next_node_idx] = true;
        is_sampled_node[next_node_num] = true;
        // printf("genes[%d] = %d\n", next_node_idx, next_node_num);

        current_node_num = next_node_num;
        current_node_idx = next_node_idx;

        delete []tmp_arr;
        delete []best_arr;
        // printf("%2d ", sampled_dist);
    }
    // === end of sample remaining nodes === //
    // printf("\n");


    // release memory
    delete[] is_sampled_idx;
    delete[] is_sampled_node;

    // print sample result chromosome
    // int total = 0;
    // for(int i = 0; i < m_problem_size; ++i){
    //     printf("%d ", genes[i]);
    //     // total += genes[i];
    // }
    // printf("\n");
    // cout <<total<< endl;
}