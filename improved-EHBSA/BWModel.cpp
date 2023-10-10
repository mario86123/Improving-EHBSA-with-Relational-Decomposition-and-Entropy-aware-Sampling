
#include "BWModel.h"
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
CBWModel::CBWModel(int problem_size, int sel_size, double b_ratio)
{
	m_problem_size=problem_size;
    m_sample_size=sel_size;
    m_samples= new int*[m_sample_size];

    m_half_relation_distance = m_problem_size / 2 + 1;
    // m_half_relation_distance = 1;
    m_half_relation_matrix = new double*[m_problem_size];

    for(int i = 0; i < m_problem_size; ++i){
        m_half_relation_matrix[i] = new double[m_problem_size]();
    }

    // b_ratio = 1 ( if don't know how to set )
    m_epsilon = double(m_problem_size * b_ratio) / (m_sample_size);

    // printf("m_epsilon: %lf\n", m_epsilon);
    // printf("m_problem_size: %d \n", m_problem_size);
    // printf("m_sample_size: %d \n", m_sample_size);
    // printf("m_half_relation_distance: %d \n", m_half_relation_distance);
}

CBWModel::~CBWModel()
{
    PrintMatrix(m_half_relation_matrix, m_problem_size);

    delete [] m_samples;
    for(int i = 0; i < m_problem_size; ++i){
        delete[] m_half_relation_matrix[i];
    }
    delete[] m_half_relation_matrix;
    
}

/*
 * Virtual learning function.
 */
// size = selection size
bool CBWModel::Learn(CPopulation * population, int size)
{   
    // population->Print();

    // clean matrix
    for(int i = 0; i < m_problem_size; ++i){
        for(int j = 0; j < m_problem_size; ++j){
            
            if (i == j) {
                // printf("i: %d, j: %d\n", i, j);
                m_half_relation_matrix[i][j] = 0.0;
            } 
            else {
                m_half_relation_matrix[i][j] = m_epsilon;
            }
        }
    }
    
    
    // selection ?????
    // selection ?????
    // selection ?????
    // printf("sample:\n");
    for(int k=0;k<size;k++){
		m_samples[k] = population->m_individuals[k]->Genes();
        
        
        // for(int i = 0; i < m_problem_size; ++i){
        //     printf("%d ", m_samples[k][i]);
        //     // total += genes[i];
        // }
        // printf("\n");

        // build matrix
        // assign value to matrix
        for(int gene_num = 0; gene_num < m_problem_size; ++gene_num){
            
            int second_gene = gene_num + m_half_relation_distance;
            if (second_gene >= m_problem_size) {
                second_gene -= m_problem_size;
            }
            m_half_relation_matrix[ m_samples[k][gene_num] ][ m_samples[k][second_gene]] += 1;
            
        }
    }
    // printMatrix(m_half_relation_matrix, m_problem_size);

    return true;
}

/*
 * Virtual sampling function.
 */
void CBWModel::Sample(int * genes)
{

    // record whether a node is sampled
    bool *is_sampled_idx= new bool[m_problem_size];
    bool *is_sampled_node= new bool[m_problem_size];
    for(int i = 0; i < m_problem_size; ++i){
        is_sampled_idx[i] = false;
        is_sampled_node[i] = false;
    }

    // initialize
    int sample_count = 0;
    int sampled_idx = 0;
    int current_idx;


    //while sample one chromosome procdure not finished yet
    while (sample_count < m_problem_size) {
        // printf("sample_count: %d\n", sample_count);

        double *sample_arr = new double[m_problem_size]();

        // two sample method
        // 1. have previous node
        // 2. do not have previous node


        // 2. do not have previous node
        // if that index is already sampled
        // or is sampling first node
        if (is_sampled_idx[sampled_idx] == true || sample_count == 0) {

            // find not sampled yet and lowest entropy column
            sampled_idx = -1;
            double min_entropy = 100000;

            // scan every column
            for (int i = 0; i < m_problem_size; ++i) {

                // count every node each column
                if (is_sampled_idx[i] == false) {
                    int *column_arr = new int[m_problem_size]();
                    for (int j = 0; j < m_sample_size; ++j) {
                        column_arr[m_samples[j][i]]++;
                    }

                    // printf("%d th column: ", i);
                    // for (int j = 0; j < m_problem_size; ++j) {
                    //     printf("%d ", column_arr[j]);
                    // }

                    double current_entropy = calculate_entropy(column_arr, m_problem_size);
                    // cout << "entropy: " << current_entropy << endl;

                    if (current_entropy < min_entropy) {
                        min_entropy = current_entropy;
                        sampled_idx = i;
                        // cout << "min_entropy: " << min_entropy << endl;
                    }
                    delete [] column_arr;
                }
            }

            // sample first node
            for (int i = 0; i < m_sample_size; ++i){
                
                if (is_sampled_node[m_samples[i][sampled_idx]] == true) {
                    continue;
                }
                else {
                    sample_arr[m_samples[i][sampled_idx]] ++;
                }
            }
        }

        // 1. have previous node
        // that index is not sampled yet
        else {
            // can not sampled already sampled node QQ
            // printf("%d: ", genes[current_idx]);
            for (int i = 0; i < m_problem_size; ++i){
                
                if (is_sampled_node[i] == true) {
                    continue;
                } 
                else {
                    sample_arr[i] = m_half_relation_matrix[ genes[current_idx] ][i];
                }
                
            }

            // sample_arr maybe all zero
            // for (int i = 0; i < m_problem_size; ++i){
            //     printf("%d ", sample_arr[i]);
            // }
            // printf("\n");

        }
        genes[sampled_idx] = sample_from_array(sample_arr, m_problem_size);
        // printf("sampled_idx: %d\n", sampled_idx);  
        // printf("genes[sampled_idx]: %d\n", genes[sampled_idx]);        


        // if add m_eta, not go through here, i.e. sample_arr always not all zero
        // sample_arr maybe all zero
        // while (genes[sampled_idx] == -1) {
        //     for (int i = 0; i < m_problem_size; ++i){
        //         printf("%lf ", sample_arr[i]);  

        //         if (is_sampled_node[i] == false) {
        //             sample_arr[i] += 1;
        //         }
        //     }
        //     // printf("\n");  
        //     genes[sampled_idx] = sample_from_array(sample_arr);
        // }



        // printf("%d\n", genes[sampled_idx]);
        is_sampled_node[genes[sampled_idx]] = true;
        is_sampled_idx[sampled_idx] = true;
        sample_count += 1;

        current_idx = sampled_idx;
        sampled_idx += m_half_relation_distance;
        
        if (sampled_idx >= m_problem_size) {
            sampled_idx -= m_problem_size;
        }
        delete[] sample_arr;
    }
    // free memory
    // printf("\n");
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