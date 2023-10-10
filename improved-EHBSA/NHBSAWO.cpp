
#include "NHBSAWO.h"
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
CNHBSAWO::CNHBSAWO(int problem_size, int sel_size, double b_ratio)
{
	m_problem_size=problem_size;
    m_sample_size=sel_size;
    m_samples= new int*[m_sample_size];

    m_NHM = new double*[problem_size];
    AllocateMatrixMemory(m_NHM, problem_size);

    m_epsilon = b_ratio * m_sample_size / (m_problem_size -1);
}

CNHBSAWO::~CNHBSAWO()
{
    // PrintMatrix(m_multiple_relation_matrix, m_problem_size);

    delete [] m_samples;

    // free matrix memory
	freeMatrixMemory(m_NHM, m_problem_size);
	delete [] m_NHM;
}

/*
 * Virtual learning function.
 */
// size: selection size
bool CNHBSAWO::Learn(CPopulation * population, int size)
{   
    // population->Print();

    // clean matrix
    // reset model
    resetNHM(m_NHM, m_problem_size, m_epsilon);
    
    // selection ?????
    // selection ?????
    // selection ?????
    // printf("sample:\n");
    for(int k=0;k<size;k++) {
		m_samples[k] = population->m_individuals[k]->Genes();
        // for(int i = 0; i < m_problem_size; ++i){
        //     printf("%d ", m_samples[k][i]);
        //     // total += genes[i];
        // }
        // printf("\n");
    }
    // printf("m_sample_size: %d\n", m_sample_size);
    buildNHM(m_samples, m_sample_size, m_NHM, m_problem_size);
    // PrintMatrix(m_NHM, m_problem_size);
    return true;
}

/*
 * Virtual sampling function.
 */
void CNHBSAWO::Sample(int * genes)
{

    // random sample
    int sample_order_idx_arr[m_problem_size];
	for (int i = 0; i < m_problem_size; ++i) {
		sample_order_idx_arr[i] = i;
	}
	random_shuffle(sample_order_idx_arr, sample_order_idx_arr + m_problem_size);


    bool is_sample_node_arr[m_problem_size];
    // set isSampledNode to default
    for(int i = 0; i < m_problem_size; ++i){
        is_sample_node_arr[i] = false;
    }




    for (int sample_count = 0; sample_count < m_problem_size; sample_count++) {
        
        double rw_vector[m_problem_size];

        for(int node_num = 0; node_num < m_problem_size; ++node_num) {

            if (is_sample_node_arr[node_num] == true) {
                rw_vector[node_num] = 0;
            }
            else { // isSampledNode[i] == false
                rw_vector[node_num] = m_NHM [sample_order_idx_arr[sample_count]] [node_num];
            }
        }

        genes[sample_order_idx_arr[sample_count]] = sample_from_array(rw_vector, m_problem_size);
        
        if (genes[sample_order_idx_arr[sample_count]] == -1) {
            cout << sample_order_idx_arr[sample_count] << endl;
            for (int i =0; i < m_problem_size; i++) {
                cout << is_sample_node_arr[i] << " ";
            }
            cout <<endl;
            PrintArray(rw_vector, m_problem_size);
        }

        is_sample_node_arr[genes[sample_order_idx_arr[sample_count]]] = true;
    }


    // print sample result chromosome
    // int total = 0;
    // for(int i = 0; i < m_problem_size; ++i){
    //     // printf("%d ", genes[i]);
    //     total += genes[i];
    // }
    // printf("\n");
    // cout <<total<< endl;
}