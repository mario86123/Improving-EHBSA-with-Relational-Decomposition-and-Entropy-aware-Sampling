
#include "NNHBSAWT.h"
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
CNNHBSAWT::CNNHBSAWT(int problem_size, int sel_size, double b_ratio, int cut_point_count, CPopulation * population)
{
	m_problem_size=problem_size;
    m_sample_size=sel_size;
    m_samples= new int*[m_sample_size];

    m_NHM = new double*[problem_size];
    AllocateMatrixMemory(m_NHM, problem_size);

    m_epsilon = b_ratio * m_sample_size / (m_problem_size);
    m_cut_point_count = cut_point_count;


    // clean matrix
    // reset model
    resetNHM(m_NHM, m_problem_size, m_epsilon);

    for(int k=0;k<m_sample_size;k++) {
		m_samples[k] = population->m_individuals[k]->Genes();
    }

    buildNHM(m_samples, m_sample_size, m_NHM, m_problem_size);


}

CNNHBSAWT::~CNNHBSAWT()
{

    delete [] m_samples;

    // free matrix memory
	freeMatrixMemory(m_NHM, m_problem_size);
	delete [] m_NHM;
}

/*
 * Virtual learning function.
 */


// update NHM
bool CNNHBSAWT::Learn(int * new_chromosome, int * templateee) {


    for (int idx = 0; idx < m_problem_size; idx++) {
        m_NHM [idx] [ new_chromosome[idx] ] ++;
        m_NHM [idx] [ templateee[idx] ] --;
    }

    // PrintMatrix(m_NHM, m_problem_size);
    return true;
}


/*
 * Virtual sampling function.
 */
void CNNHBSAWT::Sample(int * genes, int *templateee)
{
    // for example, ell = 20
    // cut point count = 1 --> need_to_sample_gene_count = 20
    // cut point count = 2 --> need_to_sample_gene_count = 10
    // cut point count = 3 --> need_to_sample_gene_count = 
    // cut point count = 4 --> need_to_sample_gene_count = 5
    // cut point count = 5 --> need_to_sample_gene_count = 4

    // so NHBSAWO == NNHBSAWT, cut point count = 1

    int need_to_sample_gene_count = m_problem_size / m_cut_point_count;
    

    // random sample
    int sample_order_idx_arr[m_problem_size];
	for (int i = 0; i < m_problem_size; ++i) {
		sample_order_idx_arr[i] = i;
	}
	random_shuffle(sample_order_idx_arr, sample_order_idx_arr + m_problem_size);
    
    
    bool * is_sample_node_arr= new bool[m_problem_size]();


    // set do not need to resample node to default
    for(int i = 0; i < m_problem_size; ++i){
        
        if (i < m_problem_size - need_to_sample_gene_count) {
            genes[sample_order_idx_arr[i]] = templateee[sample_order_idx_arr[i]];
            is_sample_node_arr[genes[sample_order_idx_arr[i]]] = true;
        }
    }

    for (int sample_count = m_problem_size - need_to_sample_gene_count; sample_count < m_problem_size; sample_count++) {
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

        is_sample_node_arr[genes[sample_order_idx_arr[sample_count]]] = true;
    }
    
    delete [] is_sample_node_arr;
}