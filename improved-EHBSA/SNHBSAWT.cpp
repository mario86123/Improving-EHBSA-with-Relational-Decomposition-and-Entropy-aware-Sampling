
#include "SNHBSAWT.h"
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
CSNHBSAWT::CSNHBSAWT(int problem_size, int sel_size, double b_ratio, int cut_point_count)
{
	m_problem_size=problem_size;
    m_sample_size=sel_size;
    m_samples= new int*[m_sample_size];

    m_NHM = new double*[problem_size];
    AllocateMatrixMemory(m_NHM, problem_size);

    m_epsilon = b_ratio * m_sample_size / (m_problem_size -1);
    m_cut_point_count = cut_point_count;
}

CSNHBSAWT::~CSNHBSAWT()
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
bool CSNHBSAWT::Learn(CPopulation * population, int size)
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
void CSNHBSAWT::Sample(int * genes, int *templateee)
{

    bool * is_sample_node_arr= new bool[m_problem_size]();

    // --- print is_sample_node_arr 
    // printf("is_sample_node_arr: ");
    // for (int i =0; i < m_problem_size; i++) {
    //     cout << is_sample_node_arr[i] << " ";
    // }
    // cout << endl;
    // --- end of print is_sample_node_arr 


    // for example, ell = 20
    // cut point count = 1 --> need_to_sample_gene_count = 20
    // cut point count = 2 --> need_to_sample_gene_count = 10
    // cut point count = 3 --> need_to_sample_gene_count = 
    // cut point count = 4 --> need_to_sample_gene_count = 5
    // cut point count = 5 --> need_to_sample_gene_count = 4

    // so NHBSAWO == SNHBSAWT, cut point count = 1

    if (m_cut_point_count == 1) // NHBSA/WO
    {
        for (int idx = 0; idx < m_problem_size; idx++) {
            // printf("%d ", idx);
            
            double *rw_vector = new double[m_problem_size];
            
            // assign rw_vector value
            for(int j = 0; j < m_problem_size; ++j) {

                if (is_sample_node_arr[j] == true) {
                    rw_vector[j] = 0;
                }
                else { // is_sample_node_arr[i] == false
                    rw_vector[j] = m_NHM [idx] [j];
                }
            }

            genes[idx] = sample_from_array(rw_vector, m_problem_size);
            
            is_sample_node_arr[genes[idx]] = true;

            delete [] rw_vector;

        }
    } 

    else  // NHBSA/WT
    {
        // cut point generator
        // cut_points[0] -> start, cut_points[1] -> end
        int cut_points[2]; 
        
        CutPointGenerator(m_problem_size, m_cut_point_count, cut_points);
        // printf("%d %d dist: %d\n", cut_points[0], cut_points[1], ((cut_points[1] + m_problem_size) - cut_points[0]) % m_problem_size);



        // initiailize new_sample wiith template
        // template = population[template_idx]
        // printf("not change genes: ");
        for (int i = cut_points[1]; i != cut_points[0];) {
            genes[i] = templateee[i];
            is_sample_node_arr[templateee[i]] = true;
            // printf("%d ", i);
            i++;
            if (i >= m_problem_size) {
                i -= m_problem_size;
            }
        }
        // printf("\n");
        

        // printf("sample part: ");
        // sample new individual
        for (int idx = cut_points[0]; idx != cut_points[1];) {
            // printf("%d ", idx);
            
            double *rw_vector = new double[m_problem_size];
            
            // assign rw_vector value
            for(int j = 0; j < m_problem_size; ++j) {

                if (is_sample_node_arr[j] == true) {
                    rw_vector[j] = 0;
                }
                else { // is_sample_node_arr[i] == false
                    rw_vector[j] = m_NHM [idx] [j];
                }
            }

            genes[idx] = sample_from_array(rw_vector, m_problem_size);
            
            is_sample_node_arr[genes[idx]] = true;

            delete [] rw_vector;

            idx++;
            if (idx >= m_problem_size) {
                idx -= m_problem_size;
            }

        }
    }

    // free the memory
    delete [] is_sample_node_arr;
}
