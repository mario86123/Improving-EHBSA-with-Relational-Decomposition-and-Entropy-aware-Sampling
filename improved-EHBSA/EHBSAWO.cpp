
#include "EHBSAWO.h"
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
CEHBSAWO::CEHBSAWO(int problem_size, int sel_size, double b_ratio)
{
	m_problem_size=problem_size;
    m_sample_size=sel_size;
    m_samples= new int*[m_sample_size];

    m_EHM = new double*[problem_size];
    AllocateMatrixMemory(m_EHM, problem_size);

    m_epsilon = b_ratio * 2 * m_sample_size / (m_problem_size -1);
}

CEHBSAWO::~CEHBSAWO()
{
    // PrintMatrix(m_multiple_relation_matrix, m_problem_size);

    delete [] m_samples;

    // free matrix memory
	freeMatrixMemory(m_EHM, m_problem_size);
	delete [] m_EHM;
}

/*
 * Virtual learning function.
 */
// size: selection size
bool CEHBSAWO::Learn(CPopulation * population, int size)
{   
    // population->Print();

    // clean matrix
    // reset model
    resetMatrix(m_EHM, m_problem_size, m_epsilon);
    
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
    buildModel(m_samples, m_sample_size, m_EHM, m_problem_size);
    // PrintMatrix(m_EHM, m_problem_size);
    return true;
}

/*
 * Virtual sampling function.
 */
void CEHBSAWO::Sample(int * genes)
{

    bool * isSampledNode = new bool[m_problem_size];

    // set isSampledNode to default
    for(int i = 0; i < m_problem_size; ++i){
        isSampledNode[i] = false;
    }

    // sample new individual
    genes[0] = rand() % m_problem_size;
    isSampledNode[genes[0]] = true;

    for (int idx = 1; idx < m_problem_size; ++idx) {
        
        double *rw_vector = new double[m_problem_size];
        
        // assign rw_vector value
        for(int j = 0; j < m_problem_size; ++j) {

            if (isSampledNode[j] == true) {
                rw_vector[j] = 0;
            }
            else { // isSampledNode[i] == false
                rw_vector[j] = m_EHM [genes [idx-1] ] [j];
            }
        }


        genes[idx] = sample_from_array(rw_vector, m_problem_size);
        isSampledNode[genes[idx]] = true;

        delete [] rw_vector;
    }
        
    // free the memory
    delete [] isSampledNode;

    // print sample result chromosome
    // int total = 0;
    // for(int i = 0; i < m_problem_size; ++i){
    //     // printf("%d ", genes[i]);
    //     total += genes[i];
    // }
    // printf("\n");
    // cout <<total<< endl;
}