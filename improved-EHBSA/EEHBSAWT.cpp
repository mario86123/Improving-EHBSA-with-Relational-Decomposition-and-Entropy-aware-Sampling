
#include "EEHBSAWT.h"
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
CEEHBSAWT::CEEHBSAWT(int problem_size, int sel_size, double b_ratio, int cut_point_count, CPopulation * population)
{
	m_problem_size=problem_size;
    m_sample_size=sel_size;
    m_samples= new int*[m_sample_size];

    m_EHM = new double*[problem_size];
    AllocateMatrixMemory(m_EHM, problem_size);

    m_epsilon = b_ratio * 2 * m_sample_size / (m_problem_size -1);
    // printf("eps: %lf\n", m_epsilon);
    m_cut_point_count = int(cut_point_count);


    resetMatrix(m_EHM, m_problem_size, m_epsilon);
    
    // selection ?????
    // selection ?????
    // selection ?????
    // printf("sample:\n");
    for(int k=0;k<m_sample_size;k++) {
		m_samples[k] = population->m_individuals[k]->Genes();\
    }
    buildModel(m_samples, m_sample_size, m_EHM, m_problem_size);
}

CEEHBSAWT::~CEEHBSAWT()
{
    // PrintMatrix(m_EHM, m_problem_size);

    delete [] m_samples;

    // free matrix memory
	freeMatrixMemory(m_EHM, m_problem_size);
	delete [] m_EHM;
}

/*
 * Virtual learning function.
 */
// size: selection size
bool CEEHBSAWT::Learn(int * new_chromosome, int * templateee)
{   
    for (int gene_num = 0; gene_num < m_problem_size; ++gene_num) {

        int front_idx = gene_num - 1, back_idx = gene_num + 1;

        // first gene
        if (gene_num == 0) {
            front_idx += m_problem_size;
        }

        // last gene
        else if (gene_num == m_problem_size - 1) {
            back_idx -= m_problem_size;
        }

        m_EHM[ new_chromosome[gene_num] ][ new_chromosome[front_idx] ] += 1;
        m_EHM[ new_chromosome[gene_num] ][ new_chromosome[back_idx] ] += 1;

        m_EHM[ templateee[gene_num] ][ templateee[front_idx] ] -= 1;
        m_EHM[ templateee[gene_num] ][ templateee[back_idx] ] -= 1;


    }
    return true;
}

/*
 * Virtual sampling function.
 */
void CEEHBSAWT::Sample(int * genes, int *templateee)
{
    // printf("\n");
    // PrintArray(templateee, m_problem_size);
    bool * isSampledNode = new bool[m_problem_size];

    // set isSampledNode to default
    for(int i = 0; i < m_problem_size; ++i){
        isSampledNode[i] = false;
    }

    // sample new individual

    // EHBSA / WO
    if (m_cut_point_count == 1) {
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
            
    }
    // EHBSA / WT
    else {

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
            isSampledNode[templateee[i]] = true;
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

                if (isSampledNode[j] == true) {
                    rw_vector[j] = 0;
                }
                else { // isSampledNode[i] == false
                    if (idx == 0){
                        rw_vector[j] = m_EHM [ genes[m_problem_size-1] ] [j];
                    }
                    else {
                        rw_vector[j] = m_EHM [ genes[idx-1] ] [j];
                    }
                }
            }

            genes[idx] = sample_from_array(rw_vector, m_problem_size);
            
            isSampledNode[genes[idx]] = true;

            delete [] rw_vector;

            idx++;
            if (idx >= m_problem_size) {
                idx -= m_problem_size;
            }

        }
    }
    // free the memory
    delete [] isSampledNode;
    // printf("\n\n");


}