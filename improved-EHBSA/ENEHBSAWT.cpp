// entropy EHBSA
#include "ENEHBSAWT.h"
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
CENEHBSAWT::CENEHBSAWT(int problem_size, int sel_size, double b_ratio, int cut_point_count)
{
	m_problem_size=problem_size;
    m_sample_size=sel_size;
    m_samples= new int*[m_sample_size];

    m_EHM = new double*[problem_size];
    AllocateMatrixMemory(m_EHM, problem_size);

    m_epsilon = b_ratio * 2 * m_sample_size / (m_problem_size -1);
    // printf("eps: %lf\n", m_epsilon);
    m_cut_point_count = int(cut_point_count);
}

CENEHBSAWT::~CENEHBSAWT()
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
bool CENEHBSAWT::Learn(CPopulation * population, int size)
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
    buildModel(m_samples, m_sample_size, m_EHM, m_problem_size);
    // PrintMatrix(m_EHM, m_problem_size);
    return true;
}

/*
 * Virtual sampling function.
 */
void CENEHBSAWT::Sample(int * genes, int *templateee)
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
        // printf("not change genes index: ");
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
        

        // === sample new individual === //
        int start_idx = cut_points[0];
        int end_idx = cut_points[1] - 1;
        if (end_idx == -1) end_idx = m_problem_size - 1;
        // int end_idx = cut_points[1];

        double *rw_vector_start = new double[m_problem_size];
        double *rw_vector_end = new double[m_problem_size];

        double start_entropy, end_entropy;

        // sample for loop
        for (; start_idx != (end_idx + 1) % m_problem_size;) {
            // printf("%d ", start_idx);
            // printf("%d\n", end_idx);
        
            
            // assign rw_vector_start value
            for(int j = 0; j < m_problem_size; ++j) {

                if (isSampledNode[j] == true) {
                    rw_vector_start[j] = 0;
                }
                else { // isSampledNode[i] == false
                    if (start_idx == 0){
                        rw_vector_start[j] = m_EHM [ genes[m_problem_size-1] ] [j];
                    }
                    else {
                        rw_vector_start[j] = m_EHM [ genes[start_idx-1] ] [j];
                    }
                }
            }
            start_entropy = calculate_entropy(rw_vector_start, m_problem_size);
            // end of assign rw_vector_start value
            // printf("start_entropy: %lf\n", start_entropy);


            // assign rw_vector_end value
            for(int j = 0; j < m_problem_size; ++j) {

                if (isSampledNode[j] == true) {
                    rw_vector_end[j] = 0;
                }
                else { // isSampledNode[i] == false
                    if (end_idx == m_problem_size - 1){
                        rw_vector_end[j] = m_EHM [ genes[0] ] [j];
                    }
                    else {
                        rw_vector_end[j] = m_EHM [ genes[end_idx+1] ] [j];
                    }
                }
            }
            end_entropy = calculate_entropy(rw_vector_end, m_problem_size);
            // end of assign rw_vector_end value


            // printf("end_entropy: %lf\n", end_entropy);

            if (start_entropy < end_entropy) { // use rw_vector_start to sample next gene
                genes[start_idx] = sample_from_array(rw_vector_start, m_problem_size);
                
                // cout << "isSampledNode 1: ";
                // PrintArray(isSampledNode, m_problem_size);

                isSampledNode[genes[start_idx]] = true;
                
                // cout << "genes[" << start_idx << "]: " << genes[start_idx] << endl;
                // cout << "isSampledNode 2: ";
                // PrintArray(isSampledNode, m_problem_size);
                start_idx ++;

                if (start_idx >= m_problem_size) {
                    start_idx -= m_problem_size;
                }
            }
        
            else { // (start_entropy > end_entropy) ==> use rw_vector_end to sample next gene
                genes[end_idx] = sample_from_array(rw_vector_end, m_problem_size);

                // cout << "isSampledNode 1: ";
                // PrintArray(isSampledNode, m_problem_size);

                isSampledNode[genes[end_idx]] = true;

                // cout << "genes[" << end_idx << "]: " << genes[end_idx] << endl;
                // cout << "isSampledNode 2: ";
                // PrintArray(isSampledNode, m_problem_size);

                end_idx --;

                if (end_idx < 0) {
                    // printf("before: %d\n", end_idx);
                    end_idx += m_problem_size;
                    // printf("after: %d\n", end_idx);
                }
            }
            // cout << endl;

        } // end of sample for loop


        delete [] rw_vector_start;
        delete [] rw_vector_end;
        // === end of sample new individual === //
    }
    // end of EHBSA / WT else


    // free the memory
    delete [] isSampledNode;
    // printf("\n\n");

    // cout << "tmemplate: ";
    // PrintArray(templateee, m_problem_size);
    // cout << "new      : ";
    // PrintArray(genes, m_problem_size);
    // cout << endl;

}