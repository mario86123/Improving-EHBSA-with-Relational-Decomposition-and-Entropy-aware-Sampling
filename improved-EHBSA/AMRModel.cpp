// Accumulate MRModel
#include "AMRModel.h"
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
CAMRModel::CAMRModel(int problem_size, int sel_size, double b_ratio, int previous_sampled_reference_count, int seed)
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


    // epsilon = b * ell / pop_size
    m_epsilon = double(m_problem_size * b_ratio) / (m_sample_size); // pop_size = sample_size



    // ***** for entropy experiment ***** //

    m_sample_genes_entropy_arr = new double[m_problem_size]();
    m_already_sampled_chromosome_count = 0;


    int iteration = m_problem_size * m_problem_size / m_sample_size * 1000;
    // printf("iteration: %d\n", iteration);

    m_entropy_change = new double*[m_problem_size];
    for (int i = 0; i < m_problem_size; ++i) {
        m_entropy_change[i] = new double[iteration]();
    }

    m_iteration = 0;

    m_previous_sampled_reference_count = previous_sampled_reference_count;
    m_seed = seed;

    // ***** for entropy experiment ***** //

}













CAMRModel::~CAMRModel()
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



    // ***** for entropy experiment ***** //

    delete[] m_sample_genes_entropy_arr;


    int iteration = m_problem_size * m_problem_size / m_sample_size * 1000;

    // --- write file ---//
    ofstream output_file;
    // output_file.open("./genes_entropy/entropy/" + std::to_string(m_problem_size) + "_" + std::to_string(m_sample_size) + "_" + std::to_string(m_previous_sampled_reference_count) + "_" + std::to_string(m_seed) + ".txt");
    output_file.open("./genes_entropy.txt");

    for (int i = 0; i < m_problem_size; ++i) {
        output_file << i << ": " << endl;
        for (int j = 0; j < iteration; ++j) {
            output_file << m_entropy_change[i][j] << " ";
        }
        output_file << endl << endl;
    }

    output_file.close();
    // --- write file end ---//


    for (int i = 0; i < m_problem_size; ++i) {
        delete[] m_entropy_change[i];
    }
    delete[] m_entropy_change;

    // ***** for entropy experiment ***** //


}










/*
 * Virtual learning function.
 */
bool CAMRModel::Learn(CPopulation * population, int size)
{   
    // population->Print();

    // ----- for entropy experiment ----- //
    m_iteration ++;
    for (int i = 0; i < m_problem_size; ++i) {
        m_sample_genes_entropy_arr[i] = 0;
    }
    m_already_sampled_chromosome_count = 0;
    // ----- end of for entropy experiment ----- //



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
    
    



    // --- first sampled node random choose --- //
    // first_node_idx = rand() % m_problem_size;



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
void CAMRModel::Sample(int * genes)
{

    // record idx number that is sampled
    int *already_sampled_idx= new int[m_problem_size];

    // record node number that is sampled
    int *already_sampled_node= new int[m_problem_size];

    for(int i = 0; i < m_problem_size; ++i){
        already_sampled_idx[i] = -1;
        already_sampled_idx[i] = -1;
    }


    /*** === sample first node start === ***/
    int first_node_num = sample_from_array(first_node_idx_arr, m_problem_size);
    genes[first_node_idx] = first_node_num;

    // printf("first sampled node: genes[%d] = %d\n", first_node_idx, first_node_num);

    already_sampled_idx[0] = first_node_idx;
    already_sampled_node[0] = first_node_num;


    m_sample_genes_entropy_arr[0] += calculate_entropy(first_node_idx_arr, m_problem_size);
    /*** === end of sample first node  === ***/

    

    /*** === sample remaining nodes start === ***/





    // ----- sample all nodes except for first node ----- //
    // "sample_count = 1" because first node is already sampled
    for (int sample_count = 1; sample_count < m_problem_size; ++sample_count) {
        
        // --- initialize entropy_matrix (use "only one" already sampled node) --- //
        double **entropy_matrix = new double*[m_problem_size];
        for (int i = 0; i < m_problem_size; ++i) {
            entropy_matrix[i] = new double[m_problem_size]();
        }


        // kkk is number of nodes to reference by previos sampled nodes
        int kkk = 1;
        while (kkk <= m_previous_sampled_reference_count && (sample_count - kkk >= 0) ) {
            // accumulate entropy using [ two dimensional matrix ]
            //                         (problem size * problem size)
            //                                  (idx * node_num)
            for (int distance = 1; distance < m_problem_size; ++distance) {

                int tmp_idx = already_sampled_idx[sample_count - kkk] + distance;
                if (tmp_idx >= m_problem_size) {
                    tmp_idx -= m_problem_size;
                }

                for (int node = 0; node < m_problem_size; node++){
                    entropy_matrix[tmp_idx][node] += m_multiple_relation_matrix[already_sampled_node[sample_count - kkk]][distance][node];
                }
            }
        
            kkk ++;
        } // while end






        // --- clear already_sampled_node and already_sampled_idx in entropy_matrix to zero
        for (int rs = 0; rs < sample_count; rs++) { // rs: already sampled

            for (int i = 0; i < m_problem_size; ++i) {
                entropy_matrix[already_sampled_idx[rs]][i] = 0;
            }

            for (int i = 0; i < m_problem_size; ++i) {
                entropy_matrix[i][already_sampled_node[rs]] = 0;
            }
        }



        // sample new gene by lowest entropy idx
        double min_entropy = 100000;
        int sampled_idx = -1;

        for (int idx = 0; idx < m_problem_size; idx ++) {

            // check if idx is in already_sampled_idx ??
            bool is_in_already_sampled_idx = false;

            for (int i = 0; i < sample_count; ++i) {
                if (already_sampled_idx[i] == idx) {
                    is_in_already_sampled_idx = true;
                    break;
                }
            }


            // if idx is not in already_sampled_idx

            if (!is_in_already_sampled_idx) {
                double current_entropy = calculate_entropy(entropy_matrix[idx], m_problem_size);

                if (current_entropy == -INFINITY) {
                    PrintArray(entropy_matrix[idx], m_problem_size);
                }
                if (current_entropy < min_entropy ) {
                // if (current_entropy != 0 && current_entropy < min_entropy ) {
                    min_entropy = current_entropy;
                    sampled_idx = idx;
                }
            }
        }

        // printf("sample count: %d, min_entropy: %lf\n", sample_count, min_entropy);
        m_sample_genes_entropy_arr[sample_count] += min_entropy;
        // printf("sampled_idx: %d\n", sampled_idx);
        // PrintMatrix(entropy_matrix, m_problem_size);

        // PrintArray(entropy_matrix[sampled_idx], m_problem_size);
        genes[sampled_idx] = sample_from_array(entropy_matrix[sampled_idx], m_problem_size);
        // genes[sampled_idx] = sample_from_array_pow(entropy_matrix[sampled_idx], m_problem_size, 1);
        already_sampled_idx[sample_count] = sampled_idx;
        already_sampled_node[sample_count] = genes[sampled_idx];


        // if (m_already_sampled_chromosome_count == 0) {
        //     PrintArray(entropy_matrix[sampled_idx], m_problem_size);
        //     printf("sample count: %d, genes[%d] = %d\n", sample_count, sampled_idx,  genes[sampled_idx]);
        // }



        // --- initialize entropy_matrix (use "only one" already sampled node) --- //
        for(int i = 0; i < m_problem_size; ++i){
            delete[] entropy_matrix[i];
        }
        delete[] entropy_matrix;
        
    }
    // PrintMRMatrix(m_multiple_relation_matrix, m_problem_size);

    // printf("\n");


    /*** === end of sample remaining nodes === ***/

    // ----- print sample distance ----- //
    // for (int i = 0; i < m_problem_size-1; i ++) {
    //     printf("%d ", already_sampled_idx[i+1] - already_sampled_idx[i]);
    // }
    // printf("\n");






    // ----- for entropy experiment ----- //

    m_already_sampled_chromosome_count ++;
    // printf("m_already_sampled_chromosome_count: %d\n", m_already_sampled_chromosome_count);

    if (m_already_sampled_chromosome_count >= m_sample_size) {
        for (int sample_order = 0; sample_order < m_problem_size; sample_order ++) {

            // printf("sample order: %d (%d), min_entropy: %lf\n",\
                        sample_order, m_problem_size - sample_order,\
                        m_sample_genes_entropy_arr[sample_order] / m_sample_size);

            // printf("m_iteration: %d\n", m_iteration);
            m_entropy_change[sample_order][m_iteration] = m_sample_genes_entropy_arr[sample_order] / m_sample_size;
        }
        // printf("\n");
    }

    // ----- end of for entropy experiment ----- //








    // ----- release memory ----- //

    delete[] already_sampled_idx;
    delete[] already_sampled_node;

    // print sample result chromosome
    // int total = 0;
    // for(int i = 0; i < m_problem_size; ++i){
    //     printf("%d ", genes[i]);
    //     // total += genes[i];
    // }
    // printf("\n");
    // cout <<total<< endl;
}