// MST model with multiple edge (MSTME)
#include "MSTME.h"
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
CMSTME::CMSTME(int problem_size, int sel_size, double b_ratio)
{
	m_problem_size=problem_size;
    m_sample_size=sel_size;
    m_samples= new int*[m_sample_size];


    m_edge_matrix = new double**[m_problem_size];
    for (int i = 0; i < m_problem_size; ++i) {
        m_edge_matrix[i] = new double*[m_problem_size];
        // m_edge_matrix[i][j][1] = node i is before node j for 1 position
        for (int j = 0; j < m_problem_size; ++j) {
            m_edge_matrix[i][j] = new double[m_problem_size]();
        }
    }


    // m_node_matrix[i][j] = node j position is node i
    m_node_matrix = new double*[m_problem_size];
    for (int i = 0; i < m_problem_size; ++i) {
        m_node_matrix[i] = new double[m_problem_size]();
    }



    // epsilon = b * ell / pop_size
    m_epsilon = double(m_problem_size * b_ratio) / (m_sample_size); // pop_size = sample_size

    m_mst_edge_arr = new edge[m_problem_size - 1];


}









CMSTME::~CMSTME()
{
    // PrintMRMatrix(m_multiple_relation_matrix, m_problem_size);

    delete [] m_samples;

    // Release edge matrix Memory
    for(int i = 0; i < m_problem_size; ++i){
        for(int j = 0; j < m_problem_size; ++j){
            delete[] m_edge_matrix[i][j];
        }
        delete[] m_edge_matrix[i];
    }
    delete[] m_edge_matrix;
    


    // Release node matrix Memory
    for(int i = 0; i < m_problem_size; ++i){
        delete[] m_node_matrix[i];
    }
    delete[] m_node_matrix;

    delete[] m_mst_edge_arr;
}







/*
 * Virtual learning function.
 */
bool CMSTME::Learn(CPopulation * population, int size)
{   

    // clean matrix
    // m_edge_matrix[i][j][1] = node i is before node j for 1 position
    for(int ref_node = 0; ref_node < m_problem_size; ++ref_node){
        for(int sample_node = 0; sample_node < m_problem_size; ++sample_node){
            for(int dist = 0; dist < m_problem_size; ++dist){
                
                if (ref_node != sample_node && dist != 0) {
                    m_edge_matrix[ref_node][sample_node][dist] = m_epsilon;
                }
                else {
                    m_edge_matrix[ref_node][sample_node][dist] = 0.0;
                }
            }
        }
    }

    for(int i = 0; i < m_problem_size; ++i) {
        for(int j = 0; j < m_problem_size; ++j) {
            m_node_matrix[i][j] = m_epsilon;
        }    
    }    
    



    // selection ?????
    // selection ?????
    // selection ?????
    // printf("sample:\n");
    for(int k=0;k<size;k++) {
		m_samples[k] = population->m_individuals[k]->Genes();
        

        // build edge matrix
        for (int reference_node_idx = 0; reference_node_idx < m_problem_size; ++reference_node_idx) {
            for (int distance = 0; distance < m_problem_size; ++distance) { // next node
                int next_node_idx = reference_node_idx + distance;
                if (next_node_idx >= m_problem_size) {
                    next_node_idx -= m_problem_size;
                }
                // m_edge_matrix[i][j][1] = node i is before node j for 1 position
                m_edge_matrix[ m_samples[k][reference_node_idx] ][ m_samples[k][next_node_idx] ][distance] ++;
            }
        }

        // build node matrix
        // m_node_matrix[i][j] = node i position is idx j
        for (int pos = 0; pos < m_problem_size; ++pos) {
            m_node_matrix[m_samples[k][pos]][pos] ++;
        }


    }
    // PrintMSTEdgeMatrix(m_edge_matrix, m_problem_size);

    // printf("pos / node:");
    // PrintMatrix(m_node_matrix, m_problem_size);


    // === build MST part === //
    
    int total_edge_num = m_problem_size * (m_problem_size - 1) / 2;

    // find every node and node biggest element edge
    vector<edge> all_edge(total_edge_num);
    int count = 0;

    for (int ref_node_num = 0; ref_node_num < m_problem_size; ++ref_node_num) {
        for (int other_node_num = ref_node_num + 1; other_node_num < m_problem_size; ++other_node_num) {

            double max_element = 0;
            int max_element_distance = 0;

            for (int distance = 1; distance < m_problem_size; ++distance) {

                if (m_edge_matrix[ref_node_num][other_node_num][distance] > max_element) {
                    max_element = m_edge_matrix[ref_node_num][other_node_num][distance];
                    max_element_distance = distance;
                }

            }
            
            // printf("%d: max element edge: (%d, %d), distance: %d, max element: %lf\n", count, ref_node_num, other_node_num, max_element_distance, max_element);
            // all_edge_element_arr[count] = max_element;
            all_edge[count].node_a = ref_node_num;
            all_edge[count].node_b = other_node_num;
            all_edge[count].length = max_element;

            count ++;
        }
    }
    // remain ell x (ell-1) / 2 edges




    // --- build MST --- //

    // sort remaining ell x (ell-1) / 2 edges
    sort( all_edge.begin(),all_edge.end(), [&](edge i,edge j) { return i.length>j.length; } );
    
    // print sorted edge array
    // for (edge i: all_edge) {
    //     cout << i.node_a << ' ' << i.node_b << endl;
    // }


    int mst_edge_count = 0;
    element node_set[m_problem_size];

    for(int i = 0; i < m_problem_size; ++i) {
        node_set[i].parent = i;
        node_set[i].rank = 0;
    }

    // cout << "\nset:" << endl;
    // while MST is not ready (for each edge in all edge)
    for (edge e: all_edge) {

        // if add that edge (AB) does not forming a cycle ( find-set(A) != find-set(B) )
        if (find_set(e.node_a, node_set) != find_set(e.node_b, node_set)) {
            
            // add that edge to graph
            union_set(e.node_a, e.node_b, node_set);
            m_mst_edge_arr[mst_edge_count] = e;
            mst_edge_count ++;

            // cout << e.node_a << ' ' << e.node_b << endl;
            if (mst_edge_count >= m_problem_size - 1) break;
        }
    }
    // --- MST finished --- //
    return true;
}








/*
 * Virtual sampling function.
 */
void CMSTME::Sample(int * genes)
{
    // two stack

    // 1. MST_edge_not_used_yet (max_length = ell - 1 ) 
    //       (m_mst_edge_arr)



    int already_sampled_idx[m_problem_size];
    int already_sampled_node[m_problem_size];
    bool is_sampled_node[m_problem_size];
    bool is_sampled_idx[m_problem_size];

    for(int i = 0; i < m_problem_size; ++i){
        already_sampled_idx[i] = -1;
        already_sampled_node[i] = -1;
        is_sampled_node[i] = false;
        is_sampled_idx[i] = false;
    }



    int sample_count = 0;
    // add node histogram to stack

    // edge edge_stack[m_problem_size - 1];
    // int edge_stack_length = 0;

    // loop until already sampled all genes (while true ?? )
    while (sample_count < m_problem_size) {

        // find min entropy element in stack
        
            
        // --- find min entropy element in node stack --- //
        int min_entropy_node_num = -1;
        
        double min_entropy_node_arr[m_problem_size];
        double tmp_entropy_node_arr[m_problem_size];

        double node_min_entropy = 100;
        double tmp_entropy;

        // for all not sampled yet node
        for(int node_num = 0; node_num < m_problem_size; ++node_num) {

            // if that node is not sampled yet
            if ( !is_sampled_node[node_num] ) {
                
                // clear already sampled idx to 0
                for(int idx = 0; idx < m_problem_size; ++idx) {
                    
                    if (is_sampled_idx[idx]) { // that idx is already sampled
                        tmp_entropy_node_arr[idx] = 0;
                    } 
                    else { // that idx is not sampled yet
                        tmp_entropy_node_arr[idx] = m_node_matrix[node_num][idx];
                    }
                }



                // calculate entropy
                tmp_entropy = calculate_entropy(tmp_entropy_node_arr, m_problem_size);



                // if that entropy smaller than node_min_entropy, update information
                if (node_min_entropy > tmp_entropy) {
                    node_min_entropy = tmp_entropy;
                    min_entropy_node_num = node_num;

                    for(int idx = 0; idx < m_problem_size; ++idx) {
                        min_entropy_node_arr[idx] = tmp_entropy_node_arr[idx];
                    }
                }

            }
        }
        // --- end of find min entropy element in node stack --- //



        // --- find min entropy element by edges --- //

        int min_entropy_edge_num = -1;
        
        double min_entropy_edge_arr[m_problem_size];
        double tmp_entropy_edge_arr[m_problem_size];

        double edge_min_entropy = 100;


        // candidate node: not sampled yet, and also have edge with already sampled node
        int candidate_node_list[m_problem_size];
        int candidate_node_count = 0;


        // candidate_node_list_arr[i][j]: the ith node in candidate_node_list's at j position's probability
        double **candidate_node_list_arr = new double*[m_problem_size];
        for (int i = 0; i < m_problem_size; ++i) {
            candidate_node_list_arr[i] = new double[m_problem_size]();
        }




        // --- find min entropy element in edge stack --- //
        for(int edge_num = 0; edge_num < m_problem_size - 1; ++edge_num) {
            
            // if can use that edge (only one node is sampled)
            if (is_sampled_node[m_mst_edge_arr[edge_num].node_a] != is_sampled_node[m_mst_edge_arr[edge_num].node_b]) {
                
                
                // === assign ref_node, need_to_sampled_node part === //
                int ref_node, need_to_sampled_node;

                // if node a is sampled, then ref_node is node a, sampled node is node b
                if (is_sampled_node[m_mst_edge_arr[edge_num].node_a]) {
                    ref_node = m_mst_edge_arr[edge_num].node_a;
                    need_to_sampled_node = m_mst_edge_arr[edge_num].node_b;
                }
                // else if node b is sampled, then ref_node is node b, sampled node is node b
                else {
                    ref_node = m_mst_edge_arr[edge_num].node_b;
                    need_to_sampled_node = m_mst_edge_arr[edge_num].node_a;
                }
                // === end of assign ref_node, need_to_sampled_node part === //


                int candidate_node_idx;
                bool already_in_candidate_node_list = false;
                for (int k = 0; k < candidate_node_count; k++) {
                    if (candidate_node_list[k] == need_to_sampled_node) {
                        already_in_candidate_node_list = true;
                        candidate_node_idx = k;
                        break;
                    }
                }

                // the need_to_sampled_node is not in candidate_node_list
                if (!already_in_candidate_node_list) { 
                    candidate_node_list[candidate_node_count] = need_to_sampled_node;
                    candidate_node_idx = candidate_node_count;
                    candidate_node_count ++;
                }


                // find ref node index
                int ref_node_idx;
                for (int s = 0; s < sample_count; s++) {
                    if (already_sampled_node[s] == ref_node) {
                        ref_node_idx = already_sampled_idx[s];
                        break;
                    }
                }
                

                // -- assign tmp arr value
                // m_edge_matrix[i][j][1] = node i is before node j for 1 position
                for(int dis = 0; dis < m_problem_size; ++dis) {
                    
                    int tmp_idx = (ref_node_idx + dis) % m_problem_size;

                    // candidate_node_list_arr[i][j]: the ith node in candidate_node_list's at j position's probability

                    if ( !is_sampled_idx[tmp_idx]) { // if that idx is not sampled yet
                        candidate_node_list_arr[candidate_node_idx][tmp_idx] += m_edge_matrix[ref_node][need_to_sampled_node][dis];
                    }
                    else { // if that idx is already sampled
                        candidate_node_list_arr[candidate_node_idx][tmp_idx] += 0;
                    }

                }
                // -- end of assign tmp arr value
            }
            // end of if can use that edge (only one node is sampled)
        }

        for (int cn_idx = 0; cn_idx < candidate_node_count; cn_idx++) {

            double tmp_entropy = calculate_entropy(candidate_node_list_arr[cn_idx], m_problem_size);

            // == if tmp entropy smaller than edge_min_entropy, update information
            if (edge_min_entropy > tmp_entropy) {
                edge_min_entropy = tmp_entropy;
                min_entropy_edge_num = candidate_node_list[cn_idx];

                for(int idx = 0; idx < m_problem_size; ++idx) {
                    min_entropy_edge_arr[idx] = candidate_node_list_arr[cn_idx][idx];
                }
            }
        }
        // == end of update information
        // --- end of find min entropy element in edge stack --- //

        for(int i = 0; i < m_problem_size; ++i){
            delete[] candidate_node_list_arr[i];
        }
        delete[] candidate_node_list_arr;



        int sampled_idx = -1;
        // --- min (node_min_entropy, edge_min_entropy) --- //
        if (node_min_entropy < edge_min_entropy) {
            sampled_idx = sample_from_array(min_entropy_node_arr, m_problem_size);
            genes[sampled_idx] = min_entropy_node_num;
        }
        else {
            sampled_idx = sample_from_array(min_entropy_edge_arr, m_problem_size);
            genes[sampled_idx] = min_entropy_edge_num;
            // if (sampled_idx == -1) {
            //     cout<< "edge" << endl;
            //     PrintArray(min_entropy_edge_arr, m_problem_size);
            //     cout << sample_count << endl;
            // }
        }
        // sample using that element


        // cout << "genes[" << sampled_idx << "]: "<< genes[sampled_idx] << endl;


        // set sampled node and idx to already sampled (already_sampled_idx and already_sampled_node)
        already_sampled_idx[sample_count] = sampled_idx;
        already_sampled_node[sample_count] = genes[sampled_idx];
        is_sampled_idx[sampled_idx] = true;
        is_sampled_node[genes[sampled_idx]] = true;


        // delete related element in stack (no need, because write in above)
        // add related element to edge stack (no need, because write in above)
        // edge_stack_length++; (no need, because write in above)

        sample_count ++;
        // if (sample count >= problem size) , then break the loop
    }
    // cout<< endl;

}



// -- union find functions part -- //



int CMSTME::find_set(int x, element * arr) {
    
    if (x != arr[x].parent) {
        arr[x].parent = find_set(arr[x].parent, arr);
    }
    return arr[x].parent;
}

void CMSTME::link(int x, int y, element * arr) {
    if (arr[x].rank > arr[y].rank) {
        arr[y].parent = x;
    }
    else {
        arr[x].parent = y;
        if (arr[x].rank == arr[y].rank) {
            arr[y].rank += 1;
        }
    }
}

void CMSTME::union_set(int x, int y, element * arr) {
    link(find_set(x, arr), find_set(y, arr), arr);
}

// -- union find functions part end -- //