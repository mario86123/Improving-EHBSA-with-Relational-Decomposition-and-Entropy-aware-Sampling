#include "FMSTModel.h"
#include "Tools.h"
#include <errno.h>
#include <fcntl.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <vector>
using std::cerr;
using std::cout;
using std::endl;

struct edge {
    int node_a;
    int node_b;
    double length;
    int already_sampled_node_count;

    int max_ele;
};

// edge * m_mst_edge_arr;

struct node {
    int node_num;
    double length;
};


/*
 * Class constructor.
 */
CFMSTModel::CFMSTModel(int problem_size, int sel_size, double b_ratio)
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

    // m_mst_edge_arr = new edge[m_problem_size - 1];
    // m_mst_node_arr = new node[m_problem_size];



    node_entropy_arr = new double[m_problem_size];
    edge_entropy_arr = new double*[m_problem_size];
    for (int i = 0; i < m_problem_size; ++i) {
        edge_entropy_arr[i] = new double[m_problem_size]();
    }

    m_sample_order = new int[m_problem_size]();

}









CFMSTModel::~CFMSTModel()
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



    delete[] node_entropy_arr;

    for(int i = 0; i < m_problem_size; ++i){
        delete[] edge_entropy_arr[i];
    }

    delete[] m_sample_order;
}


void del_e_from_list(vector<edge>& mst_edge_list, int last_already_sampled_gene) {
    for (edge & eee: mst_edge_list) {
        if (eee.node_a == last_already_sampled_gene) {
            
            if (eee.already_sampled_node_count != 0) {
                eee.already_sampled_node_count = 3;
            }

            else {
                // 1: node_a is sampled, 2: node_b is sampled
                eee.already_sampled_node_count = 1;
            }
        }
        else if (eee.node_b == last_already_sampled_gene) {
            
            if (eee.already_sampled_node_count != 0) {
                eee.already_sampled_node_count = 3;
            }

            else {
                // 1: node_a is sampled, 2: node_b is sampled
                eee.already_sampled_node_count = 2;
            }
        }
    }
}

void del_n_from_list(vector<node>& mst_node_list, int last_already_sampled_gene) {
    for (auto nit = mst_node_list.cbegin(); nit != mst_node_list.cend(); ) {
        if (nit->node_num == last_already_sampled_gene) {
            nit = mst_node_list.erase(nit);
            break;
        }
        else {
            ++nit;
        }
    }
}





/*
 * Virtual learning function.
 */
bool CFMSTModel::Learn(CPopulation * population, int size)
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
    


    // ----- build binary relation model ----- // 
    // selection ?????
    // selection ?????
    // selection ?????
    // printf("sample:\n");
    for(int k=0;k<size;k++) {
		m_samples[k] = population->m_individuals[k]->Genes();
        

        // ----- build edge matrix ----- // 
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

        // ----- build node matrix ----- // 
        // m_node_matrix[i][j] = node i position is idx j
        for (int pos = 0; pos < m_problem_size; ++pos) {
            m_node_matrix[m_samples[k][pos]][pos] ++;
        }


    }
    // ----- end of build binary relation model ----- // 

    // ----- print edge matrix ----- // 
    // PrintMSTEdgeMatrix(m_edge_matrix, m_problem_size);

    // ----- print node matrix ----- // 
    // printf("pos / node:");
    // PrintMatrix(m_node_matrix, m_problem_size);









    
    // calculate all node entropy (node num: n, time complexity: n^2)
    for (int node_num = 0; node_num < m_problem_size; ++node_num) {
        node_entropy_arr[node_num] = calculate_entropy(m_node_matrix[node_num], m_problem_size);
    }

    // calculate all edge entropy (edge num: n^2, time complexity: n^3)
    for (int i = 0; i < m_problem_size; ++i) {
        for (int j = 0; j < m_problem_size; ++j) {
            edge_entropy_arr[i][j] = calculate_entropy(m_edge_matrix[i][j], m_problem_size);
        }
    }


    // ----- normalize MST edge entropy and node entropy ----- //
        // edge: n - 1, node: n
    double max_entropy_arr[m_problem_size];
    for (int i = 0; i < m_problem_size; i++) {
        max_entropy_arr[i] = 1.0;
    }

    double node_max_entropy, edge_max_entropy;
    node_max_entropy = calculate_entropy(max_entropy_arr, m_problem_size);
    edge_max_entropy = calculate_entropy(max_entropy_arr, m_problem_size - 1);


    vector<node> mst_node_list;
    // for all node, divide another constant
    for (int node_num = 0; node_num < m_problem_size; ++node_num) {
        node tmp;
        tmp.node_num = node_num;
        tmp.length = node_entropy_arr[node_num] / node_max_entropy;
        mst_node_list.push_back(tmp);
    }
    // ----- finish normalizing MST edge entropy and node entropy ----- //

    // --- sort mst_node
    sort( mst_node_list.begin(),mst_node_list.end(), [&](node i,node j) { return i.length < j.length; } );

    // cout << "mst_node_entropy:" << endl;
    // int k = 0;
    // for (node nnn: mst_node_list) {
    //     k++;
    //     cout << k << "th node num: " << nnn.node_num << ",  normalized entopy: " << nnn.length << endl;
    // }

    // for (edge i: all_edge) {
    //     cout << i.node_a << ' ' << i.node_b << endl;
    // }



    // === find out MST edge part === //
    // --- select MST edge, delete non MST edge (build MST) --- //
    int total_edge_num = m_problem_size * (m_problem_size - 1) / 2;
    vector<edge> all_edge(total_edge_num);
    int count = 0;

    // build edge data
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
            

            all_edge[count].node_a = ref_node_num;
            all_edge[count].node_b = other_node_num;
            all_edge[count].length = edge_entropy_arr[ref_node_num][other_node_num];

            all_edge[count].max_ele = max_element_distance;

            count ++;
        }
    }

    // sort remaining ell x (ell-1) / 2 edges
    sort( all_edge.begin(),all_edge.end(), [&](edge i,edge j) { return i.length < j.length; } );
    // sort( all_edge.begin( ), all_edge.end( ), []( const edge& a, const edge&b ){return i.length < j.length;} );
    

    
    // print sorted edge array
    // for (edge i: all_edge) {
    //     cout << i.node_a << ' ' << i.node_b << endl;
    // }


    // ----- build MST using uniion find operator ----- //
    int mst_edge_count = 0;
    element node_set[m_problem_size];

    for(int i = 0; i < m_problem_size; ++i) {
        node_set[i].parent = i;
        node_set[i].rank = 0;
    }

    vector<edge> mst_edge_list;
    // cout << "\nset:" << endl;
    // while MST is not ready (for each edge in all edge)
    for (edge e: all_edge) {

        // if add that edge (AB) does not forming a cycle ( find-set(A) != find-set(B) )
        if (find_set(e.node_a, node_set) != find_set(e.node_b, node_set)) {
            
            // add that edge to graph
            union_set(e.node_a, e.node_b, node_set);

            e.length = e.length / edge_max_entropy;
            e.already_sampled_node_count = 0;
            mst_edge_list.push_back(e);
            mst_edge_count ++;

            // cout << e.node_a << ' ' << e.node_b << endl;
            if (mst_edge_count >= m_problem_size - 1) break;
        }
    }
    // ----- finish building MST ----- //
    // ----- MST information is in ""mst_edge_list"" ----- //

    // cout << "mst_edge_entropy:" << endl;
    // k = 0;
    // for (edge eee: mst_edge_list) {
    //     k++;
    //     cout << k << "th edge num: " << eee.node_a << ", " << eee.node_b <<  ",  max element: " << eee.max_ele <<  ",  normalized entopy: " << eee.length << endl;
    // }
    // cout << endl;



    // need some data structure to record sample order
    // ----- sort all nodes and MST edges by entropy to determine sample order ----- //
    
    // sample order: 
    // -------------------------------------------------------------------------------------------
    // | node: node_num         (range: 0 ~ problem_size)                                        |
    // | edge: not_sampled_yet_node_num * problem_size + already_sampled_node_num + problem size |
    // -------------------------------------------------------------------------------------------

        // first one must be node (lowest entropy node)
        
    // int m_sample_order[m_problem_size];    
    int already_sampled_gene_count = 0;

    // delete edge from m_mst_edge_arr if necessary (all two nodes are sampled) (for loop) edge --> list ??

        // add first node edges in stack and compare with another nodes entropy, until sample order is determined
            // mst_node_list: list
            // mst_edge_list: list
            // edge attribute: already sampled node count (0 ~ 2)
    while (already_sampled_gene_count < m_problem_size) {
        // cout << "\nsample order: ";
        // for (int i = 0; i < already_sampled_gene_count; i ++) {
        //     cout << m_sample_order[i] << " ";
        // }
        // cout << endl;
        // --- if is the first sampled gene
        if (already_sampled_gene_count == 0) {
            m_sample_order[already_sampled_gene_count] = mst_node_list.front().node_num;
            mst_node_list.erase (mst_node_list.begin());

            // delete edge from m_mst_edge_arr if necessary (all two nodes are sampled) (for loop) edge --> list ??
            int count = 0;
            for (edge & e: mst_edge_list) {
                // cout << "e.node_a: " << e.node_a << ", e.node_b: " << e.node_b << endl;

                // 1: node_a is sampled, 2: node_b is sampled                
                if (e.node_a == m_sample_order[already_sampled_gene_count]) {
                    e.already_sampled_node_count = 1;
                    count ++;
                }

                else if (e.node_b == m_sample_order[already_sampled_gene_count]) {
                    e.already_sampled_node_count = 2;
                    count ++;
                }

                // already discover two edge in list (speed up)
                if (count >= 2) break;
            }
            // cout << "mst_edge_entropy:" << endl;
            // k = 0;
            // for (edge eee: mst_edge_list) {
            //     k++;
            //     cout << k << "th edge num: " << eee.node_a << ", " << eee.node_b <<  ",  already_sampled_node_count: " << eee.already_sampled_node_count << endl;
            // }
        }
        // --- end of if is the first sampled gene
        

        //  ---  if not the first sampled gene
        else {

        // find edge min entropy and eligible (exactly one node is sampled) element
            double min_edge_length = 100.0;

            // int count = 0;
            bool has_available_edge = false;

            // for (edge & e: mst_edge_list) {
            for (vector<edge>::iterator eit = mst_edge_list.begin(); eit != mst_edge_list.end();) {
                // cout << "eit->node_a: " << eit->node_a << ", eit->node_b: " << eit->node_b << ", eit->count: " << eit->already_sampled_node_count << endl;
                // cout << "count: " << count << ", e.already_sampled_node_count: " << e.already_sampled_node_count << endl;
                
                // 3: two node are both sampled
                if (eit->already_sampled_node_count >= 3) {
                    // cannot use count to remove vector element
                    // because of multiple already_sampled_node_count == 3
                    // cout << "erased edge: " << endl;
                    // cout << "eit->node_a: " << eit->node_a << ", eit->node_b: " << eit->node_b << ", eit->count: " << eit->already_sampled_node_count << endl<< endl;

                    mst_edge_list.erase(eit);
                    // cout << "original continue\n";
                    // continue;
                } 

                // 1: node_a is sampled, 2: node_b is sampled
                else if (eit->already_sampled_node_count == 1 || eit->already_sampled_node_count == 2) {
                    
                    has_available_edge = true;
                    
                    // ========== pick ""node"" ========== /

                    // if entropy(edge_min) > entropy(node_min), pick ""node""
                    if (eit->length > mst_node_list.front().length) {
                        // cout << "pick node !! "<< endl;
                        // cout << "eit->node_a: " << eit->node_a<< ", eit->node_b: " << eit->node_b << endl;
                        // cout << "eit->length: "<< eit->length << ", mst_node_list.front().length: " << mst_node_list.front().length << endl;
                        // cout << "mst_node_list.front().node_num: " << mst_node_list.front().node_num << endl << endl;
                        m_sample_order[already_sampled_gene_count] = mst_node_list.front().node_num;
                        mst_node_list.erase (mst_node_list.begin());


                        // delete edge from m_mst_edge_arr if necessary (all two nodes are sampled) (for loop) edge --> list ??
                        del_e_from_list(mst_edge_list, m_sample_order[already_sampled_gene_count]);
                        // end of for loop of delete edge

                    }
                    
                    // ========== pick ""edge"" ========== /

                    // else if entropy(edge_min) < entropy(node_min), pick ""edge""
                    else {
                        // cout << "pick edge !! "<< endl;
                        // cout << "eit->node_a: " << eit->node_a<< ", eit->node_b: " << eit->node_b << endl;
                        // cout << "eit->length: "<< eit->length << ", mst_node_list.front().length: " << mst_node_list.front().length << endl << endl;


                        // check which node is already sampled, which is not sampled yet
                            // sample order: 
                            // not_sampled_yet_node_num * problem_size + already_sampled_node_num + problem size
                        
                    // node_a is already sampled
                        if (eit->already_sampled_node_count == 1) { 
                            m_sample_order[already_sampled_gene_count] = eit->node_b * m_problem_size + eit->node_a + m_problem_size;


                            // delete "node and edge" from m_mst_node_arr and m_mst_edge_arr (node --> list ??)
                            del_e_from_list(mst_edge_list, eit->node_b); // node_a is already sampled
                            del_n_from_list(mst_node_list, eit->node_b); // node_a is already sampled
                            
                        }
                    // end of node_a is already sampled
                        
                    // node_b is already sampled
                        else {
                            // cout << "node a: " << eit->node_a << ", node b: " << eit->node_b << endl;
                            m_sample_order[already_sampled_gene_count] = eit->node_a * m_problem_size + eit->node_b + m_problem_size;

                            // delete "node and edge" from m_mst_node_arr and m_mst_edge_arr (node --> list ??)
                            del_e_from_list(mst_edge_list, eit->node_a); // node_b is already sampled
                            del_n_from_list(mst_node_list, eit->node_a); // node_b is already sampled
                        }
                    // end of node_b is already sampled
                    }
                    
                    // ========== end of pick ""edge"" ========== /
                    break;                
                }
                // end of if 1: node_a is sampled, 2: node_b is sampled
                
                else {
                    ++eit;
                }
            }
            // end of iterate over mst_sdge_list
            // cout << "has_available_edge: " << has_available_edge << endl;
            
            if (has_available_edge == false) {
                // cout << "pick node !! "<< endl;
                m_sample_order[already_sampled_gene_count] = mst_node_list.front().node_num;
                mst_node_list.erase (mst_node_list.begin());


                // delete edge from m_mst_edge_arr if necessary (all two nodes are sampled) (for loop) edge --> list ??
                del_e_from_list(mst_edge_list, m_sample_order[already_sampled_gene_count]);
                // end of for loop of delete edge
            }

        }
        // end of if not the first sampled gene

        already_sampled_gene_count ++;
    }
    // end of while loop

    // cout << "sample order: " << endl;
    // for (int i = 0; i < m_problem_size; i ++) {
    //     if (m_sample_order[i] < m_problem_size) {
    //         cout << "node: " << m_sample_order[i] << endl;
    //     }
    //     else {
    //         cout << "edge: " << (m_sample_order[i] - m_problem_size) / m_problem_size << " " << (m_sample_order[i] - m_problem_size) % m_problem_size << endl;
    //     }
    // }
    // cout << endl;
    
    
    // --- MST finished --- //
    return true;
}








/*
 * Virtual sampling function.
 */
void CFMSTModel::Sample(int * genes)
{
    int already_used_idx[m_problem_size];

    double sample_arr[m_problem_size];
    int result_idx, not_sampled_yet_node_num;

    for (int already_sampled_gene_count = 0; already_sampled_gene_count < m_problem_size; already_sampled_gene_count++) {

        // check sample order is node or edge

            // m_node_matrix[i][j] = j position is node i
        
        // if it is "node"
        if (m_sample_order[already_sampled_gene_count] < m_problem_size) {
            
            not_sampled_yet_node_num = m_sample_order[already_sampled_gene_count];

            // sample using m_node_matrix
            for (int idx = 0; idx < m_problem_size; idx++) {
                sample_arr[idx] = m_node_matrix[m_sample_order[already_sampled_gene_count]][idx];
            }   
        }

        // else if it is "edge"
        else {

                // m_edge_matrix[i][j][1] = node i is before node j for 1 position

            // translate sample order to "already sampled node" and "need to sample node"
                
                // sample order: 
                // -------------------------------------------------------------------------------------------
                // | node: node_num         (range: 0 ~ problem_size)                                        |
                // | edge: not_sampled_yet_node_num * problem_size + already_sampled_node_num + problem size |
                // -------------------------------------------------------------------------------------------
            
            int already_sampled_node_num = (m_sample_order[already_sampled_gene_count] - m_problem_size) % m_problem_size;
            not_sampled_yet_node_num = (m_sample_order[already_sampled_gene_count] - m_problem_size) / m_problem_size;

            
            int already_sampled_node_num_idx;
            // find already_sampled_node_num idx
                // (can use memory to reduce time)
            for (int i = 0; i < already_sampled_gene_count; i ++) {
                if (genes[already_used_idx[i]] == already_sampled_node_num) {
                    already_sampled_node_num_idx = already_used_idx[i];
                    break;
                }
            }


            // sample using m_edge_matrix

                // change sample arr distance to idx

            for (int dis = 0; dis < m_problem_size; dis++) {
                
                int not_sample_yet_idx = already_sampled_node_num_idx + dis;
                if (not_sample_yet_idx >= m_problem_size) {
                    not_sample_yet_idx -= m_problem_size;
                }

                sample_arr[not_sample_yet_idx] = m_edge_matrix[already_sampled_node_num][not_sampled_yet_node_num][dis];
            }   

        }

            // clean already sampled index to zero
        for (int kkk = 0; kkk < already_sampled_gene_count; kkk++) {
            sample_arr[already_used_idx[kkk]] = 0.0;
        }


        result_idx = sample_from_array(sample_arr, m_problem_size);
        genes[result_idx] = not_sampled_yet_node_num;

        // record already sampled "index"
        already_used_idx[already_sampled_gene_count] = result_idx;


    } // end of sample for loop

}