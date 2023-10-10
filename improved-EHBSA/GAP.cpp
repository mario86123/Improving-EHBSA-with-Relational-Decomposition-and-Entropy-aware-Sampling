// GA project fitness function

#include "GAP.h"

/*
 *Class constructor.
 */
GAP::GAP()
{
	
}

/*
 * Class destructor.
 */
GAP::~GAP()
{

}



int GAP::Read(string filename)
{
    // if filename c++ string to int (do not need problem instance)
    // problem_type = ...
    if (filename.compare(0, 3, "pom") == 0) {
        // printf("pom\n");
        problem_type = 0;
        m_size = atoi(filename.substr(3, 4).c_str());    
        // printf("ell: %d\n", m_size);
    }

    // else if file name order problem
    // problem_type = ...
    else if (filename.compare(0, 5, "order") == 0) {
        // printf("order\n");
        problem_type = 1;
        m_size = atoi(filename.substr(5, 6).c_str());    
        // printf("ell: %d\n", m_size);
    }
    else if (filename.compare(0, 3, "aaa") == 0) {
        // printf("order\n");
        problem_type = 2;
        m_size = 10;    
        // printf("ell: %d\n", m_size);
    }
    else if (filename.compare(0, 4, "same") == 0) {
        // printf("order\n");
        problem_type = 3;
        m_size = atoi(filename.substr(4, 5).c_str());  
        // printf("ell: %d\n", m_size);
    }
	return m_size;
}

/*
 * This function evaluates the individuals for the QAP problem.
 */
double GAP::Evaluate(int * genes)
{
	double fitness=0;
    
    // permutation one max
    if (problem_type == 0) { 
        for (int i = 0; i < m_size; i++) {
            if (genes[i] == i) {
                fitness += 1;
            }
        }
    }

    // order problem
    else if (problem_type == 1) {
        for (int i = 0; i < m_size; i++) {
            
            int i_add_three = i+3;
            int i_minus_three = i-3;

            if (i_add_three >= m_size) {
                i_add_three -= m_size;
            }

            if (i_minus_three < 0) {
                i_minus_three += m_size;
            }
            if (genes[i_add_three] == genes[i] + 3 || genes[i_minus_three] == genes[i] + 3) {
                fitness += 1;
            }
        }
    }
    // aaa
    else if (problem_type == 2) {
        int pos[10] = {0};

        if (genes[2] == 2) {
            fitness += 1;
        }
        if (genes[9] == 9) {
            fitness += 1;
        }

        for (int i = 0; i < m_size; i++) {
            pos[genes[i]] = i;
        }

        if ((pos[0]+1)%10 == pos[1]) {
            fitness += 1;
        }
        if ((pos[3]+1)%10 == pos[4]) {
            fitness += 1;
        }
        if ((pos[4]+1)%10 == pos[5]) {
            fitness += 1;
        }
        if ((pos[5]+1)%10 == pos[6]) {
            fitness += 1;
        }
        if ((pos[7]+1)%10 == pos[8]) {
            fitness += 1;
        }
        
    }
    // same
    else if (problem_type == 3) {
        fitness = 0;
    }
	return fitness;
}

/*
 * This function evaluates the inverted solution of the given individual for the QAP problem.
 */
double GAP::EvaluateInv(int * genes)
{
    return 0;
}

/*
 * Returns the size of the problem.
 */
int GAP::GetProblemSize()
{
    return m_size;
}