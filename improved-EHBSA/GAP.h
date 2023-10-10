// GA project fitness function

#ifndef _GAP_H__
#define _GAP_H__

#include "PBP.h"
#include "Tools.h"
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string.h>
#include <stdio.h>

using std::ifstream;
using std::ofstream;
using std::istream;
using std::ostream;
using namespace std;
using std::cerr;
using std::cout;
using std::endl;
using std::stringstream;
using std::string;

class GAP : public PBP
{
	
public:

    int m_size; // problem_size
    
    // 0: permutation one max
    // 1: order problem
    int problem_type;

	GAP();
	
    virtual ~GAP();
	
	/*
	 * Read QAP file.
	 */
	int Read(string filename);
	
	/*
	 * This function evaluates the individuals for the QAP problem.
	 */
	double Evaluate(int * genes);
    
    /*
     * This function evaluates the inverted solution of the given individual for the QAP problem.
     */
    double EvaluateInv(int * genes);
	
    /*
     * Returns the size of the problem.
     */
    int GetProblemSize();

private:
	
};
#endif