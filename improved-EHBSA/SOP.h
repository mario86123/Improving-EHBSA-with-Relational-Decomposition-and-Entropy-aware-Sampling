/*
 *  SOP.h
 *  RankingEDAsCEC
 *
 *  Created by Josu Ceberio Uribe on 7/11/13.
 *  Copyright 2013 University of the Basque Country. All rights reserved.
 *
 */

#ifndef _SOP_H__
#define _SOP_H__

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
using std::ifstream;
using std::stringstream;
using std::string;

class SOP : public PBP
{
	
public:
	
    /*
     * Matrix of distances between the cities.
     */
	double ** m_distance_matrix;

	// vector<int[2]> m_constraint_pair;
	int ** m_constraint_pair;
	int num_of_constraint_pair;

	/*
	 * The number of cities.
	 */
	int m_size;
	
	/*
     * The constructor.
     */
	SOP();
	
    /*
     * The destructor.
     */
    virtual ~SOP();
	
	double CalculateGEODistance(double latitudeX, double latitudeY, double longitudeX, double longitudeY);

	/*
	 * Read SOP instance file that belongs to the SOPLIB library.
	 */
	int Read2(string filename);

    /*
	 * Read SOP instance file.
	 */
	int Read(string filename);
    
	/*
	 * This function evaluates the fitness of the solution for the SOP problem.
	 */
	double Evaluate(int * genes);

    /*
	 * This function evaluates the inverted solution of the given solution for the SOP problem.
	 */
	double EvaluateInv(int * genes);
    
    /*
     * Returns the size of the problem.
     */
    int GetProblemSize();
private:
	
};
#endif
