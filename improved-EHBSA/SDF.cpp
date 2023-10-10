/*
 *  QAP.cpp
 *  RankingEDAsCEC
 *
 *  Created by Josu Ceberio Uribe on 7/11/13.
 *  Copyright 2013 University of the Basque Country. All rights reserved.
 *
 */

#include "SDF.h"

/*
 *Class constructor.
 */
SDF::SDF()
{
	
}

/*
 * Class destructor.
 */
SDF::~SDF()
{
    for(int i = 0; i < m_size; ++i){
        delete[] m_constraint_matrix[i];
    }
    delete[] m_constraint_matrix;
}



int SDF::Read(string filename)
{
    char line[2048];
    ifstream indata;
    
    indata.open(filename,ios::in);

    int read_line_num=0;

    m_size = 0;
	while (read_line_num < m_size + 1)
	{
        indata.getline(line, 2048);
        stringstream ss;
		string sline;
		ss << line;
		ss >> sline;

        if (read_line_num == 0)
		{   
            m_size = atoi(sline.c_str());
            m_constraint_matrix = new int*[m_size];
            
            for (int i =0 ; i < m_size; i ++) {
                m_constraint_matrix[i] = new int[3]; // front, back, distance
            }

        } 

        else {
            char * pch;

            //front
			pch = strtok (line," ");
			m_constraint_matrix[read_line_num-1][0]=atoi(pch);
			
            //back
            pch = strtok (NULL, " ,.");
            m_constraint_matrix[read_line_num-1][1]=atoi(pch);

            //distance
            pch = strtok (NULL, " ,.");
            m_constraint_matrix[read_line_num-1][2]=atoi(pch);

        }
        read_line_num++;
    }
    indata.close();

	return m_size;
}

/*
 * This function evaluates the individuals for the QAP problem.
 */
double SDF::Evaluate(int * genes)
{
	double fitness=0;
    
    // absolute postition

    // 1965374028, max fitness = 19.5
    // if (genes[4] == 3 && genes[8] == 2 && genes[2] == 6) {
    //     fitness += 10;
    // }
    // if (genes[2] == 6 && genes[5] == 7) {
    //     fitness += 5;
    // }
    // if (genes[0] == 1 && genes[9] == 8) {
    //     fitness += 3;
    // }
    // if (genes[1] == 9 && genes[3] == 5) {
    //     fitness += 1;
    // }
    // if (genes[6] == 4 && genes[7] == 0) {
    //     fitness += 0.5;
    // }


    // relative postition
    // problem size (genes) loop
    for (int i = 0; i < m_size; i++) {
        
        // constraint loop
        for (int j = 0; j < m_size; j++) {

            // 5 is after 3 for 4 pos
            if (genes[i] == m_constraint_matrix[j][0]) {

                int distance = m_constraint_matrix[j][2];

                if (i + distance >= m_size) {
                    distance -= m_size;
                }

                if (*(&(genes[i]) + distance) == m_constraint_matrix[j][1]) {
                    fitness += 1;
                }
            }
        }
    }

	return fitness;
}

/*
 * This function evaluates the inverted solution of the given individual for the QAP problem.
 */
double SDF::EvaluateInv(int * genes)
{
    return 0;
}

/*
 * Returns the size of the problem.
 */
int SDF::GetProblemSize()
{
    return m_size;
}