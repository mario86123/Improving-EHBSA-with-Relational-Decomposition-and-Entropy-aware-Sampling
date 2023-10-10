/*
 *  SOP.cpp
 *  RankingEDAsCEC
 *
 *  Created by Josu Ceberio Uribe on 7/11/13.
 *  Copyright 2013 University of the Basque Country. All rights reserved.
 *
 */

#include "SOP.h"
#include <cmath>


/*
 *Class constructor.
 */
SOP::SOP()
{
	
}

/*
 * Class destructor.
 */
SOP::~SOP()
{
	for (int i=0;i<m_size;i++)
        delete [] m_distance_matrix[i];
    delete [] m_distance_matrix;

	for (int i=0;i<100*m_size;i++)
        delete [] m_distance_matrix[i];
    delete [] m_distance_matrix;

    delete [] m_aux;
}


double SOP::CalculateGEODistance(double latitudeX, double latitudeY, double longitudeX, double longitudeY)
{

	// printf("laX: %.2lf\n", latitudeX);
	// printf("laY: %.2lf\n", latitudeY);
	// printf("loX: %.2lf\n", longitudeX);
	// printf("loY: %.2lf\n", longitudeY);

	double PI = 3.141592;
	double RRR = 6378.388;
	
	double deg = (double)((int)latitudeX);
	double min = latitudeX - deg;
	double latitudeI = PI * (deg + 5.0 * min / 3.0 ) / 180.0;
	
	deg = (double)((int)latitudeY);
	min = latitudeY - deg;
	double longitudeI = PI * (deg + 5.0 * min / 3.0 ) / 180.0;
	
	deg = (double)((int)longitudeX);
	min = longitudeX - deg;
	double latitudeJ = PI * (deg + 5.0 * min / 3.0 ) / 180.0;
	
	deg = (double)((int)longitudeY);
	min = longitudeY - deg;
	double longitudeJ = PI * (deg + 5.0 * min / 3.0 ) / 180.0;
	
	double q1 = cos(longitudeI - longitudeJ);
	double q2 = cos(latitudeI - latitudeJ);
	double q3 = cos(latitudeI + latitudeJ);
	
	return (int)(RRR * acos(0.5 * ((1.0+q1) * q2 - (1.0 - q1) * q3)) + 1.0);
}

/*
 * Read SOP instance file that belongs to the SOPLIB library.
 */
int SOP::Read(string filename)
{
    //declaration and initialization of variables.
	bool readm_distance_matrix=false;
	bool coordinatesData=false;
	string DISTANCE_TYPE;
	double **coordinates;
	char line[2048]; // variable for input value
	string data="";
	ifstream indata;
	indata.open(filename.c_str(),ios::in);
	
	while (!indata.eof())
	{
		//LEER LA LINEA DEL FICHERO
		indata.getline(line, 2048);
		stringstream ss;
		string sline;
		ss << line;
		ss >> sline;
		// cout << "sline: " << sline << endl;
		if (sline=="EOF")
		{
			break;
		}
		if (readm_distance_matrix && coordinatesData==false)
		{
			// cout << "reading distance m_distance_matrix "<<line<< endl;
			if (data=="")
				data = line;
			else
				data = data+' '+line;
		}
		else if (readm_distance_matrix && coordinatesData==true)
		{
			//FILL DISTANCE m_distance_matrix
			char * coordPieces;
			coordPieces = strtok (line," ");
			if (coordPieces==" ")
			{
				coordPieces = strtok (NULL," ");
			}
			int coordNum = atoi(coordPieces);
			coordPieces = strtok (NULL, " ");
			double latitud = atof(coordPieces);
			coordPieces = strtok (NULL, " ");
			double longitud = atof(coordPieces);
			double *coordinate= new double[2];
			coordinate[0]=latitud;
			coordinate[1]=longitud;
			
			coordinates[coordNum-1]= coordinate;
			//cout<<"coordNum "<<coordNum-1<<" latit: "<<latitud<<" long: "<<longitud<<endl;
		}
		
		if(strContains(sline,"DIMENSION"))
		{
			char * pch;
			pch = strtok (line," ");
			pch = strtok (NULL, " ");
			if (strcmp(pch,":")==0)
			{
				pch = strtok (NULL, " ");
			}
			m_size = atoi(pch);
			// cout << "m_size: " << m_size << endl;
		}
		else if (strContains(sline,"EDGE_WEIGHT_TYPE"))
		{
			char * pch;
			pch = strtok (line," ");
			pch = strtok (NULL, " ");
			if (strcmp(pch,":")==0)
			{
				pch = strtok (NULL, " ");
			}
			stringstream s;
			string type;
			s << pch;
			s >> type;
			DISTANCE_TYPE = type;
		}
        else if (sline =="EDGE_WEIGHT_SECTION")
		{
			readm_distance_matrix=true;
			coordinatesData=false;
		}
		else if (sline=="NODE_COORD_SECTION"){
			readm_distance_matrix=true;
			coordinatesData=true;
			coordinates= new double*[m_size];
		}
		
	}
	indata.close();
	
	//BUILD DISTANCE m_distance_matrix
	m_distance_matrix = new double*[m_size];
	for (int i=0;i<m_size;i++)
	{
		m_distance_matrix[i]= new double[m_size];
	}
	
	//BUILD constarint pair m_constraint_pair
	m_constraint_pair = new int*[100*m_size];
	for (int i=0;i<100*m_size;i++)
	{
		m_constraint_pair[i]= new int[2];
	}
	num_of_constraint_pair = 0;

	//FILL DISTANCE m_distance_matrix
	if (coordinatesData==true)
	{
		// cout << "177" << endl;
		//CALCULATE EUCLIDEAN DISTANCES
		for (int i=0;i<m_size;i++)
		{
			//get coordinate A
			double *coordA=coordinates[i];
			double coordAx = coordA[0];
			double coordAy = coordA[1];
			for (int j=i;j<m_size;j++)
			{
				//get coordinate B.
				double *coordB=coordinates[j];
				double coordBx=coordB[0];
				double coordBy=coordB[1];
				double euclidean;
				if (DISTANCE_TYPE=="GEO")
				{
					//calculate geographic distance between A and B.
					euclidean=CalculateGEODistance(coordAx, coordAy, coordBx, coordBy);
				}
				else
				{
					//calculate euclidean distance between A and B.
					double absolute= fabs(pow((coordAx-coordBx),2) + pow((coordAy-coordBy),2));
					euclidean= sqrt(absolute);
				}
				// //convert double to string
				// std::ostringstream stream;
				// stream << euclidean;
				// std::string euclideanString = stream.str();

				m_distance_matrix[i][j]=  euclidean;
				m_distance_matrix[j][i]= euclidean;//<-symmetric m_distance_matrix
			}
		}
	}
	else
	{
		// cout << "full distance" << endl;
		//FILL DISTANCE m_distance_matrix
		istringstream iss(data);
		// cout << "data: "<< data <<endl;
		int i=0;
		int j=0;
		do
		{
			string sub;
			// cout << "sub: " << sub << endl;
		    iss >> sub;
			// printf("sub: %s ", sub.c_str());//)  << sub << endl;
			// cout << "m[" << i << "][" << j << "]: ";
		    
			//old
			if (sub!=""){
				//save distance in distances m_distance_matrix. Save negative distance in order to minimize fitness instead of
				//maximize.
				// cout << "sub: " << sub << endl;
		    	m_distance_matrix[i][j]= atoi(sub.c_str());

				if (m_distance_matrix[i][j] == -1) {
					m_constraint_pair[num_of_constraint_pair][0] = j;
					m_constraint_pair[num_of_constraint_pair][1] = i;
					num_of_constraint_pair ++;
					// cout << "num_of_constraint_pair: " << num_of_constraint_pair << endl;
				}
		    	// m_distance_matrix[j][i]= atoi(sub.c_str());//<-symmetric m_distance_matrix
				// printf("%lf\n", m_distance_matrix[i][j]);
				// cout << sub << endl;
		    	// if (sub=="\n")
		    	if (j >= m_size - 1)
		    	{
		    		i++;
		    		j=0;
		    	}
		    	else
		    	{
		    		j++;
		    	}
		    }
		    else
		    {
		    	break;
		    }
			// old end

			// if (sub!=""){
			// 	//save distance in distances m_distance_matrix. Save negative distance in order to minimize fitness instead of
			// 	//maximize.
			// 	printf("sub: %s ", sub.c_str());//)  << sub << endl;
			// 	cout << "m[" << i << "][" << j << "]: ";
			// 	if (i == j) {m_distance_matrix[i][j]= 0;}
		    // 	else {m_distance_matrix[i][j]= atoi(sub.c_str());}

		    // 	// m_distance_matrix[j][i]= atoi(sub.c_str());//<-symmetric m_distance_matrix
			// 	printf("%d \n", m_distance_matrix[i][j]);
			// 	// cout << sub << endl;
		    // 	if (i<m_size-1)
		    // 	{
		    // 		i++;
		    // 		// j=0;
		    // 	}
		    // 	else
		    // 	{
		    // 		j++;
			// 		i=0;
		    // 	}
		    // }
		    // else
		    // {
			// 	break;
		    // }
			// cout << "here!!" << endl;

		} while (iss);
		// cout << "end" << endl;

		// for (int i = 0; i < num_of_constraint_pair; i++) {
		// 	cout << i << ": " << m_constraint_pair[i][0] << ", " <<  m_constraint_pair[i][1] << endl;
		// }
		// cout << endl;
	}
    m_aux= new int[m_size];
    return (m_size - 2);
}



/*
 * Read SOP instance file.
 */
int SOP::Read2(string filename)
{
	char line[2048]; // variable for input value
	string data="";
	ifstream indata;
	indata.open(filename.c_str(),ios::in);
	int num=0;
	while (!indata.eof())
	{
        
		indata.getline(line, 2048);
		stringstream ss;
		string sline;
		ss << line;
		ss >> sline;
		if (sline=="")
		{
			break;
		}
		if (num==0)
		{
			m_size = atoi(line);
		}
		else
		{
			if (data=="")
				data = line;
			else
				data = data+' '+line;
		}
		num++;
	}
	indata.close();
    
	//BUILD MATRIX
	m_distance_matrix = new double*[m_size];
	for (int i=0;i<m_size;i++)
	{
		m_distance_matrix[i]= new double[m_size];
	}
    m_aux= new int[m_size];
    
	istringstream iss(data);
	int i=0;
	int j=0;
	do
	{
		string sub;
	    iss >> sub;
	    if (sub!=""){
			//save distance in distances matrix.
	    	m_distance_matrix[i][j]= atoi(sub.c_str());
	    	if (j==(m_size-1))
	    	{
	    		i++;
	    		j=0;
	    	}
	    	else
	    	{
	    		j++;
	    	}
	    }
	    else
	    {
	    	break;
	    }
	} while (iss);
	return (m_size);
}

 /*
 * This function evaluates the fitness of the solution for the SOP problem.
 */
double SOP::Evaluate(int * genes)
{
	double distanceSum=0;
	double distAB=0;
	int IDCityA, IDCityB;
    
	int new_genes[m_size];
	int new_genes_index[m_size];



	new_genes[0] = 0;
	new_genes_index[0] = 0;

	new_genes[m_size - 1] = m_size - 1;
	new_genes_index[m_size - 1] = m_size - 1;

	for (int i=1; i < m_size - 1; i++) {
		new_genes[i] = genes[i - 1] + 1;
		new_genes_index[new_genes[i]] = i;
	}

	
	// cout << "permutation: ";
	// for (int i=0;i<m_size;i++) {
	// 	cout << new_genes[i] << " ";
	// }
	// cout << endl;

	// cout << "permutation_index: ";
	// for (int i=0;i<m_size;i++) {
	// 	cout << new_genes_index[i] << " ";
	// }
	// cout << endl;

	for (int i=0;i<num_of_constraint_pair;i++) {
		// cout << m_constraint_pair[i][0] << " index: " << new_genes_index[m_constraint_pair[i][0]] << ", ";
		// cout << m_constraint_pair[i][1] << " index: " << new_genes_index[m_constraint_pair[i][1]] << endl;
		if (new_genes_index[m_constraint_pair[i][0]] > new_genes_index[m_constraint_pair[i][1]]) {
			distanceSum += 10001;
		}

	}



	for (int i=0;i<m_size-1;i++)
	{
		IDCityA = new_genes[i];
		IDCityB = new_genes[0];
		if (i+1<m_size)
		{
			IDCityB = new_genes[i+1];
		}
		
		distAB = round(m_distance_matrix[IDCityA][IDCityB]);
		// cout << "A: " << IDCityA<< ", B: " << IDCityB <<" distAB: " << distAB << endl;
		// if (distAB == -1) {
		// 	return -100000;
		// }

		distanceSum=distanceSum+distAB;
		
		//security condition
		if (IDCityA==m_size || IDCityB==m_size){
			distanceSum=0;
			break;
		}
	}
	// cout << endl;
	return -distanceSum;
}

/*
 * This function evaluates the inverted solution of the given solution for the SOP problem.
 */
double SOP::EvaluateInv(int * genes){
    Invert(genes, m_size, m_aux);
    
    double distanceSum=0;
	int distAB=0;
	int IDCityA, IDCityB;
	for (int i=0;i<m_size;i++)
	{
		IDCityA = m_aux[i];
		IDCityB = m_aux[0];
		if (i+1<m_size)
		{
			IDCityB = m_aux[i+1];
		}
		
		distAB = m_distance_matrix[IDCityA][IDCityB];
		distanceSum=distanceSum+distAB;
		
		//security condition
		if (IDCityA==m_size || IDCityB==m_size){
			distanceSum=0;
			break;
		}
	}
	return -distanceSum;
}

/*
 * Returns the size of the problem.
 */
int SOP::GetProblemSize()
{
    return m_size;
}

