/*
 *  Tools.cpp
 *  RankingEDAsCEC
 *
 *  Created by Josu Ceberio Uribe on 11/21/11.
 *  Copyright 2011 University of the Basque Country. All rights reserved.
 *
 */

#include "Tools.h"
#include <math.h>
#include <sys/time.h>
#include "Variables.h"
#include <algorithm>



/*
 * It determines if the given int sequecen if it is indeed a permutation or not.
 */
bool isPermutation(int * permutation, int size)
{
	int flags[size];
	//int * flags=new int[size];
	for (int i=0;i<size;i++) flags[i]=1;
	
	for (int i=0;i<size;i++)
	{
		int value=permutation[i];
		flags[value]=0;
	}
	
	int result,sum=0;
	for(int i=0;i<size;i++)
		sum+=flags[i];
	if (sum==0) result=true;
	else {
        result=false;
        cout <<endl;
        PrintArray(flags, size);
    }
	//delete [] flags;
	return result;
}

/*
 * Generates a random permutation of size 'n' in the given array.
 */
void GenerateRandomPermutation(int * permutation, int n)
{
	for (int i = 0; i < n; ++i)
	{
		int j = rand() % (i + 1);
		permutation[i] = permutation[j];
		permutation[j] = i;
	}
}

/*
 * Determines if a given string contains a certain substring.
 */
bool strContains(const string inputStr, const string searchStr)
{
	size_t contains;
	
	contains = inputStr.find(searchStr);
	
	if(contains != string::npos)
		return true;
	else
		return false;
}

/*
 * Prints in standard output 'length' integer elements of a given array.
 */
void PrintArray(int* array, int length, string text)
{
	cout<<text;
	for (int i=0;i<length;i++){
		cout<<array[i]<<" ";
	}
	cout<<" "<<endl;
}


/*
 * Prints in standard output 'length' long double elements of a given array.
 */
void PrintArray(long double* array, int length, string text)
{
	cout<<text;
	for (int i=0;i<length;i++){
		cout<<array[i]<<" ";
	}
	cout<<" "<<endl;
}

/*
 * Prints in the standard output the 'size' elements of the given array..
 */
void PrintArray(int * array, int size)
{
	for (int i=0;i<size;i++){
		cout<<array[i]<<" ";
	}
	cout<<" "<<endl;
}

/*
 * Prints in the standard output the 'size' elements of the given array..
 */
void PrintArray(double * array, int size)
{
	for (int i=0;i<size;i++){
		cout<<array[i]<<" ";
	}
	cout<<" "<<endl;
}

/*
 * Prints in the standard output the 'size' elements of the given array..
 */
void PrintArray(bool * array, int size)
{
	for (int i=0;i<size;i++){
		cout<<array[i]<<" ";
	}
	cout<<" "<<endl;
}

/*
 * Prints in standard output 'length' double elements of a given array.
 */
void PrintArray(float* array, int length, string text)
{
	cout<<text<<": ";
	for (int i=0;i<length;i++){
        printf(" %3.4f,",array[i]);
	}
	printf("\n ");

}

/*
 * Prints in standard output 'length' double elements of a given array.
 */
void PrintArray(double* array, int length, string text)
{
    int i;
	cout<<text;
	for (i=0;i<length;i++)
		printf(" %3.4f ",array[i]);
	printf("\n ");
}

/*
 * Prints in the standard output given matrix.
 */
void PrintMatrix(int** matrix, int length, int length2, string text)
{
    int i,j;
	cout<<text;
	for (i=0;i<length;i++)
	{
		cout<<""<<endl;
		for (j=0;j<length2;j++)
		{
			cout<<matrix[i][j]<<" ";
		}
	}
	cout<<" "<<endl;
}

void PrintMatrix(float** matrix, int length, int length2, string text)
{
    int i,j;
	cout<<text;
	for (i=0;i<length;i++)
	{
		cout<<""<<endl;
		for (j=0;j<length2;j++)
		{
			cout<<matrix[i][j]<<" ";
		}
	}
	cout<<" "<<endl;
}


/*
 * Prints in standard output 'lengthxlength' double elements of a given matrix.
 */
void PrintMatrixDouble(double** matrix, int length, int length2, string text)
{
    int i,j;
	cout<<text<<endl;
	for (i=0;i<length;i++)
	{
		for (j=0;j<length2;j++)
			printf("%3.9f, ",matrix[i][j]);
		printf("\n");
	}
}

/*
 * Calculates the tau Kendall distance between 2 permutations.
 */
int Kendall(int* permutationA, int*permutationB, int size)
{
    int i,dist;
    int * v,*composition,*invertedB,*m_aux;
    
    dist=0;
    v=new int[size-1];
    composition=new int[size];
    invertedB=new int[size];
    m_aux=new int[size];
    
    Invert(permutationB, size, invertedB);
    Compose(permutationA,invertedB,composition,size);
    vVector_Fast(v,composition,size,m_aux);
    
    for (i = 0; i < size-1; i++)
        dist += v[i];
    
    delete [] composition;
    delete [] invertedB;
    delete [] v;
    delete [] m_aux;
    
    return dist;
}

/*
 * Calculates the tau Kendall distance between 2 permutations.
 */
int Kendall(int* permutationA, int*permutationB, int size, int * m_aux)
{
    int i,dist;
    int * v,*composition,*invertedB;
    
    dist=0;
    v=new int[size-1];
    composition=new int[size];
    invertedB=new int[size];
    
    Invert(permutationB, size, invertedB);
    Compose(permutationA,invertedB,composition,size);
    vVector_Fast(v,composition,size,m_aux);

    for (i = 0; i < size-1; i++)
        dist += v[i];
    
    delete [] composition;
    delete [] invertedB;
    delete [] v;
    
    return dist;
}

/*
 * Calculates the Kendall tau distance between 2 permutations.
 * Auxiliary parameters are used for multiple continuous executions.
 */
int Kendall(int* permutationA, int*permutationB, int size, int * m_aux, int * invertedB, int *  composition, int * v)
{
    int i,dist;
    dist=0;
    Invert(permutationB, size, invertedB);
    Compose(permutationA,invertedB,composition,size);
    vVector_Fast(v,composition,size,m_aux);
    for (i = 0; i < size-1; i++)
        dist += v[i];

    return dist;
}

/*
 * Calculates the Cayley distance between 2 permutations.
 */
int Cayley(int * permutationA, int * permutationB, int size)
{
    int * invertedB= new int[size];
    int * composition= new int[size];
    int * elemsToCycles= new int[size];
    int * maxPosInCycle= new int[size];
    int * freeCycle= new int[size];

    Invert(permutationB, size, invertedB);
    Compose(permutationA,invertedB,composition,size);
    
    int index,cycle,distance;
    
    for(int i=0;i<size;i++)
    {
        elemsToCycles[i]=-1;
        maxPosInCycle[i]=-1;
        freeCycle[i]=1;
    }
    
    while((index=NextUnasignedElem(elemsToCycles,size))!=-1)
    {
        cycle=FindNewCycle(freeCycle,size);
        freeCycle[cycle]=0;
        do
        {
            elemsToCycles[index]=cycle;
            index = composition[index];//para permus de 1..n =>index = sigma[index]-1;
        }
        while(elemsToCycles[index] == -1);
    }
    distance=size-FindNewCycle(freeCycle,size);
    
    delete [] invertedB;
    delete [] composition;
    delete [] elemsToCycles;
    delete [] maxPosInCycle;
    delete [] freeCycle;
    
    return distance;
}

/*
 * Calculates the Cayley distance between 2 permutations.
 */
int Cayley(int * permutationA, int * permutationB, int size,int * invertedB, int *  composition, int * elemsToCycles, int * maxPosInCycle, int * freeCycle)
{
    
    Invert(permutationB, size, invertedB);
    Compose(permutationA,invertedB,composition,size);
    
    int index,cycle,distance;
    
    for(int i=0;i<size;i++)
    {
        elemsToCycles[i]=-1;
        maxPosInCycle[i]=-1;
        freeCycle[i]=1;
    }
    
    while((index=NextUnasignedElem(elemsToCycles,size))!=-1)
    {
        cycle=FindNewCycle(freeCycle,size);
        freeCycle[cycle]=0;
        do
        {
            elemsToCycles[index]=cycle;
            index = composition[index];//para permus de 1..n =>index = sigma[index]-1;
        }
        while(elemsToCycles[index] == -1);
    }
    distance=size-FindNewCycle(freeCycle,size);
    
    return distance;
}


int FindNewCycle(int * freeCycle, int size){
    int i;
    
    
    for(i=0;i<size;i++)
        if(freeCycle[i])
            return i;
    return size;
}

int NextUnasignedElem(int * elemsToCycles, int size)
{
    int i;
    for(i=0;i<size;i++)
        if(elemsToCycles[i]==-1)
            return i;
    return -1;
}

/*
 * Calculates the Ulam distance between 2 permutations.
 */
int Ulam(int * permutationA, int * permutationB, int size){

    int * composition=new int[size];
    int * invertedB=new int[size];
    
    Invert(permutationB, size, invertedB);
    Compose(permutationA,invertedB,composition,size);

    int dist=(size- getLISLength(composition,size));
    
    delete [] composition;
    delete [] invertedB;
    
    return dist;
}

/*
 * Calculates the length of the longest increasing subsequence in the given array of ints.
 */
int getLISLength(int*sigma,int size){
   
    // O(n log k)
    
    int i;
    vector<int> vc(1,sigma[0]);
    vector<int>::iterator vk;
    
    for (i=1;i<size;i++){
        for (vk=vc.begin(); vk != vc.end(); vk++)
            if (*vk>=sigma[i]) break;
        if (vk==vc.end())
            vc.push_back(sigma[i]);
        else *vk=sigma[i];
    }
    
    return (int)vc.size();
    
}
/*
 * Implements the compose of 2 permutations of size n.
 */
void Compose(int*s1, int*s2, int*res, int n)
{
    int i;
    for(i=0;i<n;i++)
        res[i]=s1[s2[i]];
}


/*
 * Calculates V_j-s vector.
 */
void vVector(int*v, int*permutation, int n)
{
    
    int i,j;
    for(i=0;i<n-1;i++)
        v[i]=0;
    
    for (i = n-2; i >= 0; i--)
        for (j = i+1; j < n; j++)
            if(permutation[i] > permutation[j])
                v[i]++;
}

/*
 *  Optimized version proposed by Leti for the calculation of the V_j-s vector.
 */
void vVector_Fast(int*v, int*permutation, int n, int * m_aux)
{
    int i,j, index;
    
    for(i=0;i<n-1;i++){
        v[i]=0;
        m_aux[i]=0;
    }
    m_aux[n-1]=0;
    for (j=0; j<n-1; j++){
        index=permutation[j];
        v[j]=index-m_aux[index];
        for (i=index; i<n; i++)
            m_aux[i]++;
    }
}

/*
 * Inverts a permutation.
 */
void Invert(int*permu, int n, int* inverted)
{
    int i;
    for(i=0; i<n; i++)
        inverted[permu[i]]=i;
}

/*
 * Applies the random keys sorting strategy to the vector of doubles
 */
void RandomKeys( int * a, float * criteriaValues, int size)
{
	bool * fixedValues= new bool[size];
	float criteria, min;
	int i, j;
	for (i=0;i<size;i++){
		fixedValues[i]=false;
		a[i]=0;
	}
	int minPos=0;
	for (i=0;i<size;i++)
	{
		min=MAX_INTEGER;
		for (j=0;j<size;j++)
		{
			criteria=criteriaValues[j];
			if (!fixedValues[j] && min>criteria )
			{
				min=criteria;
				minPos=j;
			}
		}
        
		fixedValues[minPos]=true;
		//a[i]=minPos;// modification por el asunto ordering /ranking
    	a[minPos]=i;// original.
	}
	delete [] fixedValues;
}



/*
 * Calculates the factorial of a solution.
 */
long double factorial(int val) {
    if(val <= 0) return 1;
    //long  N, b, c, p; // use int for fast calculation and small range of calculation..
    long   b, c;
    long double p, N;
    N=(long double)val;
    c = (long)N - 1;
    p = 1;
    while (c > 0) {
        p = 0;
        b = c;
        while (b > 0) {
            if (b & 1) {
                p += N; // p = p + N;
            }
            // if you would like to use double choose the alternative forms instead shifts
            // the code is fast even!
            // you can use the same tips on double or 64 bit int etc.... but you must... ;-)
            //b >>= 1; // b/=2; (b = b / 2;) ( b >> 1; a.s.r. is more efficent for int or long..!)
            b/=2;
            //N <<= 1; // N += N; N = N + N; N = N * 2; (N <<=1; a.s.l. is more efficent for int or long..!)
            N += N;
        } // end of: while(b>0)
        N = p;
        c--; // c = c - 1;
    } // end of: while(c > 0)
    //printf("[%d] is the factorial! \n", p);
    return p;
}

/*
 * This method applies a swap of the given i,j positions in the array.
 */
void Swap(int * array, int i, int j)
{
	int aux=array[i];
	array[i]=array[j];
	array[j]=aux;
}

// integer array
double calculate_entropy(int *arr, int arr_size) {
    
    int sum = 0;
    for (int i = 0; i < arr_size; i ++) {
        if (arr[i] != 0) {
        // if (arr[i]) {
            sum += arr[i];
        }
    }

    double entropy = 0;
    for (int i = 0; i < arr_size; i ++) {
        
        if (arr[i] != 0) {
        // if (arr[i]) {
            double p = double(arr[i]) / sum;
            entropy += p * log(p); // nature log
        }
    }

    return -entropy;
}

// floating point array
double calculate_entropy(double *arr, int arr_size) {
    
    double sum = 0;
    for (int i = 0; i < arr_size; i ++) {
        sum += arr[i];
    }
    // cout << "\n\nsum: " << sum << endl;
    double entropy = 0;
    for (int i = 0; i < arr_size; i ++) {
        
        if (arr[i] != 0) {
            double p = arr[i] / sum;
            // cout << "p: " << p << endl;
            entropy += p * log(p); // nature log
            // cout << "p * log(p): " << p * log(p) << endl;

        }
    }
    // cout << endl;

    return -entropy;
}

// floating point array with known summation
double calculate_entropy_with_known_sum(double *arr, int arr_size, double sum) {
    
    // double sum = 0;
    // for (int i = 0; i < arr_size; i ++) {
    //     sum += arr[i];
    // }

    double entropy = 0;
    for (int i = 0; i < arr_size; i ++) {
        
        if (arr[i] != 0) {
            double p = arr[i] / sum;
            entropy += p * log(p); // nature log
            // cout << "p * log(p): " << p * log(p) << endl;
        }
    }
    // cout << endl;

    return -entropy;
}

double update_entropy(double original_entropy, double deleted_element, double old_sum) {

    // original_entropy = original_entropy * (-1);
    // double new_sum = old_sum - deleted_element;
    // double deleted_element_entropy =  (deleted_element / old_sum) * log(deleted_element / old_sum);
    // return (-1) * (old_sum / new_sum) * (  original_entropy - deleted_element_entropy + (new_sum / old_sum) * log(old_sum / new_sum));
    if (deleted_element == 0.0) {
        return original_entropy;
    }
    else if (deleted_element == old_sum) {
        return 100.0;
    }
    original_entropy = original_entropy * (-1);
    double new_sum = old_sum - deleted_element;
    double deleted_element_entropy =  (deleted_element / old_sum) * log(deleted_element / old_sum);
    return (-1) * (old_sum / new_sum) * (  original_entropy - deleted_element_entropy + (new_sum / old_sum) * log(old_sum / new_sum));

}

double update_entropy(double original_entropy, int deleted_element, int old_sum) {

    // original_entropy = original_entropy * (-1);
    // double new_sum = old_sum - deleted_element;
    // double deleted_element_entropy =  (deleted_element / old_sum) * log(deleted_element / old_sum);
    // return (-1) * (old_sum / new_sum) * (  original_entropy - deleted_element_entropy + (new_sum / old_sum) * log(old_sum / new_sum));
    if (deleted_element == 0) {
        return original_entropy;
    }
    original_entropy = original_entropy * (-1);
    int new_sum = old_sum - deleted_element;
    double deleted_element_entropy =  ( double(deleted_element) / old_sum) * log( double(deleted_element) / old_sum);
    return (-1) * ( double(old_sum) / new_sum) * (  original_entropy - deleted_element_entropy + ( double(new_sum) / old_sum) * log( double(old_sum) / new_sum));

}

// sample from integer array
int sample_from_array(int *arr, int ell) {
    
    int sum = 0;
    for (int i = 0; i < ell; i++) {
        sum += arr[i];
    }

    int rand_num;
    int idx;
    
    if (sum == 0) {
        return -1;
    }

    else {
        rand_num = rand() % sum;
        // printf("rand: %d\n", rand_num);

        idx = -1;
        while (rand_num >= 0) {
            // printf("idx: %d\n", idx);
            idx += 1;
            rand_num -= arr[idx];
        }
        return idx;
    }
}

// sample from floating point array
int sample_from_array(double *arr, int ell) {
    
    double sum = 0;
    for (int i = 0; i < ell; i++) {
        sum += arr[i];
    }

    double rand_num;
    int idx;
    
    if (sum == 0) {
        return -1;
    }

    else {
        // rand_num = rand() % sum;
        rand_num = fmod(rand(), sum);
        idx = -1;
        while (rand_num >= 0) {
            idx += 1;
            rand_num -= arr[idx];
        }
        return idx;
    }
}

// sample from floating point array square
int sample_from_array_pow(double *arr, int ell, double power) {
    
    double sum = 0;
    for (int i = 0; i < ell; i++) {
        sum += pow(arr[i], power);
    }

    double rand_num;
    int idx;
    
    if (sum == 0) {
        return -1;
    }

    else {
        // rand_num = rand() % sum;
        rand_num = fmod(rand(), sum);
        idx = -1;
        while (rand_num >= 0) {
            idx += 1;
            rand_num -= pow(arr[idx], power);
        }
        return idx;
    }
}

// sample from floating point array softmax
int sample_from_array_softmax(double *arr, int ell) {
    
    double sum = 0;
    for (int i = 0; i < ell; i++) {
        sum += exp(arr[i]);
    }

    double rand_num;
    int idx;
    
    if (sum == 0) {
        return -1;
    }

    else {
        // rand_num = rand() % sum;
        rand_num = fmod(rand(), sum);
        idx = -1;
        while (rand_num >= 0) {
            idx += 1;
            rand_num -= exp(arr[idx]);
        }
        return idx;
    }
}

void AllocateMatrixMemory(double **matrix, int ell) {
    for (int i = 0; i < ell; ++i) {
        matrix[i] = new double[ell];
    }
}
void freeMatrixMemory(double **matrix, int ell) {
    for (int i = 0; i < ell; ++i) {
        delete [] matrix[i];
    }
}
void resetMatrix(double **matrix, int ell, double eps) { // eps: epsilon

	// todo: B_ratio
	// todo: B_ratio
	// todo: B_ratio

    for (int i = 0; i < ell; ++i) {
        for (int j = 0; j < ell; ++j) {
			if (i == j) {
                matrix[i][j] = 0.0;
            } 
            else {
                matrix[i][j] = eps;
            }
		}
    }
}
void resetNHM(double **matrix, int ell, double eps) { // eps: epsilon

	// todo: B_ratio
	// todo: B_ratio
	// todo: B_ratio

    for (int i = 0; i < ell; ++i) {
        for (int j = 0; j < ell; ++j) {
            matrix[i][j] = eps;
		}
    }
}
void PrintMatrix(double** matrix_ptr, int size) {

    // cout<<"half_relation_distance: " << half_relation_distance << endl;
    printf("\ni/j|");
    for(int i = 0; i < size; ++i){
        printf("    %d", i);
    }
    printf("\n----");
    for(int i = 0; i < size; ++i){
        printf("-----");
    }
    printf("\n");

    for(int i = 0; i < size; ++i){
        printf("%2d |", i);

        for(int j = 0; j < size; ++j){
            printf("%5.1lf", matrix_ptr[i][j]);
        }
        printf("   => entropy: %3lf", calculate_entropy(matrix_ptr[i], size));
        printf("\n");
    }
}

void PrintMatrix(int** matrix_ptr, int size) {

    // cout<<"half_relation_distance: " << half_relation_distance << endl;
    printf("\ni/j|");
    for(int i = 0; i < size; ++i){
        printf("    %d", i);
    }
    printf("\n----");
    for(int i = 0; i < size; ++i){
        printf("-----");
    }
    printf("\n");

    for(int i = 0; i < size; ++i){
        printf("%2d |", i);

        for(int j = 0; j < size; ++j){
            printf("%5d", matrix_ptr[i][j]);
        }
        // printf("   => entropy: %3lf", calculate_entropy(matrix_ptr[i], size));
        printf("\n");
    }
}

void PrintMRMatrix(double ***matrix, int ell) {
    printf("\n=============\n\n");
    for (int reference_node_num = 0; reference_node_num < ell; ++reference_node_num) {
        printf("reference_node: %d\n", reference_node_num);
        printf("D/N|");

        for (int next_node_num = 0; next_node_num < ell; ++next_node_num) {
             printf("%5d", next_node_num);
        }
        printf("\n");
        for (int next_node_num = 0; next_node_num <= ell; ++next_node_num) { 
            printf("-----");
        }
        printf("\n");
        for (int distance = 0; distance < ell; ++distance) { 
            // int distance = ell / 2 + 1;
            printf("%3d|", distance);
            for (int next_node_num = 0; next_node_num < ell; ++next_node_num) { 
                printf("%5.1lf", matrix[reference_node_num][distance][next_node_num]);
            }
            printf("   => entropy: %3lf", calculate_entropy(matrix[reference_node_num][distance], ell));
            printf("\n");
        }
        printf("\n");
    }
}

void PrintMRMatrix(int ***matrix, int ell) {
    printf("\n=============\n\n");
    for (int reference_node_num = 0; reference_node_num < ell; ++reference_node_num) {
        printf("reference_node: %d\n", reference_node_num);
        printf("D/N|");

        for (int next_node_num = 0; next_node_num < ell; ++next_node_num) {
             printf("%5d", next_node_num);
        }
        printf("\n");
        for (int next_node_num = 0; next_node_num <= ell; ++next_node_num) { 
            printf("-----");
        }
        printf("\n");
        for (int distance = 0; distance < ell; ++distance) { 
            // int distance = ell / 2 + 1;
            printf("%3d|", distance);
            for (int next_node_num = 0; next_node_num < ell; ++next_node_num) { 
                // printf("%5.1lf", matrix[reference_node_num][distance][next_node_num]);
                printf("%6d", matrix[reference_node_num][distance][next_node_num]);
            }
            printf("   => entropy: %3lf", calculate_entropy(matrix[reference_node_num][distance], ell));
            printf("\n");
        }
        printf("\n");
    }
}

void PrintMSTEdgeMatrix(double ***matrix, int ell) {
    printf("\n=============\n\n");
    for (int reference_node_num = 0; reference_node_num < ell; ++reference_node_num) {
        printf("reference_node: %d\n", reference_node_num);
        printf("N/D|");

        for (int distance = 0; distance < ell; ++distance) {
             printf("%5d", distance);
        }
        printf("\n");
        for (int distance = 0; distance <= ell; ++distance) { 
            printf("-----");
        }
        printf("\n");
        for (int next_node_num = 0; next_node_num < ell; ++next_node_num) { 
            printf("%3d|", next_node_num);
            for (int distance = 0; distance < ell; ++distance) { 
                printf("%5.1lf", matrix[reference_node_num][next_node_num][distance]);
            }
            printf("   => entropy: %3lf", calculate_entropy(matrix[reference_node_num][next_node_num], ell));
            printf("\n");
        }
        printf("\n");
    }
}

// build EHM
void buildModel(int **pop, int pop_size, double **ehm, int ell) {

	for (int individual_num = 0; individual_num < pop_size; ++individual_num) {
		for (int gene_num = 0; gene_num < ell; ++gene_num) {

			int front_idx = gene_num - 1, back_idx = gene_num + 1;

			// first gene
			if (gene_num == 0) {
				front_idx += ell;
			}

			// last gene
			else if (gene_num == ell - 1) {
				back_idx -= ell;
			}

			ehm[ pop[individual_num][gene_num] ][ pop[individual_num][front_idx] ] += 1;
			ehm[ pop[individual_num][gene_num] ][ pop[individual_num][back_idx] ] += 1;

		}
	}
}

// build NHM
void buildNHM(int **pop, int pop_size, double **nhm, int ell) {

    for (int individual_num = 0; individual_num < pop_size; ++individual_num) {
        for (int idx = 0; idx < ell; idx++) {
            nhm [idx] [ pop[individual_num][idx] ] ++;
        }
    }
}


double sum_arr(double *arr, int ell) {
    
    double sum = 0;
    for (int i = 0; i < ell; i++) {
        sum += arr[i];
    }
    return sum;
}

int sum_arr(int *arr, int ell) {
    
    int sum = 0;
    for (int i = 0; i < ell; i++) {
        sum += arr[i];
    }
    return sum;
}



// summation start from start_idx, but not include end_idx
double sum_arr_interval(double *arr, int start_idx, int end_idx) {
    
    double sum = 0;
    for (int i = start_idx; i < end_idx; i++) {
        sum += arr[i];
    }
    return sum;
}


// -- union find functions part -- //
int find_set(int x, element * arr) {
    
    if (x != arr[x].parent) {
        arr[x].parent = find_set(arr[x].parent, arr);
    }
    return arr[x].parent;
}

void link(int x, int y, element * arr) {
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

void union_set(int x, int y, element * arr) {
    link(find_set(x, arr), find_set(y, arr), arr);
}

// -- union find functions part end -- //

void CutPointGenerator(int problem_size, int cut_point_size, int* chosen_cut_point) {
	
	int * cut_point_arr = new int[cut_point_size];
	int * zeroToEllMinusOneArr = new int[problem_size];
	
	// initialize
	for (int i = 0; i < problem_size; ++i) {
		zeroToEllMinusOneArr[i] = i;
	}

	// random shuffle and assiign cut point value
	random_shuffle(zeroToEllMinusOneArr, zeroToEllMinusOneArr + problem_size);
	for (int i = 0; i < cut_point_size; ++i) {
		cut_point_arr[i] = zeroToEllMinusOneArr[i];
	}


	sort(cut_point_arr, cut_point_arr + cut_point_size); 
	int cut_point_idx = rand() % cut_point_size;
	chosen_cut_point[0] = cut_point_arr[cut_point_idx];
	
	if (cut_point_idx >= cut_point_size - 1) {
		chosen_cut_point[1] = cut_point_arr[0];
	} 
	else {
		chosen_cut_point[1] = cut_point_arr[cut_point_idx + 1];
	}
	delete [] cut_point_arr;
	delete [] zeroToEllMinusOneArr;
}



bool is_same_permutation(int *permu1, int *permu2, int ell) {

    for (int i = 0; i < ell; i ++) {
        if (permu1[i] != permu2[i]) {
            return false;
        }
    }
    return true;
}


int find_value_index_in_permu(int value, int *permu, int ell) {

    for (int i = 0; i < ell; i ++) {
        if (permu[i] == value) {
            return i;
        }
    }
    return -1;
}

// copy permu2 to permu1
void copy_array(int *permu1, int *permu2, int ell) {

    for (int i = 0; i < ell; i ++) {
        permu1[i] = permu2[i];
    }
}


void copy_array(double *permu1, double *permu2, int ell) {

    for (int i = 0; i < ell; i ++) {
        permu1[i] = permu2[i];
    }
}
