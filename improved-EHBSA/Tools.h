/*
 *  Tools.h
 *  RankingEDAsCEC
 *
 *  Created by Josu Ceberio Uribe on 11/21/11.
 *  Copyright 2011 University of the Basque Country. All rights reserved.
 *
 */

#ifndef TOOLS_H
#define TOOLS_H

#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string.h>
#include <stdio.h>
#include <vector>

using std::istream;
using std::ostream;
using namespace std;
using std::cerr;
using std::cout;
using std::endl;
using std::ifstream;
using std::stringstream;
using std::string;

/*
 * It determines if the given int sequecen if it is indeed a permutation or not.
 */
bool isPermutation(int * permutation, int size);

/*
 * Generates a random permutation of size 'n' in the given array.
 */
void GenerateRandomPermutation(int * permutation, int n);

/*
 * Determines if a given string contains a certain substring.
 */
bool strContains(const string inputStr, const string searchStr);

/*
 * Prints in standard output 'length' integer elements of a given array.
 */
void PrintArray(int* array, int length, string text);

/*
 * Prints in standard output 'length' float elements of a given array.
 */
void PrintArray(float* array, int length, string text);

/*
 * Prints in standard output 'length' long double elements of a given array.
 */
void PrintArray(long double* array, int length, string text);

/*
 * Prints the given doubles array in the standard output.
 */
void PrintArray(double* array, int length, string text);
void PrintArray(int * array, int size);
void PrintArray(double * array, int size);
void PrintArray(bool * array, int size);

/*
 * Prints the given integers matrix in the standard output.
 */
void PrintMatrix(int** matrix, int length, int length2, string text);

/*
 * Prints the given floats matrix in the standard output.
 */
void PrintMatrix(float** matrix, int length, int length2, string text);

/*
 * Applies the random keys sorting strategy to the vector of doubles
 */
void RandomKeys( int * a, float * criteriaValues, int size);

/*
 * Calculates the tau Kendall distance between 2 permutations.
 */
int Kendall(int* permutationA, int*permutationB, int size);

/*
 * Calculates the Kendall tau distance between 2 permutations.
 */
int Kendall(int* permutationA, int*permutationB, int size, int * m_aux);

/*
 * Calculates the Kendall tau distance between 2 permutations.
 * Auxiliary parameters are used for multiple executions.
 */
int Kendall(int* permutationA, int*permutationB, int size, int * m_aux, int * invertedB, int *  composition, int * v);

/*
 * Calculates the Cayley distance between 2 permutations.
 */
int Cayley(int * permutationA, int * permutationB, int size);

/*
 * Calculates the Cayley distance between 2 permutations.
 */
int Cayley(int * permutationA, int * permutationB, int size,int * invertedB, int *  composition, int * elemsToCycles, int * maxPosInCycle, int * freeCycle);
int FindNewCycle(int * freeCycle, int size);
int NextUnasignedElem(int * elemsToCycles, int size);
int CalculateDistance(int*sigma, int size);

/*
 * Calculates the Ulam distance between 2 permutations.
 */
int Ulam(int * permutationA, int * permutationB, int size);

/*
 * Calculates the length of the longest increasing subsequence in the given array of ints.
 */
int getLISLength(int*sigma,int size);

/*
 * Implements the compose of 2 permutations of size n.
 */
void Compose(int*s1, int*s2, int*res, int n);

/*
* Calculates V_j-s vector.
*/
void vVector(int*v, int*permutation, int n);

/*
 *  Optimized version by Leti of the V_j-s vector calculation.
 */
void vVector_Fast(int*v, int*permutation, int n, int * m_aux);

/*
 * Inverts a permutation.
 */
void Invert(int*permu, int n, int* inverted);

/*
 * Calculates the factorial of a solution.
 */
long double factorial(int val);

/*
 * This method applies a swap of the given i,j positions in the array.
 */
void Swap(int * array, int i, int j);

double calculate_entropy(int *arr, int arr_size);
double calculate_entropy(double *arr, int arr_size);
double calculate_entropy_with_known_sum(double *arr, int arr_size, double sum);
double update_entropy(double original_entropy, double deleted_element, double old_sum);
double update_entropy(double original_entropy, int deleted_element, int old_sum);
int sample_from_array(int *arr, int ell);
int sample_from_array(double *arr, int ell);
int sample_from_array_pow(double *arr, int ell, double power);
int sample_from_array_softmax(double *arr, int ell);
void AllocateMatrixMemory(double **matrix, int ell);
void freeMatrixMemory(double **matrix, int ell);
void resetMatrix(double **matrix, int ell, double eps);
void resetNHM(double **matrix, int ell, double eps);
void PrintMatrix(double** matrix_ptr, int size);
void PrintMatrix(int** matrix_ptr, int size);
void PrintMRMatrix(double ***matrix, int ell);
void PrintMRMatrix(int ***matrix, int ell);
void PrintMSTEdgeMatrix(double ***matrix, int ell);
void buildModel(int **pop, int pop_size, double **ehm, int ell);
void buildNHM(int **pop, int pop_size, double **ehm, int ell);
double sum_arr(double *arr, int ell);
int sum_arr(int *arr, int ell);
double sum_arr_interval(double *arr, int start_idx, int end_idx);

struct element {
    int parent;
    int rank;
};

int find_set(int x, element * arr);
void link(int x, int y, element * arr);
void union_set(int x, int y, element * arr);
void CutPointGenerator(int problem_size, int cut_point_size, int* chosen_cut_point);
bool is_same_permutation(int *permu1, int *permu2, int ell);
int find_value_index_in_permu(int value, int *permu, int ell);
void copy_array(int *permu1, int *permu2, int ell);
void copy_array(double *permu1, double *permu2, int ell);




struct edgeee {
    int node_a;
    int node_a_index;
    int distance;

    double element_sum;
    double entropy;
    // double entropy_after_update;

    bool can_be_used; 
        // true: this edge can be used
        // false: this edge cannot be used
};

#endif