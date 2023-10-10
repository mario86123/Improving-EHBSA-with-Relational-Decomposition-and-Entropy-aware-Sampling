//
//  Population.cc
//  RankingEDAsCEC
//
//  Created by Josu Ceberio Uribe on 11/19/13.
//  Copyright (c) 2013 Josu Ceberio Uribe. All rights reserved.
//

#include "Population.h"
#include "Tools.h"


using std::istream;
using std::ostream;
using std::cerr;
using std::cout;
using std::endl;

/*
 * Constructor function.
 */
CPopulation::CPopulation(int pop_size, int offspring_size, int individual_size)
{
    m_pop_size = pop_size;
    m_offspring_size = offspring_size;
    m_size = pop_size + offspring_size;
    m_individuals.resize(m_size);
    
    //Initialize population with empty solutions
    for (int i=0;i<m_size;i++)
    {
        m_individuals[i]= new CIndividual(individual_size);
    }
    m_individual_size = individual_size;
}

/*
 * Destructor function.
 */
CPopulation::~CPopulation()
{
    for (int i=0;i<m_size;i++)
    {
        delete m_individuals[i];
    }
    m_individuals.clear();
}

/*
 * Function to add an individual to the population (back half)
 */
void CPopulation::AddToPopulation(int * genes, int index, long int fitness)
{
    m_individuals[m_pop_size+index]->SetGenes(genes);
    m_individuals[m_pop_size+index]->SetValue(fitness);
}

/*
 * Function to set an individual in a specific position of the population (front half)
 */
void CPopulation::SetToPopulation(int * genes, int index, long int fitness)
{
    m_individuals[index]->SetGenes(genes);
    m_individuals[index]->SetValue(fitness);
}


/*
 * Prints the current population.
 */
// print m_pop_size, not m_size = pop_size + offspring_size;
void CPopulation::Print()
{
    printf("\n   ");
    for(int i = 0; i < m_individual_size; ++i){
        printf("%1d ", i);
    }
    printf("\n---");
    for(int i = 0; i < m_individual_size; ++i){
        printf("---");
    }
    printf("\n");
    // printf("m_pop_size: %d\n", m_pop_size);
    // printf("m_offspring_size: %d\n", m_offspring_size);
    // printf("m_size: %d\n", m_size);
	// for(int i=0;i<m_pop_size; i++)
	for(int i=0;i<m_size; i++) {

		cout<<i<<": "<<m_individuals[i]<<", fitness: " << m_individuals[i]->Value();
        
        if (isPermutation(m_individuals[i]->Genes(), m_individual_size) == false) {
            cout << " xxx";
        }
        cout <<endl;
    }
}

/*
 * Determines if the individuals in the population are equal.
 */
bool CPopulation::Same(int size)
{
	int best=m_individuals[0]->Value();
	for(int i=1; i<size; i++)
	{
		if (m_individuals[i]->Value()!=best)
			return false;
	}
	return true;
}

/*
 * Calculates the average fitness of the first 'size' solutions in the population
 */
float CPopulation::AverageFitnessPopulation(int size)
{
	float result=0;
	for(int i=0;i<size;i++)
    {
        result+=m_individuals[i]->Value();
    }
    return result/size;
    
}

/*
 * Calculates the average fitness from "idx" to "idx + size"  solutions in the population
 */
float CPopulation::AverageFitnessPopulationSegment(int first_idx, int size)
{
	float result=0;
	for(int i = first_idx; i < first_idx + size; i++)
    {
        // cout << m_individuals[i]->Value() << endl;
        result+=m_individuals[i]->Value();
    }
    // cout << endl;
    return result/size;
    
}

float CPopulation::FindBestFitness(int size)
{
	float best_fitness = -100000000000;
	for(int i=0;i<size;i++)
    {
        if (best_fitness < m_individuals[i]->Value()){
            best_fitness = m_individuals[i]->Value();
        }
    }
    return best_fitness;
    
}

/*
 * Sorts the individuals in the population in decreasing order of the fitness value.
 * mode = 0 means sort only the population.
 * mode = 1 means sort population and offsprings population together.
 * mode = 2 means seperate population into two half pop, and sort seperately.
 * mode = 3 means seperate population into 1/4, 1/4, 1/4, 1/4 pop, and sort seperately.
 */
void CPopulation::SortPopulation(int mode)
{
    if (mode == 0)
        sort(m_individuals.begin(), m_individuals.begin()+m_pop_size, Better);
    else if (mode == 1)
        sort(m_individuals.begin(), m_individuals.end(), Better);
    else if (mode == 2) {
        sort(m_individuals.begin(), m_individuals.begin()+ m_pop_size/2, Better);
        sort(m_individuals.begin() + m_pop_size/2, m_individuals.begin() + m_pop_size, Better);    
    }
    else if (mode == 3) {
        sort(m_individuals.begin(),                      m_individuals.begin() + m_pop_size / 4,     Better);
        sort(m_individuals.begin() + m_pop_size / 4,     m_individuals.begin() + m_pop_size / 2,     Better);    
        sort(m_individuals.begin() + m_pop_size / 2,     m_individuals.begin() + m_pop_size * 3 / 4, Better);    
        sort(m_individuals.begin() + m_pop_size * 3 / 4, m_individuals.begin() + m_pop_size,         Better);    
    }
}

// Count how many chromosome in the population has that fitness
int CPopulation::CountFitnessInPopulation(double fitness)
{
    int result = 0;
    for(int i=0;i<m_pop_size;i++)
    {
        if (m_individuals[i]->Value() == fitness) {
            result += 1;
        }
    }
    return result;
}


double CPopulation::CalculatePopulationEntropy()
{
    double entropy = 0;
    int different_individual_count = 0;
    
    int *unique_chromosome_idx = new int[m_pop_size]();
    int *unique_chromosome_count = new int[m_pop_size]();


    for(int idx = 0; idx < m_pop_size; idx++) {


        bool hasSameChromosomeInPopulation = false;

        for(int j = 0; j < different_individual_count; j++) {

            if (m_individuals[idx]->Same(m_individuals[unique_chromosome_idx[j]]->Genes()) == true) {
                unique_chromosome_count[j] ++;
                hasSameChromosomeInPopulation = true;
                break;
            }
        }

        if (hasSameChromosomeInPopulation == false) {
            
            // add this chromosome idx to array
            unique_chromosome_idx[different_individual_count] = idx;
            unique_chromosome_count[different_individual_count] = 1;

            different_individual_count ++;
        }


    }

    // printf("unique_chromosome_count: ");
    // PrintArray(unique_chromosome_count, m_pop_size);
    entropy = calculate_entropy(unique_chromosome_count, m_pop_size);


    delete [] unique_chromosome_idx;
    delete [] unique_chromosome_count;


    return entropy;
}

// (including offspring) random shuffle individuals in the population
void CPopulation::RandomShuffle()
{
    random_shuffle(m_individuals.begin(), m_individuals.end());
}

// random shuffle individuals "only" in the population, not including offspring
void CPopulation::RandomShufflePopulation()
{
    random_shuffle(m_individuals.begin(), m_individuals.begin()+m_pop_size);
}


void CPopulation::RandomShufflePopulationSplit(int start_idx_1, int start_idx_2)
{
    vector<CIndividual*> firstQuarter(m_individuals.begin() + start_idx_1, m_individuals.begin() + start_idx_1 + m_pop_size / 4);
    vector<CIndividual*> secondQuarter(m_individuals.begin() + start_idx_2, m_individuals.begin() + start_idx_2 + m_pop_size / 4);
    
    vector<CIndividual*> tmp_vec;
    tmp_vec.reserve(m_pop_size / 2);
    tmp_vec.insert(tmp_vec.end(), firstQuarter.begin(), firstQuarter.end());
    tmp_vec.insert(tmp_vec.end(), secondQuarter.begin(), secondQuarter.end());

    random_shuffle(tmp_vec.begin(), tmp_vec.end());

    copy(tmp_vec.begin(), tmp_vec.begin() + firstQuarter.size(), m_individuals.begin() + start_idx_1);
    copy(tmp_vec.begin() + firstQuarter.size(), tmp_vec.end(), m_individuals.begin() + start_idx_2);
}
