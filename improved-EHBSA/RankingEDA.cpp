/*
 *  RankingEDA.cpp
 *  RankingEDAsCEC
 *
 *  Created by Josu Ceberio Uribe on 9/15/11.
 *  Copyright 2011 University of the Basque Country. All rights reserved.
 *
 */

#include "RankingEDA.h"
#include "MallowsModel.h"
#include "Cayley.h"
#include "Kendall.h"
#include "GeneralizedMallowsModel.h"
#include "BWModel.h"
#include "MRModel.h"
#include "AMRModel.h"
#include "FMSTModel.h"
#include "MST2.h"
#include "MST2ME.h"
#include "MST3.h"
#include "MST4.h"
#include "SoftMST.h"
#include "LLMST.h"
#include "RMST.h"
#include "SRMST.h"
#include "AMST.h"
#include "MSTModel.h"
#include "MSTME.h"
#include "EHBSAWO.h"
#include "EHBSAWT.h"
#include "NHBSAWO.h"
#include "NHBSAWT.h"
#include "TEMP.h"
#include "RRBOP.h"
#include "FCRBOP.h"
#include "RBOP.h"
#include "PRBOP.h"
#include "ARBOP.h"
#include "NRBOP.h"
#include "CHI.h"
#include "ACHIGN.h"
#include "ACNRBOP.h"
#include "AGNRBOP.h"
#include "PREPRUNE.h"
#include "FOUR.h"
#include "ANRBOP.h"
#include "SRBOP.h"
#include "BM.h"
#include "FERBOP.h"
#include "NNHBSAWT.h"
#include "EEHBSAWT.h"
#include "ENEHBSAWT.h"
#include "SNHBSAWT.h"
#include "MRWT.h"
#include "time.h"

/*
 * The constructor.
 */
RankingEDA::RankingEDA(PBP * problem, int problem_size, int poplation_size, long int max_evaluations, double b_ratio, double cut_point_count, int previous_sampled_reference_count, int number_of_edge, char* model_type, char* metric_type, int inverse, int seed, char * result_file_name)
{
    //1. standard initializations
    m_problem=problem;
    m_problem_size=problem_size;
    m_max_evaluations=max_evaluations;
    m_evaluations=0;
    m_convergence_evaluations=0;
    m_best= new CIndividual(problem_size);
    m_best->SetValue(MIN_LONG_INTEGER);
    
    // total population size = m_pop_size + m_offspring_size = 2 * m_sel_size
    m_pop_size = poplation_size;
    m_sel_size = m_pop_size;
    m_offspring_size = m_pop_size;
    m_b_ratio = b_ratio;

    // m_pop_size == m_sel_size == m_offspring_size
    // cout << "m_pop_size: " << m_pop_size << endl;
    // cout << "m_sel_size: " << m_sel_size << endl;
    // cout << "m_offspring_size: " << m_offspring_size << endl;

    memcpy(m_metric_type,metric_type,sizeof(char)*10);
    memcpy(m_model_type,model_type,sizeof(char)*10);
    
    m_inverse=inverse;
    m_seed = seed;
    m_number_of_edge = number_of_edge;
    m_cut_point_count = cut_point_count;
    m_result_file_name = result_file_name;
    
    //2. initialize the population
    
    // where is 0 ??
    // why add 35 at tail??
    // int tmp[7] = {3, 2, 0, 4, 1, 6, 5};
    // int tmp[16] = {4, 5, 9, 11, 6, 15, 7, 2, 3, 13, 14, 8, 0, 1, 12, 10};
    // int tmp[60] = {46-1, 60-1, 54-1, 58-1, 22-1, 33-1, 52-1, 39-1, 36-1, 35-1, 37-1, 38-1, 40-1, 1-1, 41-1, 15-1, 21-1, 29-1, 24-1, 23-1, 30-1, 3-1, 44-1, 25-1, 26-1, 27-1, 17-1, 28-1, 20-1, 13-1, 18-1, 14-1, 16-1, 19-1, 12-1, 34-1, 10-1, 8-1, 43-1, 42-1, 50-1, 55-1, 49-1, 53-1, 32-1, 31-1, 6-1, 4-1, 5-1, 7-1, 2-1, 47-1, 45-1, 51-1, 57-1, 56-1, 48-1, 11-1, 9-1, 59-1};
    // int tmp[60] = {46-1, 60-1, 54-1, 58-1, 22-1, 33-1, 52-1, 39-1, 36-1, 35-1, 37-1, 38-1, 40-1, 1-1, 44-1, 55-1, 41-1, 15-1, 21-1, 29-1, 24-1, 25-1, 23-1, 26-1, 27-1, 17-1, 13-1, 30-1, 3-1, 28-1, 20-1, 59-1, 14-1, 12-1, 34-1, 18-1, 19-1, 10-1, 16-1, 8-1, 42-1, 43-1, 53-1, 32-1, 31-1, 4-1, 5-1, 7-1, 2-1, 47-1, 45-1, 51-1, 48-1, 50-1, 49-1, 57-1, 56-1, 11-1, 6-1, 9-1};
    // double tmp_fitness = m_problem->Evaluate(tmp);
    // cout <<"tmp fitness: "<< tmp_fitness << endl;
    // cout <<"tmp fitness: "<< endl;
    
//     double tmp_1_fitness = m_problem->Evaluate(tmp_1);
//     cout <<"tmp_1 fitness: "<< tmp_1_fitness << endl;

    m_population= new CPopulation(m_pop_size, m_offspring_size, m_problem_size);
    int * genes= new int[m_problem_size];
    for(int i=0; i<m_pop_size; i++)
	{
		//Create random individual
		GenerateRandomPermutation(genes,m_problem_size);
        if (m_inverse)
            m_population->SetToPopulation(genes, i, m_problem->EvaluateInv(genes));
        else
            m_population->SetToPopulation(genes, i, m_problem->Evaluate(genes));
        m_evaluations++;
    }
    // m_population->Print();

    // m_population->SortPopulation(0);
    // cout<<""<<m_population->m_individuals[0]->Value()<<" , "<<m_evaluations<<" , "<<m_max_evaluations-m_evaluations<<endl;
    delete [] genes;
    

    //3. Build model structures
    if (((string)model_type)=="M"){
        m_model=new CMallowsModel(m_problem_size, m_sel_size, metric_type);
    }
    else if (((string)model_type)=="GM")
    {
        m_model=new CGeneralizedMallowsModel(m_problem_size, m_sel_size, metric_type);
    }
    else if (((string)model_type)=="BW")
    {
        m_model=new CBWModel(m_problem_size, m_sel_size, b_ratio);
    }
    else if (((string)model_type)=="MR")
    {
        m_model=new CMRModel(m_problem_size, m_sel_size, b_ratio);
    }
    else if (((string)model_type)=="AMR") // Accumulate MR Model
    {
        m_model=new CAMRModel(m_problem_size, m_sel_size, b_ratio, previous_sampled_reference_count, seed);
    }
    else if (((string)model_type)=="MST2") // real MST
    {
        m_model=new CMST2(m_problem_size, m_sel_size, b_ratio, m_result_file_name);
    }
    else if (((string)model_type)=="MST2ME") // real MST with multiple edges
    {
        m_model=new CMST2ME(m_problem_size, m_sel_size, b_ratio, m_result_file_name);
    }
    else if (((string)model_type)=="MST3") // 90% edge relation ( problem_size / 10 trees )
    {
        m_model=new CMST3(m_problem_size, m_sel_size, b_ratio, m_result_file_name);
    }
    else if (((string)model_type)=="MST4") // 80% edge relation ( problem_size / 5 trees )
    {
        m_model=new CMST4(m_problem_size, m_sel_size, b_ratio, m_result_file_name);
    }
    else if (((string)model_type)=="SMST") // Softmax MST
    {
        m_model=new CSoftMST(m_problem_size, m_sel_size, b_ratio, m_result_file_name);
    }
    else if (((string)model_type)=="LLMST") // Loss Less MST Model
    {
        m_model=new CLLMST(m_problem_size, m_sel_size, b_ratio, number_of_edge, m_result_file_name);
    }
    else if (((string)model_type)=="RMST") // Random fake MST Model
    {
        m_model=new CRMST(m_problem_size, m_sel_size, b_ratio, number_of_edge, m_result_file_name);
    }
    else if (((string)model_type)=="SRMST") // Softmax Random fake MST Model
    {
        m_model=new CSRMST(m_problem_size, m_sel_size, b_ratio, number_of_edge, m_result_file_name);
    }
    else if (((string)model_type)=="AMST") // Adaptive fake MST Model
    {
        m_model=new CAMST(m_problem_size, m_sel_size, b_ratio, m_result_file_name);
    }
    else if (((string)model_type)=="FMST") // fast MST Model
    {
        m_model=new CFMSTModel(m_problem_size, m_sel_size, b_ratio);
    }
    else if (((string)model_type)=="MST") // fake MST Model
    {
        m_model=new CMSTModel(m_problem_size, m_sel_size, b_ratio);
    }
    else if (((string)model_type)=="MSTME") // MST Model with multiople edges
    {
        m_model=new CMSTME(m_problem_size, m_sel_size, b_ratio);
    }
    else if (((string)model_type)=="EO") // EHBSA WO
    {
        m_model=new CEHBSAWO(m_problem_size, m_sel_size, b_ratio);
    }
    else if (((string)model_type)=="ET") // EHBSA WT
    {
        m_model=new CEHBSAWT(m_problem_size, m_sel_size, b_ratio, cut_point_count);
    }
    else if (((string)model_type)=="ENET") // entropy EHBSA WT
    {
        m_model=new CENEHBSAWT(m_problem_size, m_sel_size, b_ratio, cut_point_count);
    }
    else if (((string)model_type)=="NO") // NHBSA WO
    {
        m_model=new CNHBSAWO(m_problem_size, m_sel_size, b_ratio);
    }
    else if (((string)model_type)=="NT") // NHBSA WT
    {
        m_model=new CNHBSAWT(m_problem_size, m_sel_size, b_ratio, cut_point_count);
    }
    else if (((string)model_type)=="TEMP") // NHBSA WT
    {
        m_model=new CTEMP(m_problem_size, m_sel_size, b_ratio, cut_point_count);
    }
    else if (((string)model_type)=="RBOP") // NHBSA WT
    {
        m_model=new CRBOP(m_problem_size, m_sel_size, b_ratio, cut_point_count);
    }
    else if (((string)model_type)=="FCRBOP") // NHBSA WT
    {
        m_model=new CFCRBOP(m_problem_size, m_sel_size, b_ratio, cut_point_count);
    }
    else if (((string)model_type)=="SRBOP") // NHBSA WT
    {
        m_model=new CSRBOP(m_problem_size, m_sel_size, b_ratio, cut_point_count);
    }
    else if (((string)model_type)=="PRBOP") // NHBSA WT
    {
        m_model=new CPRBOP(m_problem_size, m_sel_size, b_ratio, cut_point_count);
    }
    else if (((string)model_type)=="NRBOP") // NHBSA WT
    {
        m_model=new CNRBOP(m_problem_size, m_sel_size, b_ratio, cut_point_count);
    }
    else if (((string)model_type)=="ACNRBOP") // Adaptive Cut point New RBOP
    {
        // m_model=new CACNRBOP(m_problem_size, m_sel_size, b_ratio);
        m_model_1=new CACNRBOP(m_problem_size, m_sel_size, b_ratio);
        m_model_2=new CACNRBOP(m_problem_size, m_sel_size, b_ratio);
    }
    else if (((string)model_type)=="CHI") // Adaptive Cut point New RBOP
    {
        // m_model=new CCHI(m_problem_size, m_sel_size, b_ratio);
        m_model_1=new CCHI(m_problem_size, m_sel_size, b_ratio);
        m_model_2=new CCHI(m_problem_size, m_sel_size, b_ratio);
    }
    else if (((string)model_type)=="ACHIGN") // Adaptive Cut point New RBOP
    {
        // m_model=new CACHIGN(m_problem_size, m_sel_size, b_ratio);
        m_model_1=new CACHIGN(m_problem_size, m_sel_size, b_ratio);
        m_model_2=new CACHIGN(m_problem_size, m_sel_size, b_ratio);
    }
    else if (((string)model_type)=="AGNRBOP") // Adaptive Cut point New RBOP
    {
        // m_model=new CAGNRBOP(m_problem_size, m_sel_size, b_ratio);
        m_model_1=new CAGNRBOP(m_problem_size, m_sel_size, b_ratio);
        m_model_2=new CAGNRBOP(m_problem_size, m_sel_size, b_ratio);
    }
    else if (((string)model_type)=="PREPRUNE") // Adaptive Cut point New RBOP
    {
        // m_model=new CPREPRUNE(m_problem_size, m_sel_size, b_ratio);
        m_model_1=new CPREPRUNE(m_problem_size, m_sel_size, b_ratio);
        m_model_2=new CPREPRUNE(m_problem_size, m_sel_size, b_ratio);
    }
    else if (((string)model_type)=="FOUR") // Adaptive Cut point New RBOP
    {
        // m_model=new CACNRBOP(m_problem_size, m_sel_size, b_ratio);
        m_model_1=new CFOUR(m_problem_size, m_sel_size, b_ratio);
        m_model_2=new CFOUR(m_problem_size, m_sel_size, b_ratio);
    }
    else if (((string)model_type)=="ANRBOP") // NHBSA WT
    {
        m_model=new CANRBOP(m_problem_size, m_sel_size, b_ratio);
    }
    else if (((string)model_type)=="ARBOP") // NHBSA WT
    {
        int node_edge_ratio_percent = number_of_edge;
        m_model=new CARBOP(m_problem_size, m_sel_size, b_ratio, cut_point_count, node_edge_ratio_percent);
    }
    else if (((string)model_type)=="BM") // NHBSA WT
    {
        m_model=new CBM(m_problem_size, m_sel_size, b_ratio, cut_point_count);
    }
    else if (((string)model_type)=="FERBOP") // NHBSA WT
    {
        m_model=new CFERBOP(m_problem_size, m_sel_size, b_ratio, cut_point_count);
    }
    else if (((string)model_type)=="RRBOP") // NHBSA WT
    {
        m_model=new CRRBOP(m_problem_size, m_sel_size, b_ratio, cut_point_count);
    }
    else if (((string)model_type)=="NNT") // new NHBSA WT
    {
        m_model=new CNNHBSAWT(m_problem_size, m_sel_size, b_ratio, cut_point_count, m_population);
    }
    else if (((string)model_type)=="EET") // new EHBSA WT
    {
        m_model=new CEEHBSAWT(m_problem_size, m_sel_size, b_ratio, cut_point_count, m_population);
    }
    else if (((string)model_type)=="SNT") // Sequentail NHBSA WT
    {
        m_model=new CSNHBSAWT(m_problem_size, m_sel_size, b_ratio, cut_point_count);
    }

    else if (((string)model_type)=="MRWT") // MRModel With Template
    {
        m_model=new CMRWT(m_problem_size, m_sel_size, b_ratio, cut_point_count);
    }
}

/*
 * The destructor. It frees the memory allocated..
 */
RankingEDA::~RankingEDA()
{
    // printf("m_pop_size: %d\n", m_pop_size);
    delete m_best;
    delete m_population;

    if (((string)m_model_type)=="ACNRBOP" || ((string)m_model_type)=="CHI" || ((string)m_model_type)=="FOUR" || ((string)m_model_type)=="AGNRBOP" || ((string)m_model_type)=="PREPRUNE" || ((string)m_model_type)=="ACHIGN") {
        delete m_model_1;
        delete m_model_2;
    }
    else {
        delete m_model;
    }
}

/*
 * Running function
 */
int RankingEDA::Run(char* problem_type){
        // printf("Running!!\n");

    //cout<<"Running..."<<endl;
    //variable initializations.
    int i;
    long newScore = m_population->m_individuals[0]->Value();
    int * genes= new int[m_problem_size];
    int iterations=1;
    float rate=0.001;
    float population_avg_fitness = -100000000;
    int better_cut_point_count = -1;


    ofstream best_fitness_file;
    ofstream avg_fitness_file;


    string result_name("result");
    string best_fitness_file_name = m_result_file_name;
    best_fitness_file_name.replace(best_fitness_file_name.find(result_name), result_name.length(), "best_fitness");


    string avg_fitness_file_name = m_result_file_name;
    avg_fitness_file_name.replace(avg_fitness_file_name.find(result_name), result_name.length(), "avg_fitness");

    
    best_fitness_file.open(best_fitness_file_name);
    avg_fitness_file.open(avg_fitness_file_name);


    int edge_count = m_problem_size / 2;
    int edge_count_update_constant = m_problem_size / 10;


    // === cut point initialization === //
    // cut_point_percentage = 10 ---> means that there are (ell * 10 %) genes that need to resample
    int need_to_sample_genes_count = m_problem_size * 30 / 100;
    int need_to_sample_genes_count_update_constant = -m_problem_size * 15 / 100;
    int change_gene_num_NFE = 10 * m_problem_size * m_problem_size;
    int old_gene_num_evaluations = m_evaluations;


    // === end of cut point initialization === //


    // === cut point initialization === //
    // cut_point_percentage = 10 ---> means that there are (ell * 10 %) genes that need to resample
    double chi_std = 0.0;
    double chi_std_update_constant = 0.5;
    int change_chi_square_NFE = 100 * m_problem_size * m_problem_size;
    int old_chi_square_evaluations = m_evaluations;

    // int chi_std_change_NFE = 50 * m_problem_size * m_problem_size;

    // === end of cut point initialization === //
    int change_NFE = 10 * m_problem_size * m_problem_size;



    //EDA iteration. Stopping criterion is the maximum number of evaluations performed
    //               or reach global optimal
    // while (m_evaluations<m_max_evaluations && newScore > population_avg_fitness) {
    // while ( m_evaluations<m_max_evaluations && !(m_population->Same(m_pop_size)) ){
    // while (m_evaluations<m_max_evaluations && population_avg_fitness < m_problem_size) {
    // int relation_count_used_in_sampling = m_problem_size / m_cut_point_count;
    // int *improve_count_array;
    // improve_count_array = new int[relation_count_used_in_sampling]();

    // printf("here!!\n");

    // ===== random sample index start here ===== //
    int template_idx_arr[m_offspring_size];
	for (int i = 0; i < m_offspring_size; ++i) {
		template_idx_arr[i] = i;
	}
	// random_shuffle(template_idx_arr, template_idx_arr + m_offspring_size);
    // PrintArray(template_idx_arr, m_offspring_size);
    // ===== end of random sample index ===== //


    int used_distance_count_arr[m_problem_size];
	for (int i = 0; i < m_problem_size; ++i) {
		used_distance_count_arr[i] = 0;
	}





    while (m_evaluations<m_max_evaluations) {        
        
        // time_t start, end;
        // start = clock();


        if ( ((string)m_model_type)=="NNT" || ((string)m_model_type)=="EET") {

            // m_population->Print();

            int template_idx = rand() % m_pop_size;

            // template == m_population->m_individuals[i]->Genes()
            m_model->Sample(genes, m_population->m_individuals[template_idx]->Genes());
            
            long new_sample_fitness = m_problem->Evaluate(genes);
            long template_fitness = m_population->m_individuals[template_idx]->Value();

            // ----- if template better than new sample ----- //
            if ( template_fitness > new_sample_fitness) {
                m_population->AddToPopulation(genes, template_idx, new_sample_fitness);
            }

            // ----- else if template worse than new sample ----- //
            else {
                // update NHM:
                    // --- update NHM (original chromosome, new sample chromosome)

                m_model->Learn(genes, m_population->m_individuals[template_idx]->Genes());

                m_population->AddToPopulation(m_population->m_individuals[template_idx]->Genes(), template_idx, template_fitness);
                m_population->SetToPopulation(genes, template_idx, new_sample_fitness);

                if (newScore < new_sample_fitness) {
                    newScore = new_sample_fitness;
                }
            }
            // ----- end if ----- //

        
            m_evaluations++;

        }


        // === update models after sample one population === //
        else if (((string)m_model_type) != "ACNRBOP" && ((string)m_model_type) != "PREPRUNE" && ((string)m_model_type) != "CHI" && ((string)m_model_type) != "ACHIGN" && ((string)m_model_type) != "FOUR" && ((string)m_model_type) != "AGNRBOP") {


            // === learn model
            m_model->Learn(m_population, m_sel_size); // m_sel_size == m_pop_size

            // === sample the model.

            // edge model "without" template
            if (((string)m_model_type)=="EO" || \
                ((string)m_model_type)=="NO") {
                for (i=0;i< m_offspring_size && m_evaluations<m_max_evaluations;i++){
                    
                    // template == m_population->m_individuals[i]->Genes()
                    m_model->Sample(genes);
                    
                    long new_sample_fitness = m_problem->Evaluate(genes);
                    long template_fitness = m_population->m_individuals[i]->Value();

                    // if template better than new sample
                    if ( template_fitness > new_sample_fitness) {
                        m_population->AddToPopulation(genes, i, new_sample_fitness);
                    }
                    //if template worse than new sample
                    else {
                        m_population->AddToPopulation(m_population->m_individuals[i]->Genes(), i, template_fitness);
                        m_population->SetToPopulation(genes, i, new_sample_fitness);
                    }
                    m_evaluations++;
                }
                // sort only population
                m_population->SortPopulation(0);
            }

            // models using template
            else if ( ((string)m_model_type)=="ET" || \
                    ((string)m_model_type)=="ENET" || \
                    ((string)m_model_type)=="NT" || \
                    ((string)m_model_type)=="SNT" || \
                    ((string)m_model_type)=="MRWT" || \
                    ((string)m_model_type)=="MRWFT" || \
                    ((string)m_model_type)=="RRBOP" || \
                    ((string)m_model_type)=="RBOP" || \
                    ((string)m_model_type)=="FCRBOP" || \
                    ((string)m_model_type)=="SRBOP" || \
                    ((string)m_model_type)=="ARBOP" || \
                    ((string)m_model_type)=="NRBOP" || \
                    ((string)m_model_type)=="PRBOP" || \
                    ((string)m_model_type)=="FERBOP" || \
                    ((string)m_model_type)=="TEMP") 
                    {
                        
                // m_population->Print();
                for (i=0;i< m_offspring_size && m_evaluations<m_max_evaluations;i++){
                    // cout << "i: " << i << endl;
                    // cout << "i: " << i << endl;
                    // cout << "i: " << i << endl;
                    int template_idx = rand() % m_pop_size;

                    // template == m_population->m_individuals[i]->Genes()
                    m_model->Sample(genes, m_population->m_individuals[template_idx]->Genes());
                    
                    
                    long template_fitness = m_population->m_individuals[template_idx]->Value();
                    long new_sample_fitness = 0;

                    
                    // === new way to count NFE === //
                    // if (is_same_permutation(genes, m_population->m_individuals[template_idx]->Genes(), m_problem_size)) {
                    //     //  do not add NFE
                    //     new_sample_fitness = template_fitness;
                    // }
                    // else {
                    //     new_sample_fitness = m_problem->Evaluate(genes);
                    //     m_evaluations++;
                    // }
                    // === end of new way to count NFE === //


                    new_sample_fitness = m_problem->Evaluate(genes);
                    m_evaluations++;
                    
                    // ----- if template better than new sample ----- //
                    if ( template_fitness > new_sample_fitness) {
                        m_population->AddToPopulation(genes, template_idx, new_sample_fitness);
                    }

                    // ----- else if template worse than new sample ----- //
                    else {
                        m_population->AddToPopulation(m_population->m_individuals[template_idx]->Genes(), template_idx, template_fitness);
                        m_population->SetToPopulation(genes, template_idx, new_sample_fitness);
                    }

                }
                // sort only population
                m_population->SortPopulation(0);
            }

            // --- AMST start --- //
            else if ( ((string)m_model_type)=="AMST" ) 
            {
                // cout << "edge count: " << edge_count << ", edge count + edge_count_update_constant: " << edge_count + edge_count_update_constant << endl;
                if (edge_count > m_problem_size - 1 - edge_count_update_constant) { 
                    edge_count = m_problem_size - 1; 
                    edge_count_update_constant = - m_problem_size / 10;
                    // cout << "bigger, " << edge_count << endl;
                    // cout << "[update] edge count: " << edge_count << ", edge count + edge_count_update_constant: " << edge_count + edge_count_update_constant << endl;
                }
                else if (edge_count < - edge_count_update_constant) { 
                    edge_count = 0;
                    edge_count_update_constant = m_problem_size / 10;
                    // cout << "smaller, " << edge_count << endl;
                    // cout << "[update] edge count: " << edge_count << ", edge count + edge_count_update_constant: " << edge_count + edge_count_update_constant << endl;
                }

                for (i = 0; i < m_offspring_size / 2 && m_evaluations < m_max_evaluations; i++) {
                    
                    m_model->Sample(genes, edge_count);
                    if (m_inverse)
                        m_population->AddToPopulation(genes, i, m_problem->EvaluateInv(genes));
                    else
                        m_population->AddToPopulation(genes, i, m_problem->Evaluate(genes));
                    m_evaluations++;
                }

                //calculate avg fitness
                double original_edge_count_avg_fitness = m_population->AverageFitnessPopulationSegment(m_pop_size, m_pop_size / 2);

                for (i = m_offspring_size / 2; i < m_offspring_size && m_evaluations < m_max_evaluations; i++) {
                    
                    m_model->Sample(genes, edge_count + edge_count_update_constant);

                    if (m_inverse)
                        m_population->AddToPopulation(genes, i, m_problem->EvaluateInv(genes));
                    else
                        m_population->AddToPopulation(genes, i, m_problem->Evaluate(genes));
                    m_evaluations++;
                }
                //calculate avg fitness
                double new_edge_count_avg_fitness = m_population->AverageFitnessPopulationSegment(m_pop_size + m_pop_size / 2, m_pop_size / 2);
                // cout << "new_edge_count_avg_fitness: " << new_edge_count_avg_fitness << endl<< endl;

                // todo: compare avg fitness and update edge_count
                if (original_edge_count_avg_fitness > new_edge_count_avg_fitness) {
                    // if (edge_count != 0 && edge_count != m_problem_size - 1 ) {
                    edge_count_update_constant = edge_count_update_constant * (-1);
                    // }
                    // else { // edge_count == 0 || edge_count == m_problem_size - 1 --> do nothing }
                }
                // original_edge_count_avg_fitness < new_edge_count_avg_fitness
                else {
                    // if (edge_count != 1 && edge_count != m_problem_size - 2 ) {
                    edge_count = edge_count + edge_count_update_constant;
                    // }
                    // else { // edge_count == 1 || edge_count == m_problem_size - 2 --> do nothing }
                }
                // sort population and offspring together (elitism)
                m_population->SortPopulation(1);
            }
            // --- AMST end --- //



            // --- Adaptive New RBOP start --- //
            else if ( ((string)m_model_type)=="ANRBOP" ) 
            {

                if (need_to_sample_genes_count + need_to_sample_genes_count_update_constant > m_problem_size/2) { 
                    need_to_sample_genes_count = m_problem_size / 2; 
                    need_to_sample_genes_count_update_constant = - 2;
                }
                else if (need_to_sample_genes_count + need_to_sample_genes_count_update_constant < 2) { 
                    need_to_sample_genes_count = 2;
                    need_to_sample_genes_count_update_constant = 2;
                }

                // generate a random integer array to decide index to be used by two different cut point cut percentage
            	// === usage: template_idx[i] means random index === //
                random_shuffle(template_idx_arr, template_idx_arr + m_offspring_size);

                //calculate avg fitness
                double original_cut_point_count_fitness_sum = 0;
                double new_cut_point_count_fitness_sum = 0;
                
                int original_cut_point_count_improve_count = 0;
                int new_cut_point_count_improve_count = 0;


                // === original_cut_point_count === /
                for (i = 0; i < m_offspring_size / 2 && m_evaluations < m_max_evaluations; i++) {
                    
                    int template_idx = template_idx_arr[i];
                    double template_fitness = m_population->m_individuals[template_idx]->Value();

                    // sample using template index chromosome
                    m_model->Sample(genes, m_population->m_individuals[template_idx]->Genes(), need_to_sample_genes_count);

                    // calculate fitness
                    double new_sample_fitness = m_problem->Evaluate(genes);
                    m_evaluations++;

                    // === add fitness to average fitness

                    // determine whether replace the template chromosome in the population
                    if (template_fitness > new_sample_fitness) {
                        m_population->AddToPopulation(genes, template_idx, new_sample_fitness);
                    }
                    else {
                        m_population->AddToPopulation(m_population->m_individuals[template_idx]->Genes(), template_idx, template_fitness);
                        m_population->SetToPopulation(genes, template_idx, new_sample_fitness);
                        original_cut_point_count_improve_count++;
                    }
                    original_cut_point_count_fitness_sum += m_population->m_individuals[template_idx]->Value();
                }
                // === end of original_cut_point_count === /



                // === new_cut_point_count_fitness === /

                for (i = m_offspring_size / 2; i < m_offspring_size && m_evaluations < m_max_evaluations; i++) {
                    
                    int template_idx = template_idx_arr[i];
                    double template_fitness = m_population->m_individuals[template_idx]->Value();

                    // sample using template index chromosome
                    m_model->Sample(genes, m_population->m_individuals[template_idx]->Genes(), need_to_sample_genes_count + need_to_sample_genes_count_update_constant);

                    // calculate fitness
                    double new_sample_fitness = m_problem->Evaluate(genes);
                    m_evaluations++;

                    // === add fitness to average fitness

                    // determine whether replace the template chromosome in the population
                    if (template_fitness > new_sample_fitness) {
                        m_population->AddToPopulation(genes, template_idx, new_sample_fitness);
                    }
                    else {
                        m_population->AddToPopulation(m_population->m_individuals[template_idx]->Genes(), template_idx, template_fitness);
                        m_population->SetToPopulation(genes, template_idx, new_sample_fitness);
                        new_cut_point_count_improve_count++;
                    }
                    new_cut_point_count_fitness_sum += m_population->m_individuals[template_idx]->Value();
                }
                // === end of adjust cut point percentage using improvment count === //

                // === adjust cut point percentage using avg_fitness === //
                if (original_cut_point_count_fitness_sum > new_cut_point_count_fitness_sum) {
                    need_to_sample_genes_count_update_constant = need_to_sample_genes_count_update_constant * (-1);
                }
                else { // original_edge_count_avg_fitness < new_edge_count_avg_fitness
                    need_to_sample_genes_count = need_to_sample_genes_count + need_to_sample_genes_count_update_constant;
                }
                // === end of adjust cut point percentage using avg_fitness === //
            }
            // --- Adaptive New RBOP end --- //


            // --- others --- //
            else {
                for (i=0;i< m_offspring_size && m_evaluations<m_max_evaluations;i++){
                    
                    m_model->Sample(genes);

                    if (m_inverse)
                        m_population->AddToPopulation(genes, i, m_problem->EvaluateInv(genes));
                    else
                        m_population->AddToPopulation(genes, i, m_problem->Evaluate(genes));
                    m_evaluations++;
                }
                //update the model.
                
                // sort population and offspring together (elitism)
                m_population->SortPopulation(1);
            }
            // --- end of others --- //


            newScore=m_population->m_individuals[0]->Value();
        }
        // === end of update models after sample one population === //
        // === end of else if (((string)m_model_type) != "ACNRBOP") {


        else if (((string)m_model_type) == "ACNRBOP" || ((string)m_model_type) == "CHI") {
            
            if ( m_evaluations < (30 * m_problem_size * m_problem_size) ) {

                m_model_1->Learn(m_population, 1, used_distance_count_arr); // m_sel_size == m_pop_size
                m_model_2->Learn(m_population, 2, used_distance_count_arr); // m_sel_size == m_pop_size

                for (i = 0; i < m_offspring_size / 2 && m_evaluations < m_max_evaluations; i++) {

                    int template_idx = rand() % (m_pop_size / 2);
                    // cout << "model 1 template idx: " << template_idx << endl;

                    // === model 1 (cut point == 3) === //
                    long template_fitness = m_population->m_individuals[template_idx]->Value();
                    
                    m_model_1->Sample(genes, m_population->m_individuals[template_idx]->Genes(), 3);
                    long new_sample_fitness =  m_problem->Evaluate(genes);
                    m_evaluations++;

                    if ( template_fitness > new_sample_fitness) {
                        m_population->AddToPopulation(genes, template_idx, new_sample_fitness);
                    }
                    // ----- else if template worse than new sample ----- //
                    else {
                        m_population->AddToPopulation(m_population->m_individuals[template_idx]->Genes(), template_idx, template_fitness);
                        m_population->SetToPopulation(genes, template_idx, new_sample_fitness);
                    }
                    // === end of model 1 === //

                    template_idx += (m_pop_size / 2);
                    // cout << "model 2 template idx: " << template_idx << endl;
                    // === model 2 (cut point == 10) === //
                    template_fitness = m_population->m_individuals[template_idx]->Value();
                    m_model_2->Sample(genes, m_population->m_individuals[template_idx]->Genes(), 10);
                    new_sample_fitness =  m_problem->Evaluate(genes);
                    m_evaluations++;

                    if ( template_fitness > new_sample_fitness) {
                        m_population->AddToPopulation(genes, template_idx, new_sample_fitness);
                    }

                    // ----- else if template worse than new sample ----- //
                    else {
                        m_population->AddToPopulation(m_population->m_individuals[template_idx]->Genes(), template_idx, template_fitness);
                        m_population->SetToPopulation(genes, template_idx, new_sample_fitness);
                    }
                    // === end of model 2 === //
                }
                // end of population size for loop
                m_population->SortPopulation(2);
                // m_population->Print();

                if (m_population->m_individuals[0]->Value() > m_population->m_individuals[m_pop_size/2]->Value()) {
                    newScore=m_population->m_individuals[0]->Value();
                }
                else {
                    newScore=m_population->m_individuals[m_pop_size/2]->Value();
                }

                // better_cut_point_count
                // cout << "m_pop_size/2: " << m_pop_size/2 << ", m_pop_size: " << m_pop_size << endl;
                // cout << "c = 3 avg: " << m_population->AverageFitnessPopulationSegment(0, m_pop_size/2) << endl;
                // cout << "c = 10 avg: " << m_population->AverageFitnessPopulationSegment(m_pop_size/2, m_pop_size/2) << endl;
                if (m_population->AverageFitnessPopulationSegment(0, m_pop_size/2) > m_population->AverageFitnessPopulationSegment(m_pop_size/2, m_pop_size/2)) {
                    better_cut_point_count = 3;
                }
                else {
                    better_cut_point_count = 10;
                }
                // cout << "better_cut_point_count: " << better_cut_point_count << endl<< endl;

            } // end of if (m_evaluations < m_max_evaluations / 20)

            else { // if (m_evaluations >= m_max_evaluations / 20)
                m_model_1->Learn(m_population, 3, used_distance_count_arr); // m_sel_size == m_pop_size

                for (i = 0; i < m_offspring_size && m_evaluations < m_max_evaluations; i++) {

                    int template_idx = rand() % (m_pop_size);


                    // === model 1 (cut point == 3) === //
                    long template_fitness = m_population->m_individuals[template_idx]->Value();
                    
                    m_model_1->Sample(genes, m_population->m_individuals[template_idx]->Genes(), better_cut_point_count);
                    long new_sample_fitness =  m_problem->Evaluate(genes);
                    m_evaluations++;


                    if ( template_fitness > new_sample_fitness) {
                        m_population->AddToPopulation(genes, template_idx, new_sample_fitness);
                    }

                    // ----- else if template worse than new sample ----- //
                    else {
                        m_population->AddToPopulation(m_population->m_individuals[template_idx]->Genes(), template_idx, template_fitness);
                        m_population->SetToPopulation(genes, template_idx, new_sample_fitness);
                    }
                    // === end of model 1 === //
                }
                // end of population size for loop
                m_population->SortPopulation(0);
                newScore=m_population->m_individuals[0]->Value();
            }
            // end of if (m_evaluations < m_max_evaluations / 20)
        }
        // end of else if (((string)m_model_type) == "ACNRBOP")

        
        else if (((string)m_model_type) == "AGNRBOP" || ((string)m_model_type) == "PREPRUNE" ) {


            if (need_to_sample_genes_count + need_to_sample_genes_count_update_constant > m_problem_size/2) { 
                need_to_sample_genes_count = m_problem_size / 2; 
                need_to_sample_genes_count_update_constant = - abs(need_to_sample_genes_count_update_constant);
            }
            else if (need_to_sample_genes_count + need_to_sample_genes_count_update_constant < abs(need_to_sample_genes_count_update_constant)) { 
                need_to_sample_genes_count = abs(need_to_sample_genes_count_update_constant);
                need_to_sample_genes_count_update_constant = abs(need_to_sample_genes_count_update_constant);
            }
            
            int old_evaluations = m_evaluations;

            // while (m_evaluations - old_evaluations < 30 * m_problem_size * m_problem_size) {
            while (m_evaluations - old_evaluations < change_NFE && m_evaluations<m_max_evaluations) {

                // m_model_1->Learn(m_population, 1); // m_sel_size == m_pop_size
                // m_model_2->Learn(m_population, 2); // m_sel_size == m_pop_size
                
                m_model_1->Learn(m_population, 3); // m_sel_size == m_pop_size
                m_model_2->Learn(m_population, 3); // m_sel_size == m_pop_size
                for (i = 0; i < m_offspring_size / 2; i++) {

                    int template_idx = rand() % (m_pop_size / 2);
                    // cout << "model 1 template idx: " << template_idx << endl;

                    // === model 1 (cut point == 3) === //
                    long template_fitness = m_population->m_individuals[template_idx]->Value();
                    
                    m_model_1->Sample(genes, m_population->m_individuals[template_idx]->Genes(), need_to_sample_genes_count);
                    long new_sample_fitness =  m_problem->Evaluate(genes);
                    m_evaluations++;

                    if ( template_fitness > new_sample_fitness) {
                        m_population->AddToPopulation(genes, template_idx, new_sample_fitness);
                    }
                    // ----- else if template worse than new sample ----- //
                    else {
                        m_population->AddToPopulation(m_population->m_individuals[template_idx]->Genes(), template_idx, template_fitness);
                        m_population->SetToPopulation(genes, template_idx, new_sample_fitness);
                    }
                    // === end of model 1 === //

                    template_idx += (m_pop_size / 2);
                    // cout << "model 2 template idx: " << template_idx << endl;

                    // === model 2 (cut point == 10) === //
                    template_fitness = m_population->m_individuals[template_idx]->Value();
                    m_model_2->Sample(genes, m_population->m_individuals[template_idx]->Genes(), need_to_sample_genes_count + need_to_sample_genes_count_update_constant);
                    new_sample_fitness =  m_problem->Evaluate(genes);
                    m_evaluations++;

                    if ( template_fitness > new_sample_fitness) {
                        m_population->AddToPopulation(genes, template_idx, new_sample_fitness);
                    }

                    // ----- else if template worse than new sample ----- //
                    else {
                        m_population->AddToPopulation(m_population->m_individuals[template_idx]->Genes(), template_idx, template_fitness);
                        m_population->SetToPopulation(genes, template_idx, new_sample_fitness);
                    }
                    // === end of model 2 === //
                }
                // end of population size for loop
                m_population->SortPopulation(2);
                // m_population->Print();

                if (m_population->m_individuals[0]->Value() > m_population->m_individuals[m_pop_size/2]->Value()) {
                    newScore=m_population->m_individuals[0]->Value();
                }
                else {
                    newScore=m_population->m_individuals[m_pop_size/2]->Value();
                }

                population_avg_fitness = m_population->AverageFitnessPopulation(m_pop_size);
                
                best_fitness_file << m_evaluations << " " << newScore << " " << need_to_sample_genes_count << endl;
                avg_fitness_file << m_evaluations << " " << population_avg_fitness << endl;

            }
            // end of 30 * m_problem_size * m_problem_size loop
            change_NFE *= 2;
            need_to_sample_genes_count_update_constant /= 2;


            // update the number of need to sample genes count
            // cout << need_to_sample_genes_count<< "upper population genes count : " << need_to_sample_genes_count << ", avg_fitness: " << m_population->AverageFitnessPopulationSegment(0, m_pop_size/2) << endl;
            // cout << need_to_sample_genes_count + need_to_sample_genes_count_update_constant<< "lower population genes count : " << need_to_sample_genes_count + need_to_sample_genes_count_update_constant << ", avg_fitness: " << m_population->AverageFitnessPopulationSegment(m_pop_size/2, m_pop_size/2) << endl<< endl;
            if (m_population->AverageFitnessPopulationSegment(0, m_pop_size/2) > m_population->AverageFitnessPopulationSegment(m_pop_size/2, m_pop_size/2)) {
                need_to_sample_genes_count_update_constant = need_to_sample_genes_count_update_constant * (-1);
            }
            else { // original_edge_count_avg_fitness < new_edge_count_avg_fitness
                need_to_sample_genes_count = need_to_sample_genes_count + need_to_sample_genes_count_update_constant;
            }
        

            // random suffle population
            m_population->RandomShufflePopulation();
        }
        // end of else if (((string)m_model_type) == "AGNRBOP")

        else if (((string)m_model_type) == "ACHIGN" ) {

            // cout << "960: " << need_to_sample_genes_count << endl;
            // if (need_to_sample_genes_count + need_to_sample_genes_count_update_constant > m_problem_size/2) { 
            //     need_to_sample_genes_count = m_problem_size / 2; 
            //     need_to_sample_genes_count_update_constant = - abs(need_to_sample_genes_count_update_constant);
            //     cout << "964" << endl;
            // }
            if (need_to_sample_genes_count + need_to_sample_genes_count_update_constant <= 0) { 
                need_to_sample_genes_count = m_problem_size / 20;
                need_to_sample_genes_count_update_constant = abs(need_to_sample_genes_count_update_constant);
                // cout << "969" << endl;
            }
            // cout << "970: " << need_to_sample_genes_count << endl;

            if (chi_std + chi_std_update_constant < 0.0) { 
                chi_std = 0.0;
                chi_std_update_constant = abs(chi_std_update_constant);
            }



            // while (m_evaluations - old_evaluations < 30 * m_problem_size * m_problem_size) {
            // while (m_evaluations - old_evaluations < change_NFE && m_evaluations<m_max_evaluations) {

            // m_model_1->Learn(m_population, 1); // m_sel_size == m_pop_size
            // m_model_2->Learn(m_population, 2); // m_sel_size == m_pop_size
            m_model_1->Learn(m_population, 3); // m_sel_size == m_pop_size
            // m_model_2->Learn(m_population, 3); // m_sel_size == m_pop_size

            // change the following code to one function
            for (i = 0; i < m_offspring_size / 4; i++) {

                int random_template_idx = rand() % (m_pop_size / 4);
                int template_idx = random_template_idx;
                // cout << "model 1 template idx: " << template_idx << endl;

                // === model 1 (cut point == 3) === //
                long template_fitness = m_population->m_individuals[template_idx]->Value();
                
                m_model_1->Sample(genes, m_population->m_individuals[template_idx]->Genes(), need_to_sample_genes_count, chi_std);
                long new_sample_fitness =  m_problem->Evaluate(genes);
                m_evaluations++;

                if ( template_fitness > new_sample_fitness) {
                    m_population->AddToPopulation(genes, template_idx, new_sample_fitness);
                }
                // ----- else if template worse than new sample ----- //
                else {
                    m_population->AddToPopulation(m_population->m_individuals[template_idx]->Genes(), template_idx, template_fitness);
                    m_population->SetToPopulation(genes, template_idx, new_sample_fitness);
                }
                // === end of model 1 === //

                template_idx = random_template_idx + m_pop_size / 4;
                // cout << "model 2 template idx: " << template_idx << endl;

                // === model 22222 === //
                template_fitness = m_population->m_individuals[template_idx]->Value();
                m_model_1->Sample(genes, m_population->m_individuals[template_idx]->Genes(), need_to_sample_genes_count, chi_std + chi_std_update_constant);
                new_sample_fitness =  m_problem->Evaluate(genes);
                m_evaluations++;

                if ( template_fitness > new_sample_fitness) {
                    m_population->AddToPopulation(genes, template_idx, new_sample_fitness);
                }

                // ----- else if template worse than new sample ----- //
                else {
                    m_population->AddToPopulation(m_population->m_individuals[template_idx]->Genes(), template_idx, template_fitness);
                    m_population->SetToPopulation(genes, template_idx, new_sample_fitness);
                }
                // === end of model 22222 === //

                template_idx = random_template_idx + m_pop_size / 2;
                // cout << "model 2 template idx: " << template_idx << endl;

                // === model 33333 === //
                template_fitness = m_population->m_individuals[template_idx]->Value();
                m_model_1->Sample(genes, m_population->m_individuals[template_idx]->Genes(), need_to_sample_genes_count + need_to_sample_genes_count_update_constant, chi_std + chi_std_update_constant);
                new_sample_fitness =  m_problem->Evaluate(genes);
                m_evaluations++;

                if ( template_fitness > new_sample_fitness) {
                    m_population->AddToPopulation(genes, template_idx, new_sample_fitness);
                }

                // ----- else if template worse than new sample ----- //
                else {
                    m_population->AddToPopulation(m_population->m_individuals[template_idx]->Genes(), template_idx, template_fitness);
                    m_population->SetToPopulation(genes, template_idx, new_sample_fitness);
                }
                // === end of model 33333 === //

                template_idx = random_template_idx + 3 * m_pop_size / 4;
                // cout << "model 2 template idx: " << template_idx << endl;

                // === model 44444 === //
                template_fitness = m_population->m_individuals[template_idx]->Value();
                m_model_1->Sample(genes, m_population->m_individuals[template_idx]->Genes(), need_to_sample_genes_count + need_to_sample_genes_count_update_constant, chi_std);
                new_sample_fitness =  m_problem->Evaluate(genes);
                m_evaluations++;

                if ( template_fitness > new_sample_fitness) {
                    m_population->AddToPopulation(genes, template_idx, new_sample_fitness);
                }

                // ----- else if template worse than new sample ----- //
                else {
                    m_population->AddToPopulation(m_population->m_individuals[template_idx]->Genes(), template_idx, template_fitness);
                    m_population->SetToPopulation(genes, template_idx, new_sample_fitness);
                }
                // === end of model 44444 === //
            }
            // end of population size for loop

            newScore=m_population->FindBestFitness(m_pop_size);

            population_avg_fitness = m_population->AverageFitnessPopulation(m_pop_size);
            
            best_fitness_file << m_evaluations << " " << newScore << " " << need_to_sample_genes_count << endl;
            avg_fitness_file << m_evaluations << " " << population_avg_fitness << " " << chi_std << endl;
            // }
            // end of 30 * m_problem_size * m_problem_size loop
            


            if (m_evaluations - old_gene_num_evaluations >= change_gene_num_NFE) {
                // update the number of need to sample genes count
                cout << need_to_sample_genes_count<< "upper population genes count : " << need_to_sample_genes_count << ", avg_fitness: " << m_population->AverageFitnessPopulationSegment(0, m_pop_size/2) << endl;
                cout << need_to_sample_genes_count + need_to_sample_genes_count_update_constant<< "lower population genes count : " << need_to_sample_genes_count + need_to_sample_genes_count_update_constant << ", avg_fitness: " << m_population->AverageFitnessPopulationSegment(m_pop_size/2, m_pop_size/2) << endl<< endl;
                if (m_population->AverageFitnessPopulationSegment(0, m_pop_size/2) > m_population->AverageFitnessPopulationSegment(m_pop_size/2, m_pop_size/2)) {
                    need_to_sample_genes_count_update_constant = need_to_sample_genes_count_update_constant * (-1);
                }
                else { // original_edge_count_avg_fitness < new_edge_count_avg_fitness
                    need_to_sample_genes_count = need_to_sample_genes_count + need_to_sample_genes_count_update_constant;
                }
                change_gene_num_NFE *= 2;
                need_to_sample_genes_count_update_constant /= 2;
                old_gene_num_evaluations = m_evaluations;
                // shuffle part of population
                m_population->RandomShufflePopulationSplit(0, m_pop_size*3/4);
                m_population->RandomShufflePopulationSplit(m_pop_size/4, m_pop_size/2);
            }
        
            if (m_evaluations - old_chi_square_evaluations >= change_chi_square_NFE) {
                // update chi square
                double tmp_avg_fitness = (m_population->AverageFitnessPopulationSegment(0, m_pop_size/4) + m_population->AverageFitnessPopulationSegment(m_pop_size*3/4, m_pop_size/4)) / 2;
                cout << chi_std<< "upper population std : " << chi_std << ", avg_fitness: " << tmp_avg_fitness << endl;
                cout << chi_std + chi_std_update_constant<< "lower population std : " << chi_std + chi_std_update_constant << ", avg_fitness: " << m_population->AverageFitnessPopulationSegment(m_pop_size/4, m_pop_size/2) << endl<< endl;
                if (tmp_avg_fitness > m_population->AverageFitnessPopulationSegment(m_pop_size/4, m_pop_size/2)) {
                    chi_std_update_constant = chi_std_update_constant * (-1);
                }
                else { // original_edge_count_avg_fitness < new_edge_count_avg_fitness
                    chi_std = chi_std + chi_std_update_constant;
                }
                old_chi_square_evaluations = m_evaluations;
                m_population->RandomShufflePopulationSplit(0, m_pop_size/4);
                m_population->RandomShufflePopulationSplit(m_pop_size/2, m_pop_size*3/4);
            }
            
            // chi_std_update_constant /= 2;
            // random suffle population
            // m_population->RandomShufflePopulation();
        }
        // end of else if (((string)m_model_type) == "ACHIGN")



        // === four population === //
        else if (((string)m_model_type) == "FOUR" ) {
            
            if ( m_evaluations < (1000 * m_problem_size * m_problem_size / 1) ) {

                m_model_1->Learn(m_population, 0, m_pop_size / 2); // m_sel_size == m_pop_size
                m_model_2->Learn(m_population, m_pop_size / 2, m_pop_size); // m_sel_size == m_pop_size

                for (i = 0; i < m_offspring_size / 2 && m_evaluations < m_max_evaluations; i++) {

                    int template_idx = rand() % (m_pop_size / 2);
                    // cout << "model 1 template idx: " << template_idx << endl;

                    // === model 1 (cut point == 3) === //
                    long template_fitness = m_population->m_individuals[template_idx]->Value();
                    
                    // 0: sequentail, cut = 3, 1: random, , cut = 10
                    m_model_1->Sample(genes, m_population->m_individuals[template_idx]->Genes(), 0);
                    long new_sample_fitness =  m_problem->Evaluate(genes);
                    m_evaluations++;

                    if ( template_fitness > new_sample_fitness) {
                        m_population->AddToPopulation(genes, template_idx, new_sample_fitness);
                    }
                    // ----- else if template worse than new sample ----- //
                    else {
                        m_population->AddToPopulation(m_population->m_individuals[template_idx]->Genes(), template_idx, template_fitness);
                        m_population->SetToPopulation(genes, template_idx, new_sample_fitness);
                    }
                    // === end of model 1 === //
                }
                for (i = m_offspring_size / 2; i < m_offspring_size && m_evaluations < m_max_evaluations; i++) {

                    int template_idx = m_pop_size / 2 + rand() % (m_pop_size - m_pop_size / 2);

                    // cout << "model 2 template idx: " << template_idx << endl;
                    // === model 2 (cut point == 10) === //
                    long template_fitness = m_population->m_individuals[template_idx]->Value();

                    // 0: sequentail, cut = 3, 1: random, , cut = 10
                    m_model_2->Sample(genes, m_population->m_individuals[template_idx]->Genes(), 1);
                    long new_sample_fitness =  m_problem->Evaluate(genes);
                    m_evaluations++;

                    if ( template_fitness > new_sample_fitness) {
                        m_population->AddToPopulation(genes, template_idx, new_sample_fitness);
                    }

                    // ----- else if template worse than new sample ----- //
                    else {
                        m_population->AddToPopulation(m_population->m_individuals[template_idx]->Genes(), template_idx, template_fitness);
                        m_population->SetToPopulation(genes, template_idx, new_sample_fitness);
                    }
                    // === end of model 2 === //
                }
                // end of population size for loop
                m_population->SortPopulation(2);
                // m_population->Print();

                if (m_population->m_individuals[0]->Value() > m_population->m_individuals[m_pop_size/2]->Value()) {
                    newScore=m_population->m_individuals[0]->Value();
                }
                else {
                    newScore=m_population->m_individuals[m_pop_size/2]->Value();
                }

                // better_cut_point_count
                // cout << "m_pop_size/2: " << m_pop_size/2 << ", m_pop_size: " << m_pop_size << endl;
                // cout << "m_offspring_size/2: " << m_offspring_size/2 << ", m_offspring_size: " << m_offspring_size << endl;
                cout << "m_pop_size/2: " << m_pop_size/2 << ", m_pop_size: " << m_pop_size << endl;

                cout << "c = 3 best: " << m_population->m_individuals[0]->Value() << endl;
                cout << "c = 10 best: " << m_population->m_individuals[m_pop_size/2]->Value() << endl;

                cout << "c = 3 avg: " << m_population->AverageFitnessPopulationSegment(0, m_pop_size/2) << endl;
                cout << "c = 10 avg: " << m_population->AverageFitnessPopulationSegment(m_pop_size/2, m_pop_size/2) << endl;
                if (m_population->AverageFitnessPopulationSegment(0, m_pop_size/2) > m_population->AverageFitnessPopulationSegment(m_pop_size/2, m_pop_size/2)) {
                    better_cut_point_count = 0;
                    cout << "3" << endl;

                }
                else {
                    better_cut_point_count = 1;
                    cout << "10" << endl;
                }
                // cout << "better_cut_point_count: " << better_cut_point_count << endl<< endl;

            } // end of if (m_evaluations < m_max_evaluations / 20)

            else { // if (m_evaluations >= m_max_evaluations / 20)
                m_model_1->Learn(m_population, 0, m_pop_size); // m_sel_size == m_pop_size

                for (i = 0; i < m_offspring_size && m_evaluations < m_max_evaluations; i++) {

                    int template_idx = rand() % (m_pop_size);


                    // === model 1 (cut point == 3) === //
                    long template_fitness = m_population->m_individuals[template_idx]->Value();
                    
                    m_model_1->Sample(genes, m_population->m_individuals[template_idx]->Genes(), better_cut_point_count);
                    long new_sample_fitness =  m_problem->Evaluate(genes);
                    m_evaluations++;


                    if ( template_fitness > new_sample_fitness) {
                        m_population->AddToPopulation(genes, template_idx, new_sample_fitness);
                    }

                    // ----- else if template worse than new sample ----- //
                    else {
                        m_population->AddToPopulation(m_population->m_individuals[template_idx]->Genes(), template_idx, template_fitness);
                        m_population->SetToPopulation(genes, template_idx, new_sample_fitness);
                    }
                    // === end of model 1 === //
                }
                // end of population size for loop
                m_population->SortPopulation(0);
                newScore=m_population->m_individuals[0]->Value();
            }
            // end of if (m_evaluations < m_max_evaluations / 20)
        }
        // end of else if (((string)m_model_type) == "FOUR")



        // printf("iteration: %d\n", iterations);

        
        // --- write file ---//
        if (((string)m_model_type) != "AGNRBOP" && ((string)m_model_type) != "PREPRUNE" && ((string)m_model_type) != "ACHIGN") {
            population_avg_fitness = m_population->AverageFitnessPopulation(m_pop_size);
            // cout << "1259";
            best_fitness_file << m_evaluations << " " << newScore << " " << need_to_sample_genes_count << endl;
            avg_fitness_file << m_evaluations << " " << population_avg_fitness << endl;

        }
        // --- write file end ---//


        float modification=0;
        if (newScore>m_best->Value())
        {
            
            m_best->SetGenes(m_population->m_individuals[0]->Genes());
            m_best->SetValue(newScore);
            m_convergence_evaluations=m_evaluations;
          /* if (((string)m_model_type)=="M")
                cout<<""<<m_population->m_individuals[0]->Value()<<" , "<<m_evaluations<<" , "<<m_max_evaluations-m_evaluations<<"  Theta. "<<((CMallowsModel*)m_model)->m_distance_model->m_theta_parameter<<" lower: "<<((CMallowsModel*)m_model)->m_distance_model->m_lower_theta_bound<<endl;

            else
                cout<<""<<m_population->m_individuals[0]->Value()<<" , "<<m_evaluations<<" , "<<m_max_evaluations-m_evaluations<<"  Thetas. "<<((CGeneralizedMallowsModel*)m_model)->m_distance_model->m_theta_parameters[0]<<" lower: "<<((CGeneralizedMallowsModel*)m_model)->m_distance_model->m_lower_theta_bound<<endl;
        */
           modification=-rate;
        }
        else{
            modification=rate;
        }

        // cout << "best: " <<m_best->Value() << endl;
        // cout << "avg: " <<population_avg_fitness << endl;
        
        if (((string)m_model_type)=="M"){
            ((CMallowsModel*)m_model)->m_distance_model->m_lower_theta_bound= ((CMallowsModel*)m_model)->m_distance_model->m_lower_theta_bound+modification;
            if (((CMallowsModel*)m_model)->m_distance_model->m_lower_theta_bound<0.001)
                ((CMallowsModel*)m_model)->m_distance_model->m_lower_theta_bound=0.001;
        }
        else if (((string)m_model_type)=="GM"){
            ((CGeneralizedMallowsModel*)m_model)->m_distance_model->m_lower_theta_bound= ((CGeneralizedMallowsModel*)m_model)->m_distance_model->m_lower_theta_bound+modification;
            if (((CGeneralizedMallowsModel*)m_model)->m_distance_model->m_lower_theta_bound<0.001)
            ((CGeneralizedMallowsModel*)m_model)->m_distance_model->m_lower_theta_bound=0.001;
        }

        iterations++;

    } // iteration loop end


    best_fitness_file.close();
    avg_fitness_file.close();

    if (((string)m_model_type) == "ACNRBOP") {
        ofstream dsitance_file;
        string distance_file_name = m_result_file_name;
        distance_file_name.replace(distance_file_name.find(result_name), result_name.length(), "distance");
        dsitance_file.open(distance_file_name);


        for (int i = 0; i < m_problem_size; i ++) {
            dsitance_file << used_distance_count_arr[i] << " ";
        }
        dsitance_file << endl;
        // cout << "used_distance_count_arr: ";
        // PrintArray(used_distance_count_arr, m_problem_size);
        dsitance_file.close();
    }



    delete [] genes;

    return 0;
}

/*
 * Returns the number of performed evaluations.
 */
int RankingEDA::GetPerformedEvaluations(){
    return m_convergence_evaluations;
}

int RankingEDA::GetRealConvergeNFE(){
    return m_real_converge_evaluations;
}
/*
 * Returns the fitness of the best solution obtained.
 */
long int RankingEDA::GetBestSolutionFitness(){
    return m_best->Value();
}

/*
 * Returns the best solution obtained.
 */
CIndividual * RankingEDA::GetBestSolution(){
    return m_best;
}


/*
 * This method applies a swap of the given i,j positions in the array.
 */
void RankingEDA::Swap(int * array, int i, int j)
{
	int aux=array[i];
	array[i]=array[j];
	array[j]=aux;
}

/*
 * Experiment for the Metric suitability for the different problems.
 */
void RankingEDA::MetricSuitability_Experiment(char * results_file){

    int tests=1000;
    int samples=1000;
    int * consensus= new int[m_problem_size];
    double thetas_K[10]={1.48,1.8,2.08,2.33,2.6,2.9,3.28,3.75,4.5,10};
    double thetas_C[10]={2.85,3.25,3.56,3.85,4.15,4.47,4.83,5.3,6.1,12};
    double thetas_U[10]={3.3,3.75,4.1,4.4,4.7,5.02,5.4,5.9,6.7,12};
    double theta=0;
    int distance_max_K=(m_problem_size-1)*m_problem_size/2;
    int distance_max_C=m_problem_size-1;
    int distance_max_U=m_problem_size-1;
    int distance=0;
    double fitness_variation=0;
    int fitness_differential=0;
    int distance_max;
    
    int * sample= new int[m_problem_size];
    double total_tests=0;
   // ofstream output_file;
   // output_file.open(results_file);
    for (int j=0;j<10;j++){
        
        if (((string)m_metric_type)=="K")
            theta=thetas_K[j];
        else if (((string)m_metric_type)=="C")
            theta=thetas_C[j];
        else
            theta=thetas_U[j];
        total_tests=0;
        for (int i=0;i<tests;i++){
            fitness_variation=0;
            GenerateRandomPermutation(consensus, m_problem_size);
            m_model->Learn(consensus, theta);
        
            for (int z=0;z<samples;z++){
                m_model->Sample(sample);
                
                //distance of the samples solution with respect to the consensus ranking.
                if (((string)m_metric_type)=="K"){
                    distance = Kendall(consensus,sample,m_problem_size);
                    distance_max=distance_max_K;
                }
                else if (((string)m_metric_type)=="C"){
                    distance = Cayley(consensus,sample,m_problem_size);
                    distance_max=distance_max_C;
                }
                else{
                    distance = Ulam(consensus,sample,m_problem_size);
                    distance_max=distance_max_U;
                }
                
                //fitness differential
                fitness_differential=abs(m_problem->EvaluateInv(sample)-m_problem->EvaluateInv(consensus));
                fitness_variation+=fitness_differential*(1-distance/distance_max);
            }
            fitness_variation=fitness_variation/samples;
            total_tests+=fitness_variation;
            //cout<<fitness_variation<<endl;
        }
        total_tests=total_tests/tests;
        cout<<total_tests<<endl;
    }
    cout<<"--------------------------------"<<endl;
 //   output_file.close();
    delete [] consensus;
    delete [] sample;
    
}