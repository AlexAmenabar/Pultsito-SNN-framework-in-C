#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

#include <omp.h>
#include <immintrin.h> 

// public headers
#include "simulations/simulations.h"
#include "simulations/results.h"
#include "networks/snn.h"
#include "datasets/datasets.h"
#include "config/config_loader.h"
#include "training_rules/stdp.h"
#include "neuron_models/neuron_models.h"

// private headers
#include "priv_simulations.h"
#include "neuron_models/priv_neuron_models.h"


/* General simulation functions */

void simulate_batches(GPU_SNN_t *snn, GPU_dataset_t *dataset, simulation_configuration_t *conf, GPU_results_t *results){

    struct timespec start, end; 
    struct timespec start_uw, end_uw; 
    double elapsed_time;

    size_t n_batches, r_samples;
    size_t b;

    // compute number of batches // [TODO]: improve or refactorize
    n_batches = dataset->n_samples / conf->batch_size;
    r_samples = dataset->n_samples % conf->batch_size;
    
    n_batches = r_samples > 0 ? n_batches + 1 : n_batches; // one more batch if there are remaining samples

    // copy non-constant snn data for parallel batch simulation
    init_batch_snn(snn, conf);

    // initialize struct to store batch results
    GPU_results_t *batch_results = initialize_batch_results_cpu(conf, snn->n_neurons, conf->batch_size, 1, 1);
    
    
    // simulate test batch [TESTING]
    #ifdef DEBUG
        simulate_batch_CPU(snn, dataset, conf, batch_results, 0, 1);
    #endif

    
    // start timer // [TODO]: not always required
    clock_gettime(CLOCK_MONOTONIC, &start);

    // loop over batches and simulate
    for(b = 0; b<n_batches; b++){
        
        // print for feedback
        if((b+1) % 100 == 0){
            printf(" Simulating batch %zu\n", b+1);
            fflush(stdout);
        }

        // simulate batch
        simulate_batch_CPU(snn, dataset, conf, batch_results, b, 0);

        // accumulate batch execution times in general results struct // [TODO]: not required always
        acc_batch_execution_times(results, batch_results);
    }

    // deallocate memory of batch_restuls
    deallocate_results_str(batch_results);

    
    // compute execution times for each phase // [TODO]: not always required
    clock_gettime(CLOCK_MONOTONIC, &end);
    elapsed_time = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;
    printf(" > Finished in %lf!\n", elapsed_time); 


    // print information in results struct
    printf("   > Total execution time: %.3lf | Mean per batch %.3lf | Mean per time step %.3lf\n", 
                results->t, results->t / (float)n_batches, results->t / (float)n_batches / (float)conf->time_steps);
    printf("   > Total in time:        %.3lf | Mean per batch %.3lf | Mean per time step %.3lf\n", 
                results->t_in, results->t_in / (float)n_batches, results->t_in / (float)n_batches / (float)conf->time_steps);
    printf("   > Total v time:         %.3lf | Mean per batch %.3lf | Mean per time step %.3lf\n", 
                results->t_v, results->t_v / (float)n_batches, results->t_v / (float)n_batches / (float)conf->time_steps);
    printf("   > Total out time:       %.3lf | Mean per batch %.3lf | Mean per time step %.3lf\n", 
                results->t_out, results->t_out / (float)n_batches, results->t_out / (float)n_batches / (float)conf->time_steps);
    printf("   > Total reinit time:    %.3lf | Mean per batch %.3lf | Mean per time step %.3lf\n", 
                results->t_reinit, results->t_reinit / (float)n_batches, results->t_reinit / (float)n_batches / (float)conf->time_steps);
    printf("   > Total load time:      %.3lf | Mean per batch %.3lf | Mean per time step %.3lf\n", 
                results->t_load, results->t_load / (float)n_batches, results->t_load / (float)n_batches / (float)conf->time_steps);
    printf("   > Total learn time:     %.3lf | Mean per batch %.3lf | Mean per time step %.3lf\n", 
                results->t_learn, results->t_learn / (float)n_batches, results->t_learn / (float)n_batches / (float)conf->time_steps);
}

void simulate_batch_CPU(GPU_SNN_t *snn, GPU_dataset_t *dataset, simulation_configuration_t *conf, GPU_results_t *results, size_t bidx, int print_data){

    // variables to store execution times
    struct timespec start, end; 
    struct timespec start_step1, end_step1; 
    struct timespec start_step2, end_step2; 
    struct timespec start_step3, end_step3; 
    struct timespec start_step4, end_step4; 
    struct timespec start_step5, end_step5; 
    struct timespec start_load_batch, end_load_batch;
    double elapsed_time, et1 =0.0, et2=0.0, et3=0.0, et4=0.0, et5=0.0, et_load = 0.0;

    // store information in local variables
    size_t T = (size_t)conf->time_steps;
    size_t B = conf->batch_size;
    size_t LT = snn->LT;

    // indices for looping
    size_t i, t, lt;

    // reinitialize results struct
    reinitialize_batch_results_cpu(results, conf, snn->n_neurons, conf->batch_size, conf->time_steps);

    clock_gettime(CLOCK_MONOTONIC, &start);

    /* simulation step 1: (re)initialize network state */
    clock_gettime(CLOCK_MONOTONIC, &start_step1);
    reinitialize_neurons_batch(snn, conf);
    reinitialize_synapses_batch(snn, conf);
    reinitialize_spk_matrix_batch(snn, conf);
    clock_gettime(CLOCK_MONOTONIC, &end_step1);
    et1+=(end_step1.tv_sec - start_step1.tv_sec) + (end_step1.tv_nsec - start_step1.tv_nsec) / 1e9;

    
    /* simulation step 2: initialize batch matrix (struct that stores only the samples in the batch) */ 
    clock_gettime(CLOCK_MONOTONIC, &start_load_batch);
    tmp_batch_cpu_t *batch_data = initialize_batch_matrix(snn, dataset, conf, bidx);
    clock_gettime(CLOCK_MONOTONIC, &end_load_batch);
    et_load+=(end_load_batch.tv_sec - start_load_batch.tv_sec) + (end_load_batch.tv_nsec - start_load_batch.tv_nsec) / 1e9;


    /* simulationg step3: simulate over time steps */
    for(t=0; t<T; t++){

        // compute local time step for getting data in spike matrix
        lt = t % LT;

        // simulation step 3.1: load time step in SNN spike matrix
        clock_gettime(CLOCK_MONOTONIC, &start_load_batch);
        load_batch_time_step_in_SNN_batch(snn, dataset, conf, batch_data, bidx, t);
        clock_gettime(CLOCK_MONOTONIC, &end_load_batch);
        et_load+=(end_load_batch.tv_sec - start_load_batch.tv_sec) + (end_load_batch.tv_nsec - start_load_batch.tv_nsec) / 1e9;


        /* simulation step 3.2: process incoming spikes (input step, independent of the neuron model) */
        clock_gettime(CLOCK_MONOTONIC, &start_step2); 
        snn->compute_input_current(snn, conf, lt, t);
        clock_gettime(CLOCK_MONOTONIC, &end_step2);
        et2+=(end_step2.tv_sec - start_step2.tv_sec) + (end_step2.tv_nsec - start_step2.tv_nsec) / 1e9;
        

        /* simulationg step 3.3: compute neuron dynamics */
        clock_gettime(CLOCK_MONOTONIC, &start_step3);
        snn->compute_dynamics(snn, conf, lt, t);
        clock_gettime(CLOCK_MONOTONIC, &end_step3);
        et3+=(end_step3.tv_sec - start_step3.tv_sec) + (end_step3.tv_nsec - start_step3.tv_nsec) / 1e9;


        /* simulationg step 3.4: compute fired spikes (v[t] > v_thresh) */
        clock_gettime(CLOCK_MONOTONIC, &start_step4);
        snn->compute_firing(snn, conf, results, lt, t);
        clock_gettime(CLOCK_MONOTONIC, &end_step4);
        et4+=(end_step4.tv_sec - start_step4.tv_sec) + (end_step4.tv_nsec - start_step4.tv_nsec) / 1e9;


        /* simulation step 3.5: compute learning rule */
        clock_gettime(CLOCK_MONOTONIC, &start_step5);
        if(conf->learn)
            process_trace_based_STDP_batch(snn, conf); // [TODO]: not always TB-STDP
        clock_gettime(CLOCK_MONOTONIC, &end_step5);
        et5+=(end_step5.tv_sec - start_step5.tv_sec) + (end_step5.tv_nsec - start_step5.tv_nsec) / 1e9;
        
    }

    /* simulation step 4: update weights after simulating the batch */
    clock_gettime(CLOCK_MONOTONIC, &start_step5);
    if(conf->learn)
        update_weights_cpu(snn, conf);
    clock_gettime(CLOCK_MONOTONIC, &end_step5);
    et5+=(end_step5.tv_sec - start_step5.tv_sec) + (end_step5.tv_nsec - start_step5.tv_nsec) / 1e9;



    clock_gettime(CLOCK_MONOTONIC, &end);
    elapsed_time = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;


    // deallocate batch matrix structure
    deallocate_batch_matrix(batch_data);

    // store execution times and number of spikes
    results->t = elapsed_time;
    results->t_in = et2;
    results->t_v = et3;
    results->t_out = et4;
    results->t_reinit = et1;
    results->t_load = et_load;
    results->t_learn = et5;

    // print execution data if required // [TODO]: move to function??
    //if(print_data == 1){
        // print execution times
        /*printf(" > Finished in %lf! (s_rnt = %lf s_load = %lf s_in = %lf s_v = %lf s_out = %lf, s_tr = %lf)\n", 
            elapsed_time, et1, et_load, et2, et3, et4, et5);
        fflush(stdout);

        // print generated spikes
        printf(" > Generated number of spikes per neuron: ");
        for(size_t b = 0; b<conf->batch_size; b++){

            for(i = 0; i<snn->n_neurons; i++){
                printf("%d ", results->n_spks[i * B + b]);
                results->n_spks[i * B + b] = 0;
            }
            printf("\n");
        }
        printf("\n");*/
    //}
}

/* Weights updating */
void update_weights_cpu(GPU_SNN_t *snn, simulation_configuration_t *conf){

    size_t S, P, B;
    size_t i, b;

    S = snn->n_synapses;
    B = conf->batch_size;
    P = conf->n_process;

    #if defined AVX512
    {

        __m512 vZero = _mm512_set1_ps(0.0);

        if(B == 8){
            
            __mmask16 m = 0xFF; // only 256 bits used

            #pragma omp parallel for num_threads(P) private(i)
            for(i = 0; i<S; i++){

                float dw = 0.0;
                
                // sum Dw of each batch element (256 bits)
                __m512 vDw = _mm512_maskz_loadu_ps(m, &(snn->dw[i * B + 0])); // need 256 bits, 512 are loaded (if it is the final, padding is loaded)
                for(b = 8; b+7 < B; b+=8){

                    vDw = _mm512_mask_add_ps(vDw, m, vDw, _mm512_maskz_loadu_ps(m, &(snn->dw[i * B + b])));
                }
                
                // reduce to a value
                dw = _mm512_mask_reduce_add_ps(m, vDw);

                // compute mean dw
                dw = dw / (float)B;

                // update initial w and w for all copies
                vDw = _mm512_set1_ps(dw);
                for(b = 0; b+7<B; b+=8){

                    // reinitialize dw
                    _mm512_mask_storeu_ps(&(snn->dw[i * B + b]), m, vZero);

                    // init_w = init_w + dw
                    __m512 vInitW = _mm512_maskz_loadu_ps(m, &(snn->init_w[i * B + b]));
                    vInitW = _mm512_mask_add_ps(vInitW, m, vInitW, vDw);
                    _mm512_mask_storeu_ps(&(snn->init_w[i * B + b]), m, vInitW);

                    // w = init_w
                    _mm512_mask_storeu_ps(&(snn->w[i * B + b]), m, vInitW);
                }    
            }
        }
        else if(B % 16 == 0)
        {
            #pragma omp parallel for num_threads(P) private(i)
            for(i = 0; i<S; i++){

                float dw = 0.0;
                
                // sum Dw of each batch element
                __m512 vDw = _mm512_loadu_ps(&(snn->dw[i * B + 0]));
                for(b = 16; b+15 < B; b+=16){

                    vDw = _mm512_add_ps(vDw, _mm512_loadu_ps(&(snn->dw[i * B + b])));
                }
                
                // reduce to a value
                dw = _mm512_reduce_add_ps(vDw);

                // compute mean dw
                dw = dw / (float)B;

                // update initial w and w for all copies
                vDw = _mm512_set1_ps(dw);
                for(b = 0; b+15<B; b+=16){

                    // dw = 0.0
                    _mm512_storeu_ps(&(snn->dw[i * B + b]), vZero); // dw = 0.0, reinitialize

                    // init_w = init_w + dw
                    __m512 vInitW = _mm512_loadu_ps(&(snn->init_w[i * B + b]));
                    vInitW = _mm512_add_ps(vInitW, vDw);
                    _mm512_storeu_ps(&(snn->init_w[i * B + b]), vInitW);

                    // w = init_w
                    _mm512_storeu_ps(&(snn->w[i * B + b]), vInitW);
                }    
            }        
        }
        else{
            // TODO

        }
    }
    #else
    {
        
        #pragma omp parallel for num_threads(P) private(i)
        for(i = 0; i<S; i++){

            float dw = 0.0;

            // sum dw of the batch for each synapse
            for(size_t b = 0; b<B; b++){

                dw += snn->dw[i * B + b];
            }

            // compute mean dw
            dw = dw / (float)B;

            // update initial w and w for all copies
            for(size_t b = 0; b<B; b++){

                snn->dw[i * B + b] = 0.0; // reinitialize dw
                snn->init_w[i * B + b] += dw; // update initial w for incoming batches
                snn->w[i * B + b] = snn->init_w[i * B + b]; // update w
            }    
        }
    }
    #endif
}


/* Intermedaite batch bitmap management  */

tmp_batch_cpu_t* initialize_batch_matrix(GPU_SNN_t *snn, GPU_dataset_t *dataset, simulation_configuration_t *conf, size_t bidx){

    size_t i, j, n_features, b;    
    size_t N, iN, P, B, T;

    N = snn->n_neurons;
    iN = snn->n_input_neurons;
    P = conf->n_process;
    B = conf->batch_size;
    T = conf->time_steps;

    // get number of features in samples
    n_features = dataset->n_features;

    // index of the first sample in the batch
    size_t fsidx = bidx * B; 
    size_t sidx; // var to store the global index of the sample

    // allocate memory for temporal batch samples matrix // TODO: refactorize to function
    tmp_batch_cpu_t *tmp_batch = (tmp_batch_cpu_t*)calloc(1, sizeof(tmp_batch_cpu_t));
    tmp_batch->spikes = (char*)calloc(iN * B * T, sizeof(char));


    // load the samples in the batch in the matrix: loop over sample features
    #pragma omp parallel for num_threads(P) private(i, j, b, sidx)
    for(i=0; i < iN; i++){

        // loop over batch (sample = fsidx + b)
        for(b = 0; b<B; b++){

            // get sample index in batch
            sidx = fsidx + b;

            // check if the sample is in the dataset
            if(sidx < dataset->n_samples){

                // get the feature index in the sample
                size_t fidx = i; 
                
                // get the feature offset in the array of spikes
                size_t feature_offset = dataset->feature_offset[sidx * n_features + fidx];

                // get the number of spikes of the feature
                size_t n_spikes_feature = dataset->n_spikes_per_feature[sidx * n_features + fidx];
        
                // loop over spikes and store in the matrix
                for(j = 0; j<n_spikes_feature; j++){

                    // get time step of the spike
                    size_t t = dataset->spikes[feature_offset + j];

                    // store the spike in the matrix of [iN * B * T]
                    tmp_batch->spikes[(iN * B * t) + (B * i) + b] = 1; 
                }
            }
        }
    }

    return tmp_batch;
}

void deallocate_batch_matrix(tmp_batch_cpu_t* batch_matrix){
    
    if(batch_matrix->spikes) free(batch_matrix->spikes);
    if(batch_matrix) free(batch_matrix);
}


/* Simulation internal functions */

void load_batch_time_step_in_SNN_batch(GPU_SNN_t *snn, GPU_dataset_t *dataset, simulation_configuration_t *conf, tmp_batch_cpu_t *batch_data, size_t bidx, size_t t){
    
    size_t i, b;    
    size_t N, iN, P, B, L;
    size_t fsidx, sidx, lt, spk_mtr_idx, batch_mtr_idx;
    
    // get general information
    N = snn->n_neurons;
    iN = snn->n_input_neurons;
    P = conf->n_process;
    B = conf->batch_size;
    L = snn->LT;
    
    // get the index of the first sample in the batch
    fsidx = bidx * B;
    lt = t % L; // compute local time step for the matrix


    #if defined AVX512
    {
        size_t n_tasks, tsk, idx;

        // spike matrix has chars (8 bits), so in 512 bits, 512 / 8 == 64 chars

        // compute matrix indexes first time step index
        spk_mtr_idx = (iN + N) * B * lt;
        batch_mtr_idx = (iN * B * t);
        
        if(B == 8){

            #pragma omp parallel for num_threads(P) private(i, b, spk_mtr_idx, batch_mtr_idx)
            for(i = 0; i < iN; i++){
                
                // compute matrix indexes
                spk_mtr_idx = (iN + N) * B * lt + (B * i);
                batch_mtr_idx = (iN * B * t) + (B * i);

                for(b = 0; b+7<B; b += 8){
                    
                    // load spikes from batch_data
                    __m128i vSpikes = _mm_loadl_epi64((__m128i*)(&(batch_data->spikes[batch_mtr_idx + b])));

                    // store spikes in spk matrix
                    _mm_storel_epi64((__m128i*)(&(snn->spk_matrix[spk_mtr_idx + b])), vSpikes);
                }
            }

            // reinit the rest of the row
            #pragma omp parallel for num_threads(P) private(i, b, spk_mtr_idx)
            for(i = 0; i<N; i++){

                // compute matrix general index
                spk_mtr_idx = (iN + N) * B * lt + (B * (iN + i));

                for(b = 0; b+7<B; b += 8){
                
                    _mm_storel_epi64((__m128i*)(&(snn->spk_matrix[spk_mtr_idx + b])), _mm_setzero_si128());
                }
            }
        }
        else{

            #pragma omp parallel for num_threads(P) private(i, b, spk_mtr_idx, batch_mtr_idx)
            for(i = 0; i < iN; i++){
                
                // compute matrix indexes
                spk_mtr_idx = (iN + N) * B * lt + (B * i);
                batch_mtr_idx = (iN * B * t) + (B * i);

                for(b = 0; b+15<B; b+= 16){
                    
                    // load spikes from batch_data
                    __m128i vSpikes = _mm_loadu_si128((__m128i*)(&(batch_data->spikes[batch_mtr_idx + b])));

                    // store spikes in spk matrix
                    _mm_storeu_si128((__m128i*)(&(snn->spk_matrix[spk_mtr_idx + b])), vSpikes);
                }
            }

            // reinit the rest of the row
            #pragma omp parallel for num_threads(P) private(i, b, spk_mtr_idx)
            for(i = 0; i<N; i++){

                // compute matrix general index
                spk_mtr_idx = (iN + N) * B * lt + (B * (iN + i));

                for(b = 0; b+15<B; b+= 16){
                
                    _mm_storeu_si128((__m128i*)(&(snn->spk_matrix[spk_mtr_idx + b])), _mm_setzero_si128());
                }
            }
        }
    }
    #else
    {
        // load the batch samples from time step [t] to [t+L]
        #pragma omp parallel for num_threads(P) private(i, b, sidx, spk_mtr_idx, batch_mtr_idx)
        for(i=0; i < iN; i++){ // load features

            // compute matrix indexes
            spk_mtr_idx = (iN + N) * B * lt + (B * i);
            batch_mtr_idx = (iN * B * t) + (B * i);

            // loop over batch
            for(b = 0; b<B; b++){

                // get sample index 
                sidx = fsidx + b;

                // check if the sample is in the dataset
                if(sidx < dataset->n_samples){ // not really necessary, but it's okay
                    
                    //printf(" > Loading spike in neuron %zu in %zu (from %zu): %d\n", i, (iN + N) * B * (t % L) + (B * i) + b, iN * B * t + (B * i) + b, batch_data->spikes[iN * B * t + (B * i) + b]);
                    // store if the neuron fired in the spk matrix [(iN + N) * B * L] from the batch data matrix [iN * B * T]
                    snn->spk_matrix[spk_mtr_idx + b] = batch_data->spikes[batch_mtr_idx + b]; 
                }
            }
        }

        // reinit the rest of the row
        #pragma omp parallel for num_threads(P) private(i, b, spk_mtr_idx)
        for(i = 0; i<N; i++){

            // compute matrix general index
            spk_mtr_idx = (iN + N) * B * lt + (B * (iN + i));

            for(b = 0; b<B; b++){
                snn->spk_matrix[spk_mtr_idx + b] = 0; 
            }
        }
    }
    #endif
}

/* Functions to reinitialize the network(s) */

void reinitialize_neurons_batch(GPU_SNN_t *snn, simulation_configuration_t *conf){
    
    size_t i, b;
    size_t B, N, P;
    size_t g_index;

    N = snn->n_neurons;
    P = conf->n_process;
    B = conf->batch_size;
        
    // loop over neurons for reinitializing variables shared by all neuron types
    #pragma omp parallel for num_threads(P) private(i, b, g_index)
    for(i = 0; i < N; i++){

        g_index = i * B;
        
        // loop over copies in batch (serial)
        for(b = 0; b<B; b++){

            snn->post_fired[g_index + b] = 0;

            if(conf->learn == 1){
                snn->post_trace[g_index + b] = 0.0;
            }
        }
    }

    // reinitialize neuron model dependant variables
    snn->reinit_neurons(snn, conf);
}

void reinitialize_synapses_batch(GPU_SNN_t *snn, simulation_configuration_t *conf){
    
    size_t i, j, b, S, P, B, gindex, lindex;

    S = snn->n_synapses;
    P = conf->n_process;
    B = conf->batch_size;


    // reinitialize synapses 
    #pragma omp parallel for num_threads(P) private(i, b)
    for(i = 0; i<S; i++){

        size_t g_index = i * B;

        // loop over the batch
        for(b = 0; b<B; b++){

            snn->w[g_index + b] = snn->init_w[i];
            snn->pre_fired[g_index + b] = 0;
            snn->pre_trace[g_index + b] = 0.0;
        }
    }
}


void reinitialize_spk_matrix_batch(GPU_SNN_t *snn, simulation_configuration_t *conf){
    
    size_t i, j, b, lindex, gindex;
    size_t N, iN, P, LT, B;

    N = snn->n_neurons;
    iN = snn->n_input_neurons;
    P = conf->n_process;
    B = conf->batch_size;
    LT = snn->LT;

    // reinitialize spike matrix
    #pragma omp parallel for num_threads(P) private(i, j, b, gindex, lindex)
    for(i=0; i < N + iN; i++){
        
        for(j = 0; j<LT; j++){
        
            gindex = i * LT * B + j * B;

            // loop over batch
            for(b = 0; b<B; b++){
                
                lindex = gindex + b;
                snn->spk_matrix[lindex] = 0;
            }    
        }
    }
}

double get_tmp_batch_size(size_t iN, size_t nS, size_t T){

    return (

        (sizeof(tmp_batch_cpu_t) + 
        sizeof(size_t) * iN * nS * T 
        ) / 8.0
    );
}



void print_spk_matrix_2D(GPU_SNN_t *snn){
  
    size_t iN, N, LT;

    iN = snn->n_input_neurons;
    N = snn->n_neurons;

    LT = snn->LT;

    for(size_t t = 0; t<LT; t++){

        for(size_t i = 0; i<iN + N; i++){

            printf("%d ", snn->spk_matrix[(iN + N) * t + i]);
        }
        printf("\n");
    }
    printf("\n");
}

void print_spk_matrix_3D(GPU_SNN_t *snn, size_t B){
  
    size_t iN, N, LT;
    size_t t, i, b;

    iN = snn->n_input_neurons;
    N = snn->n_neurons;
    LT = snn->LT;
    
    for(t = 0; t<LT; t++){

        for(i = 0; i<iN + N - 1; i++){

            printf("[");
            for(b = 0; b < B-1; b++){
                
                printf("%d ", snn->spk_matrix[(iN + N) * B * t + i * B + b]);
            }
            printf("%d], ", snn->spk_matrix[(iN + N) * B * t + i * B + b]);
        }

        printf("[");
        for(b = 0; b < B-1; b++){
            
            printf("%d ", snn->spk_matrix[(iN + N) * B * t + i * B + b]);
        }
        printf("%d]", snn->spk_matrix[(iN + N) * B * t + i * B + b]);
        
        printf("\n");
    }
    printf("\n");
}

void print_weights_batch(GPU_SNN_t *snn, size_t B, size_t b){
  
    size_t S, LT;
    size_t t, i;

    S = snn->n_synapses;
    LT = snn->LT;
    
    for(i = 0; i < S; i++){

        //printf("[");
        //for(b = 0; b < B-1; b++){
            
        printf("%f ", snn->w[i * B + b]);
        //}
        //printf("%f], ", snn->w[i  * B + B - 1]);
    }
    printf("\n");
}

void print_weights_3D(GPU_SNN_t *snn, size_t B){
  
    size_t S, LT;
    size_t t, i, b;

    S = snn->n_synapses;
    LT = snn->LT;
    
    printf(" W = ");
    for(i = 0; i < S; i++){

        printf("[");
        for(b = 0; b < B-1; b++){
            
            printf("%f ", snn->w[i * B + b]);
        }
        printf("%f], ", snn->w[i  * B + B - 1]);
    }
    printf("\n");
}

void print_dw_3D(GPU_SNN_t *snn, size_t B){
  
    size_t S, LT;
    size_t t, i, b;

    S = snn->n_synapses;
    LT = snn->LT;
    
    printf(" dW = ");
    for(i = 0; i < S; i++){

        printf("[");
        for(b = 0; b < B-1; b++){
            
            printf("%f ", snn->dw[i * B + b]);
        }
        printf("%f], ", snn->dw[i  * B + B - 1]);
    }
    printf("\n");
}