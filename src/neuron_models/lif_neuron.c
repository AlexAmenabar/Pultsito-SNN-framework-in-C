#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <math.h>
#include <time.h>

#include <omp.h>
#include <immintrin.h> 

// libs
#include "toml_c/toml.h"

// public headers
#include "config/config_loader.h"
#include "networks/snn.h"
#include "neuron_models/neuron_models.h"
#include "networks/snn_generator.h"
#include "utils.h"

// internal headers
#include "priv_neuron_models.h"
#include "simulations/priv_simulations.h"


neurons_t load_LIF_neurons_from_file(topology_t *topology, simulation_configuration_t *conf, toml_table_t *tbl){

    toml_value_t v_thres, v_rest, t_refract, res;
    neurons_t neurons; // structure to store neurons data
    size_t i;

    FILE *f = NULL;
    open_file(&f, conf->network_neurons_file); // TOML file

    size_t N = topology->n_neurons;
    float tmp = 0.0;
    int tmp_int = 0;

    neurons.v_thresh = (float *)malloc(N * sizeof(float)); // thresholds
    neurons.v_rest = (float *)malloc(N * sizeof(float)); // resting potentials
    neurons.rft_per = (int *)malloc(N * sizeof(int)); // refractary times
    neurons.R = (float *)malloc(N * sizeof(float)); // resistances

    /* load [Neurons] section */
    // load whether different parameters are included in the input file
    v_thres = toml_table_int(tbl, "v_thres");
    v_rest = toml_table_int(tbl, "v_rest");
    t_refract = toml_table_int(tbl, "t_refract");
    res = toml_table_int(tbl, "resistance");


    // load thresholds
    if(v_thres.ok && v_thres.u.i == 1){
        for(i=0; i<N; i++)
            fscanf(f, "%f", &(neurons.v_thresh[i]));
    }
    else if(v_thres.ok && v_thres.u.i == 2){
        fscanf(f, "%f", &(tmp));
        for(i=0; i<N; i++)
            neurons.v_thresh[i] = tmp;
    }
    else{
        for(i=0; i<N; i++)
            neurons.v_thresh[i] = 1.0;
    }

    // load resting potentials
    if(v_rest.ok && v_rest.u.i == 1){
        for(i=0; i<N; i++)
            fscanf(f, "%f", &((neurons.v_rest)[i]));
    }    
    else if(v_rest.ok && v_rest.u.i == 2){
        fscanf(f, "%f", &(tmp));
        for(i=0; i<N; i++)
            neurons.v_rest[i] = tmp;
    }
    else{
        for(i=0; i<N; i++)
            neurons.v_rest[i] = 0.0;
    }

    // load refractory periods
    if(t_refract.ok && t_refract.u.i == 1){
        for(i=0; i<N; i++)
            fscanf(f, "%d", &((neurons.rft_per)[i]));
    }    
    else if(t_refract.ok && t_refract.u.i == 2){
        fscanf(f, "%d", &(tmp_int));
        for(i=0; i<N; i++)
            neurons.rft_per[i] = (int)(tmp_int);
    }
    else{
        for(i=0; i<N; i++)
            neurons.rft_per[i] = 1;
    }

    // load resistances
    if(res.ok && res.u.i == 1){
        for(i=0; i<N; i++)
            fscanf(f, "%f", &((neurons.R)[i]));
    }    
    else if(res.ok && res.u.i == 2){
        fscanf(f, "%f", &(tmp));
        for(i=0; i<N; i++)
            neurons.R[i] = tmp;
    }
    else{
        for(i=0; i<N; i++)
            neurons.R[i] = 1.0;
    }
    fclose(f);

    return neurons;
}

void set_LIF_ptrs(GPU_SNN_t *snn){

    snn->allocate_neurons = &allocate_memory_for_LIF_neurons;
    snn->init_neurons = &initialize_LIF_neurons;
    snn->cpy_neurons = &init_batch_LIF;
    snn->reinit_neurons = &reinitialize_LIF_neurons_batch;
    snn->compute_input_current = &compute_input_current_batch; 
    snn->compute_dynamics = &compute_LIF_V_batch;
    snn->compute_firing = &process_neuron_firing_batch;
    snn->deallocate_neurons = &deallocate_memory_for_LIF_neurons;
}

void allocate_memory_for_LIF_neurons(GPU_SNN_t *snn, size_t N, simulation_configuration_t *conf){

    // allocate memory for neurons and synapses
    snn->v               = (float*)malloc(N * sizeof(float));
    snn->v_thresh        = (float*)malloc(N * sizeof(float));
    snn->v_rest          = (float*)malloc(N * sizeof(float));
    snn->arrI            = (float*)malloc(N * sizeof(float));
    snn->r_period        = (int*)malloc(N * sizeof(int));
    snn->r_period_remain = (int*)malloc(N * sizeof(int));
    snn->res             = (int*)malloc(N * sizeof(int));
}

void deallocate_memory_for_LIF_neurons(GPU_SNN_t *snn){

    if(snn->v)               free(snn->v);
    if(snn->v_thresh)        free(snn->v_thresh);
    if(snn->v_rest)          free(snn->v_rest);
    if(snn->arrI)            free(snn->arrI);
    if(snn->r_period)        free(snn->r_period);
    if(snn->r_period_remain) free(snn->r_period_remain);
    if(snn->res)             free(snn->res);
}

void initialize_LIF_neurons(GPU_SNN_t *snn, topology_t *topology, simulation_configuration_t *conf){

    size_t i;
    
    // initialize neuron parameters from array
    for(i = 0; i<snn->n_neurons; i++){

        snn->v_thresh[i] = topology->neurons.v_thresh[i];
        snn->v_rest[i] = topology->neurons.v_rest[i]; // this or the next one?
        snn->v[i] = snn->v_rest[i]; 
        snn->res[i] = (int)(topology->neurons.R[i]);
        snn->r_period[i] = topology->neurons.rft_per[i];
        snn->r_period_remain[i] = -1;
        snn->arrI[i] = 0.0;
    }
}


void init_batch_LIF(GPU_SNN_t *snn, simulation_configuration_t *conf){

    size_t i, j, p, N, B;
    N = snn->n_neurons;
    B = conf->batch_size;

    size_t padding = 32;

    float *tmp_v, *tmp_arrI;
    int *tmp_rperiod_remain;

    // store old values in temporal variables
    tmp_v              = snn->v;
    tmp_rperiod_remain = snn->r_period_remain;
    tmp_arrI           = snn->arrI;

    // allocate memory for storing variables copies
    snn->v               = (float*)malloc(N * B * sizeof(float) + padding);
    snn->arrI            = (float*)malloc(N * B * sizeof(float) + padding);
    snn->r_period_remain = (int*)malloc(N * B * sizeof(int) + padding);

    // copy original data copied batch_size times
    for(i = 0; i<N; i++){

        for(j = 0; j<B; j++){

            // copy values
            snn->v[i * B + j] = tmp_v[i];
            snn->r_period_remain[i * B + j] = tmp_rperiod_remain[i];
            snn->arrI[i * B + j] = tmp_arrI[i];
        }
    }

    // initialize padding to 0
    for(p = N * B ; p<(padding / sizeof(float)) + N * B; p++){

        snn->v[p] = 0.0;                       
        snn->arrI[p] = 0.0;                         
    }
    for(p = N * B ; p<(padding / sizeof(int)) + N * B; p++){

        snn->r_period_remain[p] = 0;              
    }


    // deallocate memory allocated by original arrays
    if(tmp_v) free(tmp_v);
    if(tmp_rperiod_remain) free(tmp_rperiod_remain);
    if(tmp_arrI) free(tmp_arrI);
}

void reinitialize_LIF_neurons_batch(GPU_SNN_t *snn, simulation_configuration_t *conf){
    
    size_t i, b;
    size_t B, N, P;
    size_t g_index;

    N = snn->n_neurons;
    P = conf->n_process;
    B = conf->batch_size;
        
    // loop over neurons
    #pragma omp parallel for num_threads(P) private(i, b, g_index)
    for(i = 0; i < N; i++){

        g_index = i * B;
        
        // loop over copies in batch (serial)
        for(b = 0; b<B; b++){

            snn->v[g_index + b] = snn->v_rest[i];
            snn->r_period_remain[g_index + b] = 0;
            snn->arrI[g_index + b] = 0.0;
        }
    }
}

void compute_LIF_V_batch(GPU_SNN_t *snn, simulation_configuration_t *conf, size_t t, size_t gt){

    // get general information
    size_t N, P, B;
    size_t i, b;
    size_t neuron_index, g_neuron_index;

    N = snn->n_neurons;
    P = conf->n_process;
    B = conf->batch_size;

    // helpers
    const float alpha = 0.95f;
    const float beta = 0.05f;

    #if defined AVX512
    {
        if(B == 8){
            
            // initialize constants
            __m256 avec = _mm256_set1_ps(alpha);
            __m256 bvec = _mm256_set1_ps(beta);

            #pragma omp parallel for num_threads(P) private(i, b, neuron_index, g_neuron_index)
            for(i = 0; i<N; i++){   
                
                neuron_index = i;
                __m256 v_rest_vec = _mm256_set1_ps(snn->v_rest[neuron_index]);

                for(b = 0; b+7<B; b+=8){

                    g_neuron_index = i * B + b;

                    // load refractary period and create mask
                    __m256i refract_vec = _mm256_loadu_epi32(&(snn->r_period_remain[g_neuron_index]));
                    __mmask8 rft_mask = _mm256_cmp_epi32_mask(refract_vec, _mm256_setzero_si256(), _MM_CMPINT_LE); // if lower or eual, process

                    // load I, rest and v, if any neuron is not in refractory period
                    if(rft_mask){            

                        __m256 v_vec = _mm256_loadu_ps(&(snn->v[g_neuron_index]));
                        __m256 I_vec = _mm256_loadu_ps(&(snn->arrI[g_neuron_index]));

                        // compute dynamics
                        __m256 tmp = _mm256_add_ps(I_vec, _mm256_mul_ps(v_rest_vec, bvec));
                        __m256 rhs = _mm256_fmadd_ps(v_vec, avec, tmp);  // alpha*v + beta*v_rest + I
                        v_vec = _mm256_mask_mov_ps(v_vec, rft_mask, rhs); // only update active lanes

                        // store v
                        _mm256_storeu_ps(&(snn->v[g_neuron_index]), v_vec);
                    }
                }
            }
        }
        else if(B % 16 == 0){

            // initialize constants
            __m512 avec = _mm512_set1_ps(alpha);
            __m512 bvec = _mm512_set1_ps(beta);

            #pragma omp parallel for num_threads(P) private(i, b, neuron_index, g_neuron_index)
            for(i = 0; i<N; i++){   
                
                neuron_index = i;

                // load v_rest
                __m512 v_rest_vec = _mm512_set1_ps(snn->v_rest[neuron_index]);

                for(b = 0; b+15<B; b+=16){

                    // compute neuron index
                    g_neuron_index = i * B + b;

                    // load refractary period and create mask
                    __m512i refract_vec = _mm512_loadu_si512(&(snn->r_period_remain[g_neuron_index]));
                    __mmask16 rft_mask = _mm512_cmp_epi32_mask(refract_vec, _mm512_setzero_si512(), _MM_CMPINT_LE); // if lower or eual, process

                    // avoid computation if all neurons are in refractory period
                    if(rft_mask){

                        // load v and I
                        __m512 v_vec = _mm512_loadu_ps(&(snn->v[g_neuron_index]));
                        __m512 I_vec = _mm512_loadu_ps(&(snn->arrI[g_neuron_index]));

                        // compute dynamics v = rft ? v : compute_lif()  <-- this is done by the next lines
                        __m512 tmp = _mm512_add_ps(I_vec, _mm512_mul_ps(v_rest_vec, bvec));
                        __m512 rhs = _mm512_fmadd_ps(v_vec, avec, tmp);  // alpha*v + beta*v_rest + I
                        v_vec = _mm512_mask_mov_ps(v_vec, rft_mask, rhs); // only update active lanes

                        // store v (all values are stored, mask is not necessary since the previous operation stores the final values)
                        _mm512_storeu_ps(&(snn->v[g_neuron_index]), v_vec); 
                    }
                }
            }
        }
    }
    #else
    {
        // loop over neurons
        #pragma omp parallel for num_threads(P) private(i, b, neuron_index, g_neuron_index)
        for (i = 0; i < N; i++) {
            
            neuron_index = i;
            g_neuron_index = i * B;

            // loop over batches
            for(b = 0; b<B; b++){
                
                if (snn->r_period_remain[g_neuron_index + b] <= 0){
                    
                    snn->v[g_neuron_index + b] = alpha * snn->v[g_neuron_index + b] + beta * snn->v_rest[neuron_index] + snn->arrI[g_neuron_index + b];
                }
            }
        }
    }
    #endif
}