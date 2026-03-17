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
#include "config/config_loader.h"
#include "neuron_models/neuron_models.h"

// private headers
#include "priv_neuron_models.h"
#include "simulations/priv_simulations.h"


/* Function to compute input currents and neurons firing for spike-based neuron models */

void compute_input_current_batch(GPU_SNN_t *snn, simulation_configuration_t *conf, size_t t, size_t gt){

    size_t i, j, b;
    size_t N, iN, P, B, LT;

    // store general variables
    N = snn->n_neurons;
    iN = snn->n_input_neurons;
    P = conf->n_process;
    B = conf->batch_size;
    LT = snn->LT;
    
    // loop over neurons
    #pragma omp parallel for num_threads(P) private (i, j, b) 
    for(i=0; i<N; i++){
        
        // variables declaration
        size_t synapse_index, g_synapse_index, in_neuron_index;
        int delay, spk_time; 
        float spk;  

        size_t neuron_index = i; // for not copied variables
        size_t g_neuron_index = i * B; // for copied variables

        // get neuron information
        size_t iS = snn->n_neuron_input_synapses[neuron_index]; 
        size_t base_synapse = snn->neuron_input_synapses_offset[neuron_index]; // there are several copies of each synapse, so the offset depends

        #if defined AVX512
        {

            int preIndex = (iN + N) * B;

            // vectorize batch
            if(B == 8){ // masked operations

                size_t n_blks = B / 8;

                // initialize vector of I (for all batch elements)
                __m256 vI_arr[n_blks];
                __m256i vRefPer_arr[n_blks];
                __mmask8 mRefPer_arr[n_blks];
                __mmask8 mRefPer = 0;
                __m256i vOnes = _mm256_set1_epi8(1);
                __m256i vZeros = _mm256_set1_epi8(0);
                __m128i vZeros128 = _mm512_castsi512_si128(_mm512_set1_epi8(0));

                // the same as above but masked
                for(size_t blk = 0; blk<n_blks; blk++){

                    vI_arr[blk] = _mm256_set1_ps(0.0); // I = 0.0;
                    vRefPer_arr[blk] = _mm256_loadu_epi32(&(snn->r_period_remain[g_neuron_index + blk * 8])); 
                    mRefPer_arr[blk] = _mm256_cmp_epi32_mask(vRefPer_arr[blk], vZeros, _MM_CMPINT_LE); 

                    mRefPer |= mRefPer_arr[blk];
                }

                if(mRefPer || conf->learn){ // if at least a bit is 1

                    // loop over input synapses of the neuron
                    for(j = 0; j<iS; j++){

                        // get synapse index. In not copied variables use synapse_index, in copied g_synapse_index 
                        synapse_index = base_synapse + j;
                        g_synapse_index = synapse_index * B; // scale base synapse, since now there are B copies of each one
                        
                        // get synapse delay (not copied)
                        delay = snn->delay[synapse_index];

                        // get spike time
                        spk_time = (int)t - delay;
                    
                        // correct spk_time
                        if(spk_time < 0 && gt >= LT){ // CHECK
                            spk_time = (int)LT + spk_time;
                        }

                        // process spike 
                        if(spk_time >= 0){
                                                            
                            // get input neuron index (not copied)
                            in_neuron_index = snn->pre_neuron_index[synapse_index]; 

                            // index to load the spike from
                            size_t index = (preIndex * (size_t)spk_time) + ((in_neuron_index * B));
                            size_t b_index;

                            for(size_t blk = 0; blk < n_blks; blk++){

                                // check if this part of the batch is in refractory period
                                if(mRefPer_arr[blk] || conf->learn){
                                    
                                    b_index = index + blk * 8;
                                    
                                    // load spikes
                                    //__m128i vCharSpikes = _mm_loadl_epi64((const __m128i*)&snn->spk_matrix[b_index]); // load 8 bytes                                
                                    //__mmask16 spk_mask16 = _mm_cmpgt_epi8_mask(vCharSpikes, _mm_setzero_si128());
                                    //__mmask8 spk_mask = (__mmask8)(spk_mask16 & 0xFF); // store the neurons that received an spike
                                    
                                    // AVX512
                                    __m128i vCharSpikes = _mm_loadu_epi8((const __m128i*)&snn->spk_matrix[b_index]); // 128 bits loaded, but we only need 64 (char(8) * batch_size(8))                                
                                    __mmask16 spk_mask16 = _mm_cmpgt_epi8_mask(vCharSpikes, vZeros128);
                                    __mmask8 spk_mask = (__mmask8)(spk_mask16 & 0xFF); // store the neurons that received an spike


                                    // update I for neurons that are not in refractory period
                                    if(mRefPer_arr[blk] && spk_mask){
                                        
                                        // load weights
                                        __m256 vW = _mm256_loadu_ps(&(snn->w[g_synapse_index + blk * 8]));

                                        // update I for those neurons where r <= 0 and received a spike
                                        vI_arr[blk] = _mm256_mask_add_ps(vI_arr[blk], spk_mask & mRefPer_arr[blk], vI_arr[blk], vW);
                                    }

                                    // if TR-STDP is used, store whether the neuron fired
                                    if(conf->learn && spk_mask){
                                        
                                        // 1 if a spike is received and neuron is not in refractary period
                                        //_mm512_mask_storeu_epi8(&(snn->pre_fired[g_synapse_index + blk * 8]), spk_mask ,_mm512_set1_epi8(1));
                                        _mm_mask_storeu_epi8(&(snn->pre_fired[g_synapse_index + blk * 8]), spk_mask, _mm256_castsi256_si128(vOnes));
                                    }
                                }
                            }
                        }
                    }
                }    
                for(size_t blk = 0; blk < n_blks; blk++){
                    
                    _mm256_storeu_ps(&(snn->arrI[g_neuron_index + 8 * blk]), vI_arr[blk]);
                }                    
            }

            // [B % 16 == 0]
            else if(B % 16 == 0){

                size_t n_blks = B / 16, blk;

                // initialize vector of I (for all batch elements) and refractary periods
                __m512 vI_arr[n_blks];
                __m512i vRefPer_arr[n_blks];
                __mmask16 mRefPer_arr[n_blks];
                __mmask16 mRefPer = 0;


                for(blk = 0; blk<n_blks; blk++){

                    vI_arr[blk] = _mm512_set1_ps(0.0); // I = 0.0;

                    // refractary periods can be used to avoid some computations
                    vRefPer_arr[blk] = _mm512_loadu_si512(&(snn->r_period_remain[g_neuron_index + blk * 16])); 

                    // create mask (epi32 -- int -- 32 bits)
                    mRefPer_arr[blk] = _mm512_cmp_epi32_mask(vRefPer_arr[blk], _mm512_setzero_si512(), _MM_CMPINT_LE); 

                    mRefPer |= mRefPer_arr[blk];
                }

                // loop over input synapses only if there is a neuron in the block of 16 which is not in refractory period
                if(mRefPer || conf->learn){

                    // loop over input synapses of the neuron
                    for(j = 0; j<iS; j++){

                        // get synapse index. In not copied variables use synapse_index, in copied g_synapse_index 
                        synapse_index = base_synapse + j;
                        g_synapse_index = synapse_index * B; // scale base synapse, since now there are B copies of each one
                        
                        // get synapse delay (not copied)
                        delay = snn->delay[synapse_index];

                        // get spike time
                        spk_time = (int)t - delay;
                    
                        // correct spk_time (circular array, negative indexes the last)
                        if(spk_time < 0 && gt >= LT){ // CHECK
                            spk_time = (int)LT + spk_time;
                        }

                        // process spike if necessary
                        if(spk_time >= 0){
                                                            
                            // get input neuron index (not copied)
                            in_neuron_index = snn->pre_neuron_index[synapse_index]; 

                            // index to load the spike from in the spike matrix
                            size_t index = (preIndex * (size_t)spk_time) + ((in_neuron_index * B));
                            size_t b_index; // neuron index

                            // loop over blocks: set of 16 neurons (batched)
                            for(blk = 0; blk < n_blks; blk++){

                                // if all neurons in this mini-batch are in refractary period, avoid computation
                                if(mRefPer_arr[blk] || conf->learn){
                                    
                                    // compute first neuron index
                                    b_index = index + blk * 16;
                                    
                                    // load spikes and get mask
                                    // SSE - AVX512
                                    //__m128i vCharSpikes = _mm_loadu_si128((const __m128i*)&snn->spk_matrix[b_index]);                                
                                    //__mmask16 spk_mask = _mm_cmpgt_epi8_mask(vCharSpikes, _mm_setzero_si128()); // mask the received spikes

                                    // AVX512
                                    __m128i vCharSpikes = _mm_loadu_epi8((const __m128i*)&snn->spk_matrix[b_index]);                                
                                    __mmask16 spk_mask = _mm_cmpgt_epi8_mask(vCharSpikes, _mm_setzero_si128()); // mask the received spikes

                                    // update I only if there are neurons out of refractory period
                                    if(mRefPer_arr[blk] & spk_mask){

                                        // load weights   
                                        __m512 vW = _mm512_loadu_ps(&(snn->w[g_synapse_index + blk * 16]));

                                        // compute I: store in spk_mask & mRefPer: received an spike and not in refractory period
                                        vI_arr[blk] = _mm512_mask_add_ps(vI_arr[blk], spk_mask & mRefPer_arr[blk] , vI_arr[blk], vW);
                                    }

                                    // if TR-STDP is used, store whether the neuron fired
                                    if(conf->learn && spk_mask){
                                        
                                        // 1 if a spike is received and neuron is not in refractary period
                                        _mm512_mask_storeu_epi8(&(snn->pre_fired[g_synapse_index + blk * 16]), spk_mask ,_mm512_set1_epi8(1));
                                        //_mm_mask_storeu_epi8(&(snn->pre_fired[g_synapse_index + blk * 16]), spk_mask,_mm512_castsi512_si128(vOnes));
                                    }
                                } 
                            }
                        }
                    }
                }  

                // store I for all neurons
                for(blk = 0; blk < n_blks; blk++){
                    
                    _mm512_storeu_ps(&(snn->arrI[g_neuron_index + 16 * blk]), vI_arr[blk]);
                }                    
            }
        }   
        #else
        {
            // loop over batches
            for(b = 0; b<B; b++){

                float I = 0.0;

                // check wether the neuron is in refractary period: pre_fired flag not activated although it should be
                if(snn->r_period_remain[g_neuron_index + b] <= 0 || conf->learn){
                
                    // loop over input synapses
                    for(j=0; j<iS; j++){

                        // get synapse index. In not copied variables use synapse_index, in copied g_synapse_index 
                        synapse_index = base_synapse + j;
                        g_synapse_index = synapse_index * B; // scale base synapse, since now there are B copies of each one
                        
                        // get input neuron index (not copied)
                        in_neuron_index = snn->pre_neuron_index[synapse_index]; 

                        // get synapse delay (not copied)
                        delay = snn->delay[synapse_index];

                        // get spike time
                        spk_time = (int)t - delay;
                    
                        // correct spk_time
                        if(spk_time < 0 && gt >= LT){ // CHECK
                            spk_time = (int)LT + spk_time;
                        }

                        // process spike 
                        if(spk_time >= 0){
                            
                            // index to load the spike from
                            size_t index = ((iN + N) * B * (size_t)spk_time) + ((in_neuron_index * B));

                            // get spike
                            spk = (float)(snn->spk_matrix[index + b]); 

                            // update I if r period remain <= 0 and spike received
                            if(snn->r_period_remain[g_neuron_index + b] <= 0 && (int)spk == 1){
                                
                                I += snn->w[g_synapse_index + b] * spk; // spk = 0 / 1
                            }

                            // store if the presynaptic neuron fired for TR-STDP
                            if(conf->learn && (int)spk == 1){

                                snn->pre_fired[g_synapse_index + b] = (char)spk;
                            }
                        }
                    }
                }

                snn->arrI[g_neuron_index + b] = I;
            }
        }
        #endif
    }
}

void process_neuron_firing_batch(GPU_SNN_t *snn, simulation_configuration_t *conf, GPU_results_t *results, size_t t, size_t gt){

    // get general information
    size_t N, iN, P, B;
    size_t neuron_index, g_neuron_index;
    size_t i, b;

    N = (size_t)snn->n_neurons;
    iN = (size_t)snn->n_input_neurons;
    P = (size_t)conf->n_process;
    B = conf->batch_size;

    #if defined AVX512
    {
        if(B == 8){

            // constants
            const __m256i vOnes = _mm256_set1_epi8(1);

            #pragma omp parallel for num_threads(P) private(i, b, neuron_index, g_neuron_index)
            for(i = 0; i<N; i++){
            
                // get global and local neuron indexes
                neuron_index = i;
                g_neuron_index = i * B;
                size_t b_neuron;

                // index of the first neuron in the spike matrix in time step t
                size_t idx = (iN + N) * B * t + B * (iN + neuron_index); 
                size_t batch_idx;

                // load threshold and resting potential
                __m256 vThresh = _mm256_set1_ps(snn->v_thresh[neuron_index]);
                __m256 vRest = _mm256_set1_ps(snn->v_rest[neuron_index]);
                __m256i vRef = _mm256_set1_epi32(snn->r_period[neuron_index]);

                for(b = 0; b+7<B; b+=8){

                    batch_idx = idx + b;
                    b_neuron = g_neuron_index + b; // first neuron index to be loaded in register

                    // reduce refractary period
                    __m256i vRefPer = _mm256_loadu_epi32(&(snn->r_period_remain[b_neuron]));
                    vRefPer = _mm256_sub_epi32(vRefPer, _mm256_set1_epi32(1));

                    // check if neuron fired, and in that case, add spike to matrix
                    __m256 vV = _mm256_loadu_ps(&(snn->v[b_neuron]));

                    // mask for v >= v_thresh
                    __mmask8 fire_mask = _mm256_cmp_ps_mask(vV, vThresh, _MM_CMPINT_GE);

                    if(fire_mask){

                        // v = v_rest for those that fired
                        vV = _mm256_mask_mov_ps(vV, fire_mask, vRest);         
                        _mm256_storeu_ps(&(snn->v[b_neuron]), vV);

                        // r_period_remain = r_period for those that fired
                        vRefPer = _mm256_mask_mov_epi32(vRefPer, fire_mask, vRef); // store r_period where the neuron fired
                        
                        // load, update and store generated spikes
                        __m256i vNSpikes = _mm256_loadu_epi32(&(results->n_spks[b_neuron]));
                        vNSpikes = _mm256_mask_add_epi32(vNSpikes, fire_mask, vNSpikes, _mm256_set1_epi32(1));
                        _mm256_storeu_epi32(&(results->n_spks[b_neuron]), vNSpikes);

                        // store generated spikes in matrix (64 bits)
                        _mm256_mask_storeu_epi8(&(snn->spk_matrix[batch_idx]), fire_mask, _mm256_set1_epi8(1));


                        // if TB-STDP is used, store which neurons fired
                        if(conf->learn){
                            
                            _mm256_mask_storeu_epi8(&(snn->post_fired[b_neuron]), fire_mask, _mm256_set1_epi8(1));
                        }
                    }
                    _mm256_storeu_epi32(&(snn->r_period_remain[b_neuron]), vRefPer);
                }
            }
        }   
        else if(B % 16 == 0){     

            // constants
            const __m512i vOnes = _mm512_set1_epi8(1);

            #pragma omp parallel for num_threads(P) private(i, b, neuron_index, g_neuron_index)
            for(i = 0; i<N; i++){
            
                // get global and local neuron indexes
                neuron_index = i;
                g_neuron_index = i * B;
                size_t b_neuron;

                // index of the first neuron in the spike matrix in time step t
                size_t idx = (iN + N) * B * t + B * (iN + neuron_index); 
                size_t batch_idx;

                // load threshold, resting potential and refractory period
                __m512 vThresh = _mm512_set1_ps(snn->v_thresh[neuron_index]);
                __m512 vRest = _mm512_set1_ps(snn->v_rest[neuron_index]);
                __m512i vRef = _mm512_set1_epi32(snn->r_period[neuron_index]);

                // loop over batch
                for(b = 0; b+15<B; b+=16){

                    // get indexes
                    batch_idx = idx + b; // matrix index
                    b_neuron = g_neuron_index + b; // first neuron index to be loaded in register

                    // reduce refractary period
                    __m512i vRefPer = _mm512_loadu_si512(&(snn->r_period_remain[b_neuron]));
                    vRefPer = _mm512_sub_epi32(vRefPer, _mm512_set1_epi32(1)); // r_remain = r_remain - 1

                    // load V
                    __m512 vV = _mm512_loadu_ps(&(snn->v[b_neuron]));

                    // mask for v >= v_thresh (fired)
                    __mmask16 fire_mask = _mm512_cmp_ps_mask(vV, vThresh, _MM_CMPINT_GE);

                    // avoid computation if no one neuron fired
                    if(fire_mask){
                        
                        // v = v_rest for those that fired
                        vV = _mm512_mask_mov_ps(vV, fire_mask, vRest);         
                        _mm512_storeu_ps(&(snn->v[b_neuron]), vV); // v = v_rest

                        // r_period_remain = r_period for those that fired
                        vRefPer = _mm512_mask_mov_epi32(vRefPer, fire_mask, vRef); // store r_period where the neuron fired, r_remain = r_period

                        // load, update and store generated spikes
                        __m512i vNSpikes = _mm512_loadu_si512(&(results->n_spks[b_neuron])); 
                        vNSpikes = _mm512_mask_add_epi32(vNSpikes, fire_mask, vNSpikes, _mm512_set1_epi32(1));
                        _mm512_storeu_si512(&(results->n_spks[b_neuron]), vNSpikes); 

                        // store generated spikes in matrix
                        _mm512_mask_storeu_epi8(&(snn->spk_matrix[batch_idx]), fire_mask, _mm512_set1_epi8(1));

                        // if TB-STDP is used, store which neurons fired
                        if(conf->learn){
                            
                            _mm512_mask_storeu_epi8(&(snn->post_fired[b_neuron]), fire_mask, _mm512_set1_epi8(1));
                        }
                    }

                    _mm512_storeu_si512(&(snn->r_period_remain[b_neuron]), vRefPer); // store r_remain
                }
            }
        }
    }
    #else
    {
        
        #pragma omp parallel for num_threads(P) private(i, b, neuron_index, g_neuron_index)
        for(i = 0; i<N; i++){
                
            neuron_index = i;
            g_neuron_index = i * B;
            
            // compute index to write in the spk matrix
            // [T * B * (iN + N)]
            size_t idx = (iN + N) * B * t + B * (iN + neuron_index); // index of the first neuron in the spike matrix in time step t

            for(b = 0; b<B; b++){
            
                // reduce refractary period
                snn->r_period_remain[g_neuron_index + b] --;

                // check whether fired
                if(snn->v[g_neuron_index + b] >= snn->v_thresh[neuron_index]){
            
                    snn->spk_matrix[idx + b] = 1;

                    // reinit neuron values
                    snn->r_period_remain[g_neuron_index + b] = snn->r_period[neuron_index];
                    snn->v[g_neuron_index + b] = snn->v_rest[neuron_index]; // reinit v_rest

                    results->n_spks[g_neuron_index + b] += 1;

                    // store that the neuron fired for TB-STDP and update the trace
                    if(conf->learn){

                        snn->post_fired[g_neuron_index + b] = 1;
                    }
                }
            }
        }
    }
    #endif
}