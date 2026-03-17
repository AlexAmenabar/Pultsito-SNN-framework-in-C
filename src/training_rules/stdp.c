#include <math.h>
#include <stdlib.h>

#include <omp.h>
#include <immintrin.h> 

// public headers
#include "training_rules/stdp.h"
#include "networks/snn.h"
#include "config/config_loader.h"

void process_trace_based_STDP_batch(GPU_SNN_t *snn, simulation_configuration_t *conf){

    size_t N, S, P, B;
    size_t i, b;

    N = snn->n_neurons;
    S = snn->n_synapses;
    B = conf->batch_size;
    P = conf->n_process;

    float pA = 0.1f, mA = 0.1f, mu = 1.0f, decay = 0.6f;
    
    #if defined AVX512
    {

        if(B == 8){

            // define constants
            const __m256 vOne   = _mm256_set1_ps(1.0f);
            const __m256 vDecay = _mm256_set1_ps(decay);
            const __m256 vpA = _mm256_set1_ps(pA);
            const __m256 vmA = _mm256_set1_ps(mA);
            const __m128i vZero = _mm_set1_epi8(0); // 16 bytes of zeros

            size_t n_blks = B / 8, blk;


            #pragma omp parallel for num_threads(P) private(i, b, blk)
            for(i = 0; i<N; i++){

                // loop over neuron input synapses
                size_t base_synapse = snn->neuron_input_synapses_offset[i];
                size_t n_iS = snn->n_neuron_input_synapses[i];

                // global neuron index
                size_t g_post_neuron = i * B;

                __mmask8 mPostFired[n_blks];
                __m256 vPostTrace[n_blks];

                // loop over batch to update post synaptic trace (N * B)
                for(blk = 0; blk<n_blks; blk++){

                    // compute neuron index
                    size_t post_idx = g_post_neuron + blk * 16;
                    
                    // load if fired and trace
                    //__m128i vPostFired = _mm_loadu_si128((__m128i*)(&(snn->post_fired[post_idx]))); // char
                    __m128i vPostFired = _mm_loadu_epi8((__m128i*)(&(snn->post_fired[post_idx])));
                    mPostFired[blk] = _mm_cmpgt_epi8_mask(vPostFired, vZero); // mask
                    vPostTrace[blk] = _mm256_loadu_ps(&(snn->post_trace[post_idx])); // load trace
                    
                    // update trace and store trace
                    __m256 vPostDecay = _mm256_mul_ps(vPostTrace[blk], vDecay);
                    vPostTrace[blk] = _mm256_mask_blend_ps(mPostFired[blk], vPostDecay, vOne); // select conditionally
                    _mm256_storeu_ps(&(snn->post_trace[post_idx]), vPostTrace[blk]); // store trace
                }

                // loop over input synapses
                for(size_t j = 0; j<n_iS; j++){

                    // get synapse index
                    size_t synapse_index = base_synapse + j;
                    size_t g_synapse_index = synapse_index * B;

                    // loop over batch
                    for(blk = 0; blk<n_blks; blk++){

                        // compute post and pre indexes
                        size_t post_idx = g_post_neuron + blk * 16;
                        size_t pre_idx  = g_synapse_index + blk * 16;

                        // load wether pre fired and trace
                        //__m128i vPreFired = _mm_loadu_si128((__m128i*)(&(snn->pre_fired[pre_idx]))); // SSE
                        __m128i vPreFired = _mm_loadu_epi8(&(snn->pre_fired[pre_idx])); // AVX512
                        __mmask8 mPreFired = _mm_cmpgt_epi8_mask(vPreFired, vZero); // mask

                        // update pre trace
                        __m256 vPreTrace = _mm256_loadu_ps(&(snn->pre_trace[pre_idx])); // load trace
                        __m256 vPreDecay = _mm256_mul_ps(vPreTrace, vDecay);
                        vPreTrace = _mm256_mask_blend_ps(mPreFired, vPreDecay, vOne);
                        _mm256_storeu_ps(&(snn->pre_trace[pre_idx]), vPreTrace); // store trace

                        // avoid updating v and dw if no one neuron fired
                        if(mPreFired || mPostFired[blk]){ // time reduced to 1/2 (more or less)
                           
                            // update dw
                            __m256 vW = _mm256_loadu_ps(&(snn->w[pre_idx])); // load weight
                            __m256 vPot = _mm256_maskz_mul_ps(mPostFired[blk], vpA, _mm256_mul_ps(_mm256_sub_ps(vOne, vW), vPreTrace)); // LTP
                            __m256 vDep = _mm256_maskz_mul_ps(mPreFired, vmA, _mm256_mul_ps(vW, vPostTrace[blk])); // LTD
                            __m256 vdW = _mm256_sub_ps(vPot, vDep); // dw


                            // update and store w and dw
                            _mm256_storeu_ps(&(snn->w[pre_idx]), _mm256_add_ps(vW, vdW));
                            _mm256_storeu_ps(&(snn->dw[pre_idx]), _mm256_add_ps(_mm256_loadu_ps(&(snn->dw[pre_idx])), vdW));
                        }

                        // reinit pre fired to 0
                        if(mPreFired){
                            //_mm_storeu_si128((__m128i*)&(snn->pre_fired[pre_idx]), vZero); // SSE
                            _mm_mask_storeu_epi8(&(snn->pre_fired[pre_idx]), mPreFired, vZero); // AVX 512           
                        }         
                    }
                }

                // reinit post fired
                for(blk = 0; blk<n_blks; blk++){

                    if(mPostFired[blk]){
                        //_mm_storeu_si128((__m128i*)&(snn->post_fired[g_post_neuron + blk * 16]), vZero); // SSE  
                        _mm_mask_storeu_epi8(&(snn->post_fired[g_post_neuron + blk * 16]), mPostFired[blk], vZero); // AVX 512           
                    }
                }
            }        
        }
        else if(B % 16 == 0){

            // define constants
            const __m512 vOne   = _mm512_set1_ps(1.0f);
            const __m512 vDecay = _mm512_set1_ps(decay);
            const __m512 vpA = _mm512_set1_ps(pA);
            const __m512 vmA = _mm512_set1_ps(mA);
            const __m128i vZero = _mm_set1_epi8(0); // 16 bytes of zeros

            size_t n_blks = B / 16, blk;


            #pragma omp parallel for num_threads(P) private(i, b, blk)
            for(i = 0; i<N; i++){

                // loop over neuron input synapses
                size_t base_synapse = snn->neuron_input_synapses_offset[i];
                size_t n_iS = snn->n_neuron_input_synapses[i];

                // global neuron index
                size_t g_post_neuron = i * B;

                __mmask16 mPostFired[n_blks];
                __m512 vPostTrace[n_blks];

                // loop over batch to update post synaptic trace (N * B)
                for(blk = 0; blk<n_blks; blk++){

                    // compute neuron index
                    size_t post_idx = g_post_neuron + blk * 16;
                    
                    // load if fired and trace
                    //__m128i vPostFired = _mm_loadu_si128((__m128i*)(&(snn->post_fired[post_idx]))); // char
                    __m128i vPostFired = _mm_loadu_epi8((__m128i*)(&(snn->post_fired[post_idx])));
                    mPostFired[blk] = _mm_cmpgt_epi8_mask(vPostFired, vZero); // mask
                    vPostTrace[blk] = _mm512_loadu_ps(&(snn->post_trace[post_idx])); // load trace
                    
                    // update trace and store trace
                    __m512 vPostDecay = _mm512_mul_ps(vPostTrace[blk], vDecay);
                    vPostTrace[blk] = _mm512_mask_blend_ps(mPostFired[blk], vPostDecay, vOne); // select conditionally
                    _mm512_storeu_ps(&(snn->post_trace[post_idx]), vPostTrace[blk]); // store trace
                }

                // loop over input synapses
                for(size_t j = 0; j<n_iS; j++){

                    // get synapse index
                    size_t synapse_index = base_synapse + j;
                    size_t g_synapse_index = synapse_index * B;

                    // loop over batch
                    for(blk = 0; blk<n_blks; blk++){

                        // compute post and pre indexes
                        size_t post_idx = g_post_neuron + blk * 16;
                        size_t pre_idx  = g_synapse_index + blk * 16;

                        // load wether pre fired and trace
                        //__m128i vPreFired = _mm_loadu_si128((__m128i*)(&(snn->pre_fired[pre_idx]))); // SSE
                        __m128i vPreFired = _mm_loadu_epi8(&(snn->pre_fired[pre_idx])); // AVX512
                        __mmask16 mPreFired = _mm_cmpgt_epi8_mask(vPreFired, vZero); // mask

                        // update pre trace
                        __m512 vPreTrace = _mm512_loadu_ps(&(snn->pre_trace[pre_idx])); // load trace
                        __m512 vPreDecay = _mm512_mul_ps(vPreTrace, vDecay);
                        vPreTrace = _mm512_mask_blend_ps(mPreFired, vPreDecay, vOne);
                        _mm512_storeu_ps(&(snn->pre_trace[pre_idx]), vPreTrace); // store trace

                        // avoid updating v and dw if no one neuron fired
                        if(mPreFired || mPostFired[blk]){ // time reduced to 1/2 (more or less)
                           
                            // update dw
                            __m512 vW = _mm512_loadu_ps(&(snn->w[pre_idx])); // load weight
                            __m512 vPot = _mm512_maskz_mul_ps(mPostFired[blk], vpA, _mm512_mul_ps(_mm512_sub_ps(vOne, vW), vPreTrace)); // LTP
                            __m512 vDep = _mm512_maskz_mul_ps(mPreFired, vmA, _mm512_mul_ps(vW, vPostTrace[blk])); // LTD
                            __m512 vdW = _mm512_sub_ps(vPot, vDep); // dw


                            // update and store w and dw
                            _mm512_storeu_ps(&(snn->w[pre_idx]), _mm512_add_ps(vW, vdW));
                            _mm512_storeu_ps(&(snn->dw[pre_idx]), _mm512_add_ps(_mm512_loadu_ps(&(snn->dw[pre_idx])), vdW));
                        }

                        // reinit pre fired to 0
                        if(mPreFired){
                            //_mm_storeu_si128((__m128i*)&(snn->pre_fired[pre_idx]), vZero); // SSE
                            _mm_storeu_epi8(&(snn->pre_fired[pre_idx]), vZero); // AVX 512           
                        }         
                    }
                }

                // reinit post fired
                for(blk = 0; blk<n_blks; blk++){

                    if(mPostFired[blk]){
                        //_mm_storeu_si128((__m128i*)&(snn->post_fired[g_post_neuron + blk * 16]), vZero); // SSE  
                        _mm_storeu_epi8(&(snn->post_fired[g_post_neuron + blk * 16]), vZero); // AVX 512           
                    }
                }
            }
        }
    }
    #else
    {
        // loop over synapses, update traces and weights
        #pragma omp parallel for num_threads(P) private(i, b)
        for(i = 0; i<N; i++){

            // loop over neuron input synapses
            size_t base_synapse = snn->neuron_input_synapses_offset[i];
            size_t n_iS = snn->n_neuron_input_synapses[i];
            size_t post_idx, pre_idx;
            size_t synapse_index;
            size_t g_synapse_index;
            
            // compute neuron global index
            size_t g_post_neuron = i * B;

            // update post trace (depends only on the neuron)
            for(b = 0; b<B; b++){
                
                post_idx = g_post_neuron + b;
                snn->post_trace[post_idx] = snn->post_fired[post_idx] ? 1.0f : snn->post_trace[post_idx] * decay;
            }

            // loop over input synapses, update presynaptic traces (depends on synapses) and update weights
            for(size_t j = 0; j<n_iS; j++){

                synapse_index = base_synapse + j;
                g_synapse_index = synapse_index * B;

                for(b = 0; b<B; b++){

                    post_idx = g_post_neuron + b;
                    pre_idx  = g_synapse_index + b;

                    // update pre_trace
                    snn->pre_trace[pre_idx] = snn->pre_fired[pre_idx] ? 1.0f : snn->pre_trace[pre_idx] * decay;


                    // update w: avoid computation if no neuron fired
                    float dPot = 0.0, dDep = 0.0;
                    if(snn->post_fired[post_idx]){

                        dPot = pA * (1.0f - snn->w[pre_idx]) * snn->pre_trace[pre_idx];
                    }
                    if(snn->pre_fired[pre_idx]){

                        dDep = mA * (snn->w[pre_idx]) * snn->post_trace[post_idx];
                        
                        // reinit pre_fired
                        snn->pre_fired[pre_idx] = 0;
                    }
                    
                    if(dPot != 0 || dDep != 0){

                        float dw = dPot - dDep;
                            
                        snn->w[pre_idx] += dw;
                        snn->dw[pre_idx] += dw;
                    }
                }
            }

            // reinit whether post neuron fired for next iterations
            for(b = 0; b<B; b++){
                
                snn->post_fired[g_post_neuron + b] = 0;
            }        
        }
    }
    #endif
}