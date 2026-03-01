#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <math.h>
#include <time.h>

#include <omp.h>
#include <immintrin.h> 

#include "config/config_loader.h"
#include "networks/snn.h"
#include "neuron_models/neuron_models.h"

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
            snn->post_fired[g_index + b] = 0;

            if(conf->learn == 1){
                snn->post_trace[g_index + b] = 0.0;
            }
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