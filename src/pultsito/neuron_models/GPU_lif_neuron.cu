#include <cuda.h>
#include <cuda_runtime.h>
#include <stdio.h>
#include <omp.h>


#include "snn_library.h"
#include "neuron_models/GPU_lif_neuron.cuh"
#include "cuda/GPU_simulations.cuh"


#define cudaCheckError(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
{
    if (code != cudaSuccess)
    {
        fprintf(stderr, "CUDA Error: %s %s %d\n", cudaGetErrorString(code), file, line);
        if (abort) exit(code);
    }
}

//__constant__ size_t N;
//__constant__ size_t iN;
//__constant__ size_t S;
//__constant__ size_t T;
//__constant__ size_t LT;
//__constant__ size_t n_features; 
__constant__ int n_samples;


// initialize neurons before starting the simulation
__global__ void initialize_neurons_batch(GPU_SNN_t *snn, size_t batch_size){

    // get thread id
    size_t threadId = 
        ((size_t)blockIdx.x  + (size_t)gridDim.x  * (size_t)blockIdx.y + (size_t)gridDim.x * (size_t)gridDim.y * (size_t)blockIdx.z) *
        ((size_t)blockDim.x * (size_t)blockDim.y * (size_t)blockDim.z) +
        ((size_t)threadIdx.x + (size_t)blockDim.x * (size_t)threadIdx.y + (size_t)blockDim.x * (size_t)blockDim.y * (size_t)threadIdx.z);

    if(threadId < snn->n_neurons * batch_size){

        snn->v[threadId] = snn->v_rest[threadId / batch_size];
        snn->r_period_remain[threadId] = 0;
        snn->post_fired[threadId] = 0;
        snn->post_trace[threadId] = 0.0;
        snn->arrI[threadId] = 0.0;
    }    
}

// initialize synapses before starting the simulation
__global__ void initialize_synapses_batch(GPU_SNN_t *snn, size_t batch_size){

    // get thread id
    size_t threadId = 
        (size_t)((size_t)blockIdx.x  + (size_t)gridDim.x  * (size_t)blockIdx.y + (size_t)gridDim.x * (size_t)gridDim.y * (size_t)blockIdx.z) *
        (size_t)((size_t)blockDim.x * (size_t)blockDim.y * (size_t)blockDim.z) +
        (size_t)((size_t)threadIdx.x + (size_t)blockDim.x * (size_t)threadIdx.y + (size_t)blockDim.x * (size_t)blockDim.y * (size_t)threadIdx.z);


    //printf(" thread = %llu / %llu\n", threadId, snn->n_synapses * batch_size);

    if(threadId < snn->n_synapses * batch_size){

        size_t synapse_index = threadId / batch_size;
        snn->w[threadId] = snn->init_w[synapse_index];
        snn->dw[threadId] = 0.0;
        snn->pre_fired[threadId] = 0;
        snn->pre_trace[threadId] = 0.0;
    }    
}


__global__ void initialize_batch_matrix_GPU(tmp_batch_cpu_t *tmp_batch, GPU_dataset_t *dataset, size_t bidx, size_t n_networks, size_t iN, size_t batch_size, size_t T){

    // get thread id: iN * batch_size * T
    size_t threadId = 
        (blockIdx.x  + gridDim.x  * blockIdx.y + gridDim.x * gridDim.y * blockIdx.z) *
        (blockDim.x * blockDim.y * blockDim.z) +
        (threadIdx.x + blockDim.x * threadIdx.y + blockDim.x * blockDim.y * threadIdx.z);

    if(threadId < iN * batch_size){

        size_t fsidx = bidx * batch_size; // first sample to be processed
        size_t b = threadId % batch_size; // batch element
        size_t sidx = fsidx + b; // sample id = first sample in the batch + batch element

        size_t n_features = dataset->n_features;

        if(sidx < dataset->n_samples){
            
            size_t fidx = threadId / batch_size; // neuron index is the feature index
            
            // reinitialize batch_matrix
            for(size_t j = 0; j<T; j++){
                
                tmp_batch->spikes[(iN * batch_size * j) + (batch_size * fidx) + b] = 0; 
            }

            // get the feature offset in the array of spikes
            size_t feature_offset = dataset->feature_offset[sidx * n_features + fidx];

            // get the number of spikes of the feature
            size_t n_spikes_feature = dataset->n_spikes_per_feature[sidx * n_features + fidx];
    
            // loop over spikes and store in the matrix
            for(size_t j = 0; j<n_spikes_feature; j++){

                // get time step of the spike
                size_t t = dataset->spikes[feature_offset + j];

                // store the spike in the matrix of [iN * B * T]
                tmp_batch->spikes[(iN * batch_size * t) + (batch_size * fidx) + b] = 1; 
            }
        }
        
    }
}

__global__ void load_batch_time_step_in_SNN_GPU(GPU_SNN_t *gpu_snn, tmp_batch_cpu_t *tmp_batch, size_t iN, size_t N, size_t batch_size, size_t lt, size_t t, size_t bidx){

    // get thread id: iN * batch_size * T
    size_t threadId = 
        (blockIdx.x  + gridDim.x  * blockIdx.y + gridDim.x * gridDim.y * blockIdx.z) *
        (blockDim.x * blockDim.y * blockDim.z) +
        (threadIdx.x + blockDim.x * threadIdx.y + blockDim.x * blockDim.y * threadIdx.z);


    if(threadId < (iN + N) * batch_size){

        size_t i = threadId / batch_size; // neuron index
        size_t b = threadId % batch_size; // batch position

        size_t input = i < iN ? 1 : 0; //  

        // compute matrix indexes
        size_t spk_mtr_idx = (iN + N) * batch_size * lt + (batch_size * i) + b;

        if(input){

            size_t batch_mtr_idx = (iN * batch_size * t) + (batch_size * i) + b;
            gpu_snn->spk_matrix[spk_mtr_idx] = tmp_batch->spikes[batch_mtr_idx]; 
        }
        else{

            gpu_snn->spk_matrix[spk_mtr_idx] = 0;
        }
    }
}


__global__ void process_input_currect_batch(GPU_SNN_t *gpu_snn, size_t iN, size_t N, size_t batch_size, size_t t, size_t gt, size_t learn){

    // get thread id: iN * batch_size * T
    size_t threadId = 
        (blockIdx.x  + gridDim.x  * blockIdx.y + gridDim.x * gridDim.y * blockIdx.z) *
        (blockDim.x * blockDim.y * blockDim.z) +
        (threadIdx.x + blockDim.x * threadIdx.y + blockDim.x * blockDim.y * threadIdx.z);


    // initialize shared memory
    extern __shared__ float sharedI[];

    size_t localThreadId = threadId % 512; // TODO: tmp

    size_t i, b, p;
    i = threadId / (batch_size * thrN);
    b = (threadId / thrN) % batch_size;
    p = threadId % thrN;

    if(threadId < N * thrN * batch_size){

        sharedI[localThreadId] = 0.0;

        // variables declaration
        size_t synapse_index, g_synapse_index, in_neuron_index;
        int delay, spk_time; 
        float spk;  

        size_t neuron_index = i; // for not copied variables
        size_t g_neuron_index = i * batch_size; // for copied variables

        // get neuron information
        size_t iS = gpu_snn->n_neuron_input_synapses[neuron_index]; 
        size_t base_synapse = gpu_snn->neuron_input_synapses_offset[neuron_index]; // there are several copies of each synapse, so the offset depends
    
        float I = 0.0;

        // check wether the neuron is in refractory period: pre_fired flag not activated although it should be
        if(gpu_snn->r_period_remain[g_neuron_index + b] <= 0 || learn){
        
            // calculate synapses to be processed
            size_t first_synapse, last_synapse, n_synapses, r_synapses;
            
            n_synapses = iS / thrN;
            r_synapses = iS % thrN;
            
            // compute first and last synapses
            first_synapse = n_synapses * p;
            last_synapse = n_synapses * (p+1);

            if(p == thrN - 1){

                last_synapse += r_synapses;
            }

            // loop over input synapses
            for(size_t j = first_synapse; j<last_synapse; j++){

                // get synapse index. In not copied variables use synapse_index, in copied g_synapse_index 
                synapse_index = base_synapse + j;
                g_synapse_index = synapse_index * batch_size; // scale base synapse, since now there are B copies of each one
                
                // get input neuron index (not copied)
                in_neuron_index = gpu_snn->pre_neuron_index[synapse_index]; 

                // get synapse delay (not copied)
                delay = gpu_snn->delay[synapse_index];

                // get spike time
                spk_time = (int)t - delay;
            
                // correct spk_time
                if(spk_time < 0 && gt >= gpu_snn->LT){ // CHECK
                    spk_time = (int)(gpu_snn->LT) + spk_time;
                }

                // process spike 
                if(spk_time >= 0){
                    
                    // index to load the spike from
                    size_t index = ((iN + N) * batch_size * (size_t)spk_time) + ((in_neuron_index * batch_size));

                    // get spike
                    spk = (float)(gpu_snn->spk_matrix[index + b]); 

                    // update I if r period remain <= 0 and spike received
                    if(gpu_snn->r_period_remain[g_neuron_index + b] <= 0){
                        
                        I += gpu_snn->w[g_synapse_index + b] * spk; // spk = 0 / 1
                    }

                    // store if the presynaptic neuron fired for TR-STDP
                    if(learn && (int)spk == 1){

                        gpu_snn->pre_fired[g_synapse_index + b] = (char)spk;
                    }
                }
            }
        }

        // write I in shared memory
        sharedI[localThreadId] = I;
    }

    // sync threads
    __syncthreads();

    // the first neuron accumulates I and writes in sahred memory
    if(threadId < N * thrN * batch_size && p == 0){

        float I = 0.0;
        for(size_t j = 0; j<thrN; j++){

            I += sharedI[localThreadId + j];
        }

        gpu_snn->arrI[i * batch_size + b] = I;
    }
}


__global__ void process_V_batch(GPU_SNN_t *gpu_snn, size_t N, size_t batch_size){

    // get thread id: iN * batch_size * T
    size_t threadId = 
        (blockIdx.x  + gridDim.x  * blockIdx.y + gridDim.x * gridDim.y * blockIdx.z) *
        (blockDim.x * blockDim.y * blockDim.z) +
        (threadIdx.x + blockDim.x * threadIdx.y + blockDim.x * blockDim.y * threadIdx.z);


    if(threadId < N * batch_size){

        size_t i = threadId / batch_size; // neuron index
        size_t b = threadId % batch_size; // batch position

        size_t neuron_index = i;
        size_t g_neuron_index = i * batch_size;

        const float alpha = 0.95f;
        const float beta = 0.05f;

        // loop over batches
        if (gpu_snn->r_period_remain[g_neuron_index + b] <= 0){
            
            gpu_snn->v[g_neuron_index + b] = alpha * gpu_snn->v[g_neuron_index + b] + beta * gpu_snn->v_rest[neuron_index] + gpu_snn->arrI[g_neuron_index + b];
        }
    }
}


__global__ void process_firing_batch(GPU_SNN_t *gpu_snn, GPU_results_t *gpu_results, size_t iN, size_t N, size_t batch_size, size_t t, size_t learn){

    // get thread id: iN * batch_size * T
    size_t threadId = 
        (blockIdx.x  + gridDim.x  * blockIdx.y + gridDim.x * gridDim.y * blockIdx.z) *
        (blockDim.x * blockDim.y * blockDim.z) +
        (threadIdx.x + blockDim.x * threadIdx.y + blockDim.x * blockDim.y * threadIdx.z);


    if(threadId < N * batch_size){

        size_t i = threadId / batch_size; // neuron index
        size_t b = threadId % batch_size; // batch position

        size_t neuron_index = i;
        size_t g_neuron_index = i * batch_size;
        
        // compute index to write in the spk matrix
        // [T * B * (iN + N)]
        size_t idx = (iN + N) * batch_size * t + batch_size * (iN + neuron_index); // index of the first neuron in the spike matrix in time step t
        
        // reduce refractary period
        gpu_snn->r_period_remain[g_neuron_index + b] --;

        // check whether fired
        if(gpu_snn->v[g_neuron_index + b] >= gpu_snn->v_thresh[neuron_index]){
    
            gpu_snn->spk_matrix[idx + b] = 1;

            // reinit neuron values
            gpu_snn->r_period_remain[g_neuron_index + b] = gpu_snn->r_period[neuron_index];
            gpu_snn->v[g_neuron_index + b] = gpu_snn->v_rest[neuron_index]; // reinit v_rest

            gpu_results->n_spks[g_neuron_index + b] += 1;

            // store that the neuron fired for TB-STDP and update the trace
            if(learn){

                gpu_snn->post_fired[g_neuron_index + b] = 1;
            }
        }
    }
}

__global__ void total_LIF_batch(GPU_SNN_t *gpu_snn, GPU_results_t *gpu_results, size_t iN, size_t N, size_t batch_size, size_t t, size_t gt, size_t learn){

    // get thread id: iN * batch_size * T
    size_t threadId = 
        (blockIdx.x  + gridDim.x  * blockIdx.y + gridDim.x * gridDim.y * blockIdx.z) *
        (blockDim.x * blockDim.y * blockDim.z) +
        (threadIdx.x + blockDim.x * threadIdx.y + blockDim.x * blockDim.y * threadIdx.z);


    if(threadId < N * batch_size){

        size_t i = threadId / batch_size; // neuron index
        size_t b = threadId % batch_size; // batch position

        size_t neuron_index = i; // for not copied variables
        size_t g_neuron_index = i * batch_size; // for copied variables

        // variables declaration
        size_t synapse_index, g_synapse_index, in_neuron_index;
        int delay, spk_time; 
        float spk;  

        // get neuron information
        size_t iS = gpu_snn->n_neuron_input_synapses[neuron_index]; 
        size_t base_synapse = gpu_snn->neuron_input_synapses_offset[neuron_index]; // there are several copies of each synapse, so the offset depends
    
        
        float I = 0.0;

        // check wether the neuron is in refractory period: pre_fired flag not activated although it should be
        if(gpu_snn->r_period_remain[g_neuron_index + b] <= 0 || learn){
        
            // loop over input synapses
            for(size_t j=0; j<iS; j++){

                // get synapse index. In not copied variables use synapse_index, in copied g_synapse_index 
                synapse_index = base_synapse + j;
                g_synapse_index = synapse_index * batch_size; // scale base synapse, since now there are B copies of each one
                
                // get input neuron index (not copied)
                in_neuron_index = gpu_snn->pre_neuron_index[synapse_index]; 

                // get synapse delay (not copied)
                delay = gpu_snn->delay[synapse_index];

                // get spike time
                spk_time = (int)t - delay;
            
                // correct spk_time
                if(spk_time < 0 && gt >= gpu_snn->LT){ // CHECK
                    spk_time = (int)(gpu_snn->LT) + spk_time;
                }

                // process spike 
                if(spk_time >= 0){
                    
                    // index to load the spike from
                    size_t index = ((iN + N) * batch_size * (size_t)spk_time) + ((in_neuron_index * batch_size));

                    // get spike
                    spk = (float)(gpu_snn->spk_matrix[index + b]); 

                    // update I if r period remain <= 0 and spike received
                    if(gpu_snn->r_period_remain[g_neuron_index + b] <= 0){
                        
                        I += gpu_snn->w[g_synapse_index + b] * spk; // spk = 0 / 1
                    }

                    // store if the presynaptic neuron fired for TR-STDP
                    if(learn && (int)spk == 1){

                        gpu_snn->pre_fired[g_synapse_index + b] = (char)spk;
                    }
                }
            }
        }

        gpu_snn->arrI[g_neuron_index + b] = I;


        const float alpha = 0.95f;
        const float beta = 0.05f;

        // loop over batches
        if (gpu_snn->r_period_remain[g_neuron_index + b] <= 0){
            
            gpu_snn->v[g_neuron_index + b] = alpha * gpu_snn->v[g_neuron_index + b] + beta * gpu_snn->v_rest[neuron_index] + gpu_snn->arrI[g_neuron_index + b];
        }


        // compute index to write in the spk matrix
        // [T * B * (iN + N)]
        size_t idx = (iN + N) * batch_size * t + batch_size * (iN + neuron_index); // index of the first neuron in the spike matrix in time step t
        
        // reduce refractary period
        gpu_snn->r_period_remain[g_neuron_index + b] --;

        // check whether fired
        if(gpu_snn->v[g_neuron_index + b] >= gpu_snn->v_thresh[neuron_index]){
    
            gpu_snn->spk_matrix[idx + b] = 1;

            // reinit neuron values
            gpu_snn->r_period_remain[g_neuron_index + b] = gpu_snn->r_period[neuron_index];
            gpu_snn->v[g_neuron_index + b] = gpu_snn->v_rest[neuron_index]; // reinit v_rest

            gpu_results->n_spks[g_neuron_index + b] += 1;

            //printf(" > Neuron %zu (batch %zu): %d\n", neuron_index, b, gpu_results->n_spks[g_neuron_index + b]);

            // store that the neuron fired for TB-STDP and update the trace
            if(learn){

                gpu_snn->post_fired[g_neuron_index + b] = 1;
            }
        }
    }
}

__global__ void compute_post_traces_batch(GPU_SNN_t *gpu_snn, size_t N, size_t batch_size){

    // get thread id: iN * batch_size * T
    size_t threadId = 
        (blockIdx.x  + gridDim.x  * blockIdx.y + gridDim.x * gridDim.y * blockIdx.z) *
        (blockDim.x * blockDim.y * blockDim.z) +
        (threadIdx.x + blockDim.x * threadIdx.y + blockDim.x * blockDim.y * threadIdx.z);

    if(threadId < N * batch_size){
    
        size_t i = threadId / batch_size; // neuron index
        size_t b = threadId % batch_size; // batch position
        size_t g_neuron_index = i * batch_size; // for copied variables 
        float decay = 0.6f;

        // update post trace (depends only on the neuron)
        gpu_snn->post_trace[g_neuron_index + b] = gpu_snn->post_fired[g_neuron_index + b] ? 1.0f : gpu_snn->post_trace[g_neuron_index + b] * decay;
    }
}


__global__ void reinit_post_fired_batch(GPU_SNN_t *gpu_snn, size_t N, size_t batch_size){

    // get thread id: iN * batch_size * T
    size_t threadId = 
        (blockIdx.x  + gridDim.x  * blockIdx.y + gridDim.x * gridDim.y * blockIdx.z) *
        (blockDim.x * blockDim.y * blockDim.z) +
        (threadIdx.x + blockDim.x * threadIdx.y + blockDim.x * blockDim.y * threadIdx.z);

    if(threadId < N * batch_size){
    
        size_t i = threadId / batch_size; // neuron index
        size_t b = threadId % batch_size; // batch position
        size_t g_neuron_index = i * batch_size; // for copied variables 

        // update post trace (depends only on the neuron)
        gpu_snn->post_fired[g_neuron_index + b] = 0;
    }
}


__global__ void trace_based_STDP_batch(GPU_SNN_t *gpu_snn, size_t S, size_t batch_size){

    
    // get thread id: iN * batch_size * T
    size_t threadId = 
        (size_t)((size_t)blockIdx.x  + (size_t)gridDim.x  * (size_t)blockIdx.y + (size_t)gridDim.x * (size_t)gridDim.y * (size_t)blockIdx.z) *
        (size_t)((size_t)blockDim.x * (size_t)blockDim.y * (size_t)blockDim.z) +
        (size_t)((size_t)threadIdx.x + (size_t)blockDim.x * (size_t)threadIdx.y + (size_t)blockDim.x * (size_t)blockDim.y * (size_t)threadIdx.z);

    if(threadId < S * batch_size){

        size_t synapse_index = threadId / batch_size; // neuron index
        size_t b = threadId % batch_size; // batch position
        size_t g_synapse_index = synapse_index * batch_size; // for copied variables
        
        size_t neuron_index, g_neuron_index;
        neuron_index = gpu_snn->post_neuron_index[synapse_index] - gpu_snn->n_input_neurons; // - virtual neurons
        g_neuron_index = neuron_index * batch_size;

        // constants
        const float pA = 0.1f, mA = 0.1f, mu = 1.0f, decay = 0.6f;

        // update pre_trace
        gpu_snn->pre_trace[g_synapse_index + b] = gpu_snn->pre_fired[g_synapse_index + b] ? 1.0f : gpu_snn->pre_trace[g_synapse_index + b] * decay;


        // update w: avoid computation if no neuron fired
        float dPot = 0.0, dDep = 0.0;
        if(gpu_snn->post_fired[g_neuron_index + b]){

            dPot = pA * (1.0f - gpu_snn->w[g_synapse_index + b]) * gpu_snn->pre_trace[g_synapse_index + b];
        }
        if(gpu_snn->pre_fired[g_synapse_index + b]){

            dDep = mA * (gpu_snn->w[g_synapse_index + b]) * gpu_snn->post_trace[g_neuron_index + b];
            
            // reinit pre_fired
            gpu_snn->pre_fired[g_synapse_index + b] = 0;
        }
        
        if(dPot != 0 || dDep != 0){

            float dw = dPot - dDep;
                
            gpu_snn->w[g_synapse_index + b] += dw;
            gpu_snn->dw[g_synapse_index + b] += dw;
        }
    }
}


__global__ void trace_based_STDP_batch_old(GPU_SNN_t *gpu_snn, size_t N, size_t batch_size){

    
    // get thread id: iN * batch_size * T
    size_t threadId = 
        (blockIdx.x  + gridDim.x  * blockIdx.y + gridDim.x * gridDim.y * blockIdx.z) *
        (blockDim.x * blockDim.y * blockDim.z) +
        (threadIdx.x + blockDim.x * threadIdx.y + blockDim.x * blockDim.y * threadIdx.z);

    if(threadId < N * batch_size){

        size_t i = threadId / batch_size; // neuron index
        size_t b = threadId % batch_size; // batch position

        size_t neuron_index = i; // for not copied variables
        size_t g_neuron_index = i * batch_size; // for copied variables


        // loop over neuron input synapses
        size_t base_synapse = gpu_snn->neuron_input_synapses_offset[i];
        size_t n_iS = gpu_snn->n_neuron_input_synapses[i];
        size_t post_idx, pre_idx;
        size_t synapse_index;
        size_t g_synapse_index;

        float pA = 0.1f, mA = 0.1f, mu = 1.0f, decay = 0.6f;

        // update post trace (depends only on the neuron)
        gpu_snn->post_trace[g_neuron_index + b] = gpu_snn->post_fired[g_neuron_index + b] ? 1.0f : gpu_snn->post_trace[g_neuron_index + b] * decay;

        // loop over input synapses, update presynaptic traces (depends on synapses) and update weights
        for(size_t j = 0; j<n_iS; j++){

            synapse_index = base_synapse + j;
            g_synapse_index = synapse_index * batch_size;

            // update pre_trace
            gpu_snn->pre_trace[g_synapse_index + b] = gpu_snn->pre_fired[g_synapse_index + b] ? 1.0f : gpu_snn->pre_trace[g_synapse_index + b] * decay;


            // update w: avoid computation if no neuron fired
            float dPot = 0.0, dDep = 0.0;
            if(gpu_snn->post_fired[g_neuron_index + b]){

                dPot = pA * (1.0f - gpu_snn->w[g_synapse_index + b]) * gpu_snn->pre_trace[g_synapse_index + b];
            }
            if(gpu_snn->pre_fired[g_synapse_index + b]){

                dDep = mA * (gpu_snn->w[g_synapse_index + b]) * gpu_snn->post_trace[g_neuron_index + b];
                
                // reinit pre_fired
                gpu_snn->pre_fired[g_synapse_index + b] = 0;
            }
            
            if(dPot != 0 || dDep != 0){

                float dw = dPot - dDep;
                    
                gpu_snn->w[g_synapse_index + b] += dw;
                gpu_snn->dw[g_synapse_index + b] += dw;
            }
        }
            
        gpu_snn->post_fired[g_neuron_index + b] = 0;
    }
}


__global__ void update_weights_batch(GPU_SNN_t *gpu_snn, size_t S, size_t batch_size){

    // get thread id: iN * batch_size * T
    size_t threadId = 
        (blockIdx.x  + gridDim.x  * blockIdx.y + gridDim.x * gridDim.y * blockIdx.z) *
        (blockDim.x * blockDim.y * blockDim.z) +
        (threadIdx.x + blockDim.x * threadIdx.y + blockDim.x * blockDim.y * threadIdx.z);


    if(threadId < S){

        size_t synapse_index = threadId; // synapse index
        size_t g_synapse_index = synapse_index * batch_size; // synapse index
    
        float dw = 0.0;

        // sum dw of the batch for each synapse
        for(size_t b = 0; b<batch_size; b++){

            dw += gpu_snn->dw[g_synapse_index + b];
        }

        // compute mean dw
        dw = dw / (float)batch_size;

        // update initial w and w for all copies
        for(size_t b = 0; b<batch_size; b++){

            gpu_snn->dw[g_synapse_index + b] = 0.0; // reinitialize dw
            gpu_snn->init_w[synapse_index] += dw; // update initial w for incoming batches
            gpu_snn->w[g_synapse_index + b] = gpu_snn->init_w[synapse_index]; // update w
        }    
    }
}


__global__ void print_n_spikes(GPU_SNN_t *gpu_snn, GPU_results_t *gpu_results, size_t N, size_t batch_size){

    // get thread id: iN * batch_size * T
    size_t threadId = 
        (blockIdx.x  + gridDim.x  * blockIdx.y + gridDim.x * gridDim.y * blockIdx.z) *
        (blockDim.x * blockDim.y * blockDim.z) +
        (threadIdx.x + blockDim.x * threadIdx.y + blockDim.x * blockDim.y * threadIdx.z);


    if(threadId < N * batch_size){

        
        size_t i = threadId / batch_size; // neuron index
        size_t b = threadId % batch_size; // batch position

        printf(" ThreadId %llu (i = %llu, b = %llu): %d\n", threadId, i, b, gpu_results->n_spks[threadId]);
    }
}


extern "C" void simulate_LIF_batch_GPU(GPU_SNN_t *gpu_snn, GPU_dataset_t *gpu_dataset, GPU_results_t *gpu_results, tmp_batch_cpu_t *gpu_tmp_batch, simulation_configuration_t *conf, cuda_info_t *cuda_info, size_t bidx, GPU_SNN_t *cpu_snn, GPU_dataset_t *cpu_dataset, GPU_results_t *cpu_results){

    
    size_t N, iN, S, n_samples, n_features, n_networks, batch_size, T, LT;
    size_t t, s, snn;
    cudaError_t err;
    
    struct timespec start, end;
    struct timespec start_in, end_in, start_out, end_out, start_v, end_v;
    struct timespec start_learn, end_learn;
    struct timespec start_load, end_load, start_reinit, end_reinit;
    double et = 0.0, et_in = 0.0, et_out = 0.0, et_v = 0.0, et_learn = 0.0, et_load = 0.0, et_reinit = 0.0;


    // get general variables values
    N = cpu_snn->n_neurons;
    iN = cpu_snn->n_input_neurons;
    S = cpu_snn->n_synapses;
    T = conf->time_steps;
    LT = cpu_snn->LT;

    n_samples = cpu_dataset->n_samples;
    n_features = cpu_dataset->n_features;

    batch_size = cuda_info->batch_size;
    n_networks = cuda_info->n_networks;
    
    // set dim3 struct to launch kernels

    dim3 grid_neurons(cuda_info->n_blk_neurons_x[0], cuda_info->n_blk_neurons_y[0], cuda_info->n_blk_neurons_z[0]);
    dim3 block_neurons(cuda_info->n_thr_per_blk_neurons_x[0], cuda_info->n_thr_per_blk_neurons_y[0], cuda_info->n_thr_per_blk_neurons_z[0]);

    dim3 grid_synapses(cuda_info->n_blk_synapses_x[0], cuda_info->n_blk_synapses_y[0], cuda_info->n_blk_synapses_z[0]);
    dim3 block_synapses(cuda_info->n_threads_per_blk_synapses_x[0], cuda_info->n_threads_per_blk_synapses_y[0], cuda_info->n_threads_per_blk_synapses_z[0]);

    dim3 grid_in_neurons(cuda_info->n_blk_in_neurons_x[0], cuda_info->n_blk_in_neurons_y[0], cuda_info->n_blk_in_neurons_z[0]);
    dim3 block_in_neurons(cuda_info->n_thr_per_blk_in_neurons_x[0], cuda_info->n_thr_per_blk_in_neurons_y[0], cuda_info->n_thr_per_blk_in_neurons_z[0]);

    dim3 grid_all_neurons(cuda_info->n_blk_all_neurons_x[0], cuda_info->n_blk_all_neurons_y[0], cuda_info->n_blk_all_neurons_z[0]);
    dim3 block_all_neurons(cuda_info->n_thr_per_blk_all_neurons_x[0], cuda_info->n_thr_per_blk_all_neurons_y[0], cuda_info->n_thr_per_blk_all_neurons_z[0]);

    dim3 grid_uw(cuda_info->n_blk_uw_x[0], cuda_info->n_blk_uw_y[0], cuda_info->n_blk_uw_z[0]);
    dim3 block_uw(cuda_info->n_thr_per_blk_uw_x[0], cuda_info->n_thr_per_blk_uw_y[0], cuda_info->n_thr_per_blk_uw_z[0]);

    dim3 grid_is(cuda_info->n_blk_is_x[0], cuda_info->n_blk_is_y[0], cuda_info->n_blk_is_z[0]);
    dim3 block_is(cuda_info->n_thr_per_blk_is_x[0], cuda_info->n_thr_per_blk_is_y[0], cuda_info->n_thr_per_blk_is_z[0]);



    clock_gettime(CLOCK_MONOTONIC, &start);

    // (re)initialize neurons
    clock_gettime(CLOCK_MONOTONIC, &start_reinit);

    initialize_neurons_batch<<<grid_neurons, block_neurons>>>(gpu_snn, batch_size);
    err = cudaGetLastError();
    cudaCheckError(cudaPeekAtLastError());  // check launch errors
    cudaCheckError(cudaDeviceSynchronize());  // check runtime errors
    cudaDeviceSynchronize();

    // (re)initialize synapses
    initialize_synapses_batch<<<grid_synapses, block_synapses>>>(gpu_snn, batch_size);
    err = cudaGetLastError();
    cudaCheckError(cudaPeekAtLastError());  // check launch errors
    cudaCheckError(cudaDeviceSynchronize());  // check runtime errors
    cudaDeviceSynchronize();
    
    // load batch in tmp_batch
    initialize_batch_matrix_GPU<<<grid_in_neurons, block_in_neurons>>>(gpu_tmp_batch, gpu_dataset, bidx, batch_size, iN, batch_size, conf->max_input_spikes);
    err = cudaGetLastError();
    cudaCheckError(cudaPeekAtLastError());  // check launch errors
    cudaCheckError(cudaDeviceSynchronize());  // check runtime errors
    cudaDeviceSynchronize();

    clock_gettime(CLOCK_MONOTONIC, &end_reinit);
    et_reinit +=(end_reinit.tv_sec - start_reinit.tv_sec) + (end_reinit.tv_nsec - start_reinit.tv_nsec) / 1e9;

    // loop over networks?

    // loop over time steps
    for(t = 0; t<T; t++){

        size_t lt = t % cpu_snn->LT;

        // load sample in network
        clock_gettime(CLOCK_MONOTONIC, &start_load);

        load_batch_time_step_in_SNN_GPU<<<grid_all_neurons, block_all_neurons>>>(gpu_snn, gpu_tmp_batch, iN, N, batch_size, lt, t, bidx);
        err = cudaGetLastError();
        cudaCheckError(cudaPeekAtLastError());  // check launch errors
        cudaCheckError(cudaDeviceSynchronize());  // check runtime errors
        cudaDeviceSynchronize();
        
        clock_gettime(CLOCK_MONOTONIC, &end_load);
        et_load +=(end_load.tv_sec - start_load.tv_sec) + (end_load.tv_nsec - start_load.tv_nsec) / 1e9;

        // compute input step
        clock_gettime(CLOCK_MONOTONIC, &start_in);

        //process_input_currect_batch<<<grid_neurons, block_neurons>>>(gpu_snn, iN, N, batch_size, lt, t, (size_t)conf->learn);
        process_input_currect_batch<<<grid_is, block_is, cuda_info->n_thr_per_blk_is_x[0] * sizeof(float)>>>(gpu_snn, iN, N, batch_size, lt, t, (size_t)conf->learn);
        cudaCheckError(cudaPeekAtLastError());  // check launch errors
        cudaCheckError(cudaDeviceSynchronize());  // check runtime errors
        cudaDeviceSynchronize();
        
        clock_gettime(CLOCK_MONOTONIC, &end_in);
        et_in +=(end_in.tv_sec - start_in.tv_sec) + (end_in.tv_nsec - start_in.tv_nsec) / 1e9;

        // compute neuron dynamics
        clock_gettime(CLOCK_MONOTONIC, &start_v);

        process_V_batch<<<grid_neurons, block_neurons>>>(gpu_snn, N, batch_size);
        cudaCheckError(cudaPeekAtLastError());  // check launch errors
        cudaCheckError(cudaDeviceSynchronize());  // check runtime errors
        cudaDeviceSynchronize();

        clock_gettime(CLOCK_MONOTONIC, &end_v);
        et_v +=(end_v.tv_sec - start_v.tv_sec) + (end_v.tv_nsec - start_v.tv_nsec) / 1e9;


        // compute output step
        clock_gettime(CLOCK_MONOTONIC, &start_out);

        process_firing_batch<<<grid_neurons, block_neurons>>>(gpu_snn, gpu_results, iN, N, batch_size, lt, (size_t)conf->learn);
        cudaCheckError(cudaPeekAtLastError());  // check launch errors
        cudaCheckError(cudaDeviceSynchronize());  // check runtime errors
        cudaDeviceSynchronize();

        clock_gettime(CLOCK_MONOTONIC, &end_out);
        et_out +=(end_out.tv_sec - start_out.tv_sec) + (end_out.tv_nsec - start_out.tv_nsec) / 1e9;


        // compute learning step if training
        clock_gettime(CLOCK_MONOTONIC, &start_learn);

        if(conf->learn){

            compute_post_traces_batch<<<grid_neurons, block_neurons>>>(gpu_snn, N, batch_size);
            cudaCheckError(cudaPeekAtLastError());  // check launch errors
            cudaCheckError(cudaDeviceSynchronize());  // check runtime errors
            cudaDeviceSynchronize();

            trace_based_STDP_batch<<<grid_synapses, block_synapses>>>(gpu_snn, S, batch_size);
            cudaCheckError(cudaPeekAtLastError());  // check launch errors
            cudaCheckError(cudaDeviceSynchronize());  // check runtime errors
            cudaDeviceSynchronize();

            reinit_post_fired_batch<<<grid_neurons, block_neurons>>>(gpu_snn, N, batch_size);
            cudaCheckError(cudaPeekAtLastError());  // check launch errors
            cudaCheckError(cudaDeviceSynchronize());  // check runtime errors
            cudaDeviceSynchronize();

        }
        clock_gettime(CLOCK_MONOTONIC, &end_learn);
        et_learn +=(end_learn.tv_sec - start_learn.tv_sec) + (end_learn.tv_nsec - start_learn.tv_nsec) / 1e9;

    }

    // updaate weights if training
    clock_gettime(CLOCK_MONOTONIC, &start_learn);
    if(conf->learn){
            
        update_weights_batch<<<grid_uw, block_uw>>>(gpu_snn, S, batch_size);
        cudaCheckError(cudaPeekAtLastError());  // check launch errors
        cudaCheckError(cudaDeviceSynchronize());  // check runtime errors
        cudaDeviceSynchronize();
    }
    clock_gettime(CLOCK_MONOTONIC, &end_learn);
    et_learn +=(end_learn.tv_sec - start_learn.tv_sec) + (end_learn.tv_nsec - start_learn.tv_nsec) / 1e9;

    clock_gettime(CLOCK_MONOTONIC, &end);
    et +=(end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;


    // store in CPU results struct
    cpu_results->t = et;
    cpu_results->t_in = et_in;
    cpu_results->t_v = et_v;
    cpu_results->t_out = et_out;
    cpu_results->t_reinit = et_reinit;
    cpu_results->t_load = et_load;
    cpu_results->t_learn = et_learn;
}

extern "C" void simulate_batches_LIF_GPU(GPU_SNN_t **gpu_snn, GPU_dataset_t **gpu_dataset, GPU_results_t **gpu_results, tmp_batch_cpu_t **gpu_tmp_batch, simulation_configuration_t *conf, cuda_info_t *cuda_info, GPU_SNN_t *cpu_snn, GPU_dataset_t *cpu_dataset){
    
    // TODO: refactorize number of batches computation
    // compute number of batches
    size_t n_batches = cpu_dataset->n_samples / conf->batch_size;
    size_t r_samples = cpu_dataset->n_samples % conf->batch_size;
    // there are remaining sampels
    if(r_samples > 0){

        n_batches += 1; // last batch contains less samples
    }


    // copy constant(s) to all gpu devices
    for(size_t i = 0; i<cuda_info->nDevices; i++){
        
        cudaSetDevice(i);
        cudaError_t err = cudaMemcpyToSymbol(n_samples, &cpu_dataset->n_samples, sizeof(size_t));
    }

    
    // if it is executed in only one device
    if(cuda_info->nDevices == 1){
        
        cudaSetDevice(0);

        // initialize struct to store results of each batch
        GPU_results_t *cpu_results = initialize_batch_results_cpu(conf, cpu_snn->n_neurons, conf->batch_size, 1);
        GPU_results_t *cpu_batch_results = initialize_batch_results_cpu(conf, cpu_snn->n_neurons, conf->batch_size, 1);

        
        // simulate and print first batch
        simulate_LIF_batch_GPU(gpu_snn[0], gpu_dataset[0], gpu_results[0], gpu_tmp_batch[0], conf, cuda_info, 0, cpu_snn, cpu_dataset, cpu_batch_results);
        cpy_batch_results_GPU2CPU(cpu_results, gpu_results, cuda_info, cpu_snn->n_neurons, 1);
        
        // print results
        for(size_t batch = 0; batch<conf->batch_size; batch++){

            for(size_t i = 0; i<cpu_snn->n_neurons; i++){

                printf("%d ", cpu_results->n_spks[i * conf->batch_size + batch]);
            }
            printf("\n");
        }

        // loop over batches
        for(size_t b = 0; b<n_batches; b++){

            printf(" In batch %zu\n", b);

            // reinitialize results struct
            reinitialize_batch_results_cpu(cpu_batch_results, conf, cpu_snn->n_neurons, conf->batch_size, 1);

            // simulate batch
            simulate_LIF_batch_GPU(gpu_snn[0], gpu_dataset[0], gpu_results[0], gpu_tmp_batch[0], conf, cuda_info, b, cpu_snn, cpu_dataset, cpu_batch_results);

            // accumulate batch execution times in general results struct
            acc_batch_execution_times(cpu_results, cpu_batch_results);            
        }

        // print information in results struct
        printf("   > Total execution time: %.3lf | Mean per batch %.3lf | Mean per time step %.3lf\n", 
                    cpu_results->t, cpu_results->t / (float)n_batches, cpu_results->t / (float)n_batches / (float)conf->time_steps);
        printf("   > Total in time:        %.3lf | Mean per batch %.3lf | Mean per time step %.3lf\n", 
                    cpu_results->t_in, cpu_results->t_in / (float)n_batches, cpu_results->t_in / (float)n_batches / (float)conf->time_steps);
        printf("   > Total v time:         %.3lf | Mean per batch %.3lf | Mean per time step %.3lf\n", 
                    cpu_results->t_v, cpu_results->t_v / (float)n_batches, cpu_results->t_v / (float)n_batches / (float)conf->time_steps);
        printf("   > Total out time:       %.3lf | Mean per batch %.3lf | Mean per time step %.3lf\n", 
                    cpu_results->t_out, cpu_results->t_out / (float)n_batches, cpu_results->t_out / (float)n_batches / (float)conf->time_steps);
        printf("   > Total reinit time:    %.3lf | Mean per batch %.3lf | Mean per time step %.3lf\n", 
                    cpu_results->t_reinit, cpu_results->t_reinit / (float)n_batches, cpu_results->t_reinit / (float)n_batches / (float)conf->time_steps);
        printf("   > Total load time:      %.3lf | Mean per batch %.3lf | Mean per time step %.3lf\n", 
                    cpu_results->t_load, cpu_results->t_load / (float)n_batches, cpu_results->t_load / (float)n_batches / (float)conf->time_steps);
        printf("   > Total learn time:     %.3lf | Mean per batch %.3lf | Mean per time step %.3lf\n", 
                    cpu_results->t_learn, cpu_results->t_learn / (float)n_batches, cpu_results->t_learn / (float)n_batches / (float)conf->time_steps);

    }

    // TODO: multi GPU must be corrected
    else{

        int e, b, dev, batch_number = 0;
        int n_batches = 0;

        // allocate memory for managing dw updates
        GPU_SNN_t **tmp_snn = (GPU_SNN_t**)malloc(cuda_info->nDevices * sizeof(GPU_SNN_t*));
        float **dw = (float**)malloc(cuda_info->nDevices * sizeof(float*));
        for(int i = 0; i<cuda_info->nDevices; i++){
            tmp_snn[i] = (GPU_SNN_t*)malloc(sizeof(GPU_SNN_t));
            dw[i] = (float*)malloc(cpu_snn->n_synapses * sizeof(float));
        }


        // loop over epochs
        for(e=0; e<conf->epochs; e++){

            // compute the number of batches
            n_batches = cuda_info->n_samples / cuda_info->batch_size;
            if(cuda_info->n_samples % cuda_info->batch_size != 0)
                n_batches += 1;


            // loop over batches
            for(b=0; b<n_batches; b++){
     
                // start the simulation using several devices. Each device is managed by a openMP thread
                #pragma omp parallel num_threads(cuda_info->nDevices)
                {
                    //for(dev = 0; dev<cuda_info->nDevices; dev++){
                    
                    int dev = omp_get_thread_num();

                    //printf(" Dev = %d\n", dev);
                    fflush(stdout);

                    // set the device dev as the current one for the thread
                    cudaSetDevice(dev);

                    // simulate a part of the batch in the device
                    //simulate_LIF_in_multi_GPU(gpu_snn[dev], gpu_dataset[dev], gpu_results[dev], conf, cuda_info, cpu_snn, cpu_dataset, dev, b);

                    // wait until all GPU tasks finished
                    cudaDeviceSynchronize();

                    // if training
                    if(conf->learn == 1){
                        // sum dw values of different devices
                        // cpy weights to CPU
                        cudaMemcpy(tmp_snn[dev], gpu_snn[dev], sizeof(GPU_SNN_t), cudaMemcpyDeviceToHost); // cpy SNN structure GPU2CPU
                        cudaMemcpy(dw[dev], tmp_snn[dev]->dw, cpu_snn->n_synapses * sizeof(float), cudaMemcpyDeviceToHost); // cpy dw GPU2CPU

                        // wait until all threads obtained their weights
                        #pragma omp barrier
                        
                        // sum all dws
                        if(dev == 0){
                            for(int i = 0; i<cpu_snn->n_synapses; i++){

                                for(int j = 0; j<cuda_info->nDevices; j++){
                                        
                                    //printf(" > dev 0 (fixed) final dW[%d] = %f\n", i, dw[dev][i]);
                                    if(dev != j){
            
                                        //printf(" > dev %d final dW[%d] = %f\n", j, i, dw[j][i]);    
                                        dw[dev][i] += dw[j][i];
                                    }
                                }
                                
                                // divide by batch_size
                                dw[dev][i] = dw[dev][i] / (float)(conf->batch_size);
                                //printf(" dev %d: dW[%d] = %f\n", dev, i, dw[dev][i]);
                            }
                        }

                        #pragma omp barrier
                            
                        // cpy again to gpu
                        cudaMemcpy(tmp_snn[dev]->dw, dw[0], cpu_snn->n_synapses * sizeof(float), cudaMemcpyHostToDevice); // cpy dw CPU2GPU
                    }
                }
            }
        }
    }
}


extern "C" void simulate_LIF_in_multi_GPU(GPU_SNN_t *gpu_snn, GPU_dataset_t *gpu_dataset, GPU_results_t *gpu_results, simulation_configuration_t *conf, cuda_info_t *cuda_info, GPU_SNN_t *cpu_snn, GPU_dataset_t *cpu_dataset, int dev, int batch){
    
    // call kernels
    int time_steps, batch_size, dev_batch_size, n_networks;
    int t, s, b, snn;
    cudaError_t err;

    time_steps = cuda_info->time_steps;
    batch_size = cuda_info->batch_size;
    dev_batch_size = cuda_info->dev_batch_size[0];
    n_networks = cuda_info->n_networks_per_dev[0];

    // set dim3 struct to launch kernels
    /*dim3 grid_n_spkmatrix_reinit_blk((unsigned int)(cuda_info->n_threads_per_blk_rsm_x), cuda_info->n_threads_per_blk_rsm_y, cuda_info->n_threads_per_blk_rsm_z );
    dim3 block_n_spkmatrix_reinit_threads_per_blk((unsigned int)(cuda_info->n_threads_per_blk_rsm_x), cuda_info->n_threads_per_blk_rsm_y, cuda_info->n_threads_per_blk_rsm_z );    

    dim3 grid_n_ls_blk( (unsigned int)(cuda_info->n_blk_ls_x), cuda_info->n_blk_ls_y, cuda_info->n_blk_ls_z );
    dim3 block_n_ls_threads_per_block( (unsigned int)(cuda_info->n_threads_per_blk_ls_x), cuda_info->n_threads_per_blk_ls_y, cuda_info->n_threads_per_blk_ls_z );    

    dim3 grid_n_pn_blk( (unsigned int)(cuda_info->n_blk_nrs_x), cuda_info->n_blk_nrs_y, cuda_info->n_blk_nrs_z );
    dim3 block_n_pn_threads_per_block( (unsigned int)(cuda_info->n_threads_per_blk_nrs_x), cuda_info->n_threads_per_blk_nrs_y, cuda_info->n_threads_per_blk_nrs_z );  

    dim3 grid_n_synapses_blk((unsigned int)(cuda_info->n_blk_synapses_x), cuda_info->n_blk_synapses_y, cuda_info->n_blk_synapses_z);
    dim3 block_n_synapses_per_block((unsigned int)(cuda_info->n_threads_per_blk_synapses_x), cuda_info->n_threads_per_blk_synapses_y, cuda_info->n_threads_per_blk_synapses_z);

    dim3 grid_n_uw_blk((unsigned int)(cuda_info->n_blk_uw_x), cuda_info->n_blk_uw_y, cuda_info->n_blk_uw_z);
    dim3 block_n_uw_per_block((unsigned int)(cuda_info->n_threads_per_blk_uw_x), cuda_info->n_threads_per_blk_uw_y, cuda_info->n_threads_per_blk_uw_z);  

    
    // compute the sample index that will be simulated by this device
    int fsample = batch_size * batch + dev * dev_batch_size; 

    // update weights
    if(conf->learn == 1){
        krnl_update_weights_only<<<grid_n_uw_blk, block_n_uw_per_block>>>(gpu_snn, batch_size);
        err = cudaGetLastError();
        cudaCheckError(cudaPeekAtLastError());  // check launch errors
        cudaCheckError(cudaDeviceSynchronize());  // check runtime errors  


        // reinit synapses
        reinit_tot_synapses_kernel<<<grid_n_synapses_blk, block_n_synapses_per_block>>>(gpu_snn); // n_synapses 1440187
        err = cudaGetLastError();
        cudaCheckError(cudaPeekAtLastError());  // check launch errors
        cudaCheckError(cudaDeviceSynchronize());  // check runtime errors
    }


    // simulate a part of the batch using the available number of copies
    for(snn = 0; snn < dev_batch_size; snn = snn + n_networks){
        
        // reinit the entire matrix: number of elements = (n_neurons + n_input_neurons) * L
        // max 1024 threads per block: n_neurons * n_input_neurons threads, each L iterations
        // n_blocks = (n_neurons + n_input_neurons / 1024) + 1
        
        reinit_spk_matrix<<<grid_n_spkmatrix_reinit_blk, block_n_spkmatrix_reinit_threads_per_blk>>>(gpu_snn);
        err = cudaGetLastError();
        //printf("Launch error: %s\n", cudaGetErrorString(err));
        cudaCheckError(cudaPeekAtLastError());  // check launch errors
        cudaCheckError(cudaDeviceSynchronize());  // check runtime errors

        load_sample_kernel<<<grid_n_ls_blk, block_n_ls_threads_per_block>>>(gpu_snn, gpu_dataset, fsample + snn, fsample + dev_batch_size);
        err = cudaGetLastError();
        //printf("Launch error: %s\n", cudaGetErrorString(err));
        cudaCheckError(cudaPeekAtLastError());  // check launch errors
        cudaCheckError(cudaDeviceSynchronize());  // check runtime errors

        // reinit neurons
        reinit_neurons_kernel<<<grid_n_pn_blk, block_n_pn_threads_per_block>>>(gpu_snn); // n_neurons 8784
        err = cudaGetLastError();
        //printf("Launch error: %s\n", cudaGetErrorString(err));
        cudaCheckError(cudaPeekAtLastError());  // check launch errors
        cudaCheckError(cudaDeviceSynchronize());  // check runtime errors
        
        // reinit synapses
        reinit_synapses_kernel<<<grid_n_synapses_blk, block_n_synapses_per_block>>>(gpu_snn); // n_synapses 1440187
        err = cudaGetLastError();
        cudaCheckError(cudaPeekAtLastError());  // check launch errors
        cudaCheckError(cudaDeviceSynchronize());  // check runtime errors
        
        for(t = 0; t<time_steps; t++){

            // input step
            lif_neuron_input<<<grid_n_pn_blk, block_n_pn_threads_per_block>>>(gpu_snn, fsample + snn, t, fsample + dev_batch_size);
            err = cudaGetLastError();
            //printf("Launch error: %s\n", cudaGetErrorString(err));
            cudaCheckError(cudaPeekAtLastError());  // check launch errors
            cudaCheckError(cudaDeviceSynchronize());  // check runtime errors        
            //printf(" > Input step processed\n");
            //fflush(stdout);

            // output step
            lif_neuron_output<<<grid_n_pn_blk, block_n_pn_threads_per_block>>>(gpu_snn, fsample + snn, t, gpu_results, fsample + dev_batch_size);
            err = cudaGetLastError();
            //printf("Launch error: %s\n", cudaGetErrorString(err));
            cudaCheckError(cudaPeekAtLastError());  // check launch errors
            cudaCheckError(cudaDeviceSynchronize());  // check runtime errors          
            //printf(" > Output step processed\n");
            //fflush(stdout);
    
            // learning
            if(conf->learn == 1){
                stdp_kernel<<<grid_n_synapses_blk, block_n_synapses_per_block>>>(gpu_snn, t);
                err = cudaGetLastError();
                cudaCheckError(cudaPeekAtLastError());  // check launch errors
                cudaCheckError(cudaDeviceSynchronize());  // check runtime errors  
                //printf(" > STDP processed\n");
                //fflush(stdout);
            }
        }
    }

    // sum all copies dw in the first cpy
    if(conf->learn == 1){
        krnl_sum_dw<<<grid_n_uw_blk, block_n_uw_per_block>>>(gpu_snn);
        err = cudaGetLastError();
        cudaCheckError(cudaPeekAtLastError());  // check launch errors
        cudaCheckError(cudaDeviceSynchronize());  // check runtime errors  
    }*/
}