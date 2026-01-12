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


        //printf(" Synapse = %llu, thread = %llu\n", synapse_index, threadId);
        snn->w[threadId] = snn->init_w[synapse_index];

        //printf("   > ThreadId = %llu (s = %llu, b = %llu), w = %f, init_w = %f\n", 
        //    threadId, threadId / batch_size, threadId % batch_size, snn->w[threadId], snn->init_w[synapse_index]);

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
    
            //printf(" > Thread = %llu (i = %llu, b = %llu), loading sample %llu, number of spikes: %llu\n", threadId, fidx, b, sidx, n_spikes_feature);

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

            //if(tmp_batch->spikes[batch_mtr_idx]){

            //    printf("   > Thread = %llu (i = %llu, b = %llu), loading spike from %llu to %llu\n", threadId, i, b, batch_mtr_idx, spk_mtr_idx);
            //}
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


    if(threadId < N * batch_size){

        size_t i = threadId / batch_size; // neuron index
        size_t b = threadId % batch_size; // batch position

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
                    //spk = 0.0;
                    //printf("   > Thread = %llu (i = %llu, b = %llu), accessing to %llu, where %f\n", threadId, i, b, index, spk);

                    // update I if r period remain <= 0 and spike received
                    if(gpu_snn->r_period_remain[g_neuron_index + b] <= 0){
                        
                        //if(gpu_snn->w[g_synapse_index + b] > 0.0)
                            //printf("      > W = %f\n", gpu_snn->w[g_synapse_index + b]);
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

        //if(I > 0.0)
        //    printf(" > Neuron %llu (batch %llu) I = %f\n", neuron_index, b, I);
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
    
            //printf(" ThreadId %llu (i = %llu, b = %llu): %llu\n", threadId, i, b, idx + b);

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
                    //spk = 0.0;
                    //printf("   > Thread = %llu (i = %llu, b = %llu), accessing to %llu, where %f\n", threadId, i, b, index, spk);

                    // update I if r period remain <= 0 and spike received
                    if(gpu_snn->r_period_remain[g_neuron_index + b] <= 0){
                        
                        //if(gpu_snn->w[g_synapse_index + b] > 0.0)
                            //printf("      > W = %f\n", gpu_snn->w[g_synapse_index + b]);
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


__global__ void trace_based_STDP_batch(GPU_SNN_t *gpu_snn, size_t N, size_t batch_size){

    
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

    
    size_t N, iN, S, n_samples, n_features, n_networks, batch_size, dev_batch_size, T, LT;
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
    dev_batch_size = cuda_info->dev_batch_size;
    n_networks = cuda_info->n_networks;
    
    // set dim3 struct to launch kernels

    dim3 grid_neurons(cuda_info->n_blk_neurons_x, cuda_info->n_blk_neurons_y, cuda_info->n_blk_neurons_z);
    dim3 block_neurons(cuda_info->n_thr_per_blk_neurons_x, cuda_info->n_thr_per_blk_neurons_y, cuda_info->n_thr_per_blk_neurons_z);

    dim3 grid_synapses(cuda_info->n_blk_synapses_x, cuda_info->n_blk_synapses_y, cuda_info->n_blk_synapses_z);
    dim3 block_synapses(cuda_info->n_threads_per_blk_synapses_x, cuda_info->n_threads_per_blk_synapses_y, cuda_info->n_threads_per_blk_synapses_z);

    dim3 grid_in_neurons(cuda_info->n_blk_in_neurons_x, cuda_info->n_blk_in_neurons_y, cuda_info->n_blk_in_neurons_z);
    dim3 block_in_neurons(cuda_info->n_thr_per_blk_in_neurons_x, cuda_info->n_thr_per_blk_in_neurons_y, cuda_info->n_thr_per_blk_in_neurons_z);

    dim3 grid_all_neurons(cuda_info->n_blk_all_neurons_x, cuda_info->n_blk_all_neurons_y, cuda_info->n_blk_all_neurons_z);
    dim3 block_all_neurons(cuda_info->n_thr_per_blk_all_neurons_x, cuda_info->n_thr_per_blk_all_neurons_y, cuda_info->n_thr_per_blk_all_neurons_z);

    dim3 grid_uw(cuda_info->n_blk_uw_x, cuda_info->n_blk_uw_y, cuda_info->n_blk_uw_z);
    dim3 block_uw(cuda_info->n_thr_per_blk_uw_x, cuda_info->n_thr_per_blk_uw_y, cuda_info->n_thr_per_blk_uw_z);


    clock_gettime(CLOCK_MONOTONIC, &start);


    // loop over network copies: if we can not store batch_size neuron copies, each neuron is used several times to compute the entire batch
    //for(snn = 0; snn < batch_size; snn += n_networks){

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
        
        /*total_LIF_batch<<<grid_neurons, block_neurons>>>(gpu_snn, gpu_results, iN, N, batch_size, lt, t, (size_t)conf->learn);
        cudaCheckError(cudaPeekAtLastError());  // check launch errors
        cudaCheckError(cudaDeviceSynchronize());  // check runtime errors
        cudaDeviceSynchronize();*/
        clock_gettime(CLOCK_MONOTONIC, &start_in);

        process_input_currect_batch<<<grid_neurons, block_neurons>>>(gpu_snn, iN, N, batch_size, lt, t, (size_t)conf->learn);
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

            trace_based_STDP_batch<<<grid_neurons, block_neurons>>>(gpu_snn, N, batch_size);
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

    //print_n_spikes<<<grid_neurons, block_neurons>>>(gpu_snn, gpu_results, N, batch_size);
    //cudaCheckError(cudaDeviceSynchronize());  // check runtime errors

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
        //simulate_LIF_batch_GPU(gpu_snn[0], gpu_dataset[0], gpu_results[0], gpu_tmp_batch[0], conf, cuda_info, 0, cpu_snn, cpu_dataset, cpu_batch_results);
        //cpy_batch_results_GPU2CPU(cpu_results, gpu_results, cuda_info, cpu_snn->n_neurons, 1);
        // print results
        /*for(size_t batch = 0; batch<conf->batch_size; batch++){

            for(size_t i = 0; i<cpu_snn->n_neurons; i++){

                printf("%d ", cpu_results->n_spks[i * conf->batch_size + batch]);
            }
            printf("\n");
        }*/

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
    // if it is multi GPU
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
                    simulate_LIF_in_multi_GPU(gpu_snn[dev], gpu_dataset[dev], gpu_results[dev], conf, cuda_info, cpu_snn, cpu_dataset, dev, b);

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



//////////////////////
//////// OLD

__global__ void inititialize_spk_matrix(GPU_SNN_t *snn, size_t batch_size, size_t iN){

    // get thread id
    size_t threadId = 
        (blockIdx.x  + gridDim.x  * blockIdx.y + gridDim.x * gridDim.y * blockIdx.z) *
        (blockDim.x * blockDim.y * blockDim.z) +
        (threadIdx.x + blockDim.x * threadIdx.y + blockDim.x * blockDim.y * threadIdx.z);

    if(threadId < snn->n_synapses * batch_size){

        for(size_t t = 0; t<snn->LT; t++){
            
            //snn->spk_matrix[threadId + batch_size * t] = 0;
        }
    }    
}




__global__ void reinit_spk_matrix(GPU_SNN_t *snn){
    
    // get thread id
    size_t threadId = 
        (blockIdx.x  + gridDim.x  * blockIdx.y + gridDim.x * gridDim.y * blockIdx.z) *
        (blockDim.x * blockDim.y * blockDim.z) +
        (threadIdx.x + blockDim.x * threadIdx.y + blockDim.x * blockDim.y * threadIdx.z);


    //L = snn->max_spikes;

    // reinit in all network copies
    if(threadId < (size_t)((size_t)(snn->n_neurons + snn->n_input_neurons) * (size_t)(snn->max_spikes) * (size_t)(snn->n_networks))){
        snn->spk_matrix[threadId] = -1;
    }
}


__global__ void load_sample_kernel(GPU_SNN_t *snn, GPU_dataset_t *dataset, int b, int batch_last_sample)
{
    size_t L;

    size_t threadId = 
        (blockIdx.x  + gridDim.x  * blockIdx.y + gridDim.x * gridDim.y * blockIdx.z) *
        (blockDim.x * blockDim.y * blockDim.z) +
        (threadIdx.x + blockDim.x * threadIdx.y + blockDim.x * blockDim.y * threadIdx.z);


    L = (size_t)(snn->max_spikes);
    
    //printf(" > Neuron index = %d, input neurons = %d, sample index = %d, n_features = %d\n", 
            //neuron_index, snn->n_input_neurons, s, dataset->n_features);

    // [n_input_neurons * n_networks] threads
    size_t n_threads = (size_t)(snn->n_input_neurons) * (size_t)(snn->n_networks);

    if(threadId < n_threads){
        
        // s is the sample index and f the feature index
        // start_sample is the index in which the sample starts in the array of spikes
        // start_feature is the index in which the feature starts locally in the array of spikes
        
        // get the copy number // since there are a lot of copies, neuron_index can be bigger than max_int
        size_t network_number = (size_t)(threadId / (size_t)(snn->n_input_neurons));
            
        // local index of the neuron in the cpy
        size_t local_neuron_index = (size_t)(threadId % (size_t)(snn->n_input_neurons));

        // neuron index in the spike matrix = [(n_neurons + n_input_neurons) * n_network]
        size_t global_neuron_index = network_number * ((size_t)(snn->n_input_neurons) + (size_t)(snn->n_neurons)) + local_neuron_index;

        // get sample in batch = s + network_number
        int s = b + (int)(network_number); // the i.th network computes the s + i sample


        if(s < batch_last_sample && s < n_samples){

            // get the feature
            int f = local_neuron_index; // the neuron index is also the feature index

            //printf(" > > Feature = %d\n", f);

            size_t start_sample = dataset->sample_offset[s];
            //printf(" Thread %d, start sample = %d\n", f, start_sample);
            
            size_t start_feature = dataset->feature_offset[(size_t)(s * dataset->n_features) + size_t(f)];
            //printf(" Thread %d, start feature = %d\n", f, start_feature);

            size_t n_spikes_per_feature = size_t(dataset->n_spikes_per_feature[(size_t)(s * dataset->n_features) + size_t(f)]);

            //printf(" > Loading spike train in neuron = %d (n = %d)\n", neuron_index, n_spikes);

            for(size_t i = 0; i<n_spikes_per_feature; i++){
            
                // store a 1 value if there is a spike
                //snn->spk_matrix[global_neuron_index * L + (size_t)(dataset->spikes[start_sample + start_feature + i])] = dataset->spikes[start_sample + start_feature + i];
                snn->spk_matrix[global_neuron_index * L + (size_t)(dataset->spikes[start_feature + i])] = dataset->spikes[start_feature + i];

                //printf(" > Thread %d: Neuron %d (local %d, net %d). sample %d, feature %d, spktime %d\n", neuron_index, global_neuron_index, local_neuron_index, network_number, s, f, dataset->spikes[start_sample + start_feature + i]);
                //printf(" SPike added by thread %d to %d (value %d).\n", f, neuron_index * L + dataset->spikes[start_sample + start_feature + i], dataset->spikes[start_sample + start_feature + i]);
            }
        }
    }
}


/* OLD
__global__ void load_sample_kernel(GPU_SNN_t *snn, GPU_dataset_t *dataset, int s, int batch_size)
{
    //lortu hariaren identifikadorea

    int i, L;
    int threadsPerBlock = blockDim.x * blockDim.y;
    int threadNumInBlock = threadIdx.x + blockDim.x * threadIdx.y; //(alternatively: threadIdx.y + blockDim.y * threadIdx.x);
    int blockNumInGrid   = blockIdx.x  + gridDim.x  * blockIdx.y;  //(alternatively: blockIdx.y  + gridDim.y  * blockIdx.x);

    int globalThreadNum = blockNumInGrid * threadsPerBlock + threadNumInBlock;
    int neuron_index = globalThreadNum; // ?

    L = snn->max_spikes;

    //printf(" > Neuron index = %d, input neurons = %d, sample index = %d, n_features = %d\n", 
            //neuron_index, snn->n_input_neurons, s, dataset->n_features);

    // if it is an input neuron
    if(neuron_index < snn->n_input_neurons * snn->n_networks){
        
        // s is the sample index and f the feature index
        // start_sample is the index in which the sample starts in the array of spikes
        // start_feature is the index in which the feature starts locally in the array of spikes


        int f = neuron_index; // the neuron index is also the feature index

        //printf(" > > Feature = %d\n", f);

        int start_sample = dataset->sample_offset[s];
        //printf(" Thread %d, start sample = %d\n", f, start_sample);
        
        int start_feature = dataset->feature_offset[s * dataset->n_features + f];
        //printf(" Thread %d, start feature = %d\n", f, start_feature);

        int n_spikes_per_feature = dataset->n_spikes_per_feature[s * dataset->n_features + f];

        //printf(" > Loading spike train in neuron = %d (n = %d)\n", neuron_index, n_spikes);

        for(i = 0; i<n_spikes_per_feature; i++){
        
            // store a 1 value if there is a spike
            snn->spk_matrix[neuron_index * L + dataset->spikes[start_sample + start_feature + i]] = dataset->spikes[start_sample + start_feature + i];

            //printf(" SPike added by thread %d to %d (value %d).\n", f, neuron_index * L + dataset->spikes[start_sample + start_feature + i], dataset->spikes[start_sample + start_feature + i]);
        }
    }
}
*/

__global__ void reinit_neurons_kernel(GPU_SNN_t *snn)
{
    //lortu hariaren identifikadorea

    //size_t threadsPerBlock = blockDim.x * blockDim.y * blockDim.z;
    //size_t threadNumInBlock = threadIdx.x + blockDim.x * threadIdx.y + blockDim.x * blockDim.y * threadIdx.z; //(alternatively: threadIdx.y + blockDim.y * threadIdx.x);
    //size_t blockNumInGrid   = blockIdx.x  + gridDim.x  * blockIdx.y + gridDim.x * gridDim.y * blockIdx.z;  //(alternatively: blockIdx.y  + gridDim.y  * blockIdx.x);

    //size_t globalThreadNum = blockNumInGrid * threadsPerBlock + threadNumInBlock;

    size_t threadId = 
        (blockIdx.x  + gridDim.x  * blockIdx.y + gridDim.x * gridDim.y * blockIdx.z) *
        (blockDim.x * blockDim.y * blockDim.z) +
        (threadIdx.x + blockDim.x * threadIdx.y + blockDim.x * blockDim.y * threadIdx.z);


    if(threadId < ((size_t)(snn->n_neurons) * (size_t)(snn->n_networks))){


        // get the copy number // since there are a lot of copies, neuron_index can be bigger than max_int
        size_t network_number = (size_t)(threadId / (size_t)(snn->n_neurons));    
        size_t local_neuron_index = (size_t)(threadId % (size_t)(snn->n_neurons));
        size_t global_neuron_index = network_number * (size_t)(snn->n_neurons) + local_neuron_index;
        
        snn->v[global_neuron_index] = snn->v_rest[local_neuron_index];
        snn->r_period_remain[global_neuron_index] = 0;

        // reinit last spikes array
        //for(int i = 0; i<snn->n_last_spikes; i++)
        //    snn->t_last_spikes[global_neuron_index * snn->n_last_spikes + i] = 0;
        //snn->next_last_spike[global_neuron_index] = 0;
    }

}


__global__ void reinit_synapses_kernel(GPU_SNN_t *snn)
{
    //size_t threadsPerBlock = blockDim.x * blockDim.y * blockDim.z;
    //size_t threadNumInBlock = threadIdx.x + blockDim.x * threadIdx.y + blockDim.x * blockDim.y * threadIdx.z; //(alternatively: threadIdx.y + blockDim.y * threadIdx.x);
    //size_t blockNumInGrid   = blockIdx.x  + gridDim.x  * blockIdx.y + gridDim.x * gridDim.y * blockIdx.z;  //(alternatively: blockIdx.y  + gridDim.y  * blockIdx.x);

    //size_t globalThreadNum = blockNumInGrid * threadsPerBlock + threadNumInBlock;

    // get synapse
    size_t threadId = 
        (blockIdx.x  + gridDim.x  * blockIdx.y + gridDim.x * gridDim.y * blockIdx.z) *
        (blockDim.x * blockDim.y * blockDim.z) +
        (threadIdx.x + blockDim.x * threadIdx.y + blockDim.x * blockDim.y * threadIdx.z);

    if(threadId < (size_t)(snn->n_synapses)){

        // summ weights of all network copies
        for(size_t i = 0; i<(size_t)(snn->n_networks); i++){

            // reinitialize synapse weight
            size_t sindex  = (size_t)(snn->n_synapses) * i + threadId;

            snn->w[sindex] = snn->init_w[sindex];
            //snn->next_pre_spike[sindex] = 0;
            //snn->next_post_spike[sindex] = 0;
        }
    }
}


__global__ void reinit_tot_synapses_kernel(GPU_SNN_t *snn)
{
    //size_t threadsPerBlock = blockDim.x * blockDim.y * blockDim.z;
    //size_t threadNumInBlock = threadIdx.x + blockDim.x * threadIdx.y + blockDim.x * blockDim.y * threadIdx.z; //(alternatively: threadIdx.y + blockDim.y * threadIdx.x);
    //size_t blockNumInGrid   = blockIdx.x  + gridDim.x  * blockIdx.y + gridDim.x * gridDim.y * blockIdx.z;  //(alternatively: blockIdx.y  + gridDim.y  * blockIdx.x);

    //size_t globalThreadNum = blockNumInGrid * threadsPerBlock + threadNumInBlock;

    // get synapse
    size_t threadId = 
        (blockIdx.x  + gridDim.x  * blockIdx.y + gridDim.x * gridDim.y * blockIdx.z) *
        (blockDim.x * blockDim.y * blockDim.z) +
        (threadIdx.x + blockDim.x * threadIdx.y + blockDim.x * blockDim.y * threadIdx.z);

    if(threadId < (size_t)(snn->n_synapses)){


        // summ weights of all network copies
        for(size_t i = 0; i<(size_t)(snn->n_networks); i++){

            // reinitialize synapse weight
            size_t sindex  = (size_t)(snn->n_synapses) * i + threadId;

            snn->w[sindex] = snn->init_w[sindex];
            //snn->next_pre_spike[sindex] = 0;
            //snn->next_post_spike[sindex] = 0;
            snn->dw[sindex] = 0.0;
        }
    }
}

__global__ void lif_neuron_input(GPU_SNN_t *snn, int b, int t, int batch_last_sample)
{
    //lortu hariaren identifikadorea

    size_t L;
    
    //size_t threadsPerBlock = blockDim.x * blockDim.y * blockDim.z;
    //size_t threadNumInBlock = threadIdx.x + blockDim.x * threadIdx.y + blockDim.x * blockDim.y * threadIdx.z; //(alternatively: threadIdx.y + blockDim.y * threadIdx.x);
    //size_t blockNumInGrid   = blockIdx.x  + gridDim.x  * blockIdx.y + gridDim.x * gridDim.y * blockIdx.z;  //(alternatively: blockIdx.y  + gridDim.y  * blockIdx.x);

    //size_t globalThreadNum = blockNumInGrid * threadsPerBlock + threadNumInBlock;

    size_t threadId = 
        (blockIdx.x  + gridDim.x  * blockIdx.y + gridDim.x * gridDim.y * blockIdx.z) *
        (blockDim.x * blockDim.y * blockDim.z) +
        (threadIdx.x + blockDim.x * threadIdx.y + blockDim.x * blockDim.y * threadIdx.z);
    
    L = (size_t)(snn->max_spikes);

    if(threadId < ((size_t)(snn->n_neurons) * (size_t)(snn->n_networks))){

        float I = 0.0; // in thread local memory?

        size_t network_number = (size_t)(threadId / (size_t)(snn->n_neurons)); // network cpy number
        size_t local_n = (size_t)(threadId % (size_t)snn->n_neurons); // local neuron index
        size_t global_n = network_number * (size_t)(snn->n_neurons) + local_n; // global neuron index (in copies)

        size_t s = (size_t)(b) + network_number;

        if(s < (size_t)(batch_last_sample) && (int)(s) < n_samples){

            // loop over input synapses to check if the index should be updated in one of them
            if(snn->r_period_remain[global_n] <= 0){
                
                size_t synapse_index, in_neuron_index, spk_time_index, i;
                int latency, spk_time, spk;

                for(i = 0; i<(size_t)(snn->n_neuron_input_synapses[local_n]); i++){ // no copies

                    synapse_index = (size_t)(snn->neuron_input_synapses_offset[local_n]) + i; // no copies
                    in_neuron_index = (size_t)(snn->pre_neuron_index[synapse_index]); // no copies
                    latency = snn->delay[synapse_index]; // no copies
                    spk_time = t - latency; // actual position

                    // scale in_neuron_index in the global spike matrix
                    size_t global_in_neuron_index = 
                        network_number * (size_t)(snn->n_neurons + snn->n_input_neurons) + in_neuron_index;

                    if(spk_time >=0 ){

                        spk_time_index = (size_t)(spk_time) % L; // scale to circular array
                        spk = snn->spk_matrix[global_in_neuron_index * L + spk_time_index];
                        
                        if(spk == spk_time){

                            I += snn->w[network_number * (size_t)(snn->n_synapses) + synapse_index];
                        }
                        //if(s > b)
                        //    printf(" Neuron %d: input synapse %d (l=%d, in neuron=%d), spk_time=%d, spk_time_index=%d, spk=%d\n", global_n, synapse_index, latency, in_neuron_index, spk_time, spk_time_index, spk);
                    }    
            
                }
                
                snn->v[global_n] = (1-0.05) * snn->v[global_n] + snn->v_rest[local_n] * 0.05 + I;
                //if(s > b)
                //    printf(" > Neuron %d, V = %f, I = %f\n", global_n, snn->v[global_n], I);
            }
        }
    }
}

__global__ void lif_neuron_input_batch(GPU_SNN_t *snn, int b, int t, int batch_last_sample)
{
    //lortu hariaren identifikadorea

    size_t L;
    
    //size_t threadsPerBlock = blockDim.x * blockDim.y * blockDim.z;
    //size_t threadNumInBlock = threadIdx.x + blockDim.x * threadIdx.y + blockDim.x * blockDim.y * threadIdx.z; //(alternatively: threadIdx.y + blockDim.y * threadIdx.x);
    //size_t blockNumInGrid   = blockIdx.x  + gridDim.x  * blockIdx.y + gridDim.x * gridDim.y * blockIdx.z;  //(alternatively: blockIdx.y  + gridDim.y  * blockIdx.x);

    //size_t globalThreadNum = blockNumInGrid * threadsPerBlock + threadNumInBlock;

    size_t threadId = 
        (blockIdx.x  + gridDim.x  * blockIdx.y + gridDim.x * gridDim.y * blockIdx.z) *
        (blockDim.x * blockDim.y * blockDim.z) +
        (threadIdx.x + blockDim.x * threadIdx.y + blockDim.x * blockDim.y * threadIdx.z);
    
    L = (size_t)(snn->max_spikes);

    // threadId < N * N_networks (pass this parameters as input parameters)
    if(threadId < snn->n_neurons * snn->n_networks){

        float I = 0.0; // in thread local memory?

        size_t network_number = (size_t)(threadId / (size_t)(snn->n_neurons)); // network cpy number
        size_t local_n = (size_t)(threadId % (size_t)snn->n_neurons); // local neuron index
        size_t global_n = network_number * (size_t)(snn->n_neurons) + local_n; // global neuron index (in copies)

        size_t s = (size_t)(b) + network_number;

        if(s < (size_t)(batch_last_sample) && (int)(s) < n_samples){

            // loop over input synapses to check if the index should be updated in one of them
            if(snn->r_period_remain[global_n] <= 0){
                
                size_t synapse_index, in_neuron_index, spk_time_index, i;
                int latency, spk_time, spk;

                for(i = 0; i<(size_t)(snn->n_neuron_input_synapses[local_n]); i++){ // no copies

                    synapse_index = (size_t)(snn->neuron_input_synapses_offset[local_n]) + i; // no copies
                    in_neuron_index = (size_t)(snn->pre_neuron_index[synapse_index]); // no copies
                    latency = snn->delay[synapse_index]; // no copies
                    spk_time = t - latency; // actual position

                    // scale in_neuron_index in the global spike matrix
                    size_t global_in_neuron_index = 
                        network_number * (size_t)(snn->n_neurons + snn->n_input_neurons) + in_neuron_index;

                    if(spk_time >=0 ){

                        spk_time_index = (size_t)(spk_time) % L; // scale to circular array
                        spk = snn->spk_matrix[global_in_neuron_index * L + spk_time_index];
                        
                        if(spk == spk_time){

                            I += snn->w[network_number * (size_t)(snn->n_synapses) + synapse_index];
                        }
                        //if(s > b)
                        //    printf(" Neuron %d: input synapse %d (l=%d, in neuron=%d), spk_time=%d, spk_time_index=%d, spk=%d\n", global_n, synapse_index, latency, in_neuron_index, spk_time, spk_time_index, spk);
                    }    
            
                }
                
                snn->v[global_n] = (1-0.05) * snn->v[global_n] + snn->v_rest[local_n] * 0.05 + I;
                //if(s > b)
                //    printf(" > Neuron %d, V = %f, I = %f\n", global_n, snn->v[global_n], I);
            }
        }
    }
}


__global__ void lif_neuron_output(GPU_SNN_t *snn, int b, int t, GPU_results_t *results, int batch_last_sample)
{
    //lortu hariaren identifikadorea

    //size_t threadsPerBlock = blockDim.x * blockDim.y * blockDim.z;
    //size_t threadNumInBlock = threadIdx.x + blockDim.x * threadIdx.y + blockDim.x * blockDim.y * threadIdx.z; //(alternatively: threadIdx.y + blockDim.y * threadIdx.x);
    //size_t blockNumInGrid   = blockIdx.x  + gridDim.x  * blockIdx.y + gridDim.x * gridDim.y * blockIdx.z;  //(alternatively: blockIdx.y  + gridDim.y  * blockIdx.x);

    //size_t globalThreadNum = blockNumInGrid * threadsPerBlock + threadNumInBlock;

    size_t threadId = 
        (blockIdx.x  + gridDim.x  * blockIdx.y + gridDim.x * gridDim.y * blockIdx.z) *
        (blockDim.x * blockDim.y * blockDim.z) +
        (threadIdx.x + blockDim.x * threadIdx.y + blockDim.x * blockDim.y * threadIdx.z);
    
    
    size_t L = (size_t)snn->max_spikes;

    if(threadId < (size_t)((size_t)snn->n_neurons * (size_t)snn->n_networks)){

        size_t network_number = (size_t)(threadId / (size_t)(snn->n_neurons));
        size_t local_n = (size_t)(threadId % (size_t)snn->n_neurons);
        size_t global_n = network_number * (size_t)(snn->n_neurons) + local_n;
        size_t ind_last_spike;
        
        snn->r_period_remain[global_n] --;

        int s = b + (int)network_number;

        if(s < batch_last_sample && s < n_samples){

            // v has copies for each network, v_thresh only one original copy
            if(snn->v[global_n] >= snn->v_thresh[local_n]){

                size_t row, col, matrix_index, global_index;

                matrix_index = network_number * (size_t)(snn->n_neurons + snn->n_input_neurons) * L; // cpy starts at
                row = (local_n + (size_t)(snn->n_input_neurons)) * L; // row start at
                col = (size_t)(t) % L; // column to locate

                global_index = matrix_index + row + col;
                
                // add spike to matrix
                snn->spk_matrix[global_index] = t;

                // reinit neuron values
                snn->r_period_remain[global_n] = snn->r_period[local_n];
                
                // add last spike to the neuron
                //ind_last_spike = (size_t)(global_n * (size_t)(snn->n_last_spikes) + (size_t)snn->next_last_spike[global_n]);
                //snn->t_last_spikes[ind_last_spike] = t;
                //snn->next_last_spike[global_n] = (snn->next_last_spike[global_n] + 1) % snn->n_last_spikes;
                
                snn->v[global_n] = snn->v_rest[local_n]; // reinit v_rest

                // 
                size_t rs_index = (size_t)s * (size_t)snn->n_neurons + local_n;
                //results->nspk[rs_index] += 1;
                //printf(" > Neuron %d, t %d, SPIKE\n", n, t);
            }
        }
    }
}


__device__ int mod_gpu(int a, int m){
    return ( (a % m) + m ) % m;
}

/* STDP function */
__global__ void stdp_kernel(GPU_SNN_t *snn, int t)
{
    size_t threadId = 
        (blockIdx.x  + gridDim.x  * blockIdx.y + gridDim.x * gridDim.y * blockIdx.z) *
        (blockDim.x * blockDim.y * blockDim.z) +
        (threadIdx.x + blockDim.x * threadIdx.y + blockDim.x * blockDim.y * threadIdx.z);
    
    
    if(threadId < (size_t)((size_t)(snn->n_synapses) * (size_t)(snn->n_networks))){

        // get identifiers
        size_t network_number = (size_t)(threadId / (size_t)(snn->n_synapses));
        size_t local_n = (size_t)(threadId % (size_t)(snn->n_synapses)); // local synapse index
        size_t global_n = network_number * (size_t)(snn->n_synapses) + local_n; // global synapse index
        
        // temporal dw
        float dw = 0.0;
        size_t pre_neuron_local_index, post_neuron_local_index; 

        // get pre and post neurons local and global indexes
        pre_neuron_local_index = (size_t)(snn->pre_neuron_index[local_n]); // no copies
        post_neuron_local_index = (size_t)(snn->post_neuron_index[local_n]); // no copies


        
        int delay = snn->delay[local_n];
        size_t D = (size_t)(delay);
        size_t T = (size_t)(t);
        size_t L = (size_t)(snn->max_spikes);


        // indexes in spk matrix
        size_t matrix_index = network_number * (size_t)(snn->n_neurons + snn->n_input_neurons) * L; // cpy starts at
        size_t pre_row = pre_neuron_local_index * L; // row start at
        size_t post_row = post_neuron_local_index  * L; // row start at        
        size_t pre_neuron_spk_matrix_index, post_neuron_spk_matrix_index; // index of where the first neuron starts
        pre_neuron_spk_matrix_index = matrix_index + pre_row;
        post_neuron_spk_matrix_index = matrix_index + post_row;


        //printf("Synapse %lu, pre neuron = %lu (row = %lu), post neuron = %lu (row = %lu)\n", 
        //        threadId, pre_neuron_local_index, pre_neuron_spk_matrix_index, post_neuron_local_index, post_neuron_spk_matrix_index);

        //printf(" Pre neuron = %d, matrix index = %lu, next pre spike %d\n", pre_neuron_local_index, pre_neuron_spk_matrix_index, snn->next_pre_spike[threadId]);
        //printf(" Post neuron = %d, matrix index = %lu, next post spike %d\n", post_neuron_local_index, post_neuron_spk_matrix_index, snn->next_post_spike[threadId]);
        // initial condition: 
        // -- not input synapse (pre synapses do not spike) --> TODO
        // -- at least presynaptic spike reaches in t (LTD), or post spike occurred in t (LTP) 


        int tpost, tpre, tdiff, tmp = 0;
        float ftdiff = 0.0;
        
        if(local_n < snn->n_synapses - 1){
            // LTP: check if post neuron fired at t
            if(T > 0 && snn->spk_matrix[post_neuron_spk_matrix_index + T] == t && T >= D){

                //printf(" LTP \n");
                // loop from the next presynaptic spike to T
                size_t start = 0;//(size_t)(snn->next_pre_spike[global_n]);
                tpost = snn->spk_matrix[post_neuron_spk_matrix_index + T];

                //printf(" >  LTP--> Network[%lu], synapse[%lu], start = %lu, end = %lu\n", 
                //    network_number, local_n, start, T-D);

                for(size_t i = start; i<=T-D; i++){

                    tpre = snn->spk_matrix[pre_neuron_spk_matrix_index + i] + delay;

                    if(tpre - delay > -1){
                        tdiff = tpost - tpre;
                        ftdiff = (float)(tpost - tpre) + 0.15;                
                        
                        //printf(" > Synapse %lu, t = %d, tpost = %d, tpre = %d, delay = %d, pre_ind = %lu, i = %lu\n", threadId, t, tpost, tpre, delay, pre_neuron_spk_matrix_index, i);
                        if(ftdiff > 0.0 && ftdiff < 50.0){
                            dw += 1.0 * exp(-ftdiff / 5.0);
                        }
                        else if(ftdiff < 0.0 && ftdiff > -75.0){
                            dw -= 1.0 * exp(ftdiff / 5.0);
                        }
                        /*switch(snn->lr[local_n]){
                            case 0:
                                if(ftdiff > 0.0 && ftdiff < 50.0){
                                    dw += 1.0 * exp(-ftdiff / 5.0);
                                }
                                else if(ftdiff < 0.0 && ftdiff > -75.0){
                                    dw -= 1.0 * exp(ftdiff / 5.0);
                                }
                            break;
                            case 1:
                                if(ftdiff > 0.0 && ftdiff < 50.0){
                                    dw += 1.0 * exp(-ftdiff / 5.0);
                                }
                                else if(ftdiff < 0.0 && ftdiff > -75.0){
                                    dw -= 1.0 * exp(ftdiff / 5.0);
                                }
                            break;
                            case 2:
                                if(ftdiff > 0.0 && ftdiff < 50.0){
                                    dw -= 1.0 * exp(ftdiff / 5.0);
                                }
                                else if(ftdiff < 0.0 && ftdiff > -75.0){
                                    dw += 1.0 * exp(-ftdiff / 5.0);
                                }
                            break;
                        }*/
                        //printf(" >> [%lu in %lu] tpost = %d, tpre = %d; tdiff %lf; dw = %lf, fdw = %lf\n", threadId, T, tpost, tpre, ftdiff, dw, snn->dw[global_n]);

                    }
                }
                
                //snn->next_pre_spike[global_n] = T-D + 1;
            }


            // LTD: check if pre neuron fired at t - delay
            else if(T - D > 0 && snn->spk_matrix[pre_neuron_spk_matrix_index + T - D] == T - D){
                
                //printf(" LTD \n");
                size_t start = 0;//(size_t)(snn->next_post_spike[global_n]);
                tpre = snn->spk_matrix[pre_neuron_spk_matrix_index + T - D] + delay;

                //printf(" > LTD --> Network[%lu], synapse[%lu], start = %lu, end = %lu\n", 
                //    network_number, local_n, start, T-D);

                for(size_t i = start; i<T; i++){
                    
                    tpost = snn->spk_matrix[post_neuron_spk_matrix_index + i];
                    if(tpost > -1){
                        
                        tdiff = tpost - tpre;
                        ftdiff = (float)(tpost - tpre) - 0.15;                
                    
                        //printf(" > Synapse %lu, t = %d, tpost = %d, tpre = %d, delay = %d, post_ind = %lu, i = %lu\n", threadId, t, tpost, tpre, delay, post_neuron_spk_matrix_index, i);

                        if(ftdiff > 0.0 && ftdiff < 50.0){
                            dw += 1.0 * exp(-ftdiff / 5.0);
                        }
                        else if(ftdiff < 0.0 && ftdiff > -75.0){
                            dw -= 1.0 * exp(ftdiff / 5.0);
                        }
                        /*switch(snn->lr[local_n]){
                            case 0:
                                if(ftdiff > 0.0 && ftdiff < 50.0){
                                    dw += 1.0 * exp(-ftdiff / 5.0);
                                }
                                else if(ftdiff < 0.0 && ftdiff > -75.0){
                                    dw -= 1.0 * exp(ftdiff / 5.0);
                                }
                            break;
                            case 1:
                                if(ftdiff > 0.0 && ftdiff < 50.0){
                                    dw += 1.0 * exp(-ftdiff / 5.0);
                                }
                                else if(ftdiff < 0.0 && ftdiff > -75.0){
                                    dw -= 1.0 * exp(ftdiff / 5.0);
                                }
                            break;
                            case 2:
                                if(ftdiff > 0.0 && ftdiff < 50.0){
                                    dw -= 1.0 * exp(ftdiff / 5.0);
                                }
                                else if(ftdiff < 0.0 && ftdiff > -75.0){
                                    dw += 1.0 * exp(-ftdiff / 5.0);
                                }
                            break;
                        }*/

                        //printf(" >> [%lu in %lu] tpre = %d, tpost = %d; tdiff %lf; dw = %lf, fdw = %lf\n", threadId, T, tpre, tpost, ftdiff, dw, snn->dw[global_n]);

                    }
                }
            
                //snn->next_post_spike[global_n] = T;
            }
        }
        snn->dw[global_n] += dw;
        snn->w[global_n] += dw;

        //printf(" > Sample %lu, dw[%lu] = %f, final dw[%lu] = %f\n", network_number, local_n, dw, local_n, snn->dw[global_n]);

    }
}


/* STDP function */
__global__ void krnl_update_weights(GPU_SNN_t *snn, int batch_size)
{
    size_t threadId = 
        (blockIdx.x  + gridDim.x  * blockIdx.y + gridDim.x * gridDim.y * blockIdx.z) *
        (blockDim.x * blockDim.y * blockDim.z) +
        (threadIdx.x + blockDim.x * threadIdx.y + blockDim.x * blockDim.y * threadIdx.z);
    
    
    // update w in all copies 
    if(threadId < (size_t)((size_t)(snn->n_synapses))){

        float dw = 0.0;

        // summ weights of all network copies
        for(size_t i = 0; i<(size_t)(snn->n_networks); i++){

            size_t sindex = (size_t)(snn->n_synapses) * i + threadId;
            dw += snn->dw[sindex];
        }


        float init_w = snn->init_w[threadId] + (dw / float(batch_size)); // update init_w

        // update w and init_w in all copies
        for(size_t i = 0; i<(size_t)(snn->n_networks); i++){
            
            size_t sindex = (size_t)(snn->n_synapses) * i + threadId;
            snn->init_w[sindex] = init_w;
            snn->w[sindex] = init_w; 
        }

        // print final weights
        //printf(" dW[%ld] = %f\n", threadId, (dw / float(batch_size)));
    }
}


__global__ void krnl_sum_dw(GPU_SNN_t *snn)
{
    size_t threadId = 
        (blockIdx.x  + gridDim.x  * blockIdx.y + gridDim.x * gridDim.y * blockIdx.z) *
        (blockDim.x * blockDim.y * blockDim.z) +
        (threadIdx.x + blockDim.x * threadIdx.y + blockDim.x * blockDim.y * threadIdx.z);
    
    
    // update w in all copies 
    if(threadId < (size_t)((size_t)snn->n_synapses)){

        float dw = 0.0;

        // summ weights of all network copies
        for(size_t i = 0; i<(size_t)(snn->n_networks); i++){
            
            dw += snn->dw[(size_t)(snn->n_synapses) * i + threadId];
        }

        // store new dw in the first
        snn->dw[threadId] = dw;
    }
}

__global__ void krnl_update_weights_only(GPU_SNN_t *snn, int batch_size)
{
    size_t threadId = 
        (blockIdx.x  + gridDim.x  * blockIdx.y + gridDim.x * gridDim.y * blockIdx.z) *
        (blockDim.x * blockDim.y * blockDim.z) +
        (threadIdx.x + blockDim.x * threadIdx.y + blockDim.x * blockDim.y * threadIdx.z);
    
    
    // update w in all copies 
    if(threadId < (size_t)((size_t)snn->n_synapses)){
        // update w and init_w in all copies
        for(size_t i = 0; i<(size_t)(snn->n_networks); i++){
            
            size_t sindex = (size_t)(snn->n_synapses) * i + threadId;

            //printf(" >> init_w[%ld] = %f\n", threadId, (snn->init_w[sindex]));

            // only the first copy has the updated dw
            snn->init_w[sindex] = snn->init_w[sindex] + snn->dw[threadId];
            snn->w[sindex] = snn->init_w[sindex]; 
        }

        //printf("  dW[%ld] = %f\n", threadId, (snn->dw[threadId]));
    }
}



/*
 General simulation functions
*/

extern "C" void simulate_LIF_in_single_GPU(GPU_SNN_t *gpu_snn, GPU_dataset_t *gpu_dataset, GPU_results_t *gpu_results, simulation_configuration_t *conf, cuda_info_t *cuda_info, GPU_SNN_t *cpu_snn, GPU_dataset_t *cpu_dataset){
    
    size_t N, iN, S, n_samples, n_features, n_networks, E, batch_size, dev_batch_size, T, LT;
    size_t t, s, b, snn, e;
    double t_in_ts, t_out_ts, t_in = 0.0, t_out = 0.0;
    cudaError_t err;
    struct timespec start_in, end_in, start_out, end_out;


    // get general variables values
    N = cpu_snn->n_neurons;
    iN = cpu_snn->n_input_neurons;
    S = cpu_snn->n_synapses;
    T = conf->time_steps;
    LT = cpu_snn->LT;
    E = conf->epochs;

    n_samples = cpu_dataset->n_samples;
    n_features = cpu_dataset->n_features;

    batch_size = cuda_info->batch_size;
    dev_batch_size = cuda_info->dev_batch_size;
    n_networks = cuda_info->n_networks;
    
    // set dim3 struct to launch kernels

    dim3 grid_neurons(cuda_info->n_blk_neurons_x, cuda_info->n_blk_neurons_y, cuda_info->n_blk_neurons_z);
    dim3 block_neurons(cuda_info->n_thr_per_blk_neurons_x, cuda_info->n_thr_per_blk_neurons_y, cuda_info->n_thr_per_blk_neurons_z);

    dim3 grid_synapses(cuda_info->n_blk_synapses_x, cuda_info->n_blk_synapses_y, cuda_info->n_blk_synapses_z);
    dim3 block_synapses(cuda_info->n_thr_per_blk_synapses_x, cuda_info->n_thr_per_blk_synapses_y, cuda_info->n_thr_per_blk_synapses_z);

    dim3 grid_in_neurons(cuda_info->n_blk_in_neurons_x, cuda_info->n_blk_in_neurons_y, cuda_info->n_blk_in_neurons_z);
    dim3 block_in_neurons(cuda_info->n_thr_per_blk_in_neurons_x, cuda_info->n_thr_per_blk_in_neurons_y, cuda_info->n_thr_per_blk_in_neurons_z);

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
*/
    

    // loop over epochs
    for(e = 0; e<E; e++){

        // loop over batches //TODO: refactorize to function
        for(b = 0; b<n_samples; b+=batch_size){

            for(snn = 0; snn < batch_size; snn += n_networks){

                initialize_neurons_batch<<<grid_neurons, block_neurons>>>(gpu_snn, batch_size);
                err = cudaGetLastError();
                cudaCheckError(cudaPeekAtLastError());  // check launch errors
                cudaCheckError(cudaDeviceSynchronize());  // check runtime errors

                initialize_synapses_batch<<<grid_synapses, block_synapses>>>(gpu_snn, batch_size);
                err = cudaGetLastError();
                cudaCheckError(cudaPeekAtLastError());  // check launch errors
                cudaCheckError(cudaDeviceSynchronize());  // check runtime errors

            }
        }
    }



    /*
    for(e = 0; e < n_epochs; e++){
    
        // loop over batches
        for(b = 0; b < n_samples; b = b + batch_size){

            //printf(" In batch %d\n", b);
            //fflush(stdout);

            for(snn = 0; snn < batch_size; snn = snn + n_networks){
                
                // reinit the entire matrix: number of elements = (n_neurons + n_input_neurons) * L
                // max 1024 threads per block: n_neurons * n_input_neurons threads, each L iterations
                // n_blocks = (n_neurons + n_input_neurons / 1024) + 1
                
                reinit_spk_matrix<<<grid_n_spkmatrix_reinit_blk, block_n_spkmatrix_reinit_threads_per_blk>>>(gpu_snn);
                err = cudaGetLastError();
                cudaCheckError(cudaPeekAtLastError());  // check launch errors
                cudaCheckError(cudaDeviceSynchronize());  // check runtime errors
                cudaDeviceSynchronize();

                load_sample_kernel<<<grid_n_ls_blk, block_n_ls_threads_per_block>>>(gpu_snn, gpu_dataset, b + snn, n_samples);
                err = cudaGetLastError();
                cudaCheckError(cudaPeekAtLastError());  // check launch errors
                cudaCheckError(cudaDeviceSynchronize());  // check runtime errors
                cudaDeviceSynchronize();

                // reinit neurons
                reinit_neurons_kernel<<<grid_n_pn_blk, block_n_pn_threads_per_block>>>(gpu_snn); // n_neurons 8784
                err = cudaGetLastError();
                cudaCheckError(cudaPeekAtLastError());  // check launch errors
                cudaCheckError(cudaDeviceSynchronize());  // check runtime errors
                cudaDeviceSynchronize();

                // reinit synapses
                reinit_synapses_kernel<<<grid_n_uw_blk, block_n_uw_per_block>>>(gpu_snn); // n_synapses 1440187
                err = cudaGetLastError();
                cudaCheckError(cudaPeekAtLastError());  // check launch errors
                cudaCheckError(cudaDeviceSynchronize());  // check runtime errors
                cudaDeviceSynchronize();

                for(t = 0; t<time_steps; t++){

                    //printf(" In time step %d\n", t);
                    //fflush(stdout);
                    // input step
                    clock_gettime(CLOCK_MONOTONIC, &start_in);
                    lif_neuron_input<<<grid_n_pn_blk, block_n_pn_threads_per_block>>>(gpu_snn, b + snn, t, n_samples);
                    err = cudaGetLastError();
                    cudaCheckError(cudaPeekAtLastError());  // check launch errors
                    cudaCheckError(cudaDeviceSynchronize());  // check runtime errors        
                    cudaDeviceSynchronize();
                    clock_gettime(CLOCK_MONOTONIC, &end_in);
                    t_in_ts = (end_in.tv_sec - start_in.tv_sec) + (end_in.tv_nsec - start_in.tv_nsec) / 1e9;
                    t_in += t_in_ts;
                    printf(" >>> Elapsed time in %d input: %lf\n", t, t_in_ts);
                    fflush(stdout);
                    
                    // output step
                    clock_gettime(CLOCK_MONOTONIC, &start_out);
                    lif_neuron_output<<<grid_n_pn_blk, block_n_pn_threads_per_block>>>(gpu_snn, b + snn, t, gpu_results, n_samples);
                    err = cudaGetLastError();
                    cudaCheckError(cudaPeekAtLastError());  // check launch errors
                    cudaCheckError(cudaDeviceSynchronize());  // check runtime errors          
                    cudaDeviceSynchronize();
                    clock_gettime(CLOCK_MONOTONIC, &end_out);
                    t_out_ts = (end_out.tv_sec - start_out.tv_sec) + (end_out.tv_nsec - start_out.tv_nsec) / 1e9;
                    t_out += t_out_ts;
                    printf(" >>> Elapsed time in %d output: %lf\n", t, t_out_ts);
                    fflush(stdout);

                    //printf(" > Output step processed\n");
                    //fflush(stdout);
                    
                    // learning
                    if(conf->learn == 1){
                        stdp_kernel<<<grid_n_synapses_blk, block_n_synapses_per_block>>>(gpu_snn, t);
                        err = cudaGetLastError();
                        cudaCheckError(cudaPeekAtLastError());  // check launch errors
                        cudaCheckError(cudaDeviceSynchronize());  // check runtime errors  
                        cudaDeviceSynchronize();
                        //printf(" > STDP step processed\n");
                        //printf("Launch error: %s\n", cudaGetErrorString(err));
                    }
            //fflush(stdout);
                }


            }

            if(conf->learn == 1){
                // update weights
                krnl_update_weights<<<grid_n_uw_blk, block_n_uw_per_block>>>(gpu_snn, batch_size);
                err = cudaGetLastError();
                cudaCheckError(cudaPeekAtLastError());  // check launch errors
                cudaCheckError(cudaDeviceSynchronize());  // check runtime errors  
                cudaDeviceSynchronize();

                // reinit synapses (set dw to 0)
                reinit_tot_synapses_kernel<<<grid_n_synapses_blk, block_n_synapses_per_block>>>(gpu_snn); // n_synapses 1440187
                err = cudaGetLastError();
                cudaCheckError(cudaPeekAtLastError());  // check launch errors
                cudaCheckError(cudaDeviceSynchronize());  // check runtime errors
                cudaDeviceSynchronize();
            }
        }
    }
*/
    printf(" >>> Elapsed time in: %lf\n", t_in);
    fflush(stdout);

    printf(" >>> Elapsed time out: %lf\n", t_out);
    fflush(stdout);

}

extern "C" void simulate_LIF_in_multi_GPU(GPU_SNN_t *gpu_snn, GPU_dataset_t *gpu_dataset, GPU_results_t *gpu_results, simulation_configuration_t *conf, cuda_info_t *cuda_info, GPU_SNN_t *cpu_snn, GPU_dataset_t *cpu_dataset, int dev, int batch){
    
    // call kernels
    int time_steps, batch_size, dev_batch_size, n_networks;
    int t, s, b, snn;
    cudaError_t err;

    time_steps = cuda_info->time_steps;
    batch_size = cuda_info->batch_size;
    dev_batch_size = cuda_info->dev_batch_size;
    n_networks = cuda_info->n_networks_per_dev;

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


// general function to manage the execution in cuda
extern "C" void simulate_LIF_in_GPU(GPU_SNN_t **gpu_snn, GPU_dataset_t **gpu_dataset, GPU_results_t **gpu_results, simulation_configuration_t *conf, cuda_info_t *cuda_info, GPU_SNN_t *cpu_snn, GPU_dataset_t *cpu_dataset){


    // copy constant(s) to all gpu devices
    for(size_t i = 0; i<cuda_info->nDevices; i++){
        
        cudaSetDevice(i);
        cudaError_t err = cudaMemcpyToSymbol(n_samples, &cpu_dataset->n_samples, sizeof(size_t));
    }

    // if it is executed in only one device
    if(cuda_info->nDevices == 1){
        
        cudaSetDevice(0);
        simulate_LIF_in_single_GPU(gpu_snn[0], gpu_dataset[0], gpu_results[0], conf, cuda_info, cpu_snn, cpu_dataset);
    }

    // TODO: multi GPU must be corrected
    // if it is multi GPU
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
                    simulate_LIF_in_multi_GPU(gpu_snn[dev], gpu_dataset[dev], gpu_results[dev], conf, cuda_info, cpu_snn, cpu_dataset, dev, b);

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
