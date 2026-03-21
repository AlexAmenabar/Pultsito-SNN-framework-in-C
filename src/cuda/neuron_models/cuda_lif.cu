#include <cuda.h>
#include <cuda_runtime.h>
#include <stdio.h>
#include <omp.h>

#include "networks/snn.h"
#include "datasets/datasets.h"
#include "simulations/simulations.h"
#include "simulations/results.h"
#include "config/config_loader.h"


#include "cuda/neuron_models/cuda_lif.cuh"
#include "cuda/cuda_simulations.cuh"
#include "cuda/cuda_simulations_conf.cuh"

// private headers
#include "cuda/priv_cuda_simulations.cuh"
#include "cuda/neuron_models/priv_cuda_lif.cuh"


#define cudaCheckError(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
{
    if (code != cudaSuccess)
    {
        fprintf(stderr, "CUDA Error: %s %s %d\n", cudaGetErrorString(code), file, line);
        if (abort) exit(code);
    }
}

void set_LIF_cuda_kernels(cuda_info_t *cuda_info){

    // CPU-GPU communication
    cuda_info->allocate_neuron_memory_cuda = &allocate_LIF_memory_cuda;
    cuda_info->cpy_neurons_CPU2GPU_cuda = &cpy_LIF_memory_cuda;

    // kernel wrappers
    cuda_info->initialize_neurons_cuda = &wrap_initialize_LIF_neurons_batch_kernel;
    cuda_info->process_input_current_cuda = &wrap_process_input_current_batch_kernel;
    cuda_info->process_dynamics_cuda = &wrap_process_V_batch_kernel;
    cuda_info->process_firing_cuda = &wrap_process_firing_batch_kernel;
}

/* GPU memory management */

void allocate_LIF_memory_cuda(GPU_SNN_t *d_snn, size_t N, size_t batch){
    
    cudaError_t err;

    err = cudaMalloc(&(d_snn->v), N * batch * sizeof(float)); 
    if (err != cudaSuccess) 
        printf("v allocation failed: %s\n", cudaGetErrorString(err));
    
    err = cudaMalloc(&(d_snn->arrI), N * batch * sizeof(float));
    if (err != cudaSuccess) 
        printf("arrI allocation failed: %s\n", cudaGetErrorString(err)); 
    
    err = cudaMalloc(&(d_snn->v_thresh), N * sizeof(float)); 
    if (err != cudaSuccess) 
        printf("v_thresh allocation failed: %s\n", cudaGetErrorString(err));
    
    err = cudaMalloc(&(d_snn->v_rest), N * sizeof(float)); 
    if (err != cudaSuccess) 
        printf("v_rest allocation failed: %s\n", cudaGetErrorString(err));
    
    err = cudaMalloc(&(d_snn->r_period), N * sizeof(int)); 
    if (err != cudaSuccess) 
        printf("r_period allocation failed: %s\n", cudaGetErrorString(err));
    
    err = cudaMalloc(&(d_snn->r_period_remain), N * batch * sizeof(int));
    if (err != cudaSuccess) 
        printf("r_period_remain allocation failed: %s\n", cudaGetErrorString(err)); 
    
    err = cudaMalloc(&(d_snn->res), N * sizeof(int)); 
    if (err != cudaSuccess) 
        printf("res allocation failed: %s\n", cudaGetErrorString(err));
}


void cpy_LIF_memory_cuda(GPU_SNN_t *d_snn, GPU_SNN_t *h_snn, size_t N){

    cudaError_t err;

    err = cudaMemcpy(d_snn->v_thresh, h_snn->v_thresh, N * sizeof(float), cudaMemcpyHostToDevice); // copy neurons information
    if (err != cudaSuccess) 
        printf("res allocation failed: %s\n", cudaGetErrorString(err));

    err = cudaMemcpy(d_snn->v_rest, h_snn->v_rest, N * sizeof(float), cudaMemcpyHostToDevice); // copy neurons information
    if (err != cudaSuccess) 
        printf("res allocation failed: %s\n", cudaGetErrorString(err));

    err = cudaMemcpy(d_snn->r_period, h_snn->r_period, N * sizeof(int), cudaMemcpyHostToDevice); // copy neurons information
    if (err != cudaSuccess) 
        printf("res allocation failed: %s\n", cudaGetErrorString(err));

    err = cudaMemcpy(d_snn->res, h_snn->res, N * sizeof(int), cudaMemcpyHostToDevice); // copy neurons information
    if (err != cudaSuccess) 
        printf("res allocation failed: %s\n", cudaGetErrorString(err));
}

/* LIF kernels */

__global__ void initialize_LIF_neurons_batch(GPU_SNN_t *snn){

    // get thread id
    size_t threadId = 
        ((size_t)blockIdx.x  + (size_t)gridDim.x  * (size_t)blockIdx.y + (size_t)gridDim.x * (size_t)gridDim.y * (size_t)blockIdx.z) *
        ((size_t)blockDim.x * (size_t)blockDim.y * (size_t)blockDim.z) +
        ((size_t)threadIdx.x + (size_t)blockDim.x * (size_t)threadIdx.y + (size_t)blockDim.x * (size_t)blockDim.y * (size_t)threadIdx.z);

    if(threadId < snn->n_neurons * dev_batch_size){

        snn->v[threadId] = snn->v_rest[threadId / dev_batch_size];
        snn->r_period_remain[threadId] = 0;
        snn->arrI[threadId] = 0.0;

        if(learn){

            snn->post_fired[threadId] = 0;
            snn->post_trace[threadId] = 0.0;
        }
    }    
}

__global__ void process_input_currect_batch_new(GPU_SNN_t *gpu_snn, size_t iN, size_t N, size_t t, size_t gt, size_t N_per_launch, size_t launch){

    // get thread id: iN * batch_size * T
    size_t threadId = 
        (blockIdx.x  + gridDim.x  * blockIdx.y + gridDim.x * gridDim.y * blockIdx.z) *
        (blockDim.x * blockDim.y * blockDim.z) +
        (threadIdx.x + blockDim.x * threadIdx.y + blockDim.x * blockDim.y * threadIdx.z);


    // initialize shared memory
    extern __shared__ float sharedI[];


    size_t neuron_index, b, p, blk_b, localThreadId, init_neuron;

    // get initial neuron index: launch number * neurons launched per launch
    init_neuron = launch * N_per_launch;

    // get neuron index
    neuron_index = threadId / (dev_batch_size * thrN) + init_neuron;
    
    // get section identifier
    p = (threadId / batch_size_per_block) % thrN;

    // batch identifier in the CUDA block
    blk_b = threadId % batch_size_per_block;
    
    // compute global batch identifier 
    b = threadId / (thrN * batch_size_per_block) % blocks_per_batch; 
    b = b * batch_size_per_block + blk_b;

    // shared memory index
    localThreadId = p * batch_size_per_block + blk_b;
    
    //printf(" theadId = %llu, N per launch = %llu, thrN = %llu, dev_batch_size = %llu, neuron_index = %llu, N = %llu\n", 
    //        threadId, N_per_launch, thrN, dev_batch_size, neuron_index, N);

    if(threadId < N_per_launch * thrN * dev_batch_size && neuron_index < N){

        //printf(" Thread %llu (local thread = %llu), neuron = %llu, p = %llu, b = %llu, blk_b = %llu\n", threadId, localThreadId, neuron_index, p, b, blk_b);
        //sharedI[localThreadId] = 0.0;

        // variables declaration
        size_t synapse_index, g_synapse_index, in_neuron_index;
        int delay, spk_time; 
        float spk;  

        size_t g_neuron_index = neuron_index * dev_batch_size; // for copied variables

        // get neuron information
        size_t iS = gpu_snn->n_neuron_input_synapses[neuron_index]; 
        size_t base_synapse = gpu_snn->neuron_input_synapses_offset[neuron_index]; // there are several copies of each synapse, so the offset depends
    
        float I = 0.0;

        // check wether the neuron is in refractory period: pre_fired flag not activated although it should be
        if(gpu_snn->r_period_remain[g_neuron_index + b] <= 0){
        
            // calculate synapses to be processed
            size_t first_synapse, last_synapse, n_synapses, r_synapses;
            
            /*if(thrN == 1){

                first_synapse = 0;
                last_synapse = iS;
            }
            else{*/
                
                n_synapses = iS / thrN;

                // set remaining
                r_synapses = iS % thrN;

                if(p < r_synapses)
                    n_synapses += 1;
                
                // compute first and last synapses
                first_synapse = n_synapses * p;
                last_synapse = n_synapses * (p+1);

                if(p >= r_synapses){
                    
                    first_synapse += p;
                    last_synapse += p;
                }
            //}
            //printf(" Neuron %llu: first synapse = %llu, last synapse %llu (iS = %llu, n_synapses = %llu, r_synapses = %llu, p = %llu, thrN = %llu)\n", 
            //            neuron_index, first_synapse, last_synapse, iS, n_synapses, r_synapses, p, thrN);

            // loop over input synapses
            for(size_t j = first_synapse; j<last_synapse; j++){

                // get synapse index. In not copied variables use synapse_index, in copied g_synapse_index 
                synapse_index = base_synapse + j;
                g_synapse_index = synapse_index * dev_batch_size; // scale base synapse, since now there are B copies of each one
                
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
                    size_t index = ((iN + N) * dev_batch_size * (size_t)spk_time) + ((in_neuron_index * dev_batch_size));

                    // get spike
                    spk = (float)(gpu_snn->spk_matrix[index + b]); 

                    // update I if r period remain <= 0 and spike received
                    //if(gpu_snn->r_period_remain[g_neuron_index + b] <= 0){
                    I += gpu_snn->w[g_synapse_index + b] * spk; // spk = 0 / 1
                    //}
                }
            }
        }

        // write I in shared memory
        sharedI[localThreadId] = I;
        
        //atomicAdd(&gpu_snn->arrI[neuron_index * batch_size + b], I);
    }

    // sync threads
    __syncthreads();

    // the first neuron accumulates I and writes in sahred memory
    if(threadId < N_per_launch * thrN * dev_batch_size && neuron_index < N && p == 0){

        float I = 0.0;
        for(size_t j = 0; j<thrN; j++){

            //I += sharedI[localThreadId + j];
            I += sharedI[localThreadId + j * batch_size_per_block];
        }

        gpu_snn->arrI[neuron_index * dev_batch_size + b] = I;
    }
    

    // tree-based sum
    //if(threadId < N_per_launch * thrN * dev_batch_size && neuron_index < N){

    /*char active = (threadId < N_per_launch * thrN * dev_batch_size && neuron_index < N);


        for(int stride = thrN/2; stride > 0; stride >>= 1){
            
            __syncthreads();
            if(active && p < (size_t)stride)
                sharedI[localThreadId] += sharedI[localThreadId + (size_t)(stride) * batch_size_per_block];
        }

    __syncthreads();


        // the first thread of each section stores the final result in the memory
        if(active && p == 0)
            gpu_snn->arrI[neuron_index * dev_batch_size + b] = sharedI[localThreadId];*/
    //}
}


/*
__global__ void process_input_currect_batch(GPU_SNN_t *gpu_snn, size_t iN, size_t N, size_t t, size_t gt, size_t N_per_launch, size_t launch){

    // get thread id: iN * batch_size * T
    size_t threadId = 
        (blockIdx.x  + gridDim.x  * blockIdx.y + gridDim.x * gridDim.y * blockIdx.z) *
        (blockDim.x * blockDim.y * blockDim.z) +
        (threadIdx.x + blockDim.x * threadIdx.y + blockDim.x * blockDim.y * threadIdx.z);


    // initialize shared memory
    extern __shared__ float sharedI[];


    size_t neuron_index, b, p, blk_b, localThreadId, init_neuron;

    // get initial neuron index: launch number * neurons launched per launch
    init_neuron = launch * N_per_launch;

    // get neuron index
    neuron_index = threadId / (dev_batch_size * thrN) + init_neuron;
    
    // get section identifier
    p = (threadId / batch_size_per_block) % thrN;

    // batch identifier in the CUDA block
    blk_b = threadId % batch_size_per_block;
    
    // compute global batch identifier 
    b = threadId / (thrN * batch_size_per_block) % blocks_per_batch; 
    b = b * batch_size_per_block + blk_b;

    // shared memory index
    localThreadId = p * batch_size_per_block + blk_b;


    
    if(threadId < N_per_launch * thrN * dev_batch_size && neuron_index < N){

        //printf(" Thread %llu (local thread = %llu), neuron = %llu, p = %llu, b = %llu, blk_b = %llu\n", threadId, localThreadId, neuron_index, p, b, blk_b);
        sharedI[localThreadId] = 0.0;

        // variables declaration
        size_t synapse_index, g_synapse_index, in_neuron_index;
        int delay, spk_time; 
        float spk;  

        size_t g_neuron_index = neuron_index * dev_batch_size; // for copied variables

        // get neuron information
        size_t iS = gpu_snn->n_neuron_input_synapses[neuron_index]; 
        size_t base_synapse = gpu_snn->neuron_input_synapses_offset[neuron_index]; // there are several copies of each synapse, so the offset depends
    
        float I = 0.0;

        // check wether the neuron is in refractory period: pre_fired flag not activated although it should be
        if(gpu_snn->r_period_remain[g_neuron_index + b] <= 0){
        
            // calculate synapses to be processed
            size_t first_synapse, last_synapse, n_synapses, r_synapses;
            
            n_synapses = iS / thrN;
            r_synapses = iS % thrN;
            
            // in case n_synapses is 0
            if(n_synapses == 0 && p < iS){
                n_synapses = 1;
            }

            // compute first and last synapses

            first_synapse = n_synapses * p;
            last_synapse = n_synapses * (p+1);

            if(n_synapses > 0){

                if(p < r_synapses){

                    first_synapse = first_synapse + p;
                    last_synapse = last_synapse + p + 1;
                }
                else{

                    first_synapse = first_synapse + r_synapses + 1;
                    last_synapse = last_synapse + r_synapses + 1;
                }
            }

            // loop over input synapses
            for(size_t j = first_synapse; j<last_synapse; j++){

                // get synapse index. In not copied variables use synapse_index, in copied g_synapse_index 
                synapse_index = base_synapse + j;
                g_synapse_index = synapse_index * dev_batch_size; // scale base synapse, since now there are B copies of each one
                
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
                    size_t index = ((iN + N) * dev_batch_size * (size_t)spk_time) + ((in_neuron_index * dev_batch_size));

                    // get spike
                    spk = (float)(gpu_snn->spk_matrix[index + b]); 

                    // update I if r period remain <= 0 and spike received
                    //if(gpu_snn->r_period_remain[g_neuron_index + b] <= 0){
                        
                        I += gpu_snn->w[g_synapse_index + b] * spk; // spk = 0 / 1
                    //}
                }
            }
        }

        // write I in shared memory
        sharedI[localThreadId] = I;
        
        //atomicAdd(&gpu_snn->arrI[neuron_index * batch_size + b], I);
    }

    // sync threads
    __syncthreads();

    // the first neuron accumulates I and writes in sahred memory
    /*if(threadId < N_per_launch * thrN * dev_batch_size && neuron_index < N && p == 0){

        float I = 0.0;
        for(size_t j = 0; j<thrN; j++){

            //I += sharedI[localThreadId + j];
            I += sharedI[localThreadId + j * batch_size_per_block];
        }

        gpu_snn->arrI[neuron_index * dev_batch_size + b] = I;
    }*/
    

    // tree-based sum
    //if(threadId < N_per_launch * thrN * dev_batch_size && neuron_index < N){
/*
    char active = (threadId < N_per_launch * thrN * dev_batch_size && neuron_index < N);


        for(int stride = thrN/2; stride > 0; stride >>= 1){
            
            __syncthreads();
            if(active && p < (size_t)stride)
                sharedI[localThreadId] += sharedI[localThreadId + (size_t)(stride) * batch_size_per_block];
        }

    __syncthreads();


        // the first thread of each section stores the final result in the memory
        if(active && p == 0)
            gpu_snn->arrI[neuron_index * dev_batch_size + b] = sharedI[localThreadId];
    //}
}*/


__global__ void process_V_batch(GPU_SNN_t *gpu_snn, size_t N){

    // get thread id: iN * batch_size * T
    size_t threadId = 
        (blockIdx.x  + gridDim.x  * blockIdx.y + gridDim.x * gridDim.y * blockIdx.z) *
        (blockDim.x * blockDim.y * blockDim.z) +
        (threadIdx.x + blockDim.x * threadIdx.y + blockDim.x * blockDim.y * threadIdx.z);


    if(threadId < N * dev_batch_size){

        size_t i = threadId / dev_batch_size; // neuron index
        size_t b = threadId % dev_batch_size; // batch position

        size_t neuron_index = i;
        size_t g_neuron_index = i * dev_batch_size;

        const float alpha = 0.95f;
        const float beta = 0.05f;

        // loop over batches
        if (gpu_snn->r_period_remain[g_neuron_index + b] <= 0){
            
            gpu_snn->v[g_neuron_index + b] = alpha * gpu_snn->v[g_neuron_index + b] + beta * gpu_snn->v_rest[neuron_index] + gpu_snn->arrI[g_neuron_index + b];
        }
    }
}

__global__ void process_firing_batch(GPU_SNN_t *gpu_snn, GPU_results_t *gpu_results, size_t iN, size_t N, size_t t){

    // get thread id: iN * batch_size * T
    size_t threadId = 
        (blockIdx.x  + gridDim.x  * blockIdx.y + gridDim.x * gridDim.y * blockIdx.z) *
        (blockDim.x * blockDim.y * blockDim.z) +
        (threadIdx.x + blockDim.x * threadIdx.y + blockDim.x * blockDim.y * threadIdx.z);


    if(threadId < N * dev_batch_size){

        size_t i = threadId / dev_batch_size; // neuron index
        size_t b = threadId % dev_batch_size; // batch position

        size_t neuron_index = i;
        size_t g_neuron_index = i * dev_batch_size;
        
        // compute index to write in the spk matrix
        // [T * B * (iN + N)]
        size_t idx = (iN + N) * dev_batch_size * t + dev_batch_size * (iN + neuron_index); // index of the first neuron in the spike matrix in time step t
        
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


/* Kernel wrappers */

void wrap_initialize_LIF_neurons_batch_kernel(GPU_SNN_t *gpu_snn, cuda_info_t *cuda_info, size_t dev){
    
    //printf(" In wrapper..\n");
    //fflush(stdout);

    dim3 grid_neurons(cuda_info->n_blk_neurons_x[dev], cuda_info->n_blk_neurons_y[dev], cuda_info->n_blk_neurons_z[dev]);
    dim3 block_neurons(cuda_info->n_thr_per_blk_neurons_x[dev], cuda_info->n_thr_per_blk_neurons_y[dev], cuda_info->n_thr_per_blk_neurons_z[dev]);

    initialize_LIF_neurons_batch<<<grid_neurons, block_neurons>>>(gpu_snn);
    cudaError_t err = cudaGetLastError();
    cudaCheckError(cudaPeekAtLastError());  // check launch errors
    cudaCheckError(cudaDeviceSynchronize());  // check runtime errors
    cudaDeviceSynchronize();
}

void wrap_process_input_current_batch_kernel(GPU_SNN_t *gpu_snn, size_t iN, size_t N, size_t t, size_t gt, cuda_info_t *cuda_info, size_t dev){
    
    dim3 grid_is(cuda_info->n_blk_is_x[dev], cuda_info->n_blk_is_y[dev], cuda_info->n_blk_is_z[dev]);
    dim3 block_is(cuda_info->n_thr_per_blk_is_x[dev], cuda_info->n_thr_per_blk_is_y[dev], cuda_info->n_thr_per_blk_is_z[dev]);
    cudaError_t err;

    for(size_t l = 0; l<cuda_info->n_neuron_launchs; l++){

        //process_input_currect_batch<<<grid_is, block_is, conf->thrN * cuda_info->batch_size_per_block * sizeof(float)>>>(gpu_snn, iN, N, lt, t, cuda_info->n_neurons_per_launch, l);
        process_input_currect_batch_new<<<grid_is, block_is, cuda_info->thrN * cuda_info->batch_size_per_block * sizeof(float)>>>(gpu_snn, iN, N, t, gt, cuda_info->n_neurons_per_launch, l);
        err = cudaGetLastError();
        cudaCheckError(cudaPeekAtLastError());  // check launch errors
        cudaCheckError(cudaDeviceSynchronize());  // check runtime errors
        cudaDeviceSynchronize();   
    }
}

void wrap_process_V_batch_kernel(GPU_SNN_t *gpu_snn, size_t N, cuda_info_t *cuda_info, size_t dev){
    
    dim3 grid_neurons(cuda_info->n_blk_neurons_x[dev], cuda_info->n_blk_neurons_y[dev], cuda_info->n_blk_neurons_z[dev]);
    dim3 block_neurons(cuda_info->n_thr_per_blk_neurons_x[dev], cuda_info->n_thr_per_blk_neurons_y[dev], cuda_info->n_thr_per_blk_neurons_z[dev]);

    process_V_batch<<<grid_neurons, block_neurons>>>(gpu_snn, N);
    cudaError_t err = cudaGetLastError();
    cudaCheckError(cudaPeekAtLastError());  // check launch errors
    cudaCheckError(cudaDeviceSynchronize());  // check runtime errors
    cudaDeviceSynchronize();   
}

void wrap_process_firing_batch_kernel(GPU_SNN_t *gpu_snn, GPU_results_t *gpu_results, size_t iN, size_t N, size_t t, cuda_info_t *cuda_info, size_t dev){
    
    dim3 grid_neurons(cuda_info->n_blk_neurons_x[dev], cuda_info->n_blk_neurons_y[dev], cuda_info->n_blk_neurons_z[dev]);
    dim3 block_neurons(cuda_info->n_thr_per_blk_neurons_x[dev], cuda_info->n_thr_per_blk_neurons_y[dev], cuda_info->n_thr_per_blk_neurons_z[dev]);

    process_firing_batch<<<grid_neurons, block_neurons>>>(gpu_snn, gpu_results, iN, N, t);
    cudaError_t err = cudaGetLastError();
    cudaCheckError(cudaPeekAtLastError());  // check launch errors
    cudaCheckError(cudaDeviceSynchronize());  // check runtime errors
    cudaDeviceSynchronize();
}