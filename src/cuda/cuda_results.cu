#include <cuda.h>
#include <cuda_runtime.h>
#include <stdio.h>
#include <stdlib.h>

// 
#include "simulations/results.h"

// cuda related includes
#include "cuda/cuda_results.cuh"
#include "cuda/cuda_simulations_conf.cuh"
#include "cuda/cuda_simulations.cuh"

// private headers
#include "cuda/priv_cuda_results.cuh"
#include "cuda/priv_cuda_simulations.cuh"


extern "C" GPU_results_t** initialize_results_str_in_GPU(size_t nDevices, size_t n, cuda_info_t *cuda_info){

    size_t dev, dev_batch_size;
    GPU_results_t *tmp_r, **d_results;
    cudaError_t err;

    // allocate memory
    tmp_r = (GPU_results_t*)malloc(sizeof(GPU_results_t));
    d_results = (GPU_results_t**)malloc(nDevices * sizeof(GPU_results_t*));

    // loop over devices and allocate memory for results
    for(dev = 0; dev < nDevices; dev ++){
    
        // set dev as active one
        cudaSetDevice(dev);

        // batch size per dev
        dev_batch_size = cuda_info->dev_batch_size[dev];

        // allocate memory to store the number of spikes
        err = cudaMalloc(&(tmp_r->n_spks), n * dev_batch_size * sizeof(int));
        if (err != cudaSuccess) 
            printf("neuron_input_synapses_offset allocation failed: %s\n", cudaGetErrorString(err));


        // final structure
        err = cudaMalloc(&(d_results[dev]), sizeof(GPU_results_t)); // allocate memory for neurons
        if (err != cudaSuccess) 
            printf("neuron_input_synapses_offset allocation failed: %s\n", cudaGetErrorString(err));        
        
        // copy data to GPU
        cudaMemcpy(d_results[dev], tmp_r, sizeof(GPU_results_t), cudaMemcpyHostToDevice); // copy neurons information
    }
    
    // deallocate temporal struct memory
    free(tmp_r);

    return d_results;
}


extern "C" GPU_results_t** initialize_results_str_in_CPU(size_t nDevices, size_t n, cuda_info_t *cuda_info){

    size_t dev, dev_batch_size;
    GPU_results_t **d_results;

    // allocate memory
    GPU_results_t** h_results = (GPU_results_t**)malloc(nDevices * sizeof(GPU_results_t*));

    // loop over devices and allocate memory for results
    for(dev = 0; dev < nDevices; dev ++){

        // batch size per dev
        dev_batch_size = cuda_info->dev_batch_size[dev];

        // allocate memory to store the number of spikes
        h_results[dev] = (GPU_results_t*)calloc(1, sizeof(GPU_results_t));
        h_results[dev]->n_spks = (int*)calloc(n * dev_batch_size, sizeof(int));
    }
    
    return h_results;
}

extern "C" void cpy_batch_results_struct_GPU2CPU(GPU_results_t *cpu_results, GPU_results_t *gpu_results, cuda_info_t *cuda_info, size_t N, size_t T){

    size_t dev;
    size_t batch_size = cuda_info->batch_size;
    size_t dev_batch_size, dev_batch_offset;
    size_t nDevices = cuda_info->nDevices;

    // allocate memory for temporal results structure
    GPU_results_t *tmp_results = (GPU_results_t*)malloc(sizeof(GPU_results_t));
    

    // get batch size and offset of the device
    dev_batch_size = cuda_info->dev_batch_size[dev];
    dev_batch_offset = cuda_info->dev_batch_offset[dev];

    // copy from GPU to intermediate CPU structure
    cudaMemcpy(tmp_results, gpu_results, sizeof(GPU_results_t), cudaMemcpyDeviceToHost);

    // copy results (REVISE)
    cudaMemcpy(cpu_results->n_spks, tmp_results->n_spks, N * dev_batch_size * sizeof(int), cudaMemcpyDeviceToHost);

    free(tmp_results);
}

extern "C" void cpy_batch_results_GPU2CPU(GPU_results_t *cpu_results, GPU_results_t **gpu_results, cuda_info_t *cuda_info, size_t N, size_t T){

    size_t dev;
    size_t batch_size = cuda_info->batch_size;
    size_t dev_batch_size, dev_batch_offset;
    size_t nDevices = cuda_info->nDevices;

    // allocate memory for temporal results structure
    GPU_results_t *tmp_results = (GPU_results_t*)calloc(1, sizeof(GPU_results_t));
    GPU_results_t *cpu_results_per_dev = (GPU_results_t*)calloc(1, sizeof(GPU_results_t));
    
    // copy results on each device to the CPU
    for(dev = 0; dev < nDevices; dev ++){

        // set dev as active one
        cudaSetDevice(dev);

        // get batch size and offset of the device
        dev_batch_size = cuda_info->dev_batch_size[dev];
        dev_batch_offset = cuda_info->dev_batch_offset[dev];
        cpu_results_per_dev->n_spks = (int*)calloc(N * dev_batch_size, sizeof(int));

        // copy internal GPU pointers to the CPU
        cudaMemcpy(tmp_results, gpu_results[dev], sizeof(GPU_results_t), cudaMemcpyDeviceToHost);

        // copy results to CPU structures // n_spks [N x batch_size]
        cudaMemcpy(cpu_results_per_dev->n_spks, tmp_results->n_spks, N * dev_batch_size * sizeof(int), cudaMemcpyDeviceToHost);

        // copy to final structure
        for(size_t i = 0; i<dev_batch_size; i++){

            for(size_t j = 0; j<N; j++){

                cpu_results->n_spks[j * batch_size + i + dev_batch_offset] = cpu_results_per_dev->n_spks[j * dev_batch_size + i];
            }
        }

        free(cpu_results_per_dev->n_spks);
    }

    free(tmp_results);
}

__global__ void reinitialize_results_number_of_spikes_GPU(GPU_results_t *results, size_t N){

    // get thread id
    size_t threadId = 
        ((size_t)blockIdx.x  + (size_t)gridDim.x  * (size_t)blockIdx.y + (size_t)gridDim.x * (size_t)gridDim.y * (size_t)blockIdx.z) *
        ((size_t)blockDim.x * (size_t)blockDim.y * (size_t)blockDim.z) +
        ((size_t)threadIdx.x + (size_t)blockDim.x * (size_t)threadIdx.y + (size_t)blockDim.x * (size_t)blockDim.y * (size_t)threadIdx.z);

    if(threadId < N * dev_batch_size){

        results->n_spks[threadId] = 0;
    }   
}