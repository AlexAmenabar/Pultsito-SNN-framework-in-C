#include <cuda.h>
#include <cuda_runtime.h>
#include <stdio.h>
#include <stdlib.h>

// 
#include "datasets/datasets.h"

// cuda related includes
#include "cuda/cuda_datasets.cuh"
#include "cuda/cuda_simulations_conf.cuh"


extern "C" GPU_dataset_t** cpy_dataset2GPU(GPU_dataset_t *cpu_dataset, cuda_info_t *cuda_info){

    size_t nS, nF, nSpks;
    size_t dev;

    nS = cpu_dataset->n_samples;
    nF = cpu_dataset->n_features;
    nSpks = cpu_dataset->n_spikes;

    // allocate to temporal
    GPU_dataset_t *tmp_gpu_dataset = (GPU_dataset_t *)malloc(sizeof(GPU_dataset_t));
    GPU_dataset_t **d_gpu_dataset = (GPU_dataset_t**)malloc(cuda_info->nDevices * sizeof(GPU_dataset_t*));
    cudaError_t err;

    // loop over devices
    for(dev = 0; dev < cuda_info->nDevices; dev++){
        
        // set dev as active one
        cudaSetDevice(dev);

        
        // allocate memory in GPU
        err = cudaMalloc(&(tmp_gpu_dataset->spikes), nSpks * sizeof(size_t)); // allocate memory for neurons
        if (err != cudaSuccess) 
            printf("neuron_input_synapses_offset allocation failed: %s\n", cudaGetErrorString(err));
        
        err = cudaMalloc(&(tmp_gpu_dataset->n_spikes_per_feature), nS * nF * sizeof(size_t)); // allocate memory for neurons
        if (err != cudaSuccess) 
            printf("neuron_input_synapses_offset allocation failed: %s\n", cudaGetErrorString(err));

        err = cudaMalloc(&(tmp_gpu_dataset->sample_offset), nS * sizeof(size_t)); // allocate memory for neurons
        if (err != cudaSuccess) 
            printf("neuron_input_synapses_offset allocation failed: %s\n", cudaGetErrorString(err));

        err = cudaMalloc(&(tmp_gpu_dataset->feature_offset), nS * nF * sizeof(size_t)); // allocate memory for neurons
        if (err != cudaSuccess) 
            printf("neuron_input_synapses_offset allocation failed: %s\n", cudaGetErrorString(err));

        err = cudaMalloc(&(tmp_gpu_dataset->freq), nS * nF * sizeof(size_t)); // allocate memory for neurons
        if (err != cudaSuccess) 
            printf("neuron_input_synapses_offset allocation failed: %s\n", cudaGetErrorString(err));

        err = cudaMalloc(&(tmp_gpu_dataset->first_spk), nS * nF * sizeof(size_t)); // allocate memory for neurons
        if (err != cudaSuccess) 
            printf("neuron_input_synapses_offset allocation failed: %s\n", cudaGetErrorString(err));

        // copy data to GPU
        tmp_gpu_dataset->n_samples = cpu_dataset->n_samples;
        tmp_gpu_dataset->type = cpu_dataset->type;
        tmp_gpu_dataset->n_classes = cpu_dataset->n_classes;
        tmp_gpu_dataset->n_features = cpu_dataset->n_features;
        tmp_gpu_dataset->n_spikes = cpu_dataset->n_spikes;
        cudaMemcpy(tmp_gpu_dataset->spikes, cpu_dataset->spikes, nSpks * sizeof(size_t), cudaMemcpyHostToDevice); // copy neurons information
        cudaMemcpy(tmp_gpu_dataset->n_spikes_per_feature, cpu_dataset->n_spikes_per_feature, nS * nF * sizeof(size_t), cudaMemcpyHostToDevice); // copy neurons information
        cudaMemcpy(tmp_gpu_dataset->sample_offset, cpu_dataset->sample_offset, nS * sizeof(size_t), cudaMemcpyHostToDevice); // copy neurons information
        cudaMemcpy(tmp_gpu_dataset->feature_offset, cpu_dataset->feature_offset, nS * nF * sizeof(size_t), cudaMemcpyHostToDevice); // copy neurons information
        cudaMemcpy(tmp_gpu_dataset->freq, cpu_dataset->freq, nS * nF * sizeof(size_t), cudaMemcpyHostToDevice); // copy neurons information
        cudaMemcpy(tmp_gpu_dataset->first_spk, cpu_dataset->first_spk, nS * nF * sizeof(size_t), cudaMemcpyHostToDevice); // copy neurons information

        // copy to GPU
        err = cudaMalloc(&(d_gpu_dataset[dev]), sizeof(GPU_dataset_t)); // allocate memory for neurons
        if (err != cudaSuccess) 
            printf("neuron_input_synapses_offset allocation failed: %s\n", cudaGetErrorString(err));

        cudaMemcpy(d_gpu_dataset[dev], tmp_gpu_dataset, sizeof(GPU_dataset_t), cudaMemcpyHostToDevice); // copy neurons information
    }

    // deallocate memory of temporal dataset
    free(tmp_gpu_dataset);

    return d_gpu_dataset;
}