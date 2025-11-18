#include <cuda.h>
#include <cuda_runtime.h>
#include <stdio.h>

#include "snn_library.h"
#include "cuda/GPU_simulations.cuh"
#include "neuron_models/GPU_lif_neuron.cuh"

// functions below should be moved to a more general file

extern "C" GPU_SNN_t* cpy_SNN2GPU(GPU_SNN_t *cpu_snn, cuda_info_t *cuda_info){
    
    int i;

    // tmp structure to 
    GPU_SNN_t *tmp_gpu_snn = (GPU_SNN_t*)malloc(sizeof(GPU_SNN_t)); // input neurons???
    
    //uint64t n_neurons, n_input_neurons, n_networks, 


    // allocate memory in GPU
    cudaMalloc(&(tmp_gpu_snn->v), (size_t)((size_t)(cpu_snn->n_neurons) * (size_t)(cuda_info->n_networks) * (size_t)(sizeof(float)))); // allocate memory for neurons
    cudaMalloc(&(tmp_gpu_snn->v_thresh), (size_t)(cpu_snn->n_neurons * sizeof(float))); // allocate memory for neurons
    cudaMalloc(&(tmp_gpu_snn->v_rest), (size_t)(cpu_snn->n_neurons * sizeof(float))); // allocate memory for neurons
    cudaMalloc(&(tmp_gpu_snn->r_period), (size_t)(cpu_snn->n_neurons * sizeof(int))); // allocate memory for neurons
    cudaMalloc(&(tmp_gpu_snn->r_period_remain), (size_t)(cpu_snn->n_neurons * cuda_info->n_networks * sizeof(int))); // allocate memory for neurons
    cudaMalloc(&(tmp_gpu_snn->res), (size_t)(cpu_snn->n_neurons * sizeof(int))); // allocate memory for neurons
    cudaMalloc(&(tmp_gpu_snn->t_last_spike), (size_t)(cpu_snn->n_neurons * cuda_info->n_networks * sizeof(int))); // allocate memory for neurons
    cudaMalloc(&(tmp_gpu_snn->n_neuron_input_synapses), (size_t)(cpu_snn->n_neurons * sizeof(int)));
    cudaMalloc(&(tmp_gpu_snn->neuron_input_synapses_offset), (size_t)(cpu_snn->n_neurons * sizeof(int)));
    //cudaMalloc(&(gpu_snn->last_spk), snn->n_neurons * sizeof(int));
    //cudaMalloc(&(gpu_snn->next_spk), (snn->n_synapses - snn->n_output_synapses) * sizeof(int));
    
    cudaMalloc(&(tmp_gpu_snn->w), (size_t)(cpu_snn->n_synapses * cuda_info->n_networks * sizeof(float))); // allocate memory for neurons
    cudaMalloc(&(tmp_gpu_snn->dw), (size_t)(cpu_snn->n_synapses * cuda_info->n_networks * sizeof(float))); // allocate memory for neurons
    cudaMalloc(&(tmp_gpu_snn->delay), (size_t)(cpu_snn->n_synapses * sizeof(int))); // allocate memory for neurons
    cudaMalloc(&(tmp_gpu_snn->lr), (size_t)(cpu_snn->n_synapses * sizeof(int))); // allocate memory for neurons
    cudaMalloc(&(tmp_gpu_snn->pre_neuron_index), (size_t)(cpu_snn->n_synapses * sizeof(int))); // allocate memory for neurons
    cudaMalloc(&(tmp_gpu_snn->post_neuron_index), (size_t)(cpu_snn->n_synapses * sizeof(int))); // allocate memory for neurons

    // TODO: n_input_synapses???
    cudaMalloc(&(tmp_gpu_snn->spk_matrix), (size_t)(size_t(cpu_snn->n_neurons + cpu_snn->n_input_neurons) * (size_t)cpu_snn->max_spikes * (size_t)cuda_info->n_networks * sizeof(int))); // allocate memory for neurons

    // cpy arrays and data to GPU
    tmp_gpu_snn->n_neurons = cpu_snn->n_neurons;
    tmp_gpu_snn->n_input_neurons = cpu_snn->n_input_neurons;
    tmp_gpu_snn->n_output_neurons = cpu_snn->n_output_neurons;
    tmp_gpu_snn->n_synapses = cpu_snn->n_synapses;
    tmp_gpu_snn->n_input_synapses = cpu_snn->n_input_synapses;
    tmp_gpu_snn->n_output_synapses = cpu_snn->n_output_synapses;
    tmp_gpu_snn->max_spikes = cpu_snn->max_spikes;
    tmp_gpu_snn->n_networks = cpu_snn->n_networks;

    cudaMemcpy(tmp_gpu_snn->v_thresh, cpu_snn->v_thresh, (size_t)(cpu_snn->n_neurons * sizeof(float)), cudaMemcpyHostToDevice); // copy neurons information
    cudaMemcpy(tmp_gpu_snn->v_rest, cpu_snn->v_rest, (size_t)(cpu_snn->n_neurons * sizeof(float)), cudaMemcpyHostToDevice); // copy neurons information
    cudaMemcpy(tmp_gpu_snn->r_period, cpu_snn->r_period, (size_t)(cpu_snn->n_neurons * sizeof(int)), cudaMemcpyHostToDevice); // copy neurons information
    cudaMemcpy(tmp_gpu_snn->res, cpu_snn->res, (size_t)(cpu_snn->n_neurons * sizeof(int)), cudaMemcpyHostToDevice); // copy neurons information
    cudaMemcpy(tmp_gpu_snn->n_neuron_input_synapses, cpu_snn->n_neuron_input_synapses, (size_t)(cpu_snn->n_neurons * sizeof(int)), cudaMemcpyHostToDevice); // copy neurons information
    cudaMemcpy(tmp_gpu_snn->neuron_input_synapses_offset, cpu_snn->neuron_input_synapses_offset, (size_t)(cpu_snn->n_neurons * sizeof(int)), cudaMemcpyHostToDevice); // copy neurons information
    //cudaMemcpy(gpu_snn->last_spk, gpu_snn_mapper->last_spk, snn->n_neurons * sizeof(int), cudaMemcpyHostToDevice); // copy neurons information
    //cudaMemcpy(gpu_snn->next_spk, gpu_snn_mapper->next_spk, (snn->n_synapses - snn->n_output_synapses) * sizeof(int), cudaMemcpyHostToDevice); // copy neurons information



    cudaMemcpy(tmp_gpu_snn->delay, cpu_snn->delay, (size_t)(cpu_snn->n_synapses * sizeof(int)), cudaMemcpyHostToDevice); // copy neurons information
    cudaMemcpy(tmp_gpu_snn->lr, cpu_snn->lr, (size_t)(cpu_snn->n_synapses * sizeof(int)), cudaMemcpyHostToDevice); // copy neurons information
    cudaMemcpy(tmp_gpu_snn->pre_neuron_index, cpu_snn->pre_neuron_index, (size_t)(cpu_snn->n_synapses * sizeof(int)), cudaMemcpyHostToDevice); // copy neurons information
    cudaMemcpy(tmp_gpu_snn->post_neuron_index, cpu_snn->post_neuron_index, (size_t)(cpu_snn->n_synapses * sizeof(int)), cudaMemcpyHostToDevice); // copy neurons information


    // copy n times writible arrays
    for(i = 0; i<cuda_info->n_networks; i++){
        
        cudaMemcpy(tmp_gpu_snn->v + i * cpu_snn->n_neurons, cpu_snn->v, (size_t)(cpu_snn->n_neurons * sizeof(float)), cudaMemcpyHostToDevice); // copy neurons information
        cudaMemcpy(tmp_gpu_snn->r_period_remain +i * cpu_snn->n_neurons, cpu_snn->r_period_remain, (size_t)(cpu_snn->n_neurons * sizeof(int)), cudaMemcpyHostToDevice); // copy neurons information
        cudaMemcpy(tmp_gpu_snn->t_last_spike + i * cpu_snn->n_neurons, cpu_snn->t_last_spike, (size_t)(cpu_snn->n_neurons * sizeof(int)), cudaMemcpyHostToDevice); // copy neurons information

        cudaMemcpy(tmp_gpu_snn->w + i * cpu_snn->n_synapses, cpu_snn->w, (size_t)(cpu_snn->n_synapses * sizeof(float)), cudaMemcpyHostToDevice); // copy neurons information
        cudaMemcpy(tmp_gpu_snn->dw + i * cpu_snn->n_synapses, cpu_snn->dw, (size_t)(cpu_snn->n_synapses * sizeof(float)), cudaMemcpyHostToDevice); // copy neurons information
        
        // n_input?
        //cudaMemcpy(tmp_gpu_snn->spk_matrix + i * (cpu_snn->n_neurons + cpu_snn->n_input_neurons) * cpu_snn->max_spikes, 
        //            cpu_snn->spk_matrix, (size_t)((cpu_snn->n_neurons + cpu_snn->n_input_neurons) * cpu_snn->max_spikes * sizeof(int)), cudaMemcpyHostToDevice);
    }
    

    // final structure allocation and data copy
    GPU_SNN_t *d_GPU_SNN;
    cudaMalloc(&d_GPU_SNN, sizeof(GPU_SNN_t));
    cudaMemcpy(d_GPU_SNN, tmp_gpu_snn, sizeof(GPU_SNN_t), cudaMemcpyHostToDevice); // copy neurons information

    // free temporal GPU SNN structure
    free(tmp_gpu_snn);


    return d_GPU_SNN;
}

extern "C" GPU_dataset_t* cpy_dataset2GPU(GPU_dataset_t *cpu_dataset, cuda_info_t *cuda_info){
        
    // allocate to temporal
    GPU_dataset_t *tmp_gpu_dataset = (GPU_dataset_t *)malloc(sizeof(GPU_dataset_t));
    cudaMalloc(&(tmp_gpu_dataset->spikes), cpu_dataset->n_spikes * sizeof(int)); // allocate memory for neurons
    cudaMalloc(&(tmp_gpu_dataset->n_spikes_per_feature), cpu_dataset->n_samples * cpu_dataset->n_features * sizeof(int)); // allocate memory for neurons
    cudaMalloc(&(tmp_gpu_dataset->sample_offset), cpu_dataset->n_samples * sizeof(size_t)); // allocate memory for neurons
    cudaMalloc(&(tmp_gpu_dataset->feature_offset), cpu_dataset->n_samples * cpu_dataset->n_features * sizeof(size_t)); // allocate memory for neurons

    // memcpy
    tmp_gpu_dataset->n_samples = cpu_dataset->n_samples;
    tmp_gpu_dataset->type = cpu_dataset->type;
    tmp_gpu_dataset->n_classes = cpu_dataset->n_classes;
    tmp_gpu_dataset->n_features = cpu_dataset->n_features;
    tmp_gpu_dataset->n_spikes = cpu_dataset->n_spikes;
    cudaMemcpy(tmp_gpu_dataset->spikes, cpu_dataset->spikes, cpu_dataset->n_spikes * sizeof(int), cudaMemcpyHostToDevice); // copy neurons information
    cudaMemcpy(tmp_gpu_dataset->n_spikes_per_feature, cpu_dataset->n_spikes_per_feature, cpu_dataset->n_samples * cpu_dataset->n_features * sizeof(int), cudaMemcpyHostToDevice); // copy neurons information
    cudaMemcpy(tmp_gpu_dataset->sample_offset, cpu_dataset->sample_offset, cpu_dataset->n_samples * sizeof(size_t), cudaMemcpyHostToDevice); // copy neurons information
    cudaMemcpy(tmp_gpu_dataset->feature_offset, cpu_dataset->feature_offset, cpu_dataset->n_samples * cpu_dataset->n_features * sizeof(size_t), cudaMemcpyHostToDevice); // copy neurons information


    // copy to GPU
    GPU_dataset_t *d_gpu_dataset;
    cudaMalloc(&(d_gpu_dataset), sizeof(GPU_dataset_t)); // allocate memory for neurons
    cudaMemcpy(d_gpu_dataset, tmp_gpu_dataset, sizeof(GPU_dataset_t), cudaMemcpyHostToDevice); // copy neurons information

    // deallocate memory of temporal dataset
    free(tmp_gpu_dataset);


    return d_gpu_dataset;
}


extern "C" GPU_results_t* initialize_results_str_in_GPU(simulation_configuration_t *conf, int n, int s){
    
    GPU_results_t *tmp_r, *d_results;
    tmp_r = (GPU_results_t*)malloc(sizeof(GPU_results_t));
    
    // allocate memory in GPU
    cudaMalloc(&(tmp_r->nspk), n * s * sizeof(int)); // allocate memory for number of spikes per each sample [n_neurons * n_samples]
    cudaMalloc(&(tmp_r->gs), n * s * conf->max_input_spikes * sizeof(int)); // allocate memory for neurons

    // final structure
    cudaMalloc(&(d_results), sizeof(GPU_results_t)); // allocate memory for neurons
    cudaMemcpy(d_results, tmp_r, sizeof(GPU_results_t), cudaMemcpyHostToDevice); // copy neurons information

    // deallocate temporal struct memory
    free(tmp_r);

    return d_results;
}

extern "C" GPU_results_t* simulate_in_GPU(GPU_SNN_t *gpu_snn, GPU_dataset_t *gpu_dataset, simulation_configuration_t *conf, cuda_info_t *cuda_info, spiking_nn_t *cpu_snn, input_data_t *cpu_dataset){

    double elpt;
    struct timespec start_gpu, end_gpu;


    // simulate
    clock_gettime(CLOCK_MONOTONIC, &start_gpu);

    GPU_results_t *gpu_results = initialize_results_str_in_GPU(conf, cpu_snn->n_neurons, cpu_dataset->n_samples);

    switch (conf->neuron_type){
        // LIF neurons
        case 0:
            simulate_LIF_in_GPU(gpu_snn, gpu_dataset, gpu_results, conf, cuda_info, cpu_snn, cpu_dataset);
            break;
        default:
            simulate_LIF_in_GPU(gpu_snn, gpu_dataset, gpu_results, conf, cuda_info, cpu_snn, cpu_dataset);
            break;
    }

    // sync
    cudaDeviceSynchronize();


    // TODO: refactorize the following
    clock_gettime(CLOCK_MONOTONIC, &end_gpu);
    elpt = (end_gpu.tv_sec - start_gpu.tv_sec) + (end_gpu.tv_nsec - start_gpu.tv_nsec) / 1e9;
    printf(" Elapsed GPU time: %lf\n", elpt);

    GPU_results_t *cpu_results = (GPU_results_t*)malloc(sizeof(GPU_results_t));
    int *nspk = (int*)malloc(cpu_snn->n_neurons * cpu_dataset->n_samples * sizeof(int));
    cudaMemcpy(cpu_results, gpu_results, sizeof(GPU_results_t), cudaMemcpyDeviceToHost);
    cudaMemcpy(nspk, cpu_results->nspk, cpu_snn->n_neurons * cpu_dataset->n_samples * sizeof(int), cudaMemcpyDeviceToHost);

    printf("\n");
    for(int i = 0; i<150; i++){

        printf(" > Sample %d: ", i);
        for(int j = 0; j<cpu_snn->n_neurons; j++){

            printf("%d ", nspk[i * cpu_snn->n_neurons + j]);
        }
        printf("\n");
    }

    return (cpu_results);
}
