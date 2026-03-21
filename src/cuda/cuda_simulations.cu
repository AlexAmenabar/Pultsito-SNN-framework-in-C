#include <cuda.h>
#include <cuda_runtime.h>
#include <stdio.h>
#include <omp.h>
#include<unistd.h>

// public headers
#include "networks/snn.h"
#include "datasets/datasets.h"
#include "simulations/simulations.h"
#include "simulations/results.h"
#include "config/config_loader.h"

#include "cuda/cuda_simulations.cuh"
#include "cuda/cuda_results.cuh"
#include "cuda/cuda_simulations_conf.cuh"

// private cuda headers
#include "cuda/priv_cuda_simulations.cuh"
#include "cuda/priv_cuda_results.cuh"

#define cudaCheckError(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
{
    if (code != cudaSuccess)
    {
        fprintf(stderr, "CUDA Error: %s %s %d\n", cudaGetErrorString(code), file, line);
        if (abort) exit(code);
    }
}

// constant values for helping the simulations
__constant__ size_t learn;
__constant__ size_t batch_size_per_block; // samples processed on each cuda block per neuron
__constant__ size_t blocks_per_batch; // number of cuda blocks for processing all neuron batches
__constant__ size_t dev_batch_size; // batch size proccessed by each device
__constant__ size_t dev_batch_offset; // offset to the sample proccessed by each device
__constant__ size_t batch_size;
__constant__ size_t thrN;



extern "C" void simulate_batches_GPU(GPU_results_t **cpu_results, GPU_SNN_t **gpu_snn, GPU_dataset_t **gpu_dataset, simulation_configuration_t *conf, cuda_info_t *cuda_info){

    size_t n_batches, r_samples, b;
    
    GPU_dataset_t *cpu_dataset = cuda_info->cpu_dataset;

    // compute the number of batches
    n_batches = cpu_dataset->n_samples / conf->batch_size;
    r_samples = cpu_dataset->n_samples % conf->batch_size;
    if(r_samples > 0){

        n_batches += 1; // last batch contains less samples
    }
    
    // simulate
    for(b = 0; b<n_batches; b++){
        simulate_batch_GPU(cpu_results[b], gpu_snn, gpu_dataset, conf, cuda_info, b);
    }
}

extern "C" void simulate_batch_GPU(GPU_results_t *cpu_results, GPU_SNN_t **gpu_snn, GPU_dataset_t **gpu_dataset, simulation_configuration_t *conf, cuda_info_t *cuda_info, size_t bidx){

    size_t times_dev = 0;
    struct timespec start, end;
    double et = 0.0;
    clock_gettime(CLOCK_MONOTONIC, &start);

    // initialize tmp_batch struct in gpu
    tmp_batch_cpu_t **gpu_tmp_batch = allocate_batch_matrix_in_GPU(cuda_info->nDevices, cuda_info->cpu_snn->n_input_neurons, conf->max_input_spikes, cuda_info); // refactorize to function

    // call function depending on whether it is single GPU or multiGPU
    if(cuda_info->nDevices == 1){
        
        simulate_batch_single_GPU(gpu_snn[0], gpu_dataset[0], gpu_tmp_batch[0], conf, cuda_info, bidx);
    }
    // nDevices > 1
    else{

        simulate_batch_multi_GPU(gpu_snn, gpu_dataset, gpu_tmp_batch, conf, cuda_info, bidx);
    }

    // deallocate tmp_batch memory
    deallocate_batch_matrix_in_GPU(gpu_tmp_batch, cuda_info);

    clock_gettime(CLOCK_MONOTONIC, &end);
    et +=(end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // accumulate results (move from GPU to the CPU and then accumulate in only one structure)
    // if only 1 device, copy directly
    if(cuda_info->nDevices == 1){

        cpy_batch_results_struct_GPU2CPU(cpu_results, cuda_info->gpu_results[0], cuda_info, cuda_info->cpu_snn->n_neurons, conf->time_steps);
    }
    // if more than 1 devices, copy and accumulate
    else{

        // accumulate results from several devices into one
        cpy_batch_results_GPU2CPU(cpu_results, cuda_info->gpu_results, cuda_info, cuda_info->cpu_snn->n_neurons, conf->time_steps);

        // get the device that needed more time and store in the first results structure
        double max = 0;
        size_t max_index = 0;
        for(size_t dev = 0; dev < cuda_info->nDevices; dev++){

            if(cuda_info->intermedaite_cpu_results[dev]->t > max){

                max = cuda_info->intermedaite_cpu_results[dev]->t;
                max_index = dev;
            }
        }
        times_dev = max_index;
    }

    // store total execution time
    cpu_results->t = et;

    // store the rest of execution times
    cpu_results->t_in     = cuda_info->intermedaite_cpu_results[times_dev]->t_in;     // time processing input spikes
    cpu_results->t_v      = cuda_info->intermedaite_cpu_results[times_dev]->t_v;      // time processing neurons dynamics
    cpu_results->t_out    = cuda_info->intermedaite_cpu_results[times_dev]->t_out;    // time processing output spikes
    cpu_results->t_learn  = cuda_info->intermedaite_cpu_results[times_dev]->t_learn;  // time learning
    cpu_results->t_reinit = cuda_info->intermedaite_cpu_results[times_dev]->t_reinit; // network reinitialization
    cpu_results->t_load   = cuda_info->intermedaite_cpu_results[times_dev]->t_load;   // loading sample or batch in network
}

extern "C" void simulate_batch_single_GPU(GPU_SNN_t *gpu_snn, GPU_dataset_t *gpu_dataset, tmp_batch_cpu_t *gpu_tmp_batch, simulation_configuration_t *conf, cuda_info_t *cuda_info, size_t bidx){

    // get intermediate results structures
    GPU_results_t *cpu_results = cuda_info->intermedaite_cpu_results[0]; 

    // reinitialize intermediate results structure (CPU)
    reinitialize_batch_results_cpu(cpu_results, conf, cuda_info->cpu_snn->n_neurons, conf->batch_size, 1);

    // run kernels
    run_simulation_batch_GPU(gpu_snn, gpu_dataset, gpu_tmp_batch, conf, cuda_info, bidx, 0);

    // update weights
    update_weights_single_GPU(gpu_snn, conf, cuda_info);
}

extern "C" void update_weights_single_GPU(GPU_SNN_t *gpu_snn, simulation_configuration_t *conf, cuda_info_t *cuda_info){

    struct timespec start_learn, end_learn;
    double et_learn = 0.0;

    clock_gettime(CLOCK_MONOTONIC, &start_learn);

    // get intermediate results structures
    GPU_results_t *cpu_results = cuda_info->intermedaite_cpu_results[0]; 
    GPU_SNN_t *cpu_snn = cuda_info->cpu_snn;

    // define grid and block
    dim3 grid_uw(cuda_info->n_blk_uw_x[0], cuda_info->n_blk_uw_y[0], cuda_info->n_blk_uw_z[0]);
    dim3 block_uw(cuda_info->n_thr_per_blk_uw_x[0], cuda_info->n_thr_per_blk_uw_y[0], cuda_info->n_thr_per_blk_uw_z[0]);

    // update weights     
    if(conf->learn){
            
        update_weights_batch<<<grid_uw, block_uw>>>(gpu_snn, cpu_snn->n_synapses);
        cudaCheckError(cudaPeekAtLastError());  // check launch errors
        cudaCheckError(cudaDeviceSynchronize());  // check runtime errors
        cudaDeviceSynchronize();
    }

    clock_gettime(CLOCK_MONOTONIC, &end_learn);
    et_learn +=(end_learn.tv_sec - start_learn.tv_sec) + (end_learn.tv_nsec - start_learn.tv_nsec) / 1e9;

    // store execution time
    cpu_results->t_learn += et_learn;
}

extern "C" void simulate_batch_multi_GPU(GPU_SNN_t **gpu_snn, GPU_dataset_t **gpu_dataset, tmp_batch_cpu_t **gpu_tmp_batch, simulation_configuration_t *conf, cuda_info_t *cuda_info, size_t bidx){

    size_t dev, nDevices;

    // get intermediate results structures
    GPU_results_t **cpu_results = cuda_info->intermedaite_cpu_results; 
    GPU_SNN_t *cpu_snn = cuda_info->cpu_snn;
    nDevices = cuda_info->nDevices;

    // reinitialize results struct
    for(dev = 0; dev<nDevices; dev++)
        reinitialize_batch_results_cpu(cpu_results[dev], conf, cpu_snn->n_neurons, cuda_info->dev_batch_size[dev], 1);


    // start the simulation using several devices. Each device is managed by a openMP thread
    #pragma omp parallel num_threads(nDevices) private(dev)
    {        
        // set the device dev as the current one for the thread
        dev = (size_t)omp_get_thread_num();
        cudaSetDevice(dev);

        // run kernels
        run_simulation_batch_GPU(gpu_snn[dev], gpu_dataset[dev], gpu_tmp_batch[dev], conf, cuda_info, bidx, dev);

        // wait until device tasks finished
        cudaDeviceSynchronize();
    }
    
    // update weights
    update_weights_multi_GPU(gpu_snn, conf, cuda_info);
}

extern "C" void update_weights_multi_GPU(GPU_SNN_t **gpu_snn, simulation_configuration_t *conf, cuda_info_t *cuda_info){

    size_t dev;

    size_t nDevices = cuda_info->nDevices;
    GPU_SNN_t **tmp_snn = cuda_info->tmp_snn;
    float **dw = cuda_info->dw;
    GPU_SNN_t *cpu_snn = cuda_info->cpu_snn;

    #pragma omp parallel num_threads(nDevices) private(dev)
    {       
        
        struct timespec start_learn, end_learn;
        double et_learn = 0.0;
        clock_gettime(CLOCK_MONOTONIC, &start_learn);


        // set the device dev as the current one for the thread
        dev = (size_t)omp_get_thread_num();
        cudaSetDevice(dev);

        // get results structure to store execution time for each thread
        GPU_results_t *cpu_results = cuda_info->intermedaite_cpu_results[dev];


        // define grid and block
        dim3 grid_uw(cuda_info->n_blk_uw_x[dev], cuda_info->n_blk_uw_y[dev], cuda_info->n_blk_uw_z[dev]);
        dim3 block_uw(cuda_info->n_thr_per_blk_uw_x[dev], cuda_info->n_thr_per_blk_uw_y[dev], cuda_info->n_thr_per_blk_uw_z[dev]);

        // update weights 
        if(conf->learn == 1){

            // accummulate all dw in one network for updating
            acc_weights_batch<<<grid_uw, block_uw>>>(gpu_snn[dev], cpu_snn->n_synapses);
            cudaCheckError(cudaPeekAtLastError());  // check launch errors
            cudaCheckError(cudaDeviceSynchronize());  // check runtime errors

            // wait until device tasks finished
            cudaDeviceSynchronize();

            // cpy weights to CPU
            cudaMemcpy(tmp_snn[dev], gpu_snn[dev], sizeof(GPU_SNN_t), cudaMemcpyDeviceToHost); // cpy dw GPU2CPU
            cudaMemcpy(dw[dev], tmp_snn[dev]->acc_dw, cpu_snn->n_synapses * sizeof(float), cudaMemcpyDeviceToHost); // cpy dw GPU2CPU

            // wait until all threads obtained their weights
            #pragma omp barrier
            
            // sum dw of different devices on device 0
            if(dev == 0){
                for(size_t s = 0; s<cpu_snn->n_synapses; s++){

                    for(size_t d = 0; d<cuda_info->nDevices; d++){
                            
                        if(dev != d){

                            dw[dev][s] += dw[d][s]; // sum weights of other devices
                        }
                    }
                    
                    // compute mean
                    dw[dev][s] = dw[dev][s] / (float)(conf->batch_size);
                }
            }

            #pragma omp barrier
                
            // cpy again to gpu
            cudaMemcpy(tmp_snn[dev]->acc_dw, dw[0], cpu_snn->n_synapses * sizeof(float), cudaMemcpyHostToDevice); // cpy dw CPU2GPU

            // update weights in GPU
            update_weights_dev<<<grid_uw, block_uw>>>(gpu_snn[dev], cpu_snn->n_synapses);
            cudaCheckError(cudaPeekAtLastError());  // check launch errors
            cudaCheckError(cudaDeviceSynchronize());  // check runtime errors
            cudaDeviceSynchronize();
        }

        clock_gettime(CLOCK_MONOTONIC, &end_learn);
        et_learn +=(end_learn.tv_sec - start_learn.tv_sec) + (end_learn.tv_nsec - start_learn.tv_nsec) / 1e9;

        // store execution times
        cpu_results->t_learn += et_learn;
        cpu_results->t += et_learn;
    }
}

extern "C" void run_simulation_batch_GPU(GPU_SNN_t *gpu_snn, GPU_dataset_t *gpu_dataset, tmp_batch_cpu_t *gpu_tmp_batch, simulation_configuration_t *conf, cuda_info_t *cuda_info, size_t bidx, size_t dev){

    size_t N, iN, S, n_networks, T, LT;
    size_t t, s, snn;
    cudaError_t err;
    
    struct timespec start_in, end_in, start_out, end_out, start_v, end_v;
    struct timespec start_learn, end_learn;
    struct timespec start_load, end_load, start_reinit, end_reinit;
    double et = 0.0, et_in = 0.0, et_out = 0.0, et_v = 0.0, et_learn = 0.0, et_load = 0.0, et_reinit = 0.0;
    
    // get results structures
    GPU_results_t *cpu_results = cuda_info->intermedaite_cpu_results[dev];
    GPU_results_t *gpu_results = cuda_info->gpu_results[dev]; 

    // get network and dataset
    GPU_SNN_t *cpu_snn = cuda_info->cpu_snn;
    GPU_dataset_t *cpu_dataset = cuda_info->cpu_dataset;

    // get general variables values
    N = cpu_snn->n_neurons;
    iN = cpu_snn->n_input_neurons;
    S = cpu_snn->n_synapses;
    T = conf->time_steps;
    LT = cpu_snn->LT;

    n_networks = cuda_info->n_networks;
    
    // set dim3 struct to launch kernels

    dim3 grid_neurons(cuda_info->n_blk_neurons_x[dev], cuda_info->n_blk_neurons_y[dev], cuda_info->n_blk_neurons_z[dev]);
    dim3 block_neurons(cuda_info->n_thr_per_blk_neurons_x[dev], cuda_info->n_thr_per_blk_neurons_y[dev], cuda_info->n_thr_per_blk_neurons_z[dev]);

    dim3 grid_synapses(cuda_info->n_blk_synapses_x[dev], cuda_info->n_blk_synapses_y[dev], cuda_info->n_blk_synapses_z[dev]);
    dim3 block_synapses(cuda_info->n_threads_per_blk_synapses_x[dev], cuda_info->n_threads_per_blk_synapses_y[dev], cuda_info->n_threads_per_blk_synapses_z[dev]);

    dim3 grid_in_neurons(cuda_info->n_blk_in_neurons_x[dev], cuda_info->n_blk_in_neurons_y[dev], cuda_info->n_blk_in_neurons_z[dev]);
    dim3 block_in_neurons(cuda_info->n_thr_per_blk_in_neurons_x[dev], cuda_info->n_thr_per_blk_in_neurons_y[dev], cuda_info->n_thr_per_blk_in_neurons_z[dev]);

    dim3 grid_all_neurons(cuda_info->n_blk_all_neurons_x[dev], cuda_info->n_blk_all_neurons_y[dev], cuda_info->n_blk_all_neurons_z[dev]);
    dim3 block_all_neurons(cuda_info->n_thr_per_blk_all_neurons_x[dev], cuda_info->n_thr_per_blk_all_neurons_y[dev], cuda_info->n_thr_per_blk_all_neurons_z[dev]);

    dim3 grid_is(cuda_info->n_blk_is_x[dev], cuda_info->n_blk_is_y[dev], cuda_info->n_blk_is_z[dev]);
    dim3 block_is(cuda_info->n_thr_per_blk_is_x[dev], cuda_info->n_thr_per_blk_is_y[dev], cuda_info->n_thr_per_blk_is_z[dev]);

    // (re)initialize neurons
    clock_gettime(CLOCK_MONOTONIC, &start_reinit);
    cuda_info->initialize_neurons_cuda(gpu_snn, cuda_info, dev);

    // (re)initialize synapses
    initialize_synapses_batch<<<grid_synapses, block_synapses>>>(gpu_snn);
    err = cudaGetLastError();
    cudaCheckError(cudaPeekAtLastError());  // check launch errors
    cudaCheckError(cudaDeviceSynchronize());  // check runtime errors
    cudaDeviceSynchronize();
    

    // load data from the dataset to the intermediate bitmap structure
    initialize_batch_matrix_GPU<<<grid_in_neurons, block_in_neurons>>>(gpu_tmp_batch, gpu_dataset, bidx, iN, conf->max_input_spikes);
    err = cudaGetLastError();
    cudaCheckError(cudaPeekAtLastError());  // check launch errors
    cudaCheckError(cudaDeviceSynchronize());  // check runtime errors
    cudaDeviceSynchronize();


    // reinitialize data is results structure
    reinitialize_results_number_of_spikes_GPU<<<grid_all_neurons, block_all_neurons>>>(gpu_results, N);
    err = cudaGetLastError();
    cudaCheckError(cudaPeekAtLastError());  // check launch errors
    cudaCheckError(cudaDeviceSynchronize());  // check runtime errors
    cudaDeviceSynchronize();

    clock_gettime(CLOCK_MONOTONIC, &end_reinit);
    et_reinit +=(end_reinit.tv_sec - start_reinit.tv_sec) + (end_reinit.tv_nsec - start_reinit.tv_nsec) / 1e9;


    // loop over time steps
    for(t = 0; t<T; t++){

        // compute time step index in the spike matrix
        size_t lt = t % cpu_snn->LT;

        // load sample time step in the SNN spike matrix
        clock_gettime(CLOCK_MONOTONIC, &start_load);

        load_batch_time_step_in_SNN_GPU<<<grid_all_neurons, block_all_neurons>>>(gpu_snn, gpu_tmp_batch, iN, N, lt, t, bidx);
        err = cudaGetLastError();
        cudaCheckError(cudaPeekAtLastError()); 
        cudaCheckError(cudaDeviceSynchronize());  
        cudaDeviceSynchronize();
        
        clock_gettime(CLOCK_MONOTONIC, &end_load);
        et_load +=(end_load.tv_sec - start_load.tv_sec) + (end_load.tv_nsec - start_load.tv_nsec) / 1e9;


        // compute input currents
        clock_gettime(CLOCK_MONOTONIC, &start_in);

        cuda_info->process_input_current_cuda(gpu_snn, iN, N, lt, t, cuda_info, dev);

        clock_gettime(CLOCK_MONOTONIC, &start_learn);

        // update whether neurons fired for learning
        if(conf->learn){

            compute_pre_fired_batch<<<grid_synapses, block_synapses>>>(gpu_snn, N, iN, S, lt, t);
            cudaCheckError(cudaPeekAtLastError()); 
            cudaCheckError(cudaDeviceSynchronize());  
            cudaDeviceSynchronize();
                
        }

        clock_gettime(CLOCK_MONOTONIC, &end_learn);
        et_learn +=(end_learn.tv_sec - start_learn.tv_sec) + (end_learn.tv_nsec - start_learn.tv_nsec) / 1e9;


        clock_gettime(CLOCK_MONOTONIC, &end_in);
        et_in +=(end_in.tv_sec - start_in.tv_sec) + (end_in.tv_nsec - start_in.tv_nsec) / 1e9;


        // compute neuron dynamics
        clock_gettime(CLOCK_MONOTONIC, &start_v);

        cuda_info->process_dynamics_cuda(gpu_snn, N, cuda_info, dev);

        clock_gettime(CLOCK_MONOTONIC, &end_v);
        et_v +=(end_v.tv_sec - start_v.tv_sec) + (end_v.tv_nsec - start_v.tv_nsec) / 1e9;


        // compute firing
        clock_gettime(CLOCK_MONOTONIC, &start_out);

        cuda_info->process_firing_cuda(gpu_snn, gpu_results, iN, N, lt, cuda_info, dev);

        clock_gettime(CLOCK_MONOTONIC, &end_out);
        et_out +=(end_out.tv_sec - start_out.tv_sec) + (end_out.tv_nsec - start_out.tv_nsec) / 1e9;


        // compute learning step if training
        clock_gettime(CLOCK_MONOTONIC, &start_learn);

        if(conf->learn){

            // update traces
            compute_post_traces_batch<<<grid_neurons, block_neurons>>>(gpu_snn, N);
            cudaCheckError(cudaPeekAtLastError()); 
            cudaCheckError(cudaDeviceSynchronize());  
            cudaDeviceSynchronize();

            // compute STDP
            trace_based_STDP_batch<<<grid_synapses, block_synapses>>>(gpu_snn, S);
            cudaCheckError(cudaPeekAtLastError()); 
            cudaCheckError(cudaDeviceSynchronize());  
            cudaDeviceSynchronize();

            // reinit flags that indicates if neurons fired
            reinit_post_fired_batch<<<grid_neurons, block_neurons>>>(gpu_snn, N);
            cudaCheckError(cudaPeekAtLastError()); 
            cudaCheckError(cudaDeviceSynchronize());  
            cudaDeviceSynchronize();

        }
        clock_gettime(CLOCK_MONOTONIC, &end_learn);
        et_learn +=(end_learn.tv_sec - start_learn.tv_sec) + (end_learn.tv_nsec - start_learn.tv_nsec) / 1e9;
    }


    // store execution times in the CPU structure
    cpu_results->t_in = et_in;
    cpu_results->t_v = et_v;
    cpu_results->t_out = et_out;
    cpu_results->t_reinit = et_reinit;
    cpu_results->t_load = et_load;
    cpu_results->t_learn = et_learn;

    // store execution times of the device
    cpu_results->t = et_in + et_v + et_out + et_reinit + et_load + et_learn;
}


extern "C" tmp_batch_cpu_t** allocate_batch_matrix_in_GPU(size_t nDevices, size_t iN, size_t T, cuda_info_t *cuda_info){

    size_t dev, dev_batch_size;
    tmp_batch_cpu_t *tmp_batch, **d_tmp_batch;
    cudaError_t err;


    // allocate memory
    tmp_batch = (tmp_batch_cpu_t*)calloc(1, sizeof(tmp_batch_cpu_t));
    d_tmp_batch = (tmp_batch_cpu_t**)calloc(nDevices, sizeof(tmp_batch_cpu_t*));


    // loop over devices
    for(dev = 0; dev < nDevices; dev++){
    
        // set dev as active one
        cudaSetDevice(dev);

        // batch size per dev
        dev_batch_size = cuda_info->dev_batch_size[dev];

        printf(" In allocate batch matrix: iN = %zu, dev_batch = %zu, T = %zu\n", iN, dev_batch_size, T);
        
        // allocate memory to store the number of spikes
        err = cudaMalloc(&(tmp_batch->spikes), iN * dev_batch_size * T * sizeof(char));
        if (err != cudaSuccess) 
            printf("neuron_input_synapses_offset allocation failed: %s\n", cudaGetErrorString(err));

        // final structure
        err = cudaMalloc(&(d_tmp_batch[dev]), sizeof(tmp_batch_cpu_t)); // allocate memory for neurons
        if (err != cudaSuccess) 
            printf("neuron_input_synapses_offset allocation failed: %s\n", cudaGetErrorString(err));

        // copy data to GPU
        err = cudaMemcpy(d_tmp_batch[dev], tmp_batch, sizeof(tmp_batch_cpu_t), cudaMemcpyHostToDevice); // copy neurons information
        if (err != cudaSuccess) 
            printf("neuron_input_synapses_offset allocation failed: %s\n", cudaGetErrorString(err));
    }
    
    // deallocate temporal struct memory
    free(tmp_batch);

    return d_tmp_batch;
}

extern "C" void deallocate_batch_matrix_in_GPU(tmp_batch_cpu_t **batch_matrix, cuda_info_t *cuda_info){

    size_t dev;
    cudaError_t err;
    tmp_batch_cpu_t *tmp_batch = (tmp_batch_cpu_t*)calloc(1, sizeof(tmp_batch_cpu_t));

    // loop over devices
    for(dev = 0; dev < cuda_info->nDevices; dev ++){
    
        // set dev as active one
        cudaSetDevice(dev);

        // copy struct to the CPU to deallocate internal pointers
        cudaMemcpy(tmp_batch, batch_matrix[dev], sizeof(tmp_batch_cpu_t), cudaMemcpyDeviceToHost);

        // deallocate internal pointers
        err = cudaFree(tmp_batch->spikes);
        if (err != cudaSuccess) 
            printf("neuron_input_synapses_offset allocation failed: %s\n", cudaGetErrorString(err));

        // deallocate structure
        err = cudaFree(batch_matrix[dev]);
        if (err != cudaSuccess) 
            printf("neuron_input_synapses_offset allocation failed: %s\n", cudaGetErrorString(err));
    }
    
    // deallocate temporal struct memory
    free(tmp_batch);
}


/* Kernels */

// initialize synapses before starting the simulation
__global__ void initialize_synapses_batch(GPU_SNN_t *snn){

    // get thread id
    size_t threadId = 
        (size_t)((size_t)blockIdx.x  + (size_t)gridDim.x  * (size_t)blockIdx.y + (size_t)gridDim.x * (size_t)gridDim.y * (size_t)blockIdx.z) *
        (size_t)((size_t)blockDim.x * (size_t)blockDim.y * (size_t)blockDim.z) +
        (size_t)((size_t)threadIdx.x + (size_t)blockDim.x * (size_t)threadIdx.y + (size_t)blockDim.x * (size_t)blockDim.y * (size_t)threadIdx.z);


    //printf(" thread = %llu / %llu\n", threadId, snn->n_synapses * batch_size);

    if(threadId < snn->n_synapses * dev_batch_size){

        size_t synapse_index = threadId / dev_batch_size;
        snn->w[threadId] = snn->init_w[synapse_index];
        snn->dw[threadId] = 0.0;

        if(learn){

            snn->pre_fired[threadId] = 0;
            snn->pre_trace[threadId] = 0.0;
        }
    }    
}

__global__ void initialize_batch_matrix_GPU(tmp_batch_cpu_t *tmp_batch, GPU_dataset_t *dataset, size_t bidx, size_t iN, size_t T){

    // get thread id: iN * batch_size * T
    size_t threadId = 
        (blockIdx.x  + gridDim.x  * blockIdx.y + gridDim.x * gridDim.y * blockIdx.z) *
        (blockDim.x * blockDim.y * blockDim.z) +
        (threadIdx.x + blockDim.x * threadIdx.y + blockDim.x * blockDim.y * threadIdx.z);

    if(threadId < iN * dev_batch_size){

        size_t fsidx = bidx * batch_size + dev_batch_offset; // first sample to be processed
        size_t b = threadId % dev_batch_size; // batch element
        size_t sidx = fsidx + b; // sample id = first sample in the batch + batch element

        size_t n_features = dataset->n_features;

        if(sidx < dataset->n_samples){
            
            size_t fidx = threadId / dev_batch_size; // neuron index is the feature index
            
            // reinitialize batch_matrix
            for(size_t j = 0; j<T; j++){
                
                tmp_batch->spikes[(iN * dev_batch_size * j) + (dev_batch_size * fidx) + b] = 0; 
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
                tmp_batch->spikes[(iN * dev_batch_size * t) + (dev_batch_size * fidx) + b] = 1; 
            }
        }
        
    }
}

__global__ void load_batch_time_step_in_SNN_GPU(GPU_SNN_t *gpu_snn, tmp_batch_cpu_t *tmp_batch, size_t iN, size_t N, size_t lt, size_t t, size_t bidx){

    // get thread id: iN * batch_size * T
    size_t threadId = 
        (blockIdx.x  + gridDim.x  * blockIdx.y + gridDim.x * gridDim.y * blockIdx.z) *
        (blockDim.x * blockDim.y * blockDim.z) +
        (threadIdx.x + blockDim.x * threadIdx.y + blockDim.x * blockDim.y * threadIdx.z);


    if(threadId < (iN + N) * dev_batch_size){

        size_t i = threadId / dev_batch_size; // neuron index
        size_t b = threadId % dev_batch_size; // batch position

        size_t input = i < iN ? 1 : 0; //  

        // compute matrix indexes
        size_t spk_mtr_idx = (iN + N) * dev_batch_size * lt + (dev_batch_size * i) + b;

        if(input){

            size_t batch_mtr_idx = (iN * dev_batch_size * t) + (dev_batch_size * i) + b;
            gpu_snn->spk_matrix[spk_mtr_idx] = tmp_batch->spikes[batch_mtr_idx]; 
        }
        else{

            gpu_snn->spk_matrix[spk_mtr_idx] = 0;
        }
    }
}

__global__ void compute_pre_fired_batch(GPU_SNN_t *gpu_snn, size_t N, size_t iN, size_t S, size_t t, size_t gt){

    // get thread id: iN * batch_size * T
    size_t threadId = 
        (blockIdx.x  + gridDim.x  * blockIdx.y + gridDim.x * gridDim.y * blockIdx.z) *
        (blockDim.x * blockDim.y * blockDim.z) +
        (threadIdx.x + blockDim.x * threadIdx.y + blockDim.x * blockDim.y * threadIdx.z);


    if(threadId < S * dev_batch_size){

        // variables declaration
        size_t synapse_index, g_synapse_index, neuron_index, g_neuron_index, b;
        synapse_index = threadId / dev_batch_size;
        g_synapse_index = synapse_index * dev_batch_size;
        b = threadId % dev_batch_size;

        neuron_index = gpu_snn->pre_neuron_index[synapse_index];
        g_neuron_index = neuron_index * dev_batch_size;

        int delay, spk_time; 
        char spk;  

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
            size_t index = ((iN + N) * dev_batch_size * (size_t)spk_time) + ((g_neuron_index));

            // get spike
            spk = gpu_snn->spk_matrix[index + b]; 

            // store if the presynaptic neuron fired for TR-STDP
            //if((int)spk == 1){

            gpu_snn->pre_fired[g_synapse_index + b] = (char)spk;
            //}
        }
    }
}


__global__ void compute_post_traces_batch(GPU_SNN_t *gpu_snn, size_t N){

    // get thread id: iN * batch_size * T
    size_t threadId = 
        (blockIdx.x  + gridDim.x  * blockIdx.y + gridDim.x * gridDim.y * blockIdx.z) *
        (blockDim.x * blockDim.y * blockDim.z) +
        (threadIdx.x + blockDim.x * threadIdx.y + blockDim.x * blockDim.y * threadIdx.z);

    if(threadId < N * dev_batch_size){
    
        size_t i = threadId / dev_batch_size; // neuron index
        size_t b = threadId % dev_batch_size; // batch position
        size_t g_neuron_index = i * dev_batch_size; // for copied variables 
        float decay = 0.6f;

        // update post trace (depends only on the neuron)
        gpu_snn->post_trace[g_neuron_index + b] = gpu_snn->post_fired[g_neuron_index + b] ? 1.0f : gpu_snn->post_trace[g_neuron_index + b] * decay;
    }
}

__global__ void reinit_post_fired_batch(GPU_SNN_t *gpu_snn, size_t N){

    // get thread id: iN * batch_size * T
    size_t threadId = 
        (blockIdx.x  + gridDim.x  * blockIdx.y + gridDim.x * gridDim.y * blockIdx.z) *
        (blockDim.x * blockDim.y * blockDim.z) +
        (threadIdx.x + blockDim.x * threadIdx.y + blockDim.x * blockDim.y * threadIdx.z);

    if(threadId < N * dev_batch_size){
    
        size_t i = threadId / dev_batch_size; // neuron index
        size_t b = threadId % dev_batch_size; // batch position
        size_t g_neuron_index = i * dev_batch_size; // for copied variables 

        // update post trace (depends only on the neuron)
        gpu_snn->post_fired[g_neuron_index + b] = 0;
    }
}


__global__ void trace_based_STDP_batch(GPU_SNN_t *gpu_snn, size_t S){

    
    // get thread id: iN * batch_size * T
    size_t threadId = 
        (size_t)((size_t)blockIdx.x  + (size_t)gridDim.x  * (size_t)blockIdx.y + (size_t)gridDim.x * (size_t)gridDim.y * (size_t)blockIdx.z) *
        (size_t)((size_t)blockDim.x * (size_t)blockDim.y * (size_t)blockDim.z) +
        (size_t)((size_t)threadIdx.x + (size_t)blockDim.x * (size_t)threadIdx.y + (size_t)blockDim.x * (size_t)blockDim.y * (size_t)threadIdx.z);

    if(threadId < S * dev_batch_size){

        size_t synapse_index = threadId / dev_batch_size; // neuron index
        size_t b = threadId % dev_batch_size; // batch position
        size_t g_synapse_index = synapse_index * dev_batch_size; // for copied variables
        
        size_t neuron_index, g_neuron_index;
        neuron_index = gpu_snn->post_neuron_index[synapse_index] - gpu_snn->n_input_neurons; // - virtual neurons
        g_neuron_index = neuron_index * dev_batch_size;

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


__global__ void update_weights_batch(GPU_SNN_t *gpu_snn, size_t S){

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
        
        // update init_w
        gpu_snn->init_w[synapse_index] += dw; // update initial w for incoming batches

        // update initial w and w for all copies
        for(size_t b = 0; b<batch_size; b++){

            gpu_snn->dw[g_synapse_index + b] = 0.0; // reinitialize dw
            gpu_snn->w[g_synapse_index + b] = gpu_snn->init_w[synapse_index]; // update w
        }    
    }
}

__global__ void acc_weights_batch(GPU_SNN_t *gpu_snn, size_t S){

    // get thread id: iN * batch_size * T
    size_t threadId = 
        (blockIdx.x  + gridDim.x  * blockIdx.y + gridDim.x * gridDim.y * blockIdx.z) *
        (blockDim.x * blockDim.y * blockDim.z) +
        (threadIdx.x + blockDim.x * threadIdx.y + blockDim.x * blockDim.y * threadIdx.z);


    if(threadId < S){

        size_t synapse_index = threadId; // synapse index
        size_t g_synapse_index = synapse_index * dev_batch_size; // synapse index
    
        float dw = 0.0;

        // sum dw of the batch for each synapse
        for(size_t b = 0; b<dev_batch_size; b++){

            dw += gpu_snn->dw[g_synapse_index + b];
        }   

        // store accumulated dw
        gpu_snn->acc_dw[synapse_index] = dw;
    }
}

__global__ void update_weights_dev(GPU_SNN_t *gpu_snn, size_t S){

    // get thread id: iN * batch_size * T
    size_t threadId = 
        (blockIdx.x  + gridDim.x  * blockIdx.y + gridDim.x * gridDim.y * blockIdx.z) *
        (blockDim.x * blockDim.y * blockDim.z) +
        (threadIdx.x + blockDim.x * threadIdx.y + blockDim.x * blockDim.y * threadIdx.z);


    if(threadId < S){

        size_t synapse_index = threadId; // synapse index
        size_t g_synapse_index = synapse_index * dev_batch_size; // synapse index
    
        float dw = gpu_snn->acc_dw[synapse_index];
        gpu_snn->init_w[synapse_index] += dw;

        // sum dw of the batch for each synapse
        for(size_t b = 0; b<dev_batch_size; b++){

            gpu_snn->dw[g_synapse_index + b] = 0.0;
            gpu_snn->w[g_synapse_index + b] = gpu_snn->init_w[synapse_index];
        }   
    }
}