#include <cuda.h>
#include <cuda_runtime.h>
#include <stdio.h>
#include <omp.h>


#include "neuron_models/GPU_lif_neuron.cuh"
#include "cuda/GPU_simulations.cuh"
#include "cuda/cuda_simulations_conf.h"

#include "networks/snn.h"
#include "datasets/datasets.h"
#include "simulations/simulations.h"
#include "simulations/results.h"
#include "config/config_loader.h"


#define cudaCheckError(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
{
    if (code != cudaSuccess)
    {
        fprintf(stderr, "CUDA Error: %s %s %d\n", cudaGetErrorString(code), file, line);
        if (abort) exit(code);
    }
}

__constant__ size_t learn;
__constant__ size_t batch_size_per_block; // samples processed on each cuda block per neuron
__constant__ size_t blocks_per_batch; // number of cuda blocks for processing all neuron batches
__constant__ size_t dev_batch_size; // batch size proccessed by each device
__constant__ size_t dev_batch_offset; // offset to the sample proccessed by each device
__constant__ size_t batch_size;
__constant__ size_t thrN;


extern "C" void simulate_batches_LIF_GPU(GPU_SNN_t **gpu_snn, GPU_dataset_t **gpu_dataset, GPU_results_t **gpu_results, tmp_batch_cpu_t **gpu_tmp_batch, simulation_configuration_t *conf, cuda_info_t *cuda_info, GPU_SNN_t *cpu_snn, GPU_dataset_t *cpu_dataset){
    
    size_t n_batches, r_samples, b, i, dev;
    cudaError_t err;
    size_t nDevices = cuda_info->nDevices;
    

    // compute the number of batches
    n_batches = cpu_dataset->n_samples / conf->batch_size;
    r_samples = cpu_dataset->n_samples % conf->batch_size;
    // there are remaining sampels
    if(r_samples > 0){

        n_batches += 1; // last batch contains less samples
    }

    // copy constants to all gpu devices
    for(dev = 0; dev<nDevices; dev++){
        
        cudaSetDevice(dev);
        
        err = cudaMemcpyToSymbol(learn, &conf->learn, sizeof(size_t));
        err = cudaMemcpyToSymbol(batch_size, &conf->batch_size, sizeof(size_t));
        err = cudaMemcpyToSymbol(batch_size_per_block, &cuda_info->batch_size_per_block, sizeof(size_t));
        err = cudaMemcpyToSymbol(blocks_per_batch, &cuda_info->blocks_per_batch[dev], sizeof(size_t));
        err = cudaMemcpyToSymbol(dev_batch_size, &cuda_info->dev_batch_size[dev], sizeof(size_t)); // how much samples simulates each device
        err = cudaMemcpyToSymbol(dev_batch_offset, &cuda_info->dev_batch_offset[dev], sizeof(size_t)); // offset to the first sample simulated by the device in the batch
        err = cudaMemcpyToSymbol(thrN, &conf->thrN, sizeof(size_t)); // offset to the first sample simulated by the device in the batch
    }


    // allocate memory for results structs (a struct per device)
    GPU_results_t **cpu_results = (GPU_results_t**)calloc(nDevices, sizeof(GPU_results_t*));
    GPU_results_t **cpu_batch_results = (GPU_results_t**)calloc(nDevices, sizeof(GPU_results_t*));
        
    // copy constant(s) to all gpu devices: batch size per device and batch offset
    for(dev = 0; dev<nDevices; dev++){
        
        // initialize batch results for device
        cpu_results[dev] = initialize_batch_results_cpu(conf, cpu_snn->n_neurons, cuda_info->dev_batch_size[dev], 1, 0); 
        cpu_batch_results[dev] = initialize_batch_results_cpu(conf, cpu_snn->n_neurons, cuda_info->dev_batch_size[dev], 1, 0);
    }
    

    // if it is executed in only one device
    if(nDevices == 1){
        
        cudaSetDevice(0);

        // loop over batches
        for(b = 0; b<n_batches; b++){

            // simulate batch
            simulate_LIF_batch_single_GPU(gpu_snn[0], gpu_dataset[0], gpu_results[0], gpu_tmp_batch[0], conf, cuda_info, b, cpu_snn, cpu_dataset, cpu_batch_results[0]);

            // accumulate batch execution times in general results struct
            acc_batch_execution_times(cpu_results[0], cpu_batch_results[0]);            
        }
    }

    // TODO: multi GPU must be corrected
    else{
        
        // simulate one batch and print
        /*simulate_LIF_batch_multi_GPU(gpu_snn, gpu_dataset, gpu_results, gpu_tmp_batch, conf, cuda_info, 0, cpu_snn, cpu_dataset, cpu_batch_results);
        cpy_batch_results_GPU2CPU(cpu_results, gpu_results, cuda_info, cpu_snn->n_neurons, 1);
        for(dev = 0; dev<cuda_info->nDevices; dev++){
            
            for(size_t batch = 0; batch<cuda_info->dev_batch_size[dev]; batch++){

                for(size_t i = 0; i<cpu_snn->n_neurons; i++){

                    printf("%d ", cpu_results[dev]->n_spks[i * cuda_info->dev_batch_size[dev] + batch]);
                }
                printf("\n");
            }
        }*/


        // loop over batches
        for(b = 0; b<n_batches; b++){

            // simulate batch
            simulate_LIF_batch_multi_GPU(gpu_snn, gpu_dataset, gpu_results, gpu_tmp_batch, conf, cuda_info, b, cpu_snn, cpu_dataset, cpu_batch_results);

            // accumulate batch execution times in general results struct
            for(dev = 0; dev<nDevices; dev++)
                acc_batch_execution_times(cpu_results[dev], cpu_batch_results[dev]);            
        }


        /*simulate_LIF_batch_multi_GPU(gpu_snn, gpu_dataset, gpu_results, gpu_tmp_batch, conf, cuda_info, 0, cpu_snn, cpu_dataset, cpu_batch_results);
        cpy_batch_results_GPU2CPU(cpu_results, gpu_results, cuda_info, cpu_snn->n_neurons, 1);
        for(dev = 0; dev<cuda_info->nDevices; dev++){
            
            for(size_t batch = 0; batch<cuda_info->dev_batch_size[dev]; batch++){

                for(size_t i = 0; i<cpu_snn->n_neurons; i++){

                    printf("%d ", cpu_results[dev]->n_spks[i * cuda_info->dev_batch_size[dev] + batch]);
                }
                printf("\n");
            }
        }*/


        for(dev = 1; dev<nDevices; dev++){

            //cpu_results[0]->t += cpu_results[dev]->t;
            cpu_results[0]->t_in += cpu_results[dev]->t_in;
            cpu_results[0]->t_v += cpu_results[dev]->t_v;
            cpu_results[0]->t_out += cpu_results[dev]->t_out;
            cpu_results[0]->t_reinit += cpu_results[dev]->t_reinit;
            cpu_results[0]->t_load += cpu_results[dev]->t_load;
            cpu_results[0]->t_learn += cpu_results[dev]->t_learn;
        }

        // compute the mean
        cpu_results[0]->t_in     /= nDevices;
        cpu_results[0]->t_v      /= nDevices;
        cpu_results[0]->t_out    /= nDevices;
        cpu_results[0]->t_reinit /= nDevices;
        cpu_results[0]->t_load   /= nDevices;
        cpu_results[0]->t_learn  /= nDevices;
    }


    // print information in results struct
    printf("   > Total execution time: %.3lf | Mean per batch %.3lf | Mean per time step %.3lf\n", 
                cpu_results[0]->t, cpu_results[0]->t / (float)n_batches, cpu_results[0]->t / (float)n_batches / (float)conf->time_steps);
    printf("   > Total in time:        %.3lf | Mean per batch %.3lf | Mean per time step %.3lf\n", 
                cpu_results[0]->t_in, cpu_results[0]->t_in / (float)n_batches, cpu_results[0]->t_in / (float)n_batches / (float)conf->time_steps);
    printf("   > Total v time:         %.3lf | Mean per batch %.3lf | Mean per time step %.3lf\n", 
                cpu_results[0]->t_v, cpu_results[0]->t_v / (float)n_batches, cpu_results[0]->t_v / (float)n_batches / (float)conf->time_steps);
    printf("   > Total out time:       %.3lf | Mean per batch %.3lf | Mean per time step %.3lf\n", 
                cpu_results[0]->t_out, cpu_results[0]->t_out / (float)n_batches, cpu_results[0]->t_out / (float)n_batches / (float)conf->time_steps);
    printf("   > Total reinit time:    %.3lf | Mean per batch %.3lf | Mean per time step %.3lf\n", 
                cpu_results[0]->t_reinit, cpu_results[0]->t_reinit / (float)n_batches, cpu_results[0]->t_reinit / (float)n_batches / (float)conf->time_steps);
    printf("   > Total load time:      %.3lf | Mean per batch %.3lf | Mean per time step %.3lf\n", 
                cpu_results[0]->t_load, cpu_results[0]->t_load / (float)n_batches, cpu_results[0]->t_load / (float)n_batches / (float)conf->time_steps);
    printf("   > Total learn time:     %.3lf | Mean per batch %.3lf | Mean per time step %.3lf\n", 
                cpu_results[0]->t_learn, cpu_results[0]->t_learn / (float)n_batches, cpu_results[0]->t_learn / (float)n_batches / (float)conf->time_steps);
}

extern "C" void simulate_LIF_batch_multi_GPU(GPU_SNN_t **gpu_snn, GPU_dataset_t **gpu_dataset, GPU_results_t **gpu_results, tmp_batch_cpu_t **gpu_tmp_batch, simulation_configuration_t *conf, cuda_info_t *cuda_info, size_t bidx, GPU_SNN_t *cpu_snn, GPU_dataset_t *cpu_dataset, GPU_results_t **cpu_results){

    struct timespec start, end, start_learn, end_learn, start_gpu, end_gpu;
    double et = 0.0, et_learn = 0.0, et_gpu = 0.0;
    size_t dev, nDevices;
    GPU_SNN_t **tmp_snn;
    float **dw;

    tmp_snn = cuda_info->tmp_snn;
    dw = cuda_info->dw;
    nDevices = cuda_info->nDevices;


    // reinitialize results struct
    for(dev = 0; dev<nDevices; dev++)
        reinitialize_batch_results_cpu(cpu_results[dev], conf, cpu_snn->n_neurons, cuda_info->dev_batch_size[dev], 1);


    clock_gettime(CLOCK_MONOTONIC, &start);

    // start the simulation using several devices. Each device is managed by a openMP thread
    #pragma omp parallel num_threads(nDevices) private(dev, start_learn, end_learn, et_learn, start_gpu, end_gpu, et_gpu)
    {        
        // set the device dev as the current one for the thread
        dev = (size_t)omp_get_thread_num();
        cudaSetDevice(dev);

        // define grid and block
        dim3 grid_uw(cuda_info->n_blk_uw_x[dev], cuda_info->n_blk_uw_y[dev], cuda_info->n_blk_uw_z[dev]);
        dim3 block_uw(cuda_info->n_thr_per_blk_uw_x[dev], cuda_info->n_thr_per_blk_uw_y[dev], cuda_info->n_thr_per_blk_uw_z[dev]);

        // simulate batch using several devices
        simulate_LIF_batch_GPU(gpu_snn[dev], gpu_dataset[dev], gpu_results[dev], gpu_tmp_batch[dev], conf, cuda_info, bidx, cpu_snn, cpu_dataset, cpu_results[dev], dev);

        // wait until device tasks finished
        cudaDeviceSynchronize();
        
        // if training
        clock_gettime(CLOCK_MONOTONIC, &start_learn);
        if(conf->learn == 1){

            // acc all dw of the batch on the first
            acc_weights_batch<<<grid_uw, block_uw>>>(gpu_snn[dev], cpu_snn->n_synapses);
            cudaCheckError(cudaPeekAtLastError());  // check launch errors
            cudaCheckError(cudaDeviceSynchronize());  // check runtime errors
            cudaDeviceSynchronize();

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
            //cudaMemcpy(gpu_snn[dev], tmp_snn[dev], sizeof(GPU_SNN_t), cudaMemcpyHostToDevice); // cpy dw CPU2GPU

            // update weights in GPU
            update_weights_dev<<<grid_uw, block_uw>>>(gpu_snn[dev], cpu_snn->n_synapses);
            cudaCheckError(cudaPeekAtLastError());  // check launch errors
            cudaCheckError(cudaDeviceSynchronize());  // check runtime errors
            cudaDeviceSynchronize();
        }

        clock_gettime(CLOCK_MONOTONIC, &end_learn);
        et_learn +=(end_learn.tv_sec - start_learn.tv_sec) + (end_learn.tv_nsec - start_learn.tv_nsec) / 1e9;


        // store execution times
        cpu_results[dev]->t_learn += et_learn;
    }

    clock_gettime(CLOCK_MONOTONIC, &end);

    // global batch execution time
    et +=(end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    clock_gettime(CLOCK_MONOTONIC, &end);
    cpu_results[0]->t = et;
}

extern "C" void simulate_LIF_batch_single_GPU(GPU_SNN_t *gpu_snn, GPU_dataset_t *gpu_dataset, GPU_results_t *gpu_results, tmp_batch_cpu_t *gpu_tmp_batch, simulation_configuration_t *conf, cuda_info_t *cuda_info, size_t bidx, GPU_SNN_t *cpu_snn, GPU_dataset_t *cpu_dataset, GPU_results_t *cpu_results){

    struct timespec start, end, start_learn, end_learn;
    double et = 0.0, et_learn = 0.0;

    // define grid and block
    dim3 grid_uw(cuda_info->n_blk_uw_x[0], cuda_info->n_blk_uw_y[0], cuda_info->n_blk_uw_z[0]);
    dim3 block_uw(cuda_info->n_thr_per_blk_uw_x[0], cuda_info->n_thr_per_blk_uw_y[0], cuda_info->n_thr_per_blk_uw_z[0]);


    // reinitialize results struct
    reinitialize_batch_results_cpu(cpu_results, conf, cpu_snn->n_neurons, conf->batch_size, 1);


    clock_gettime(CLOCK_MONOTONIC, &start);

    // simulate batch
    simulate_LIF_batch_GPU(gpu_snn, gpu_dataset, gpu_results, gpu_tmp_batch, conf, cuda_info, bidx, cpu_snn, cpu_dataset, cpu_results, 0);

    // update weights 
    clock_gettime(CLOCK_MONOTONIC, &start_learn);
    if(conf->learn){
            
        update_weights_batch<<<grid_uw, block_uw>>>(gpu_snn, cpu_snn->n_synapses);
        cudaCheckError(cudaPeekAtLastError());  // check launch errors
        cudaCheckError(cudaDeviceSynchronize());  // check runtime errors
        cudaDeviceSynchronize();
    }
    clock_gettime(CLOCK_MONOTONIC, &end_learn);
    et_learn +=(end_learn.tv_sec - start_learn.tv_sec) + (end_learn.tv_nsec - start_learn.tv_nsec) / 1e9;


    clock_gettime(CLOCK_MONOTONIC, &end);
    et +=(end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // store execution times
    cpu_results->t = et;
    cpu_results->t_learn += et_learn;
}




/** GPU KERNELS */

__global__ void initialize_neurons_batch(GPU_SNN_t *snn){

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

extern "C" void simulate_LIF_batch_GPU(GPU_SNN_t *gpu_snn, GPU_dataset_t *gpu_dataset, GPU_results_t *gpu_results, tmp_batch_cpu_t *gpu_tmp_batch, simulation_configuration_t *conf, cuda_info_t *cuda_info, size_t bidx, GPU_SNN_t *cpu_snn, GPU_dataset_t *cpu_dataset, GPU_results_t *cpu_results, size_t dev){

    
    size_t N, iN, S, n_networks, T, LT;
    size_t t, s, snn;
    cudaError_t err;
    
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

    initialize_neurons_batch<<<grid_neurons, block_neurons>>>(gpu_snn);
    err = cudaGetLastError();
    cudaCheckError(cudaPeekAtLastError());  // check launch errors
    cudaCheckError(cudaDeviceSynchronize());  // check runtime errors
    cudaDeviceSynchronize();

    // (re)initialize synapses
    initialize_synapses_batch<<<grid_synapses, block_synapses>>>(gpu_snn);
    err = cudaGetLastError();
    cudaCheckError(cudaPeekAtLastError());  // check launch errors
    cudaCheckError(cudaDeviceSynchronize());  // check runtime errors
    cudaDeviceSynchronize();
    
    // load batch in tmp_batch
    initialize_batch_matrix_GPU<<<grid_in_neurons, block_in_neurons>>>(gpu_tmp_batch, gpu_dataset, bidx, iN, conf->max_input_spikes);
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

        load_batch_time_step_in_SNN_GPU<<<grid_all_neurons, block_all_neurons>>>(gpu_snn, gpu_tmp_batch, iN, N, lt, t, bidx);
        err = cudaGetLastError();
        cudaCheckError(cudaPeekAtLastError());  // check launch errors
        cudaCheckError(cudaDeviceSynchronize());  // check runtime errors
        cudaDeviceSynchronize();
        
        clock_gettime(CLOCK_MONOTONIC, &end_load);
        et_load +=(end_load.tv_sec - start_load.tv_sec) + (end_load.tv_nsec - start_load.tv_nsec) / 1e9;

        // compute input step
        clock_gettime(CLOCK_MONOTONIC, &start_in);

        //process_input_currect_batch<<<grid_neurons, block_neurons>>>(gpu_snn, iN, N, batch_size, lt, t, (size_t)conf->learn);

        for(size_t l = 0; l<cuda_info->n_neuron_launchs; l++){

            //process_input_currect_batch<<<grid_is, block_is, conf->thrN * cuda_info->batch_size_per_block * sizeof(float)>>>(gpu_snn, iN, N, lt, t, cuda_info->n_neurons_per_launch, l);
            process_input_currect_batch_new<<<grid_is, block_is, conf->thrN * cuda_info->batch_size_per_block * sizeof(float)>>>(gpu_snn, iN, N, lt, t, cuda_info->n_neurons_per_launch, l);
            
            cudaCheckError(cudaPeekAtLastError());  // check launch errors
            cudaCheckError(cudaDeviceSynchronize());  // check runtime errors
            cudaDeviceSynchronize();   
        }

        clock_gettime(CLOCK_MONOTONIC, &start_learn);

        if(conf->learn){

            compute_pre_fired_batch<<<grid_synapses, block_synapses>>>(gpu_snn, N, iN, S, lt, t);
            cudaCheckError(cudaPeekAtLastError());  // check launch errors
            cudaCheckError(cudaDeviceSynchronize());  // check runtime errors
            cudaDeviceSynchronize();
                
        }

        clock_gettime(CLOCK_MONOTONIC, &end_learn);
        et_learn +=(end_learn.tv_sec - start_learn.tv_sec) + (end_learn.tv_nsec - start_learn.tv_nsec) / 1e9;


        clock_gettime(CLOCK_MONOTONIC, &end_in);
        et_in +=(end_in.tv_sec - start_in.tv_sec) + (end_in.tv_nsec - start_in.tv_nsec) / 1e9;

        // compute neuron dynamics
        clock_gettime(CLOCK_MONOTONIC, &start_v);

        process_V_batch<<<grid_neurons, block_neurons>>>(gpu_snn, N);
        cudaCheckError(cudaPeekAtLastError());  // check launch errors
        cudaCheckError(cudaDeviceSynchronize());  // check runtime errors
        cudaDeviceSynchronize();

        clock_gettime(CLOCK_MONOTONIC, &end_v);
        et_v +=(end_v.tv_sec - start_v.tv_sec) + (end_v.tv_nsec - start_v.tv_nsec) / 1e9;


        // compute output step
        clock_gettime(CLOCK_MONOTONIC, &start_out);

        process_firing_batch<<<grid_neurons, block_neurons>>>(gpu_snn, gpu_results, iN, N, lt);
        cudaCheckError(cudaPeekAtLastError());  // check launch errors
        cudaCheckError(cudaDeviceSynchronize());  // check runtime errors
        cudaDeviceSynchronize();

        clock_gettime(CLOCK_MONOTONIC, &end_out);
        et_out +=(end_out.tv_sec - start_out.tv_sec) + (end_out.tv_nsec - start_out.tv_nsec) / 1e9;


        // compute learning step if training
        clock_gettime(CLOCK_MONOTONIC, &start_learn);

        if(conf->learn){

            compute_post_traces_batch<<<grid_neurons, block_neurons>>>(gpu_snn, N);
            cudaCheckError(cudaPeekAtLastError());  // check launch errors
            cudaCheckError(cudaDeviceSynchronize());  // check runtime errors
            cudaDeviceSynchronize();

            trace_based_STDP_batch<<<grid_synapses, block_synapses>>>(gpu_snn, S);
            cudaCheckError(cudaPeekAtLastError());  // check launch errors
            cudaCheckError(cudaDeviceSynchronize());  // check runtime errors
            cudaDeviceSynchronize();

            reinit_post_fired_batch<<<grid_neurons, block_neurons>>>(gpu_snn, N);
            cudaCheckError(cudaPeekAtLastError());  // check launch errors
            cudaCheckError(cudaDeviceSynchronize());  // check runtime errors
            cudaDeviceSynchronize();

        }
        clock_gettime(CLOCK_MONOTONIC, &end_learn);
        et_learn +=(end_learn.tv_sec - start_learn.tv_sec) + (end_learn.tv_nsec - start_learn.tv_nsec) / 1e9;
    }


    // store in CPU results struct
    cpu_results->t_in = et_in;
    cpu_results->t_v = et_v;
    cpu_results->t_out = et_out;
    cpu_results->t_reinit = et_reinit;
    cpu_results->t_load = et_load;
    cpu_results->t_learn = et_learn;
}