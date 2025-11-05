#include <cuda.h>
#include <cuda_runtime.h>
#include <stdio.h>
#include "snn_library.h"
#include "neuron_models/lif_neuron.h"
#include "neuron_models/GPU_lif_neuron.cuh"
#include "training_rules/stdp.h"

#define THR_PER_BLOCK 1024 
//#define T 500

/**
D = A * B + C kalkulua egiten duen kernela
*/
__global__ void cuda_add_dot_matrix(int rowsAC, int colsBC, int colsArowsB, float *A, float *B, float *C, float *D)
{
    //lortu hariaren identifikadorea
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j;

    //hariak kalkulatu behar duen Dko elementuaren errenkada eta zutabea
    int zutab = i%colsBC;
    int errenk = i/colsBC;

    //egiaztatu hariak kalkulua egin behar duela
    if(i<(rowsAC * colsBC))
    {
        //D kalkulatzeko Ako errenkada eta Bko zutabea prozesatu
        for(j=0; j<colsArowsB; j++)  
            D[i]+=A[errenk * colsArowsB + j]*B[j * colsBC + zutab];

        D[i] += C[i];
    }
}


__global__ void load_sample_kernel(GPU_SNN_t *snn, GPU_input_info_t *dataset, int s)
{
    //lortu hariaren identifikadorea

    int threadsPerBlock = blockDim.x * blockDim.y;
    int threadNumInBlock = threadIdx.x + blockDim.x * threadIdx.y; //(alternatively: threadIdx.y + blockDim.y * threadIdx.x);
    int blockNumInGrid   = blockIdx.x  + gridDim.x  * blockIdx.y;  //(alternatively: blockIdx.y  + gridDim.y  * blockIdx.x);

    int globalThreadNum = blockNumInGrid * threadsPerBlock + threadNumInBlock;
    int neuron_index = globalThreadNum; // ?

    printf(" > Loading spike train in neuron = %d\n", neuron_index);

    int i;

    for(i = 0; i<dataset->n_spikes[s]; i++){
        //snn.in_spk_matrix[neuron_index * ]

        int spk = ;

        snn->spk_matrix[];
    }
}

__global__ void reinit_neurons_kernel(GPU_SNN_t *snn, int s)
{
    //lortu hariaren identifikadorea

    int threadsPerBlock = blockDim.x * blockDim.y;
    int threadNumInBlock = threadIdx.x + blockDim.x * threadIdx.y; //(alternatively: threadIdx.y + blockDim.y * threadIdx.x);
    int blockNumInGrid   = blockIdx.x  + gridDim.x  * blockIdx.y;  //(alternatively: blockIdx.y  + gridDim.y  * blockIdx.x);

    int globalThreadNum = blockNumInGrid * threadsPerBlock + threadNumInBlock;
    int neuron_index = globalThreadNum;

}

__global__ void reinit_synapses_kernel(GPU_SNN_t *snn, int s)
{
    //lortu hariaren identifikadorea

    int threadsPerBlock = blockDim.x * blockDim.y;
    int threadNumInBlock = threadIdx.x + blockDim.x * threadIdx.y; //(alternatively: threadIdx.y + blockDim.y * threadIdx.x);
    int blockNumInGrid   = blockIdx.x  + gridDim.x  * blockIdx.y;  //(alternatively: blockIdx.y  + gridDim.y  * blockIdx.x);

    int globalThreadNum = blockNumInGrid * threadsPerBlock + threadNumInBlock;
    int neuron_index = globalThreadNum;

}

__global__ void lif_neuron_input(GPU_SNN_t *snn, int s, int t)
{
    //lortu hariaren identifikadorea

    int threadsPerBlock = blockDim.x * blockDim.y;
    int threadNumInBlock = threadIdx.x + blockDim.x * threadIdx.y; //(alternatively: threadIdx.y + blockDim.y * threadIdx.x);
    int blockNumInGrid   = blockIdx.x  + gridDim.x  * blockIdx.y;  //(alternatively: blockIdx.y  + gridDim.y  * blockIdx.x);

    int globalThreadNum = blockNumInGrid * threadsPerBlock + threadNumInBlock;
    int neuron_index = globalThreadNum;

}

__global__ void lif_neuron_output(GPU_SNN_t *snn, int s, int t)
{
    //lortu hariaren identifikadorea

    int threadsPerBlock = blockDim.x * blockDim.y;
    int threadNumInBlock = threadIdx.x + blockDim.x * threadIdx.y; //(alternatively: threadIdx.y + blockDim.y * threadIdx.x);
    int blockNumInGrid   = blockIdx.x  + gridDim.x  * blockIdx.y;  //(alternatively: blockIdx.y  + gridDim.y  * blockIdx.x);

    int globalThreadNum = blockNumInGrid * threadsPerBlock + threadNumInBlock;
    int neuron_index = globalThreadNum;

}



__global__ void cuda_simulation_step_lif_neuron(int rowsAC, int colsBC, int colsArowsB, float *A, float *B, float *C, float *D)
{
    //lortu hariaren identifikadorea
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j;

    //hariak kalkulatu behar duen Dko elementuaren errenkada eta zutabea
    int zutab = i%colsBC;
    int errenk = i/colsBC;

    //egiaztatu hariak kalkulua egin behar duela
    if(i<(rowsAC * colsBC))
    {
        //D kalkulatzeko Ako errenkada eta Bko zutabea prozesatu
        for(j=0; j<colsArowsB; j++)  
            D[i]+=A[errenk * colsArowsB + j]*B[j * colsBC + zutab];

        D[i] += C[i];
    }
}


extern "C" void simulate_in_GPU(spiking_nn_t *snn, simulation_configuration_t *conf, input_data_t *dataset, simulation_results_t *results){

    int n_networks = 1;
    // Copy SNN to GPU n times
    GPU_SNN_t *d_GPU_SNN = copy_snn_structure_to_GPU(snn, n_networks);
    GPU_input_info_t *d_dataset = copy_dataset_to_GPU(dataset);

    // simulate
    simulate_snn_in_GPU(d_GPU_SNN, snn, conf, d_dataset, dataset, results);
}


extern "C" void simulate_snn_in_GPU(GPU_SNN_t *gpu_snn, spiking_nn_t *snn, simulation_configuration_t *conf, GPU_input_info_t *gpu_dataset, input_data_t *dataset, simulation_results_t *results){

    // call kernels
    int t, s;

    printf(" > n_samples = %d\n", dataset->n_samples);
    for(s = 0; s<dataset->n_samples; s++){

        // load sample
        printf(" > Running kernel...\n");
        fflush(stdout);
        load_sample_kernel<<<4, 196>>>(gpu_snn, gpu_dataset, s); // number of threads = n_input // 784
        cudaDeviceSynchronize();

        // reinit neurons
        reinit_neurons_kernel<<<17, 516>>>(gpu_snn, s); // n_neurons 8784
        cudaDeviceSynchronize();
        
        // reinit synapses
        reinit_synapses_kernel<<<2792, 516>>>(gpu_snn, s); // n_synapses 1440187
        cudaDeviceSynchronize();
        
        for(t = 0; t<conf->time_steps; t++){

            printf(" T = %d\n", t);

            // input step
            lif_neuron_input<<<17, 516>>>(gpu_snn, s, t);
            cudaDeviceSynchronize();

            // output step
            lif_neuron_output<<<17, 516>>>(gpu_snn, s, t);
            cudaDeviceSynchronize();
            
            // learning
            
        }
    }

}


extern "C" GPU_SNN_t* copy_snn_structure_to_GPU(spiking_nn_t *snn, int n){
    
    // event variables
    cudaEvent_t start, stop;
    
    // control variables
    int i, j;

    // cuda things
    float milliseconds = 0;
    int thr_per_blk_neurons, blk_in_grid_neurons, thr_per_blk_synapses, blk_in_grid_synapses;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
    
    // Map CPU SNN structure to GPU SNN structure
    // Map neurons indexes too

    // allocate memory for neurons
    GPU_SNN_t *gpu_snn_mapper = (GPU_SNN_t*)malloc(sizeof(GPU_SNN_t)); // input neurons???
    gpu_snn_mapper->v = (float *)malloc((snn->n_neurons) * sizeof(float));
    gpu_snn_mapper->v_thresh = (float *)malloc((snn->n_neurons) * sizeof(float));
    gpu_snn_mapper->v_rest = (float *)malloc((snn->n_neurons) * sizeof(float));
    gpu_snn_mapper->r_period = (int *)malloc((snn->n_neurons) * sizeof(int));
    gpu_snn_mapper->r_period_remain = (int *)malloc((snn->n_neurons) * sizeof(int));
    gpu_snn_mapper->res = (int *)malloc((snn->n_neurons) * sizeof(int));


    // allocate memory for spike matrix (virtual neurons are included too here)
    gpu_snn_mapper->spk_matrix = (int *)malloc((snn->n_neurons + snn->n_input) * T * sizeof(int));
    
    // count total number of input and output synapses
    int tmp_in = 0, tmp_out = 0;
    for(i = 0; i<snn->n_neurons; i++){
        tmp_in += snn->lif_neurons[i].n_input_synapse;
        tmp_out += snn->lif_neurons[i].n_output_synapse;
    }

    // allocate memory for input and output synapse indexes and offsets
    gpu_snn_mapper->in_synapses = (int*)malloc(tmp_in * sizeof(int));
    gpu_snn_mapper->in_off = (int*)malloc(snn->n_neurons * sizeof(int));
    gpu_snn_mapper->out_synapses = (int*)malloc(tmp_out * sizeof(int));
    gpu_snn_mapper->out_off = (int*)malloc(snn->n_neurons * sizeof(int));

    // copy CPU structure neurons data to GPU structure neurons
    tmp_in = 0; tmp_out = 0;
    for(i = 0; i<snn->n_neurons; i++){
        gpu_snn_mapper->v[i] = (float)(snn->lif_neurons[i].v);
        gpu_snn_mapper->v_thresh[i] = (float)(snn->lif_neurons[i].v_tresh);
        gpu_snn_mapper->v_rest[i] = (float)(snn->lif_neurons[i].v_rest);
        gpu_snn_mapper->r_period[i] = snn->lif_neurons[i].r_time;
        gpu_snn_mapper->r_period_remain[i] = 0;
        gpu_snn_mapper->res[i] = snn->lif_neurons[i].r;

        // offsets
        gpu_snn_mapper->in_off[i] = tmp_in;
        gpu_snn_mapper->out_off[i] = tmp_out;
        
        // input synapses
        for(j = 0; j<snn->lif_neurons[i].n_input_synapse; j++){
            gpu_snn_mapper->in_synapses[tmp_in] = snn->lif_neurons[i].input_synapse_indexes[j];
            tmp_in ++;
        }

        // output synapses
        for(j = 0; j<snn->lif_neurons[i].n_output_synapse; j++){
            gpu_snn_mapper->out_synapses[tmp_out] = snn->lif_neurons[i].output_synapse_indexes[j];
            tmp_out ++;            
        }
    }

    // initialize spike matrix
    for(i = 0; i<snn->n_neurons + snn->n_input; i++)
        for(j = 0; j<T; j++)
            gpu_snn_mapper->spk_matrix[i * T + j] = -1;


    // allocate memory for synapses
    gpu_snn_mapper->w = (float*)malloc(snn->n_synapses * sizeof(float));
    gpu_snn_mapper->dw = (float*)malloc(snn->n_synapses * sizeof(float));
    gpu_snn_mapper->delay = (int*)malloc(snn->n_synapses * sizeof(int));
    gpu_snn_mapper->lr = (int*)malloc(snn->n_synapses * sizeof(int));
    gpu_snn_mapper->pre_neuron_index = (int*)malloc(snn->n_synapses * sizeof(int));
    gpu_snn_mapper->post_neuron_index = (int*)malloc(snn->n_synapses * sizeof(int));


    // copy synapses and map computation neurons
    for(i = 0; i<snn->n_synapses; i++){

        gpu_snn_mapper->w[i] = (float)(snn->synapses[i].w);
        gpu_snn_mapper->dw[i] = 0.0;
        gpu_snn_mapper->delay[i] = snn->synapses[i].delay;
        gpu_snn_mapper->lr[i] = snn->synapses[i].lr;
        gpu_snn_mapper->pre_neuron_index[i] = snn->synapses[i].pre_neuron_index;
        gpu_snn_mapper->post_neuron_index[i] = snn->synapses[i].post_neuron_index;

        // if it is an input synapse, then it's pre neuron is a virtual one
        if(i > snn->n_input_synapses)
            gpu_snn_mapper->pre_neuron_index[i] += snn->n_input;
        gpu_snn_mapper->post_neuron_index[i] += snn->n_input;
    }


    /* Allocate memory on GPU and copy SNN */
    GPU_SNN_t *gpu_snn = (GPU_SNN_t*)malloc(sizeof(GPU_SNN_t)); // input neurons???
    
    cudaMalloc(&(gpu_snn->v), snn->n_neurons * sizeof(float)); // allocate memory for neurons
    cudaMalloc(&(gpu_snn->v_thresh), snn->n_neurons * sizeof(float)); // allocate memory for neurons
    cudaMalloc(&(gpu_snn->v_rest), snn->n_neurons * sizeof(float)); // allocate memory for neurons
    cudaMalloc(&(gpu_snn->r_period), snn->n_neurons * sizeof(int)); // allocate memory for neurons
    cudaMalloc(&(gpu_snn->r_period_remain), snn->n_neurons * sizeof(int)); // allocate memory for neurons
    cudaMalloc(&(gpu_snn->res), snn->n_neurons * sizeof(int)); // allocate memory for neurons

    cudaMalloc(&(gpu_snn->in_synapses), tmp_in * sizeof(int)); // allocate memory for neurons
    cudaMalloc(&(gpu_snn->in_off), snn->n_neurons * sizeof(int)); // allocate memory for neurons
    cudaMalloc(&(gpu_snn->out_synapses), tmp_out * sizeof(int)); // allocate memory for neurons
    cudaMalloc(&(gpu_snn->out_off), snn->n_neurons * sizeof(int)); // allocate memory for neurons
    
    cudaMalloc(&(gpu_snn->w), snn->n_synapses * sizeof(float)); // allocate memory for neurons
    cudaMalloc(&(gpu_snn->dw), snn->n_synapses * sizeof(float)); // allocate memory for neurons
    cudaMalloc(&(gpu_snn->delay), snn->n_synapses * sizeof(int)); // allocate memory for neurons
    cudaMalloc(&(gpu_snn->lr), snn->n_synapses * sizeof(int)); // allocate memory for neurons
    cudaMalloc(&(gpu_snn->pre_neuron_index), snn->n_synapses * sizeof(int)); // allocate memory for neurons
    cudaMalloc(&(gpu_snn->post_neuron_index), snn->n_synapses * sizeof(int)); // allocate memory for neurons

    cudaMalloc(&(gpu_snn->spk_matrix), (snn->n_neurons + snn->n_input) * T * sizeof(int)); // allocate memory for neurons

    
    // copy
    cudaMemcpy(gpu_snn->v, gpu_snn_mapper->v, snn->n_neurons * sizeof(float), cudaMemcpyHostToDevice); // copy neurons information
    cudaMemcpy(gpu_snn->v_thresh, gpu_snn_mapper->v_thresh, snn->n_neurons * sizeof(float), cudaMemcpyHostToDevice); // copy neurons information
    cudaMemcpy(gpu_snn->v_rest, gpu_snn_mapper->v_rest, snn->n_neurons * sizeof(float), cudaMemcpyHostToDevice); // copy neurons information
    cudaMemcpy(gpu_snn->r_period, gpu_snn_mapper->r_period, snn->n_neurons * sizeof(int), cudaMemcpyHostToDevice); // copy neurons information
    cudaMemcpy(gpu_snn->r_period_remain, gpu_snn_mapper->r_period_remain, snn->n_neurons * sizeof(int), cudaMemcpyHostToDevice); // copy neurons information
    cudaMemcpy(gpu_snn->res, gpu_snn_mapper->res, snn->n_neurons * sizeof(int), cudaMemcpyHostToDevice); // copy neurons information
    cudaMemcpy(gpu_snn->in_synapses, gpu_snn_mapper->in_synapses, tmp_in * sizeof(int), cudaMemcpyHostToDevice); // copy neurons information
    cudaMemcpy(gpu_snn->out_synapses, gpu_snn_mapper->out_synapses, tmp_out * sizeof(int), cudaMemcpyHostToDevice); // copy neurons information
    cudaMemcpy(gpu_snn->in_off, gpu_snn_mapper->in_off, snn->n_neurons * sizeof(int), cudaMemcpyHostToDevice); // copy neurons information
    cudaMemcpy(gpu_snn->out_off, gpu_snn_mapper->out_off, snn->n_neurons * sizeof(int), cudaMemcpyHostToDevice); // copy neurons information

    cudaMemcpy(gpu_snn->w, gpu_snn_mapper->w, snn->n_synapses * sizeof(float), cudaMemcpyHostToDevice); // copy neurons information
    cudaMemcpy(gpu_snn->dw, gpu_snn_mapper->dw, snn->n_synapses * sizeof(float), cudaMemcpyHostToDevice); // copy neurons information
    cudaMemcpy(gpu_snn->delay, gpu_snn_mapper->delay, snn->n_synapses * sizeof(int), cudaMemcpyHostToDevice); // copy neurons information
    cudaMemcpy(gpu_snn->lr, gpu_snn_mapper->lr, snn->n_synapses * sizeof(int), cudaMemcpyHostToDevice); // copy neurons information
    cudaMemcpy(gpu_snn->pre_neuron_index, gpu_snn_mapper->pre_neuron_index, snn->n_synapses * sizeof(int), cudaMemcpyHostToDevice); // copy neurons information
    cudaMemcpy(gpu_snn->post_neuron_index, gpu_snn_mapper->post_neuron_index, snn->n_synapses * sizeof(int), cudaMemcpyHostToDevice); // copy neurons information

    cudaMemcpy(gpu_snn->spk_matrix, gpu_snn_mapper->spk_matrix, (snn->n_neurons + snn->n_input) * sizeof(int), cudaMemcpyHostToDevice);

    GPU_SNN_t *d_GPU_SNN;
    cudaMalloc(&d_GPU_SNN, sizeof(GPU_SNN_t));
    cudaMemcpy(d_GPU_SNN, gpu_snn, sizeof(GPU_SNN_t), cudaMemcpyHostToDevice); // copy neurons information

    printf(" > SNN structure succesfully copied to GPU!!\n");
    fflush(stdout);

    return d_GPU_SNN;
}

extern "C" GPU_input_info_t* copy_dataset_to_GPU(input_data_t *dataset){
    
    int i, j, l, total_n_spikes = 0;

    GPU_input_info_t *tmp_gpu_input_dataset, *gpu_input_dataset, *d_gpu_input_dataset;
    tmp_gpu_input_dataset = (GPU_input_info_t *)malloc(sizeof(GPU_input_info_t));
    
    // copy first to a temporal struct in the CPU
    tmp_gpu_input_dataset->type = dataset->type;
    tmp_gpu_input_dataset->n_classes = dataset->n_classes;
    tmp_gpu_input_dataset->n_samples = dataset->n_samples;

    // array to store the size of each sample (number of elements / features)
    tmp_gpu_input_dataset->n_features = (int*)malloc(dataset->n_samples * sizeof(int)); // number of features for each sample
    tmp_gpu_input_dataset->n_spikes = (int*)malloc(dataset->n_samples * dataset->image_size * sizeof(int)); // number of spikes for each feature (the number of features can not be constant)

    tmp_gpu_input_dataset->sample_offset = (int*)malloc(dataset->n_samples * sizeof(int)); // offset of the sample
    tmp_gpu_input_dataset->feature_offset = (int*)malloc(dataset->n_samples * dataset->image_size * sizeof(int)); // offset of each feature (the number of features can not be constnant)


    // count total number of spikes in the dataset (or the part of the dataset)
    // loop over samples
    for(i = 0; i<dataset->n_samples; i++){

        // loop over the features of each sample
        for(j = 0; j<dataset->image_size; j++){

            // count the number of spikes
            total_n_spikes += dataset->samples[i].st[j].n_spikes;
        }
    }

    // allocate memory to store all the spikes in a 1D array
    tmp_gpu_input_dataset->spikes = (int*)malloc(total_n_spikes * sizeof(int));


    // count total number of spikes
    int h_spk_ind = 0; // spike position index
    int next_feature = 0;

    // loop over all samples
    for(i = 0; i<dataset->n_samples; i++){

        // set the number of features for the sample i
        tmp_gpu_input_dataset->n_features[i] = dataset->image_size; 

        // store the offset of the sample (the index of the first spike)
        tmp_gpu_input_dataset->sample_offset[i] = h_spk_ind; 
        

        // loop over sample features
        for(j = 0; j<dataset->image_size; j++){

            // store the number of spikes in the feature
            tmp_gpu_input_dataset->n_spikes[next_feature] = dataset->samples[i].st[j].n_spikes; 

            // store the local index of the first element in the feature
            tmp_gpu_input_dataset->feature_offset[next_feature] = h_spk_ind - tmp_gpu_input_dataset->sample_offset[i];
            next_feature ++; // we moved to the next feature

            // loop over feature array of spikes
            for(l = 0; l<dataset->samples[i].st[j].n_spikes; l++){

                tmp_gpu_input_dataset->spikes[h_spk_ind] = dataset->samples[i].st[j].stimes[l]; // store the spikes
                h_spk_ind ++; // update the index for the spikes
            }
        }
    }

    printf(" > Dataset copied to new struct type!\n");
    fflush(stdout);
 
    // allocate to temporal
    gpu_input_dataset = (GPU_input_info_t *)malloc(sizeof(GPU_input_info_t));
    cudaMalloc(&(gpu_input_dataset->spikes), total_n_spikes * sizeof(int)); // allocate memory for neurons
    cudaMalloc(&(gpu_input_dataset->n_features), dataset->n_samples * sizeof(int)); // allocate memory for neurons
    cudaMalloc(&(gpu_input_dataset->n_spikes), dataset->n_samples * dataset->image_size * sizeof(int)); // allocate memory for neurons
    cudaMalloc(&(gpu_input_dataset->sample_offset), dataset->n_samples * sizeof(int)); // allocate memory for neurons
    cudaMalloc(&(gpu_input_dataset->feature_offset), dataset->n_samples * dataset->image_size * sizeof(int)); // allocate memory for neurons

    // memcpy
    cudaMemcpy(gpu_input_dataset->spikes, tmp_gpu_input_dataset->spikes, total_n_spikes * sizeof(int), cudaMemcpyHostToDevice); // copy neurons information
    cudaMemcpy(gpu_input_dataset->n_features, tmp_gpu_input_dataset->n_features, dataset->n_samples * sizeof(int), cudaMemcpyHostToDevice); // copy neurons information
    cudaMemcpy(gpu_input_dataset->n_spikes, tmp_gpu_input_dataset->n_spikes, dataset->n_samples * dataset->image_size * sizeof(int), cudaMemcpyHostToDevice); // copy neurons information
    cudaMemcpy(gpu_input_dataset->sample_offset, tmp_gpu_input_dataset->sample_offset, dataset->n_samples * sizeof(int), cudaMemcpyHostToDevice); // copy neurons information
    cudaMemcpy(gpu_input_dataset->feature_offset, tmp_gpu_input_dataset->feature_offset, dataset->n_samples * dataset->image_size * sizeof(int), cudaMemcpyHostToDevice); // copy neurons information
    gpu_input_dataset->n_samples = tmp_gpu_input_dataset->n_samples;
    gpu_input_dataset->type = tmp_gpu_input_dataset->type;
    gpu_input_dataset->n_classes = tmp_gpu_input_dataset->n_classes;
    

    // copy to GPU
    cudaMalloc(&(d_gpu_input_dataset), sizeof(GPU_input_info_t)); // allocate memory for neurons
    cudaMemcpy(d_gpu_input_dataset, gpu_input_dataset->spikes, sizeof(GPU_input_info_t), cudaMemcpyHostToDevice); // copy neurons information


    printf(" > Dataset copied to the GPU!\n");
    fflush(stdout);
 
    return d_gpu_input_dataset;
}


void getProperties(){
    int nDevices;
  cudaGetDeviceCount(&nDevices);
  
  printf("Number of devices: %d\n", nDevices);
  
  for (int i = 0; i < nDevices; i++) {
    cudaDeviceProp prop;
    cudaGetDeviceProperties(&prop, i);
    printf("Device Number: %d\n", i);
    printf("  Device name: %s\n", prop.name);
    printf("  Memory Clock Rate (MHz): %d\n",
           prop.memoryClockRate/1024);
    printf("  Memory Bus Width (bits): %d\n",
           prop.memoryBusWidth);
    printf("  Peak Memory Bandwidth (GB/s): %.1f\n",
           2.0*prop.memoryClockRate*(prop.memoryBusWidth/8)/1.0e6);
    printf("  Total global memory (Gbytes) %.1f\n",(float)(prop.totalGlobalMem)/1024.0/1024.0/1024.0);
    printf("  Shared memory per block (Kbytes) %.1f\n",(float)(prop.sharedMemPerBlock)/1024.0);
    printf("  minor-major: %d-%d\n", prop.minor, prop.major);
    printf("  Warp-size: %d\n", prop.warpSize);
    printf("  Concurrent kernels: %s\n", prop.concurrentKernels ? "yes" : "no");
    printf("  Concurrent computation/communication: %s\n\n",prop.deviceOverlap ? "yes" : "no");
  }
}
