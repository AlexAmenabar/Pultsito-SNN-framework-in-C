#include <cuda.h>
#include <cuda_runtime.h>
#include <stdio.h>
#include "snn_library.h"
#include "neuron_models/lif_neuron.h"
#include "neuron_models/GPU_lif_neuron.cuh"
#include "training_rules/stdp.h"

#define THR_PER_BLOCK 1024 


__global__ void reinit_spk_matrix(GPU_SNN_t *snn){
    
    int i;
    int threadsPerBlock = blockDim.x * blockDim.y;
    int threadNumInBlock = threadIdx.x + blockDim.x * threadIdx.y; //(alternatively: threadIdx.y + blockDim.y * threadIdx.x);
    int blockNumInGrid   = blockIdx.x  + gridDim.x  * blockIdx.y;  //(alternatively: blockIdx.y  + gridDim.y  * blockIdx.x);

    int globalThreadNum = blockNumInGrid * threadsPerBlock + threadNumInBlock;
    int matrix_index = globalThreadNum; // ?

    if(matrix_index < (snn->n_neurons + snn->n_input_neurons) * L){
        snn->spk_matrix[matrix_index] = -1;
    }
}


__global__ void load_sample_kernel(GPU_SNN_t *snn, GPU_input_info_t *dataset, int s)
{
    //lortu hariaren identifikadorea

    int i;
    int threadsPerBlock = blockDim.x * blockDim.y;
    int threadNumInBlock = threadIdx.x + blockDim.x * threadIdx.y; //(alternatively: threadIdx.y + blockDim.y * threadIdx.x);
    int blockNumInGrid   = blockIdx.x  + gridDim.x  * blockIdx.y;  //(alternatively: blockIdx.y  + gridDim.y  * blockIdx.x);

    int globalThreadNum = blockNumInGrid * threadsPerBlock + threadNumInBlock;
    int neuron_index = globalThreadNum; // ?

    //printf(" > Neuron index = %d, input neurons = %d, sample index = %d, n_features = %d\n", 
            //neuron_index, snn->n_input_neurons, s, dataset->n_features);

    // if it is an input neuron
    if(neuron_index < snn->n_input_neurons){
        
        // s is the sample index and f the feature index
        // start_sample is the index in which the sample starts in the array of spikes
        // start_feature is the index in which the feature starts locally in the array of spikes

        int f = neuron_index; // the neuron index is also the feature index

        //printf(" > > Feature = %d\n", f);

        int start_sample = dataset->sample_offset[s];
        //printf(" Thread %d, start sample = %d\n", f, start_sample);
        
        int start_feature = dataset->feature_offset[s * dataset->n_features + f];
        //printf(" Thread %d, start feature = %d\n", f, start_feature);

        int n_spikes = dataset->n_spikes[s * dataset->n_features + f];

        //printf(" > Loading spike train in neuron = %d (n = %d)\n", neuron_index, n_spikes);

        for(i = 0; i<n_spikes; i++){
        
            // store a 1 value if there is a spike
            snn->spk_matrix[neuron_index * L + dataset->spikes[start_sample + start_feature + i]] = dataset->spikes[start_sample + start_feature + i];

            //printf(" SPike added by thread %d to %d (value %d).\n", f, neuron_index * L + dataset->spikes[start_sample + start_feature + i], dataset->spikes[start_sample + start_feature + i]);
        }
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


    if(neuron_index < snn->n_neurons){

        snn->v[neuron_index] = snn->v_rest[neuron_index];
        snn->r_period_remain[neuron_index] = 0;
        snn->t_last_spike[neuron_index] = 0;
    }

}

__global__ void reinit_synapses_kernel(GPU_SNN_t *snn, int s)
{
    //lortu hariaren identifikadorea

    int threadsPerBlock = blockDim.x * blockDim.y;
    int threadNumInBlock = threadIdx.x + blockDim.x * threadIdx.y; //(alternatively: threadIdx.y + blockDim.y * threadIdx.x);
    int blockNumInGrid   = blockIdx.x  + gridDim.x  * blockIdx.y;  //(alternatively: blockIdx.y  + gridDim.y  * blockIdx.x);

    int globalThreadNum = blockNumInGrid * threadsPerBlock + threadNumInBlock;
    int synapse_index = globalThreadNum;

    if(synapse_index < snn->n_synapses){

        //continue;
        
    }
}

__global__ void lif_neuron_input(GPU_SNN_t *snn, int s, int t)
{
    //lortu hariaren identifikadorea

    int threadsPerBlock = blockDim.x * blockDim.y;
    int threadNumInBlock = threadIdx.x + blockDim.x * threadIdx.y; //(alternatively: threadIdx.y + blockDim.y * threadIdx.x);
    int blockNumInGrid   = blockIdx.x  + gridDim.x  * blockIdx.y;  //(alternatively: blockIdx.y  + gridDim.y  * blockIdx.x);

    int globalThreadNum = blockNumInGrid * threadsPerBlock + threadNumInBlock;
    int neuron_index = globalThreadNum;
    int n, i;


    if(neuron_index < snn->n_neurons){

        n = neuron_index;
        float I = 0.0;

        // loop over input synapses to check if the index should be updated in one of them
        if(snn->r_period_remain[n] <= 0){
            
            for(i = 0; i<snn->n_neuron_input_synapses[n]; i++){

                int synapse_index = snn->neuron_input_synapses_offset[n] + i;
                int in_neuron_index = snn->pre_neuron_index[synapse_index];
                int latency = snn->delay[synapse_index];

                int spk_time = t - snn->delay[synapse_index]; // actual position

                if(spk_time >=0 ){

                    int spk_time_index = spk_time % L; // scale to circular array
                    int spk = snn->spk_matrix[in_neuron_index * L + spk_time_index];
                    if(spk == spk_time){

                        I += snn->w[synapse_index];
                    }
                    printf(" Neuron %d: input synapse %d (l=%d, in neuron=%d), spk_time=%d, spk_time_index=%d, spk=%d\n", n, synapse_index, latency, in_neuron_index, spk_time, spk_time_index, spk);

                }    
          
            }
            
            snn->v[n] = (1-0.05) * snn->v[n] + snn->v_rest[n] * 0.05 + I;

            printf(" > Neuron %d, V = %f, I = %f\n", n, snn->v[n], I);

        }
    }
}

__global__ void lif_neuron_output(GPU_SNN_t *snn, int s, int t, GPU_results_t *results)
{
    //lortu hariaren identifikadorea

    int threadsPerBlock = blockDim.x * blockDim.y;
    int threadNumInBlock = threadIdx.x + blockDim.x * threadIdx.y; //(alternatively: threadIdx.y + blockDim.y * threadIdx.x);
    int blockNumInGrid   = blockIdx.x  + gridDim.x  * blockIdx.y;  //(alternatively: blockIdx.y  + gridDim.y  * blockIdx.x);

    int globalThreadNum = blockNumInGrid * threadsPerBlock + threadNumInBlock;
    int neuron_index = globalThreadNum;
    int n = neuron_index;


    if(n < snn->n_neurons){

        snn->r_period_remain[n] --;

        if(snn->v[n] >= snn->v_thresh[n]){

            int index = t % L;
            index = ((n + snn->n_input_neurons) * L + index);
            
            snn->spk_matrix[index] = t;
            snn->r_period_remain[n] = snn->r_period[n];
            snn->t_last_spike[n] = t;
            snn->v[n] = snn->v_rest[n];
            results->nspk[s * snn->n_neurons + n] += 1;
            printf(" > Neuron %d, t %d, SPIKE\n", n, t);
        }
    }

}



extern "C" void simulate_in_GPU(spiking_nn_t *snn, simulation_configuration_t *conf, input_data_t *dataset, simulation_results_t *results){

    int n_networks = 1;
    // Copy SNN to GPU n times
    GPU_SNN_t *d_GPU_SNN = copy_snn_structure_to_GPU(snn, n_networks);
    GPU_input_info_t *d_dataset = copy_dataset_to_GPU(dataset);
    GPU_results_t *d_results = copy_results_to_GPU(snn, dataset, results);

    // simulate
    simulate_snn_in_GPU(d_GPU_SNN, snn, conf, d_dataset, dataset, d_results, results);

    cudaDeviceSynchronize();

    GPU_results_t *tmp_results = (GPU_results_t*)malloc(sizeof(GPU_results_t));
    int *nspk = (int*)malloc(snn->n_neurons * dataset->n_samples * sizeof(int));
    cudaMemcpy(tmp_results, d_results, sizeof(GPU_results_t), cudaMemcpyDeviceToHost);
    cudaMemcpy(nspk, tmp_results->nspk, snn->n_neurons * dataset->n_samples * sizeof(int), cudaMemcpyDeviceToHost);

    printf("\n");
    for(int i = 0; i<dataset->n_samples; i++){

        printf(" > Sample %d: ", i);
        for(int j = 0; j<snn->n_neurons; j++){

            printf("%d ", nspk[i * snn->n_neurons + j]);
        }
        printf("\n");
    }
}


extern "C" void simulate_snn_in_GPU(GPU_SNN_t *gpu_snn, spiking_nn_t *snn, simulation_configuration_t *conf, GPU_input_info_t *gpu_dataset, input_data_t *dataset, GPU_results_t *gpu_results, simulation_results_t *results){

    // call kernels
    int t, s;

    printf(" > n_samples = %d\n", dataset->n_samples);
    for(s = 0; s<dataset->n_samples; s++){

        // load sample
        printf(" > Running kernel...\n");
        fflush(stdout);
        //load_sample_kernel<<<4, 196>>>(gpu_snn, gpu_dataset, s); // number of threads = n_input // 784
        
        reinit_spk_matrix<<<(snn->n_neurons + snn->n_input), L>>>(gpu_snn);
        cudaDeviceSynchronize();

        load_sample_kernel<<<1,128>>>(gpu_snn, gpu_dataset, s);
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
            lif_neuron_output<<<17, 516>>>(gpu_snn, s, t, gpu_results);
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
    int offset;

    // cuda things
    float milliseconds = 0;
    int thr_per_blk_neurons, blk_in_grid_neurons, thr_per_blk_synapses, blk_in_grid_synapses;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
    
    // Map CPU SNN structure to GPU SNN structure
    // Map neurons indexes too


    // STEP 1: create structure to map from CPU old structure to CPU new structure

    // allocate memory
    GPU_SNN_t *gpu_snn_mapper = (GPU_SNN_t*)malloc(sizeof(GPU_SNN_t)); // input neurons???

    // neurons
    gpu_snn_mapper->v = (float *)malloc((snn->n_neurons) * sizeof(float));
    gpu_snn_mapper->v_thresh = (float *)malloc((snn->n_neurons) * sizeof(float));
    gpu_snn_mapper->v_rest = (float *)malloc((snn->n_neurons) * sizeof(float));
    gpu_snn_mapper->r_period = (int *)malloc((snn->n_neurons) * sizeof(int));
    gpu_snn_mapper->r_period_remain = (int *)calloc(snn->n_neurons, sizeof(int)); // 0
    gpu_snn_mapper->res = (int *)malloc((snn->n_neurons) * sizeof(int));
    gpu_snn_mapper->t_last_spike = (int*)calloc(snn->n_neurons, sizeof(int)); // 0
    gpu_snn_mapper->n_neuron_input_synapses = (int*)malloc(snn->n_neurons * sizeof(int));
    gpu_snn_mapper->neuron_input_synapses_offset = (int*)malloc(snn->n_neurons * sizeof(int));
    gpu_snn_mapper->last_spk = (int*)calloc(snn->n_neurons, sizeof(int));
    gpu_snn_mapper->next_spk = (int*)calloc(snn->n_synapses - snn->n_output_synapses, sizeof(int));

    // synapses
    gpu_snn_mapper->w = (float*)malloc(snn->n_synapses * sizeof(float));
    gpu_snn_mapper->dw = (float*)calloc(snn->n_synapses, sizeof(float)); // 0
    gpu_snn_mapper->delay = (int*)malloc(snn->n_synapses * sizeof(int));
    gpu_snn_mapper->lr = (int*)malloc(snn->n_synapses * sizeof(int));
    gpu_snn_mapper->pre_neuron_index = (int*)malloc(snn->n_synapses * sizeof(int));
    gpu_snn_mapper->post_neuron_index = (int*)malloc(snn->n_synapses * sizeof(int));


    // allocate memory for spike matrix (virtual neurons are included too here) // it should be n_input_synapses?
    gpu_snn_mapper->spk_matrix = (int *)malloc((snn->n_neurons + snn->n_input) * L * sizeof(int));

    // count total number of input and output synapses
    int tmp_in = 0, tmp_out = 0;
    for(i = 0; i<snn->n_neurons; i++){
        tmp_in += snn->lif_neurons[i].n_input_synapse;
        tmp_out += snn->lif_neurons[i].n_output_synapse;
    }


    // STEP 2: copy CPU struct information to GPU struct in CPU
    tmp_in = 0; tmp_out = 0;

    // cpy neurons
    offset = 0;
    for(i = 0; i<snn->n_neurons; i++){
        
        // copy information from CPU
        gpu_snn_mapper->v[i] = (float)(snn->lif_neurons[i].v);
        gpu_snn_mapper->v_thresh[i] = (float)(snn->lif_neurons[i].v_tresh);
        gpu_snn_mapper->v_rest[i] = (float)(snn->lif_neurons[i].v_rest);
        gpu_snn_mapper->r_period[i] = snn->lif_neurons[i].r_time;
        gpu_snn_mapper->r_period_remain[i] = 0;
        gpu_snn_mapper->res[i] = snn->lif_neurons[i].r;
        gpu_snn_mapper->t_last_spike[i] = 0;

        // input synapses and offsets
        gpu_snn_mapper->n_neuron_input_synapses[i] = snn->lif_neurons[i].n_input_synapse;
        gpu_snn_mapper->neuron_input_synapses_offset[i] = offset;
        
        for(j = offset; j<snn->lif_neurons[i].n_input_synapse; j++){

            gpu_snn_mapper->next_spk[j] = 0;
        }
        gpu_snn_mapper->last_spk[i] = 0; // no spikes yet

        // update offset for the next iteration
        offset += gpu_snn_mapper->n_neuron_input_synapses[i];

        // offsets
        /*gpu_snn_mapper->in_off[i] = tmp_in;
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
        }*/
    }


    // cpy synapses
    for(i = 0; i<snn->n_synapses; i++){

        gpu_snn_mapper->w[i] = (float)(snn->synapses[i].w);
        gpu_snn_mapper->dw[i] = 0.0;
        gpu_snn_mapper->delay[i] = snn->synapses[i].delay;
        gpu_snn_mapper->lr[i] = snn->synapses[i].lr;
        gpu_snn_mapper->pre_neuron_index[i] = snn->synapses[i].pre_neuron_index;
        gpu_snn_mapper->post_neuron_index[i] = snn->synapses[i].post_neuron_index;

        // if it is an input synapse, then it's pre neuron is a virtual one. Map indexes
        if(i >= snn->n_input_synapses)
            gpu_snn_mapper->pre_neuron_index[i] += snn->n_input; // n_input_synapses??
        gpu_snn_mapper->post_neuron_index[i] += snn->n_input;

        //printf(" Printing Synapse %d: w (%f), ")
    }


    // initialize spike matrix // TODO: n_input_synapses???
    for(i = 0; i<snn->n_neurons + snn->n_input; i++)
        for(j = 0; j<L; j++)
            gpu_snn_mapper->spk_matrix[i * L + j] = -1; // no spikes



    /* STEP 3: Allocate memory on GPU and copy internal arrays */
    
    // allocate
    GPU_SNN_t *gpu_snn = (GPU_SNN_t*)malloc(sizeof(GPU_SNN_t)); // input neurons???
    
    gpu_snn->n_neurons = snn->n_neurons;
    gpu_snn->n_input_neurons = snn->n_input;
    gpu_snn->n_output_neurons = snn->n_output;
    gpu_snn->n_synapses = snn->n_synapses;
    gpu_snn->n_input_synapses = snn->n_input_synapses;
    gpu_snn->n_output_synapses = snn->n_output_synapses;
    gpu_snn->T = L;

    cudaMalloc(&(gpu_snn->v), snn->n_neurons * sizeof(float)); // allocate memory for neurons
    cudaMalloc(&(gpu_snn->v_thresh), snn->n_neurons * sizeof(float)); // allocate memory for neurons
    cudaMalloc(&(gpu_snn->v_rest), snn->n_neurons * sizeof(float)); // allocate memory for neurons
    cudaMalloc(&(gpu_snn->r_period), snn->n_neurons * sizeof(int)); // allocate memory for neurons
    cudaMalloc(&(gpu_snn->r_period_remain), snn->n_neurons * sizeof(int)); // allocate memory for neurons
    cudaMalloc(&(gpu_snn->res), snn->n_neurons * sizeof(int)); // allocate memory for neurons
    cudaMalloc(&(gpu_snn->t_last_spike), snn->n_neurons * sizeof(int)); // allocate memory for neurons
    cudaMalloc(&(gpu_snn->n_neuron_input_synapses), snn->n_neurons * sizeof(int));
    cudaMalloc(&(gpu_snn->neuron_input_synapses_offset), snn->n_neurons * sizeof(int));
    cudaMalloc(&(gpu_snn->last_spk), snn->n_neurons * sizeof(int));
    cudaMalloc(&(gpu_snn->next_spk), (snn->n_synapses - snn->n_output_synapses) * sizeof(int));


    /*cudaMalloc(&(gpu_snn->in_synapses), tmp_in * sizeof(int)); // allocate memory for neurons
    cudaMalloc(&(gpu_snn->in_off), snn->n_neurons * sizeof(int)); // allocate memory for neurons
    cudaMalloc(&(gpu_snn->out_synapses), tmp_out * sizeof(int)); // allocate memory for neurons
    cudaMalloc(&(gpu_snn->out_off), snn->n_neurons * sizeof(int)); // allocate memory for neurons*/
    
    cudaMalloc(&(gpu_snn->w), snn->n_synapses * sizeof(float)); // allocate memory for neurons
    cudaMalloc(&(gpu_snn->dw), snn->n_synapses * sizeof(float)); // allocate memory for neurons
    cudaMalloc(&(gpu_snn->delay), snn->n_synapses * sizeof(int)); // allocate memory for neurons
    cudaMalloc(&(gpu_snn->lr), snn->n_synapses * sizeof(int)); // allocate memory for neurons
    cudaMalloc(&(gpu_snn->pre_neuron_index), snn->n_synapses * sizeof(int)); // allocate memory for neurons
    cudaMalloc(&(gpu_snn->post_neuron_index), snn->n_synapses * sizeof(int)); // allocate memory for neurons

    // TODO: n_input_synapses???
    cudaMalloc(&(gpu_snn->spk_matrix), (snn->n_neurons + snn->n_input) * L * sizeof(int)); // allocate memory for neurons

    // cpy
    cudaMemcpy(gpu_snn->v, gpu_snn_mapper->v, snn->n_neurons * sizeof(float), cudaMemcpyHostToDevice); // copy neurons information
    cudaMemcpy(gpu_snn->v_thresh, gpu_snn_mapper->v_thresh, snn->n_neurons * sizeof(float), cudaMemcpyHostToDevice); // copy neurons information
    cudaMemcpy(gpu_snn->v_rest, gpu_snn_mapper->v_rest, snn->n_neurons * sizeof(float), cudaMemcpyHostToDevice); // copy neurons information
    cudaMemcpy(gpu_snn->r_period, gpu_snn_mapper->r_period, snn->n_neurons * sizeof(int), cudaMemcpyHostToDevice); // copy neurons information
    cudaMemcpy(gpu_snn->r_period_remain, gpu_snn_mapper->r_period_remain, snn->n_neurons * sizeof(int), cudaMemcpyHostToDevice); // copy neurons information
    cudaMemcpy(gpu_snn->res, gpu_snn_mapper->res, snn->n_neurons * sizeof(int), cudaMemcpyHostToDevice); // copy neurons information
    cudaMemcpy(gpu_snn->t_last_spike, gpu_snn_mapper->t_last_spike, snn->n_neurons * sizeof(int), cudaMemcpyHostToDevice); // copy neurons information
    cudaMemcpy(gpu_snn->n_neuron_input_synapses, gpu_snn_mapper->n_neuron_input_synapses, snn->n_neurons * sizeof(int), cudaMemcpyHostToDevice); // copy neurons information
    cudaMemcpy(gpu_snn->neuron_input_synapses_offset, gpu_snn_mapper->neuron_input_synapses_offset, snn->n_neurons * sizeof(int), cudaMemcpyHostToDevice); // copy neurons information
    cudaMemcpy(gpu_snn->last_spk, gpu_snn_mapper->last_spk, snn->n_neurons * sizeof(int), cudaMemcpyHostToDevice); // copy neurons information
    cudaMemcpy(gpu_snn->next_spk, gpu_snn_mapper->next_spk, (snn->n_synapses - snn->n_output_synapses) * sizeof(int), cudaMemcpyHostToDevice); // copy neurons information

    
    /*cudaMemcpy(gpu_snn->in_synapses, gpu_snn_mapper->in_synapses, tmp_in * sizeof(int), cudaMemcpyHostToDevice); // copy neurons information
    cudaMemcpy(gpu_snn->out_synapses, gpu_snn_mapper->out_synapses, tmp_out * sizeof(int), cudaMemcpyHostToDevice); // copy neurons information
    cudaMemcpy(gpu_snn->in_off, gpu_snn_mapper->in_off, snn->n_neurons * sizeof(int), cudaMemcpyHostToDevice); // copy neurons information
    cudaMemcpy(gpu_snn->out_off, gpu_snn_mapper->out_off, snn->n_neurons * sizeof(int), cudaMemcpyHostToDevice); // copy neurons information*/

    cudaMemcpy(gpu_snn->w, gpu_snn_mapper->w, snn->n_synapses * sizeof(float), cudaMemcpyHostToDevice); // copy neurons information
    cudaMemcpy(gpu_snn->dw, gpu_snn_mapper->dw, snn->n_synapses * sizeof(float), cudaMemcpyHostToDevice); // copy neurons information
    cudaMemcpy(gpu_snn->delay, gpu_snn_mapper->delay, snn->n_synapses * sizeof(int), cudaMemcpyHostToDevice); // copy neurons information
    cudaMemcpy(gpu_snn->lr, gpu_snn_mapper->lr, snn->n_synapses * sizeof(int), cudaMemcpyHostToDevice); // copy neurons information
    cudaMemcpy(gpu_snn->pre_neuron_index, gpu_snn_mapper->pre_neuron_index, snn->n_synapses * sizeof(int), cudaMemcpyHostToDevice); // copy neurons information
    cudaMemcpy(gpu_snn->post_neuron_index, gpu_snn_mapper->post_neuron_index, snn->n_synapses * sizeof(int), cudaMemcpyHostToDevice); // copy neurons information

    // n_input?
    cudaMemcpy(gpu_snn->spk_matrix, gpu_snn_mapper->spk_matrix, (snn->n_neurons + snn->n_input) * L * sizeof(int), cudaMemcpyHostToDevice);

    // STEP 5: copy entire structure to GPU
    GPU_SNN_t *d_GPU_SNN;
    cudaMalloc(&d_GPU_SNN, sizeof(GPU_SNN_t));
    cudaMemcpy(d_GPU_SNN, gpu_snn, sizeof(GPU_SNN_t), cudaMemcpyHostToDevice); // copy neurons information


    // TODO: deallocate memory

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
    tmp_gpu_input_dataset->n_features = dataset->image_size; // TODO: generalize
    //tmp_gpu_input_dataset->n_features = (int*)malloc(dataset->n_samples * sizeof(int)); // number of features for each sample
    tmp_gpu_input_dataset->n_spikes = (int*)malloc(dataset->n_samples * dataset->image_size * sizeof(int)); 

    tmp_gpu_input_dataset->sample_offset = (int*)malloc(tmp_gpu_input_dataset->n_samples * sizeof(int));
    tmp_gpu_input_dataset->feature_offset = (int*)malloc(tmp_gpu_input_dataset->n_samples * tmp_gpu_input_dataset->n_features * sizeof(int));


    // loop over samples to count the total number of spikes in the dataset
    for(i = 0; i<dataset->n_samples; i++){
        // loop over the features of each sample
        for(j = 0; j<dataset->image_size; j++){
            // count the number of spikes
            total_n_spikes += dataset->samples[i].st[j].n_spikes;
        }
    }

    // allocate memory to store all the spikes in a 1D array
    tmp_gpu_input_dataset->spikes = (int*)malloc(total_n_spikes * sizeof(int));


    // get offsets for samples and features
    int next_spike = 0;
    int next_feature = 0; // local to sample
    int n_samples, n_features;
    n_samples = tmp_gpu_input_dataset->n_samples;
    n_features = tmp_gpu_input_dataset->n_features;

    for(i = 0; i<n_samples; i++){

        tmp_gpu_input_dataset->sample_offset[i] = next_spike;
        next_feature = 0; // local for sample

        // loop over features
        for(j = 0; j<n_features; j++){

            tmp_gpu_input_dataset->feature_offset[i * n_features + j] = next_feature;

            // store the number of spikes of the feature
            tmp_gpu_input_dataset->n_spikes[i * n_features + j] = dataset->samples[i].st[j].n_spikes;
            printf(" > Sample %d, feature %d, n. spikes = %d\n", i, j, dataset->samples[i].st[j].n_spikes);

            // loop over spikes
            for(l = 0; l<dataset->samples[i].st[j].n_spikes; l++){

                tmp_gpu_input_dataset->spikes[next_spike] = dataset->samples[i].st[j].stimes[l];
                next_spike ++;
                next_feature ++;
            }
        }
    }


    printf(" >>> Printing sample offsets: ");
    for(i = 0; i<n_samples; i++){
        printf("%d ", tmp_gpu_input_dataset->sample_offset[i]);
    }
    printf("\n");

    printf(" >>> Printing feature offsets: ");
    for(i = 0; i<n_samples * n_features; i++){
        printf("%d ", tmp_gpu_input_dataset->feature_offset[i]);
    }
    printf("\n");

    printf(" >>> Printing dataset: ");
    for(i = 0; i<next_spike; i++){
        printf("%d ", tmp_gpu_input_dataset->spikes[i]);
    }
    printf("\n");
    // count total number of spikes
    /*int h_spk_ind = 0; // spike position index
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
    }*/

    printf(" > Dataset copied to new struct type!\n");
    fflush(stdout);
 
    // allocate to temporal
    gpu_input_dataset = (GPU_input_info_t *)malloc(sizeof(GPU_input_info_t));
    cudaMalloc(&(gpu_input_dataset->spikes), total_n_spikes * sizeof(int)); // allocate memory for neurons
    cudaMalloc(&(gpu_input_dataset->n_spikes), dataset->n_samples * dataset->image_size * sizeof(int)); // allocate memory for neurons
    cudaMalloc(&(gpu_input_dataset->sample_offset), dataset->n_samples * sizeof(int)); // allocate memory for neurons
    cudaMalloc(&(gpu_input_dataset->feature_offset), dataset->n_samples * dataset->image_size * sizeof(int)); // allocate memory for neurons

    // memcpy
    cudaMemcpy(gpu_input_dataset->spikes, tmp_gpu_input_dataset->spikes, total_n_spikes * sizeof(int), cudaMemcpyHostToDevice); // copy neurons information
    cudaMemcpy(gpu_input_dataset->n_spikes, tmp_gpu_input_dataset->n_spikes, dataset->n_samples * dataset->image_size * sizeof(int), cudaMemcpyHostToDevice); // copy neurons information
    cudaMemcpy(gpu_input_dataset->sample_offset, tmp_gpu_input_dataset->sample_offset, dataset->n_samples * sizeof(int), cudaMemcpyHostToDevice); // copy neurons information
    cudaMemcpy(gpu_input_dataset->feature_offset, tmp_gpu_input_dataset->feature_offset, dataset->n_samples * dataset->image_size * sizeof(int), cudaMemcpyHostToDevice); // copy neurons information
    gpu_input_dataset->n_samples = tmp_gpu_input_dataset->n_samples;
    gpu_input_dataset->type = tmp_gpu_input_dataset->type;
    gpu_input_dataset->n_classes = tmp_gpu_input_dataset->n_classes;
    gpu_input_dataset->n_features = tmp_gpu_input_dataset->n_features;

    // copy to GPU
    cudaMalloc(&(d_gpu_input_dataset), sizeof(GPU_input_info_t)); // allocate memory for neurons
    cudaMemcpy(d_gpu_input_dataset, gpu_input_dataset, sizeof(GPU_input_info_t), cudaMemcpyHostToDevice); // copy neurons information


    printf(" > Dataset copied to the GPU (feat %d)!\n", gpu_input_dataset->n_features);
    fflush(stdout);
 
    return d_gpu_input_dataset;
}


extern "C" GPU_results_t* copy_results_to_GPU(spiking_nn_t *snn, input_data_t *dataset, simulation_results_t *results){
    
    GPU_results_t *tmp_r;
    tmp_r = (GPU_results_t*)malloc(sizeof(GPU_results_t));
    
    cudaMalloc(&(tmp_r->nspk), snn->n_neurons * dataset->n_samples * sizeof(int)); // allocate memory for neurons
    cudaMalloc(&(tmp_r->nspk), snn->n_neurons * dataset->n_samples * 8 * sizeof(int)); // allocate memory for neurons

    GPU_results_t *d_results;
    cudaMalloc(&(d_results), sizeof(GPU_results_t)); // allocate memory for neurons
    cudaMemcpy(d_results, tmp_r, sizeof(GPU_results_t), cudaMemcpyHostToDevice); // copy neurons information

    return d_results;
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
