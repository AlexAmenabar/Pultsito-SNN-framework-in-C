#include <cuda.h>
#include <cuda_runtime.h>
#include <stdio.h>

#include "snn_library.h"
#include "neuron_models/GPU_lif_neuron.cuh"

#define cudaCheckError(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
{
    if (code != cudaSuccess)
    {
        fprintf(stderr, "CUDA Error: %s %s %d\n", cudaGetErrorString(code), file, line);
        if (abort) exit(code);
    }
}


__global__ void reinit_spk_matrix(GPU_SNN_t *snn){
    
    int i, L;
    int threadsPerBlock = blockDim.x * blockDim.y;
    int threadNumInBlock = threadIdx.x + blockDim.x * threadIdx.y; //(alternatively: threadIdx.y + blockDim.y * threadIdx.x);
    int blockNumInGrid   = blockIdx.x  + gridDim.x  * blockIdx.y;  //(alternatively: blockIdx.y  + gridDim.y  * blockIdx.x);

    int globalThreadNum = blockNumInGrid * threadsPerBlock + threadNumInBlock;
    int matrix_index = globalThreadNum; // ?

    L = snn->max_spikes;

    // reinit in all network copies
    if(matrix_index < (snn->n_neurons + snn->n_input_neurons) * L * snn->n_networks){
        snn->spk_matrix[matrix_index] = -1;
    }
}


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

        // get the copy number
        int network_number = (int)(neuron_index / snn->n_input_neurons);

        // scale neuron identifier === network number * first network neuron index + neuron local index
        int local_neuron_index = (neuron_index % snn->n_input_neurons);
        int global_neuron_index = network_number * (snn->n_input_neurons + snn->n_neurons) + local_neuron_index;

        // get sample in batch
        s = s + network_number; // the i.th network computes the s + i sample

        // get the feature
        int f = local_neuron_index; // the neuron index is also the feature index

        //printf(" > > Feature = %d\n", f);

        int start_sample = dataset->sample_offset[s];
        //printf(" Thread %d, start sample = %d\n", f, start_sample);
        
        int start_feature = dataset->feature_offset[s * dataset->n_features + f];
        //printf(" Thread %d, start feature = %d\n", f, start_feature);

        int n_spikes_per_feature = dataset->n_spikes_per_feature[s * dataset->n_features + f];

        //printf(" > Loading spike train in neuron = %d (n = %d)\n", neuron_index, n_spikes);

        for(i = 0; i<n_spikes_per_feature; i++){
        
            // store a 1 value if there is a spike
            snn->spk_matrix[global_neuron_index * L + dataset->spikes[start_sample + start_feature + i]] = dataset->spikes[start_sample + start_feature + i];

            //printf(" SPike added by thread %d to %d (value %d).\n", f, neuron_index * L + dataset->spikes[start_sample + start_feature + i], dataset->spikes[start_sample + start_feature + i]);
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

    int L;

    int threadsPerBlock = blockDim.x * blockDim.y;
    int threadNumInBlock = threadIdx.x + blockDim.x * threadIdx.y; //(alternatively: threadIdx.y + blockDim.y * threadIdx.x);
    int blockNumInGrid   = blockIdx.x  + gridDim.x  * blockIdx.y;  //(alternatively: blockIdx.y  + gridDim.y  * blockIdx.x);

    int n = blockNumInGrid * threadsPerBlock + threadNumInBlock;
    
    L = snn->max_spikes;

    if(n < snn->n_neurons){

        float I = 0.0; // in thread local memory?

        // loop over input synapses to check if the index should be updated in one of them
        if(snn->r_period_remain[n] <= 0){
            
            int synapse_index, in_neuron_index, latency, spk_time, spk_time_index, spk, i;


            for(i = 0; i<snn->n_neuron_input_synapses[n]; i++){

                synapse_index = snn->neuron_input_synapses_offset[n] + i;
                in_neuron_index = snn->pre_neuron_index[synapse_index];
                latency = snn->delay[synapse_index];

                spk_time = t - latency; // actual position

                if(spk_time >=0 ){

                    spk_time_index = spk_time % L; // scale to circular array
                    spk = snn->spk_matrix[in_neuron_index * L + spk_time_index];
                    
                    if(spk == spk_time){

                        I += snn->w[synapse_index];
                    }
                    //if(n < snn->n_input_neurons)
                        //printf(" Neuron %d: input synapse %d (l=%d, in neuron=%d), spk_time=%d, spk_time_index=%d, spk=%d\n", n, synapse_index, latency, in_neuron_index, spk_time, spk_time_index, spk);

                }    
          
            }
            
            snn->v[n] = (1-0.05) * snn->v[n] + snn->v_rest[n] * 0.05 + I;
            //if(n < snn->n_input_neurons)
                //printf(" > Neuron %d, V = %f, I = %f\n", n, snn->v[n], I);
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

    int L = snn->max_spikes;

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
            //printf(" > Neuron %d, t %d, SPIKE\n", n, t);
        }
    }

}


extern "C" void simulate_LIF_in_GPU(GPU_SNN_t *gpu_snn, GPU_dataset_t *gpu_dataset, GPU_results_t *gpu_results, simulation_configuration_t *conf, cuda_info_t *cuda_info, spiking_nn_t *cpu_snn, input_data_t *cpu_dataset){
    

    printf(" > Simulating LIF neurons network in GPU...\n");
    fflush(stdout);

    // call kernels
    int t, s, n_samples, time_steps;

    // these numbers are only temporal
    int n_spkmatrix_reinit_blk = ((cpu_snn->n_neurons + cpu_snn->n_input) * conf->max_input_spikes * cuda_info->n_networks) / 1024 + 1;
    int n_spkmatrix_reinit_threads_per_blk = 1024;
    
    int n_ls_threads_per_block = 1024;
    int n_ls_blk = (cpu_snn->n_input * cuda_info->n_networks) / 1024 + 1;

    int n_pn_threads_per_block = 516;
    int n_pn_blk = (cpu_snn->n_neurons * cuda_info->n_networks) / 516 + 1;

    n_samples = cpu_dataset->n_samples;
    time_steps = conf->time_steps;


    // loop over batches


    for(s = 0; s<n_samples; s++){

        // reinit the entire matrix: number of elements = (n_neurons + n_input_neurons) * L
        // max 1024 threads per block: n_neurons * n_input_neurons threads, each L iterations
        // n_blocks = (n_neurons + n_input_neurons / 1024) + 1
        reinit_spk_matrix<<<n_spkmatrix_reinit_blk, n_spkmatrix_reinit_threads_per_blk>>>(gpu_snn);
        cudaCheckError(cudaPeekAtLastError());  // check launch errors
        cudaCheckError(cudaDeviceSynchronize());  // check runtime errors

        load_sample_kernel<<<n_ls_blk, n_ls_threads_per_block>>>(gpu_snn, gpu_dataset, s);
        cudaCheckError(cudaPeekAtLastError());  // check launch errors
        cudaCheckError(cudaDeviceSynchronize());  // check runtime errors

        // reinit neurons
        reinit_neurons_kernel<<<n_pn_blk, n_pn_threads_per_block>>>(gpu_snn, s); // n_neurons 8784
        cudaCheckError(cudaPeekAtLastError());  // check launch errors
        cudaCheckError(cudaDeviceSynchronize());  // check runtime errors
        
        // reinit synapses
        //reinit_synapses_kernel<<<2792, 516>>>(gpu_snn, s); // n_synapses 1440187
        

        for(t = 0; t<time_steps; t++){


            // input step
            lif_neuron_input<<<n_pn_blk, n_pn_threads_per_block>>>(gpu_snn, s, t);
            cudaCheckError(cudaPeekAtLastError());  // check launch errors
            cudaCheckError(cudaDeviceSynchronize());  // check runtime errors        
            
            // output step
            lif_neuron_output<<<n_pn_blk, n_pn_threads_per_block>>>(gpu_snn, s, t, gpu_results);
            cudaCheckError(cudaPeekAtLastError());  // check launch errors
            cudaCheckError(cudaDeviceSynchronize());  // check runtime errors          
            // learning
            
        }
    }
}