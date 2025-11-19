#include <cuda.h>
#include <cuda_runtime.h>
#include <stdio.h>

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

__constant__ int c_n_samples;

//__constant__ size_t L;

__global__ void reinit_spk_matrix(GPU_SNN_t *snn){
    
    //int L;
    
    //size_t threadsPerBlock = blockDim.x * blockDim.y * blockDim.z;
    //size_t threadNumInBlock = threadIdx.x + blockDim.x * threadIdx.y + blockDim.x * blockDim.y * threadIdx.z; //(alternatively: threadIdx.y + blockDim.y * threadIdx.x);
    //size_t blockNumInGrid   = blockIdx.x  + gridDim.x  * blockIdx.y + gridDim.x * gridDim.y * blockIdx.z;  //(alternatively: blockIdx.y  + gridDim.y  * blockIdx.x);

    //size_t globalThreadNum = blockNumInGrid * threadsPerBlock + threadNumInBlock;

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


        if(s < batch_last_sample && s < c_n_samples){

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
        snn->t_last_spike[global_neuron_index] = 0;
    }

}

/* OLD
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

*/

__global__ void reinit_synapses_kernel(GPU_SNN_t *snn, int s)
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

    if(threadId < snn->n_synapses){

        //continue;
        
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

        if(s < (size_t)(batch_last_sample) && (int)(s) < c_n_samples){

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

        snn->r_period_remain[global_n] --;

        int s = b + (int)network_number;

        if(s < batch_last_sample && s < c_n_samples){

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
                snn->t_last_spike[global_n] = t;
                snn->v[global_n] = snn->v_rest[local_n];

                // 
                size_t rs_index = (size_t)s * (size_t)snn->n_neurons + local_n;
                results->nspk[rs_index] += 1;
                //printf(" > Neuron %d, t %d, SPIKE\n", n, t);
            }
        }
    }
}


extern "C" void simulate_LIF_in_single_GPU(GPU_SNN_t *gpu_snn, GPU_dataset_t *gpu_dataset, GPU_results_t *gpu_results, simulation_configuration_t *conf, cuda_info_t *cuda_info, spiking_nn_t *cpu_snn, input_data_t *cpu_dataset){
    

    printf(" > Simulating LIF neurons network in GPU...\n");
    fflush(stdout);

    // call kernels
    int n_samples, time_steps, batch_size, n_networks, n_epochs;
    int t, s, b, snn, e;
    cudaError_t err;

    n_samples = cuda_info->n_samples;
    time_steps = cuda_info->time_steps;
    batch_size = cuda_info->batch_size;
    n_networks = cuda_info->n_networks;
    n_epochs = conf->epochs;

    // set dim3 struct to launch kernels
    dim3 grid_n_spkmatrix_reinit_blk((unsigned int)(cuda_info->n_threads_per_blk_rsm_x), cuda_info->n_threads_per_blk_rsm_y, cuda_info->n_threads_per_blk_rsm_z );
    dim3 block_n_spkmatrix_reinit_threads_per_blk((unsigned int)(cuda_info->n_threads_per_blk_rsm_x), cuda_info->n_threads_per_blk_rsm_y, cuda_info->n_threads_per_blk_rsm_z );    

    dim3 grid_n_ls_blk( (unsigned int)(cuda_info->n_blk_ls_x), cuda_info->n_blk_ls_y, cuda_info->n_blk_ls_z );
    dim3 block_n_ls_threads_per_block( (unsigned int)(cuda_info->n_threads_per_blk_ls_x), cuda_info->n_threads_per_blk_ls_y, cuda_info->n_threads_per_blk_ls_z );    

    dim3 grid_n_pn_blk( (unsigned int)(cuda_info->n_blk_nrs_x), cuda_info->n_blk_nrs_y, cuda_info->n_blk_nrs_z );
    dim3 block_n_pn_threads_per_block( (unsigned int)(cuda_info->n_threads_per_blk_nrs_x), cuda_info->n_threads_per_blk_nrs_y, cuda_info->n_threads_per_blk_nrs_z );  



    // compute epochs
    for(e = 0; e < n_epochs; e++){
    
        // loop over batches
        for(b = 0; b < n_samples; b = b + batch_size){

            for(snn = 0; snn < batch_size; snn = snn + n_networks){
                
                // reinit the entire matrix: number of elements = (n_neurons + n_input_neurons) * L
                // max 1024 threads per block: n_neurons * n_input_neurons threads, each L iterations
                // n_blocks = (n_neurons + n_input_neurons / 1024) + 1
                
                reinit_spk_matrix<<<grid_n_spkmatrix_reinit_blk, block_n_spkmatrix_reinit_threads_per_blk>>>(gpu_snn);
                err = cudaGetLastError();
                //printf("Launch error: %s\n", cudaGetErrorString(err));
                cudaCheckError(cudaPeekAtLastError());  // check launch errors
                cudaCheckError(cudaDeviceSynchronize());  // check runtime errors

                load_sample_kernel<<<grid_n_ls_blk, block_n_ls_threads_per_block>>>(gpu_snn, gpu_dataset, b + snn, n_samples);
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
                //reinit_synapses_kernel<<<2792, 516>>>(gpu_snn, s); // n_synapses 1440187
                
                for(t = 0; t<time_steps; t++){

                    // input step
                    lif_neuron_input<<<grid_n_pn_blk, block_n_pn_threads_per_block>>>(gpu_snn, b + snn, t, n_samples);
                    err = cudaGetLastError();
                    //printf("Launch error: %s\n", cudaGetErrorString(err));
                    cudaCheckError(cudaPeekAtLastError());  // check launch errors
                    cudaCheckError(cudaDeviceSynchronize());  // check runtime errors        
                    
                    // output step
                    lif_neuron_output<<<grid_n_pn_blk, block_n_pn_threads_per_block>>>(gpu_snn, b + snn, t, gpu_results, n_samples);
                    err = cudaGetLastError();
                    //printf("Launch error: %s\n", cudaGetErrorString(err));
                    cudaCheckError(cudaPeekAtLastError());  // check launch errors
                    cudaCheckError(cudaDeviceSynchronize());  // check runtime errors          
                    // learning
                    
                }
            }

            // update weights
        }
    }
}

extern "C" void simulate_LIF_in_multi_GPU(GPU_SNN_t *gpu_snn, GPU_dataset_t *gpu_dataset, GPU_results_t *gpu_results, simulation_configuration_t *conf, cuda_info_t *cuda_info, spiking_nn_t *cpu_snn, input_data_t *cpu_dataset, int dev, int batch){
    
    // call kernels
    int time_steps, batch_size, dev_batch_size, n_networks;
    int t, s, b, snn;
    cudaError_t err;

    time_steps = cuda_info->time_steps;
    batch_size = cuda_info->batch_size;
    dev_batch_size = cuda_info->n_batch_per_dev;
    n_networks = cuda_info->n_networks_per_dev;

    // set dim3 struct to launch kernels
    dim3 grid_n_spkmatrix_reinit_blk((unsigned int)(cuda_info->n_threads_per_blk_rsm_x), cuda_info->n_threads_per_blk_rsm_y, cuda_info->n_threads_per_blk_rsm_z );
    dim3 block_n_spkmatrix_reinit_threads_per_blk((unsigned int)(cuda_info->n_threads_per_blk_rsm_x), cuda_info->n_threads_per_blk_rsm_y, cuda_info->n_threads_per_blk_rsm_z );    

    dim3 grid_n_ls_blk( (unsigned int)(cuda_info->n_blk_ls_x), cuda_info->n_blk_ls_y, cuda_info->n_blk_ls_z );
    dim3 block_n_ls_threads_per_block( (unsigned int)(cuda_info->n_threads_per_blk_ls_x), cuda_info->n_threads_per_blk_ls_y, cuda_info->n_threads_per_blk_ls_z );    

    dim3 grid_n_pn_blk( (unsigned int)(cuda_info->n_blk_nrs_x), cuda_info->n_blk_nrs_y, cuda_info->n_blk_nrs_z );
    dim3 block_n_pn_threads_per_block( (unsigned int)(cuda_info->n_threads_per_blk_nrs_x), cuda_info->n_threads_per_blk_nrs_y, cuda_info->n_threads_per_blk_nrs_z );  

    
    // compute the sample index that will be simulated by this device
    int fsample = batch_size * batch + dev * dev_batch_size; 


    // simulate the batch using the available number of copies
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
        //reinit_synapses_kernel<<<2792, 516>>>(gpu_snn, s); // n_synapses 1440187
        
        for(t = 0; t<time_steps; t++){

            // input step
            lif_neuron_input<<<grid_n_pn_blk, block_n_pn_threads_per_block>>>(gpu_snn, fsample + snn, t, fsample + dev_batch_size);
            err = cudaGetLastError();
            //printf("Launch error: %s\n", cudaGetErrorString(err));
            cudaCheckError(cudaPeekAtLastError());  // check launch errors
            cudaCheckError(cudaDeviceSynchronize());  // check runtime errors        
            
            // output step
            lif_neuron_output<<<grid_n_pn_blk, block_n_pn_threads_per_block>>>(gpu_snn, fsample + snn, t, gpu_results, fsample + dev_batch_size);
            err = cudaGetLastError();
            //printf("Launch error: %s\n", cudaGetErrorString(err));
            cudaCheckError(cudaPeekAtLastError());  // check launch errors
            cudaCheckError(cudaDeviceSynchronize());  // check runtime errors          
            // learning
            
        }
    }
}


extern "C" void simulate_LIF_in_GPU(GPU_SNN_t **gpu_snn, GPU_dataset_t **gpu_dataset, GPU_results_t **gpu_results, simulation_configuration_t *conf, cuda_info_t *cuda_info, spiking_nn_t *cpu_snn, input_data_t *cpu_dataset){


    // copy constant(s) to gpu
    for(int i = 0; i<cuda_info->nDevices; i++){
        
        cudaSetDevice(i);
        cudaError_t err = cudaMemcpyToSymbol(c_n_samples,
                                            &cpu_dataset->n_samples,
                                            sizeof(int));
    }

    if(cuda_info->nDevices == 1){
        
        cudaSetDevice(0);
        simulate_LIF_in_single_GPU(gpu_snn[0], gpu_dataset[0], gpu_results[0], conf, cuda_info, cpu_snn, cpu_dataset);
    }
    // multi GPU
    else{

        int e, b, dev, batch_number = 0;

        // simulate epochs
        for(e=0; e<conf->epochs; e++){

            // simulate batches in GPUs
            int n_batches = cuda_info->n_samples / cuda_info->batch_size;

            if(cuda_info->n_samples % cuda_info->batch_size != 0)
                n_batches += 1;

            for(b=0; b<n_batches; b++){
     
                printf(" > In batch %d\n", b);
                fflush(stdout);
     
                // loop over devices (openMP threads)
                #pragma omp parallel for private(dev)
                for(dev = 0; dev<cuda_info->nDevices; dev++){

                    cudaSetDevice(dev);
                    simulate_LIF_in_multi_GPU(gpu_snn[dev], gpu_dataset[dev], gpu_results[dev], conf, cuda_info, cpu_snn, cpu_dataset, dev, b);
                }


                // wait until all finished and receive weights // TODO
            }
        }
    }
}
