#include <cuda.h>
#include <cuda_runtime.h>
#include <stdio.h>
#include "snn_library.h"
#include "neuron_models/GPU_lif_neuron.cuh"

#define THR_PER_BLOCK 1024 

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


void simulate_in_GPU(spiking_nn_t *snn, simulation_configuration_t *conf, simulation_results_t *results){

    // copy data to GPU
    copy_snn_structure_to_GPU(snn);

    // simulate
}


double copy_snn_structure_to_GPU(spiking_nn_t *snn){
    
    // event variables
    cudaEvent_t start, stop;
    
    // structs for snn, neurons and synapses
    spiking_nn_t *d_snn;
    lif_neuron_t *d_lif_neurons; 
    synapse_t *d_synapses;

    // control variables
    int i, j;

    // cuda things
    float milliseconds = 0;
    int thr_per_blk_neurons, blk_in_grid_neurons, thr_per_blk_synapses, blk_in_grid_synapses;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);

    // reserve memory for neurons and synapses lists
    //gpuErrchk(cudaMalloc(&d_lif_neurons, snn->n_neurons * sizeof(lif_neuron_t)));
    //gpuErrchk(cudaMalloc(&d_synapses, snn->n_synapses * sizeof(synapse_t)));



    /* Memory allocation and transfer */

    cudaMalloc(&d_snn, sizeof(spiking_nn_t)); // allocate memory for neurons
    cudaMemcpy(d_snn, snn, sizeof(spiking_nn_t), cudaMemcpyHostToDevice); // copy neurons information

    cudaMalloc(&d_lif_neurons, snn->n_neurons * sizeof(lif_neuron_t)); // allocate memory for neurons
    cudaMemcpy(d_lif_neurons, snn->lif_neurons, snn->n_neurons * sizeof(lif_neuron_t), cudaMemcpyHostToDevice); // copy neurons information

    cudaMalloc(&d_synapses, snn->n_synapses * sizeof(synapse_t)); // allocate memory for synapses
    cudaMemcpy(d_synapses, snn->synapses, snn->n_synapses * sizeof(synapse_t), cudaMemcpyHostToDevice); // copy synapses information



    d_snn->lif_neurons = d_lif_neurons;
    d_snn->synapses = d_synapses;
    // TODO: connect neuron initializer... if it is necessary


    // allocate memory for each neuron elements
    for(i = 0; i<snn->n_neurons; i++){
        //gpuErrchk(cudaMalloc(&d_lif_neurons[i].input_synapse_indexes, snn->lif_neurons[i].n_input_synapse * sizeof(int)));
        
        // allocate memory for input and outpyt synapse indexes, and copy
        cudaMalloc(&d_lif_neurons[i].input_synapse_indexes, snn->lif_neurons[i].n_input_synapse * sizeof(int));
        cudaMalloc(&d_lif_neurons[i].output_synapse_indexes, snn->lif_neurons[i].n_output_synapse * sizeof(int));
        cudaMemcpy(&d_lif_neurons[i].input_synapse, snn->lif_neurons[i].input_synapse, snn->lif_neurons[i].n_input_synapse * sizeof(int), cudaMemcpyHostToDevice);
        cudaMemcpy(&d_lif_neurons[i].output_synapse, snn->lif_neurons[i].output_synapse, snn->lif_neurons[i].n_output_synapse * sizeof(int), cudaMemcpyHostToDevice);

        // allocate memory for next spike indexes and initialize
        cudaMalloc(&d_lif_neurons[i].next_synapse_index, snn->lif_neurons[i].n_input_synapse * sizeof(int));     
        for(j = 0; j<lif_neurons[i].n_input_synapse; j++){
            d_lif_neurons[i].next_synapse_index[j] = 0;
        }  

        // allocate memory for spike times and initialize
        cudaMalloc(&d_lif_neurons[i].spike_times_arr, snn->lif_neurons[i].max_spikes * sizeof(int));
        for(j = 0; j<lif_neurons[i].max_spikes; j++){
            d_lif_neurons[i].spike_times_arr[j] = -1;
        }  
    }


    /* initialize references (pointers) */
    
    // connect neurons and synapses
    lif_neuron_t *d_neuron;
    synapse_t *d_synapse;
    for(i=0; i<snn->n_neurons; i++){

        d_neuron = &(d_neurons[i]);
        for(j=0; j<snn->lif_neurons[i].n_input_synapse){

            d_synapse = &(d_synapses[snn->lif_neurons[i].input_synapses[j]]);
            d_synapse->post_synaptic_lif_neuron = d_neuron;
        }

        for(j=0; j<snn->lif_neurons[i].n_output_synapse){

            d_synapse = &(d_synapses[snn->lif_neurons[i].output_synapses[j]]);
            d_synapse->pre_synaptic_lif_neuron = d_neuron;
        }
    }
    
    // set training rule for synapse // TODO: this must be refactorized to a function
    synapse->lr = lists->training_zones[synapse_id];

    for(i=0; i<snn.n_synapses; i++){
        synapse_t *d_synapse = &(d_synapses[i]);

        switch (d_synapse->lr) // get synapse training zone from list
        {
            case 0:
                d_synapse->learning_rule = &add_stdp;//(void (*)())&add_stdp;
                break;
            case 1:
                d_synapse->learning_rule = &mult_stdp;//(void (*)())&mult_stdp;
                break;
            case 2:
                d_synapse->learning_rule = &anti_stdp;//(void (*)())&anti_stdp;
                break;
            //case 3:
            //    synapse->learning_rule = &triplet_stdp;//(void (*)())&triplet_stdp;
            //    break;*/
            default:
                d_synapse->learning_rule = &add_stdp;//(void (*)())&add_stdp;
                break;
        }
    }
    


    /* Copy information from CPU to GPU */

    // copy information to gpu
    //gpuErrchk(cudaMemcpy(d_lif_neurons, snn->lif_neurons, snn->n_neurons * sizeof(lif_neuron_t), cudaMemcpyHostToDevice));
    //gpuErrchk(cudaMemcpy(d_synapses, snn->synapses, snn->n_synapses * sizeof(synapse_t), cudaMemcpyHostToDevice));

    // copy neurons information

    // copy synapses information

    // copy info of neurons (ONLY POINTERS; HOW IS THE REST OF INFORMATION PASSED?)
    for(int i = 0; i<snn->n_neurons; i++){
        cudaMemcpy(d_lif_neurons[i].input_synapse_indexes, snn->lif_neurons[i].input_synapse_indexes, snn->lif_neurons[i].n_input_synapse * sizeof(int), cudaMemcpyHostToDevice);
        cudaMemcpy(d_lif_neurons[i].output_synapse_indexes, snn->lif_neurons[i].output_synapse_indexes, snn->lif_neurons[i].n_output_synapse * sizeof(int), cudaMemcpyHostToDevice);
   }


    return 0.0;    
}


/**
GPUko memoriara mugitu matrizeak eta jaurti kernela
*/
double process_simulation_lif_neuron(spiking_nn_t *snn, int n, int m, int time_steps){
    printf("Running CUDA code, YUJUUUU\n");
    
    cudaEvent_t start, stop;
    
    // list of neurons and synapses
    lif_neuron_t *d_lif_neurons; 
    synapse_t *d_synapses;

    // cuda things
    float milliseconds = 0;
    int thr_per_blk_neurons, blk_in_grid_neurons, thr_per_blk_synapses, blk_in_grid_synapses;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);

    // reserve memory for neurons and synapses lists
    //gpuErrchk(cudaMalloc(&d_lif_neurons, snn->n_neurons * sizeof(lif_neuron_t)));
    //gpuErrchk(cudaMalloc(&d_synapses, snn->n_synapses * sizeof(synapse_t)));
    cudaMalloc(&d_lif_neurons, snn->n_neurons * sizeof(lif_neuron_t));
    cudaMalloc(&d_synapses, snn->n_synapses * sizeof(synapse_t));

    // reserve memory for each neuron synapse list indexes
    for(int i = 0; i<snn->n_neurons; i++){
        //gpuErrchk(cudaMalloc(&d_lif_neurons[i].input_synapse_indexes, snn->lif_neurons[i].n_input_synapse * sizeof(int)));
        cudaMalloc(&d_lif_neurons[i].input_synapse_indexes, snn->lif_neurons[i].n_input_synapse * sizeof(int));
        //gpuErrchk(cudaMalloc(&d_lif_neurons[i].output_synapse_indexes, snn->lif_neurons[i].n_output_synapse * sizeof(int)));
        cudaMalloc(&d_lif_neurons[i].output_synapse_indexes, snn->lif_neurons[i].n_output_synapse * sizeof(int));
    }

    // reserve memory for synapse pointers
    for(int i = 0; i<snn->n_synapses; i++){
        //gpuErrchk(cudaMalloc(&d_synapses[i].l_spike_times, snn->synapses[i].max_spikes * sizeof(int)));
        cudaMalloc(&d_synapses[i].l_spike_times, snn->synapses[i].max_spikes * sizeof(int));
        //gpuErrchk(cudaMalloc(&d_synapses[i].pre_synaptic_lif_neuron, sizeof(lif_neuron_t)));
        cudaMalloc(&d_synapses[i].pre_synaptic_lif_neuron, sizeof(lif_neuron_t));
        //gpuErrchk(cudaMalloc(&d_synapses[i].post_synaptic_lif_neuron, sizeof(lif_neuron_t)));
        cudaMalloc(&d_synapses[i].post_synaptic_lif_neuron, sizeof(lif_neuron_t));
        //d_synapse[i].learning_rule = snn->synapses[i].learning_rule;
    }


    // copy information to gpu
    //gpuErrchk(cudaMemcpy(d_lif_neurons, snn->lif_neurons, snn->n_neurons * sizeof(lif_neuron_t), cudaMemcpyHostToDevice));
    //gpuErrchk(cudaMemcpy(d_synapses, snn->synapses, snn->n_synapses * sizeof(synapse_t), cudaMemcpyHostToDevice));
    cudaMemcpy(d_lif_neurons, snn->lif_neurons, snn->n_neurons * sizeof(lif_neuron_t), cudaMemcpyHostToDevice);
    cudaMemcpy(d_synapses, snn->synapses, snn->n_synapses * sizeof(synapse_t), cudaMemcpyHostToDevice);

    // copy info of neurons (ONLY POINTERS; HOW IS THE REST OF INFORMATION PASSED?)
    for(int i = 0; i<snn->n_neurons; i++){
        cudaMemcpy(d_lif_neurons[i].input_synapse_indexes, snn->lif_neurons[i].input_synapse_indexes, snn->lif_neurons[i].n_input_synapse * sizeof(int), cudaMemcpyHostToDevice);
        cudaMemcpy(d_lif_neurons[i].output_synapse_indexes, snn->lif_neurons[i].output_synapse_indexes, snn->lif_neurons[i].n_output_synapse * sizeof(int), cudaMemcpyHostToDevice);
   }

    // reserve memory for synapse pointers
    for(int i = 0; i<snn->n_synapses; i++){
        cudaMemcpy(&d_synapses[i].l_spike_times, snn->synapses[i].l_spike_times, snn->synapses[i].max_spikes * sizeof(int), cudaMemcpyHostToDevice);
        cudaMemcpy(&d_synapses[i].pre_synaptic_lif_neuron, snn->synapses[i].pre_synaptic_lif_neuron, sizeof(lif_neuron_t), cudaMemcpyHostToDevice);
        cudaMemcpy(&d_synapses[i].post_synaptic_lif_neuron, snn->synapses[i].post_synaptic_lif_neuron, sizeof(lif_neuron_t), cudaMemcpyHostToDevice);
    }

    return 0.0;

    // grid for neurons kernel launching
    //thr_per_blk_neurons = colsBC; //hari bakoitzean emaitzeko matrizearen osagai bat kalkulatzen da (errenkada x zutabea)
    //blk_in_grid_neurons = rowsAC; //bloke bakoitzean A-ren errenkada bat

    // grid for synapses kernel launching
    //thr_per_blk_synapses = colsBC; //hari bakoitzean emaitzeko matrizearen osagai bat kalkulatzen da (errenkada x zutabea)
    //blk_in_grid_synapses = rowsAC;

    //blokeko gehienez 1024 hari
    /*if(colsBC > 1024)
    {
        int total = rowsAC * colsBC;
        blk_in_grid = total / 1024;
        thr_per_blk = 1024;
    }


    // simulation loop
    //for()
        // launch neuron input synapse kernel
        // launch neuron output synapse kernel
        // launch synapse learning

        // store information?
        
    // launch kernel
    gpuErrchk(cudaEventRecord(start));
    cuda_add_dot_matrix<<<blk_in_grid, thr_per_blk>>>(rowsAC, colsBC, colsArowsB, d_A, d_B, d_C, d_D);
    gpuErrchk(cudaEventRecord(stop));

   //Kopiatu D GPUko memoriatik CPUra
    cudaMemcpy(D, d_D, rowsAC * colsBC * sizeof(float), cudaMemcpyDeviceToHost);

    cudaEventSynchronize(stop);
    cudaEventElapsedTime(&milliseconds, start, stop);

    //Askatu GPUko memoria
    cudaFree(d_A);
    cudaFree(d_B);
    cudaFree(d_C);
    cudaFree(d_D);

    return(milliseconds);*/
}


/*double add_dot_matrix_GPU(int rowsAC, int colsBC, int colsArowsB, float* A, float* B, float* C, float* D)
{   
    cudaEvent_t start, stop;
    float *d_A, *d_B, *d_C, *d_D;
    float milliseconds = 0;
    int thr_per_blk, blk_in_grid;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);

    //erreserbatu memorian GPUan
    gpuErrchk(cudaMalloc(&d_A, rowsAC * colsArowsB * sizeof(float)));
    gpuErrchk(cudaMalloc(&d_B, colsArowsB * colsBC * sizeof(float)));
    gpuErrchk(cudaMalloc(&d_C, rowsAC * colsBC * sizeof(float)));
    gpuErrchk(cudaMalloc(&d_D, rowsAC * colsBC * sizeof(float)));

    //kopiatu A, B eta C matrizeak GPUko memorian
    gpuErrchk(cudaMemcpy(d_A, A, rowsAC * colsArowsB * sizeof(float), cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(d_B, B, colsArowsB * colsBC * sizeof(float), cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(d_C, C, rowsAC * colsBC * sizeof(float), cudaMemcpyHostToDevice));

    //Sareta eta blokeen egitura zehaztu
    thr_per_blk = colsBC; //hari bakoitzean emaitzeko matrizearen osagai bat kalkulatzen da (errenkada x zutabea)
    blk_in_grid = rowsAC; //bloke bakoitzean A-ren errenkada bat

    //blokeko gehienez 1024 hari
    if(colsBC > 1024)
    {
        int total = rowsAC * colsBC;
        blk_in_grid = total / 1024;
        thr_per_blk = 1024;
    }

    //jaurti kernela
    gpuErrchk(cudaEventRecord(start));
    cuda_add_dot_matrix<<<blk_in_grid, thr_per_blk>>>(rowsAC, colsBC, colsArowsB, d_A, d_B, d_C, d_D);
    gpuErrchk(cudaEventRecord(stop));

   //Kopiatu D GPUko memoriatik CPUra
    cudaMemcpy(D, d_D, rowsAC * colsBC * sizeof(float), cudaMemcpyDeviceToHost);

    cudaEventSynchronize(stop);
    cudaEventElapsedTime(&milliseconds, start, stop);

    //Askatu GPUko memoria
    cudaFree(d_A);
    cudaFree(d_B);
    cudaFree(d_C);
    cudaFree(d_D);

    return(milliseconds);
}*/

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
