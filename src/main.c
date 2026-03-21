#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <unistd.h>

#include "arceus.h"
/*#include "config/config_loader.h"
#include "networks/snn.h"
#include "simulations/simulations.h"
#include "simulations/results.h"
#include "datasets/datasets.h"


// include cuda files if defined
#ifdef CUDA
    #include "cuda/cuda_utils.cuh"
    #include "cuda/cuda_simulations_conf.h"
    #include "cuda/GPU_simulations.cuh"
#endif*/



/* main.c */
int main(int argc, char *argv[]) {

    // randomize execution
    srand(time(NULL));

    /*
    Load and initialize
    */
    printf(" ============================= \n Loading and initializing data \n ============================= \n");

    // load configuration parameters from input file
    printf(" > Loading configuration file...\n");
    simulation_configuration_t *conf = load_configuration_params_from_toml(argv[1]);
    printf(" > Configuration file loaded!\n\n");
    fflush(stdout);

    // initialize network
    printf(" > Initializing network...\n");
    GPU_SNN_t *cpu_snn = initialize_network_cpu(conf);
    printf(" > Network initialized!\n");
    fflush(stdout);

    //print_network(cpu_snn);

    // load dataset
    printf(" > Loading dataset... \n");
    GPU_dataset_t *cpu_dataset = load_dataset_from_file_cpu(conf->dataset, conf->labels, conf->n_samples, conf);
    printf(" > Dataset loaded!\n");
    fflush(stdout);

    // initialize results struct
    /*printf(" > Initializing results struct...\n");
    GPU_results_t *cpu_results = initialize_batch_results_cpu(conf, cpu_snn->n_neurons, conf->batch_size, 1, 1);
    printf(" > Results struct initialized!\n");
    fflush(stdout);


    // simulate in CPU
    if(conf->cuda == 0){

        printf("\n ============================= \n ==== Starting simulation ==== \n ============================= \n");
        simulate_batches(cpu_snn, cpu_dataset, conf, cpu_results);
    }
    // simulate in GPU if device is founded and CUDA is defined
    else{

        #if defined CUDA
        {
            // get GPU information
            cuda_info_t *cuda_info = getGPUProperties();

            // configure GPU simulation
            configure_cuda_simulation(cuda_info, cpu_snn, cpu_dataset, conf);

            // simulate
            cpu_results = simulate_batches_GPU(cuda_info, cpu_snn, cpu_dataset, conf);
        }
        // no device
        #else
        {
            printf(" > No cuda defice founded, simulating in CPU!\n");
            printf("\n ============================= \n ==== Starting simulation ==== \n ============================= \n");
            simulate_batches(cpu_snn, cpu_dataset, conf, cpu_results);
        }
        #endif
    }*/

    //deallocate_snn_str(cpu_snn);
    //deallocate_dataset_str(cpu_dataset);
    //deallocate_results_str(cpu_results);

    // define 
    //deallocate_memory();
    // called when simulation finishes


    // [CPU]

    size_t n_batches, r_samples;
    size_t b;

    // compute number of batches // [TODO]: improve or refactorize
    n_batches = cpu_dataset->n_samples / conf->batch_size;
    r_samples = cpu_dataset->n_samples % conf->batch_size;
    
    n_batches = r_samples > 0 ? n_batches + 1 : n_batches; // one more batch if there are remaining samples

    // copy non-constant snn data for parallel batch simulation
    /*init_batch_snn(cpu_snn, conf);

    // initialize struct to store batch results
    GPU_results_t **results = initialize_batch_results_array(conf, cpu_snn->n_neurons, conf->batch_size, 1, 1, n_batches);

    // loop over batches and simulate
    for(b = 0; b<n_batches; b++){
        
        // print for feedback
        if((b+1) % 100 == 0){
            printf(" Simulating batch %zu\n", b+1);
            fflush(stdout);
        }

        // simulate batch
        simulate_batch_CPU(cpu_snn, cpu_dataset, conf, results[b], b, 0);
    }

    store_number_of_spikes_array(results, conf, cpu_snn->n_neurons, conf->batch_size, n_batches);*/


    // [GPU]
    // init cuda_info
    cuda_info_t *cuda_info = configure_cuda_simulation(cpu_snn, cpu_dataset, conf);
    printf(" Cuda simulation configured\n");
    fflush(stdout);

    // move data to the GPU
    GPU_SNN_t **gpu_snn = cpy_SNN2GPU(cpu_snn, cuda_info, conf);
    printf(" SNN copied to the GPU\n");
    fflush(stdout);

    // move dataset to the GPU
    GPU_dataset_t **gpu_dataset = cpy_dataset2GPU(cpu_dataset, cuda_info);
    printf(" Dataset copied to the GPU\n");
    fflush(stdout);

    // initialize results structure
    GPU_results_t **results = initialize_batch_results_array(conf, cpu_snn->n_neurons, conf->batch_size, 1, 1, n_batches);
    printf(" Results array initialized\n");
    fflush(stdout);

    // call simulation
    for(b = 0; b<n_batches; b++){
        
        // print for feedback
        if((b+1) % 100 == 0){
            printf(" Simulating batch %zu\n", b+1);
            fflush(stdout);
        }

        // simulate batch
        simulate_batch_GPU(results[b], gpu_snn, gpu_dataset, conf, cuda_info, b);
    }

    printf(" Simulation finished\n");
    fflush(stdout);

    store_number_of_spikes_array(results, conf, cpu_snn->n_neurons, conf->batch_size, n_batches);


    return 0;
}
