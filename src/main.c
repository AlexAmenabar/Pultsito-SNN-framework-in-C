#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <unistd.h>

#include "config/config_loader.h"
#include "networks/snn.h"
#include "simulations/simulations.h"
#include "simulations/results.h"
#include "datasets/datasets.h"


// include cuda files if defined
#ifdef CUDA
    #include "cuda/cuda_utils.cuh"
    #include "cuda/cuda_simulations_conf.h"
    #include "cuda/GPU_simulations.cuh"
#endif



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
    printf(" > Initializing results struct...\n");
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
    }

    return 0;
}
