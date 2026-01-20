#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

#include "snn_library.h"
#include "load_data.h"

#include "simulations/simulations.h"

#include "neuron_models/lif_neuron.h"

#include "training_rules/stdp.h"

#include "helpers.h"

#ifdef CUDA
    #include "cuda/GPU_simulations.cuh"

    #include "neuron_models/GPU_lif_neuron.cuh"
    #include "cuda/cuda_helpers.cuh"
#endif

#include <immintrin.h>  // AVX intrinsics
#include <omp.h>


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

    // load network initialization arrays
    printf(" > Loading network information in lists...\n");
    network_construction_lists_t *data = load_network_information_in_lists(conf);
    printf(" > Network information loaded!\n\n");
    fflush(stdout);

    // initialize network
    printf(" > Initializing network...\n");
    GPU_SNN_t *cpu_snn = initialize_network_cpu(conf, data);
    printf(" > Network initialized!\n");
    fflush(stdout);

    // load dataset
    printf(" > Loading dataset... \n");
    GPU_dataset_t *cpu_dataset = load_dataset_from_file_cpu(conf->train_set, conf->train_labels, conf->n_train, conf);
    printf(" > Dataset loaded!\n");
    fflush(stdout);

    // initialize results struct
    printf(" > Initializing results struct...\n");
    GPU_results_t *cpu_results = initialize_batch_results_cpu(conf, cpu_snn->n_neurons, conf->batch_size, 1);
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
        #else
        {
            printf(" > No cuda defice founded, simulating in CPU!\n");
            simulate_batches(cpu_snn, cpu_dataset, conf, cpu_results);
        }
        #endif
    }

    return 0;
}
