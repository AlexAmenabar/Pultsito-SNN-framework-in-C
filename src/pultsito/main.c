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

/* main.c */
int main(int argc, char *argv[]) {

    // I think that too much structures are used, probably this should be refactorized
    spiking_nn_t snn; // base SNN structure and copies to process samples in parallel
    simulation_configuration_t conf; // simulation configuration data
    simulation_results_t train_results, test_results; // simulation results
    network_construction_lists_t lists; // structures to 
    input_data_t train_dataset, test_dataset; // train and test datasets


    // randomize execution
    srand(time(NULL));

    /*
    Load and initialize
    */

    printf(" ============================= \n Loading and initializing data \n ============================= \n");

    // load configuration parameters from input file
    printf(" > Loading configuration file...\n");
    load_configuration_params_from_toml(argv[1], &conf);
    printf(" > Configuration file loaded!\n\n");
    fflush(stdout);

    // load information about the snn from the network file //
    printf(" > Loading network data from file...\n");
    load_network_information(conf.network_file, &snn, &lists, &conf); // I don't like that this function loads some data into the SNN structure directly, it's confusing
    printf(" > Network data loaded!\n\n");
    fflush(stdout);

    // initialize the network
    printf(" > Initializing network...\n"); 
    initialize_network(&snn, &conf, &lists);
    printf(" > Network initialized!\n\n");
    fflush(stdout);

    // load input spike train from file (different depending on execution type) // ESTO DEBERÍA CAMBIARLO; NO ME TERMINA DE GUSTAR COMO ESTÁ PASANDO DIRECTAMENTE EL PARÁMETRO DE ENTRADA
    printf(" > Loading dataset...\n");
    load_dataset_from_file(&train_dataset, conf.train_set, conf.train_labels, conf.n_train, &conf);
    if(conf.test_provided == 1)
        load_dataset_from_file(&test_dataset, conf.test_set, conf.test_labels, conf.n_test, &conf);
    printf(" > Dataset loaded!\n");


    // initialize struct to store results // REVISE THIS
    printf(" > Initializing results struct...\n");
    initialize_results_struct(&train_results, &conf, train_dataset.n_samples, snn.n_neurons);
    if(conf.test_provided == 1)
        initialize_results_struct(&test_results, &conf, test_dataset.n_samples, snn.n_neurons);

    printf(" > Results strcut initialized!\n");
    fflush(stdout);


#ifdef REORDER

    printf(" > Reordering array of synapses...\n");

    #ifdef CUDA
    // if cuda is used, initialized delays to 0
    if(conf.cuda != 0)
        for(int i = 0; i<snn.n_input_synapses; i++){
            snn.synapses[i].delay = 0;
        }
    #endif

    // reorder array
    reorder_synapse_list(&snn);

    printf(" > Array of synapses reordered!\n");
    fflush(stdout);

#endif
 
    printf("\n ============================= \n ==== Starting simulation ==== \n ============================= \n");

    //printf(" > Simulation properties:\n");
    //printf(" > > ")

    fflush(stdout);

    // train the network
    if(conf.cuda == 0){

        // biological simulation
        if(conf.simulation_obj == 0){

            simulate_samples(&snn, &conf, &train_results, &train_dataset);
        }
        // ML
        else{
            if(conf.train_provided == 1)
                train_network(&snn, &conf, &train_results, &train_dataset);

            // test the network
            if(conf.test_provided == 1)
                test_network(&snn, &conf, &test_results, &test_dataset);

            // another function to integrate train and test?
        }
    }
    else{
        #ifdef CUDA
        
        double elpt;
        struct timespec start, end; // to measure simulation complete time
        struct timespec start_cpy, end_cpy;

        // start timing
        clock_gettime(CLOCK_MONOTONIC, &start);

        // get gpu properties
        cuda_info_t *cuda_info = getProperties();
        // compute_cuda_sim_conf(cuda_info);


        // map data structures to GPU structures in CPU
        GPU_SNN_t *gpu_snn = SNN_CPU2GPU_mapping(&snn, &conf);
        GPU_dataset_t *gpu_dataset = dataset_CPU2GPU_mapping(&train_dataset, &conf);

        // configure cuda simulation
        configure_cuda_simulation(cuda_info, gpu_snn, gpu_dataset, &conf);

        // copy data structures to GPU memory
        clock_gettime(CLOCK_MONOTONIC, &start_cpy);        
        GPU_SNN_t **d_gpu_snn = cpy_SNN2GPU(gpu_snn, cuda_info); // structure in GPU
        GPU_dataset_t **d_gpu_dataset = cpy_dataset2GPU(gpu_dataset, cuda_info);
        clock_gettime(CLOCK_MONOTONIC, &end_cpy);
        elpt = (end_cpy.tv_sec - start_cpy.tv_sec) + (end_cpy.tv_nsec - start_cpy.tv_nsec) / 1e9;
        printf(" Elapsed time copying data: %lf\n", elpt);

        // free memory of GPU structs allocated in CPU memory
        free_gpu_snn_in_CPU(gpu_snn);
        free_gpu_dataset_in_CPU(gpu_dataset);
     

        // perform the simulation
        simulate_in_GPU(d_gpu_snn, d_gpu_dataset, &conf, cuda_info, &snn, &train_dataset); // refactorize


        // end timing
        clock_gettime(CLOCK_MONOTONIC, &end);
        elpt = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;
    
        printf(" Elapsed total time: %lf\n", elpt);

        #endif
    }

    
    // free memory

    return 0;
}
