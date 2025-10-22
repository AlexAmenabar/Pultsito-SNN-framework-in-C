#include "snn_library.h"
#include "load_data.h"
#include "helpers.h"
#include "training_rules/stdp.h"

#include "neuron_models/lif_neuron.h"
#include "neuron_models/GPU_lif_neuron.cuh"

#include "simulations/simulations.h"


#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>



/* main.c */
int main(int argc, char *argv[]) {
    // variables definition
    int i, j;

    // I think that too much structures are used, probably this should be refactorized
    spiking_nn_t snn; // base SNN structure and copies to process samples in parallel
    simulation_configuration_t conf; // simulation configuration data
    simulation_results_t train_results, test_results; // simulation results
    network_construction_lists_t lists; // structures to 
    input_data_t train_dataset, test_dataset; // train and test datasets

    // randomize execution
    srand(time(NULL));

    // load configuration parameters from input file
    //load_configuration_params(argv[1], &conf);
    printf(" > Loading simulation configuration file information...\n");
    load_configuration_params_from_toml(argv[1], &conf);
    printf(" > Simulation configuration file loaded!\n\n");

    // load information about the snn from the network file //
    printf(" > Loading network data...\n");
    load_network_information(conf.network_file, &snn, &lists, &conf); // I don't like that this function loads some data into the SNN structure directly, it's confusing
    printf(" > Network data loaded!\n\n");


    // initialize the network
    printf(" > Initializing network...\n"); 
    initialize_network(&snn, &conf, &lists);
    printf(" > Network initialized!\n\n");


    // load input spike train from file (different depending on execution type) // ESTO DEBERÍA CAMBIARLO; NO ME TERMINA DE GUSTAR COMO ESTÁ PASANDO DIRECTAMENTE EL PARÁMETRO DE ENTRADA
    printf(" > Loading datasets...\n");
    load_dataset_from_file(&train_dataset, conf.train_set, conf.train_labels, conf.n_train, &conf);
    if(conf.test_provided == 1)
        load_dataset_from_file(&test_dataset, conf.test_set, conf.test_labels, conf.n_test, &conf);
    

    printf(" > Datasets loaded!\n");
    /*printf(" > Loading input spike trains...\n");
    load_input_spike_trains_on_snn(conf.input_spikes_file, &snn);
    printf(" > Spike trains loaded!\n");*/

    // initialize struct to store results
    printf(" > Initializing results struct...\n");
    
    // TODO: this is too much, probably one struct for each process and a file for each sample? Epoch?
    initialize_results_struct(&train_results, &conf, train_dataset.n_samples, snn.n_neurons);
    if(conf.test_provided == 1)
        initialize_results_struct(&test_results, &conf, test_dataset.n_samples, snn.n_neurons);

    printf(" > Results strcut initialized!\n");


#ifdef REORDER
    printf(" > Reordering synapses list...\n");
    reorder_synapse_list(&snn);
    printf(" > List of synapses reordered!\n");
#endif
 

    // copy the network n times // TODO: this can be paralelized
/*#ifdef PAR_SAMPLES 
    printf(" > Copying network %d times...\n", conf.n_process);
    pr_snns = (spiking_nn_t *)malloc(conf.n_process * sizeof(spiking_nn_t));
    for(i = 0; i<conf.n_process; i++){
        cp_network(&(pr_snns[i]), &snn, &conf);
    }
    printf(" > Network copied!\n");
    fflush(stdout);
#endif*/

    printf("Initializing training / simulation\n");



    // train the network
    train_network(&snn, &conf, &train_results, &train_dataset);

    // test the network
    if(conf.test_provided == 1)
        test_network(&snn, &conf, &test_results, &test_dataset);

    // store results
    store_results(&train_results, &conf, &snn, &train_dataset);

    // free memory


    return 0;
}
