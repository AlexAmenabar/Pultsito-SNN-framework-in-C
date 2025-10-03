#include "snn_library.h"
#include "load_data.h"
#include "helpers.h"
#include "training_rules/stdp.h"

#include "neuron_models/lif_neuron.h"
//#include "neuron_models/GPU_lif_neuron.cuh"

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
    spiking_nn_t snn; // SNN structure
    simulation_configuration_t conf; // simulation configuration data
    results_configuration_t results_conf; // configuration for results after simulation
    simulation_results_t results; // simulation results
    network_construction_lists_t lists; // structures to 

    // randomize execution
    srand(time(NULL));

    // load configuration parameters from input file
    //load_configuration_params(argv[1], &conf);
    printf(" > Loading simulation configuration...\n");
    load_configuration_params_from_toml(argv[1], &conf);
    printf(" > Simulation configuration loaded!\n\n");

    // load information about the snn from the network file
    printf(" > Loading network data...\n");
    load_network_information(conf.network_file, &snn, &lists, &conf); // I don't like that this function loads some data into the SNN structure directly, it's confusing
    printf(" > Network data loaded!\n\n");


    // initialize the network
    printf(" > Initializing network...\n"); 
    initialize_network(&snn, &conf, &lists);
    printf(" > Network initialized!\n\n");


    // load input spike train from file (different depending on execution type) // ESTO DEBERÍA CAMBIARLO; NO ME TERMINA DE GUSTAR COMO ESTÁ PASANDO DIRECTAMENTE EL PARÁMETRO DE ENTRADA
    printf(" > Loading input spike trains...\n");
    load_input_spike_trains_on_snn(conf.input_spikes_file, &snn);
    printf(" > Spike trains loaded!\n");

    // initialize struct to store results
    printf(" > Initializing results struct...\n");
    results_conf.n_neurons = snn.n_neurons;
    results_conf.n_samples = conf.n_samples;
    results_conf.time_steps = conf.time_steps;
    initialize_results_struct(&results, &results_conf);
    printf(" > Results strcut initialized!\n");


#ifdef REORDER
    printf(" > Reordering synapses list...\n");
    reorder_synapse_list(&snn);
    printf(" > List of synapses reordered!\n");
#endif
 

    // free memory TODO: must be corrected
    //free_lists_memory(&lists, &snn);

    printf("Initializing training / simulation\n");

    // Run the simulation

#ifndef BY_SAMPLE
    int reps = 1;
    for(i=0; i<reps; i++){
        simulate(&snn, &conf, &results);
    }
    // compute means of execution times
    results.results_per_sample[0].elapsed_time = results.results_per_sample[0].elapsed_time / reps;
    results.results_per_sample[0].elapsed_time_neurons = results.results_per_sample[0].elapsed_time_neurons / reps;
    results.results_per_sample[0].elapsed_time_neurons_input = results.results_per_sample[0].elapsed_time_neurons_input / reps;
    results.results_per_sample[0].elapsed_time_neurons_output = results.results_per_sample[0].elapsed_time_neurons_output / reps;
    results.results_per_sample[0].elapsed_time_synapses = results.results_per_sample[0].elapsed_time_synapses / reps;
    results.results_per_sample[0].elapsed_time_synapses_input = results.results_per_sample[0].elapsed_time_synapses_input / reps;
    results.results_per_sample[0].elapsed_time_synapses_output = results.results_per_sample[0].elapsed_time_synapses_output / reps;
    results.results_per_sample[0].elapsed_time_learning = results.results_per_sample[0].elapsed_time_learning / reps; 
#else
    simulate_by_samples();
#endif


    // store results (fnal network, execution times...)
    store_results(&results, &conf, &snn);

    // free memory
    // TODO


    return 0;
}
