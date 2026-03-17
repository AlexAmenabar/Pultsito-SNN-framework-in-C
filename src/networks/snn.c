#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>

#include "toml_c/toml.h"

#include "config/config_loader.h"
#include "networks/snn.h"
#include "networks/snn_generator.h"
#include "utils.h"

// internal
#include "priv_snn.h"
#include "neuron_models/priv_neuron_models.h"


/* FILE private functions */

/// @brief Compute the input synapses offsets for each neuron
/// @param snn SNN structure
/// @param topology Structure containing helper values
void count_neurons_input_synapses_and_conpute_offsets(GPU_SNN_t *snn, topology_t *topology){
    
    size_t i, j, off = 0;

    // loop over neurons and count number of input synapses
    for(i = 0; i<snn->n_neurons; i++){

        // loop over input synapses of the neuron
        for(j = 0; j<(size_t)topology->input_neurons_per_neuron[i][0]; j++){

            // set the number of input synapses
            // [i * 2 + 1] is the index of the neuron
            // [i * 2 + 2] the number of synapses from that neuron to the current 
            snn->n_neuron_input_synapses[i] += (size_t)topology->input_neurons_per_neuron[i][j*2+2];
        }
    }

    // loop over neurons and compute offsets
    for(i = 0; i<snn->n_neurons; i++){

        // set offset
        snn->neuron_input_synapses_offset[i] = off;

        // update offset
        off += snn->n_neuron_input_synapses[i];
    }
}

/// @brief Set the input and output neurons of the synapses
/// @param snn SNN structure
/// @param topology Structure containing helper values
void connect_synapses_to_network(GPU_SNN_t *snn, topology_t *topology){

    size_t i, j, l, next_syn = 0;

    // loop over the neurons
    for(i = 0; i<snn->n_neurons; i++){

        // loop over the input synapses of the neuron
        for(j = 0; j<topology->input_neurons_per_neuron[i][0]; j++){

            // set the input and output neurons of the synapse
            for(l = 0; l<topology->input_neurons_per_neuron[i][j * 2 + 2]; l++){
                
                snn->post_neuron_index[next_syn] = i + snn->n_input_neurons; // first [n_input_neurons] neurons are virtual input
                snn->pre_neuron_index[next_syn] = topology->input_neurons_per_neuron[i][j * 2 + 1];
                next_syn ++;
            }
        }
    }
}



/* FUNCTIONS */

/// @brief Load the information of the neuron from the input files to a intermediate structure.
/// @param conf Configuration structure with information about the network
/// @return Structure with the arrays storing the data of the neurons, synapses and connectivity to initialize the network
topology_t* load_network_information_in_topology_from_file(simulation_configuration_t *conf){

    FILE *f = NULL, *f_neurons, *f_synapses;
    char errbuf[100];
    size_t i, j;

    // struct to store arrays
    //network_construction_lists_t *lists = (network_construction_lists_t*)calloc(1, sizeof(network_construction_lists_t));
    topology_t *topology = (topology_t*)calloc(1, sizeof(topology_t)); // structure to store all the network data
    neurons_t neurons; // structure to store synapses data
    synapses_t synapses; // structure to store synapses data

    // define table and parameters variables
    toml_table_t *tbl, *tbl_general, *tbl_neurons, *tbl_synapses;

    toml_array_t *latency_lst, *weights_lst, *training_zones_lst, *connection_lst_lst, *connection_lst;

    toml_value_t n_neurons, n_input_neurons, n_output_neurons, n_synapses, latency, 
                    training_zone, weight, n_connections, neuron_id, n_synapses_to_neuron;


    // open file and load in TOML struct
    open_file(&f, conf->network_file);
    tbl = toml_parse_file(f, errbuf, 100);
    fclose(f);
    

    /* get sections from file */
    tbl_general = toml_table_table(tbl, "general");
    tbl_neurons = toml_table_table(tbl, "neurons");
    tbl_synapses = toml_table_table(tbl, "synapsis");
    
    /* load [General] section */
    n_neurons = toml_table_int(tbl_general, "neurons");
    n_input_neurons = toml_table_int(tbl_general, "input_neurons");
    n_output_neurons = toml_table_int(tbl_general, "output_neurons");
    n_synapses = toml_table_int(tbl_general, "synapsis");

    size_t N, iN, oN, S;
    N = (size_t)n_neurons.u.i;
    iN = (size_t)n_input_neurons.u.i;
    oN = (size_t)n_output_neurons.u.i;
    S = (size_t)n_synapses.u.i;
    
    // load information in the topology structure
    topology->n_neurons = N;
    topology->n_input = iN;
    topology->n_output_neurons = oN;
    topology->n_synapses = S;

    // check that all the information has been loaded correctly
    if(!(n_neurons.ok && n_input_neurons.ok && n_output_neurons.ok && n_synapses.ok)){
        printf("The number of neurons, input neurons, output neurons, synapses, input synapses and output synapses must be provided in the network file!");
        fflush(stdout);
        exit(1);
    }
    
    // allocate memory for arrays related to neurons data
    switch (conf->neuron_type){  
        case 1:
            neurons = load_LIF_neurons_from_file(topology, conf, tbl_neurons);
            break;
        default:
            neurons = load_LIF_neurons_from_file(topology, conf, tbl_neurons);
            break;
    }
        
    
    /* Synapses section */
    open_file(&f_synapses, conf->network_synapses_file); // TOML file

    // allocate memory for arrays related to synapses data
    synapses.w = (float *)malloc(S * sizeof(float)); // weights
    synapses.delay = (int *)malloc(S * sizeof(int)); // delays
    synapses.lr = (int *)malloc(S * sizeof(int)); // learnin rules

    // load whether different parameters are included in the input file
    latency = toml_table_int(tbl_synapses, "delay");
    weight = toml_table_int(tbl_synapses, "weight");
    training_zone = toml_table_int(tbl_synapses, "training_zone");

    int tmp_int;
    float tmp;

    // load delays and compute the maximum value
    if(latency.ok && latency.u.i == 1){
    
        for(i=0; i<S; i++)
            fscanf(f_synapses, "%d", &(synapses.delay[i]));
    }
    else if(latency.ok && latency.u.i == 2){

        fscanf(f_synapses, "%d", &(tmp_int));
        for(i = 0; i<S; i++)
            synapses.delay[i] = tmp_int;
    }
    else{
        for(i=0; i<S; i++)
            synapses.delay[i] = 1;
    }

    // load weights
    if(weight.ok && weight.u.i == 1){
        for(i=0; i<S; i++)
            fscanf(f_synapses, "%f", &(synapses.w[i]));
    }
    else if(weight.ok && weight.u.i == 2){
        fscanf(f_synapses, "%f", &(tmp));
        for(i = 0; i<S; i++)
            synapses.w[i] = tmp;
    }
    else{
        for(i = 0; i<S; i++)
            synapses.w[i] = 0.25;
    }

    // load training zones
    if(training_zone.ok && training_zone.u.i == 1){
        for(i=0; i<S; i++){
            fscanf(f_synapses, "%d", &(synapses.lr[i]));
        }
    }
    else if(training_zone.ok && training_zone.u.i == 2){
        fscanf(f_synapses, "%d", &(tmp_int));
        for(i = 0; i<S; i++)
            synapses.lr[i] = tmp_int;
    }
    else{
        for(i = 0; i<S; i++)
            synapses.lr[i] = 0;
    }

    printf(" >> Synapses loaded from file!\n");
    fflush(stdout);


    /* Connectivity */

    // allocate memory for connectibity
    topology->input_neurons_per_neuron = (size_t **)malloc(N * sizeof(size_t *)); // +1, since the input layer is stored too

    int number_connections;
    for(i=0; i<N; i++){ // network input synapses are loaded first and output synapses last
        fscanf(f_synapses, "%d", &number_connections);

        // alloc memory
        (topology->input_neurons_per_neuron)[i] = (size_t*)malloc((number_connections * 2 + 1) * sizeof(size_t)); // for each connection the neuron id and the number of synapses must be stored
        (topology->input_neurons_per_neuron)[i][0] = (size_t)(number_connections);

        for(j = 0; j<number_connections; j++){
            fscanf(f_synapses, "%zu", &((topology->input_neurons_per_neuron)[i][j * 2 + 1])); // number of synapses connected to that neuron
            fscanf(f_synapses, "%zu", &((topology->input_neurons_per_neuron)[i][j * 2 + 2])); // number of synapses connected to that neuron
        }
    }
    fclose(f_synapses);

    printf(" >> Synapses section loaded\n");
    fflush(stdout);


    // link neurons and synapses data for the topology
    topology->neurons = neurons;
    topology->synapses = synapses;

    return topology;
}

GPU_SNN_t* initialize_network_cpu(simulation_configuration_t *conf){

    // load network information into intermedaite arrays
    topology_t *topology; 

    // store information of network, neurons and synapses in an intermediate structure
    if(conf->load_network == 0){ // load from file
        topology = load_network_information_in_topology_from_file(conf);
    }
    else if(conf->load_network == 1){ // generate
        
        // TODO
    }
    else if(conf->load_network == 2){ // other?
        
        // TODO
    }

    // get network maximum delay value
    int max_delay = get_max_value(topology->synapses.delay, topology->n_synapses);

    // allocate memory for the SNN structure arrays
    GPU_SNN_t* snn = allocate_memory_for_SNN(topology->n_neurons, topology->n_input, topology->n_synapses, max_delay + 1, conf);
    printf(" Memory allocated");
    fflush(stdout);

    // load information from intermediate structure to the SNN

    // initialize general info
    snn->neuron_type = conf->neuron_type;
    snn->n_neurons = topology->n_neurons;
    snn->n_input_neurons = topology->n_input;
    snn->n_output_neurons = topology->n_output_neurons;
    snn->n_synapses = topology->n_synapses;
    snn->LT = max_delay + 1;

    // set function pointers for managing neurons reinitialization and dynamics simualtion
    set_network_neuron_ptr(snn, conf);

    // initialize the network neurons
    initialize_neurons_CPU(snn, topology, conf);

    // initialize synapses
    initialize_synapses_CPU(snn, topology, conf);

    // connect synapses and neurons
    connect_network_input_criteria(snn, topology, conf);

    // deallocate intermeadite structure used for initializing the network
    deallocate_topology_str(topology);

    // return the initialized SNN structure
    return snn;
}

GPU_SNN_t* allocate_memory_for_SNN(size_t N, size_t iN, size_t S, size_t mD, simulation_configuration_t *conf){

    GPU_SNN_t *snn = (GPU_SNN_t *)calloc(1, sizeof(GPU_SNN_t));

    // allocate memory for neurons and synapses
    switch (conf->neuron_type){
        case 1:
            allocate_memory_for_LIF_neurons(snn, N, conf);
            break;
        default:
            allocate_memory_for_LIF_neurons(snn, N, conf);
            break;
    }
    
    // general neuron properties
    snn->post_fired                    = (char*)malloc(N * sizeof(char));
    snn->post_trace                    = (float*)malloc(N * sizeof(float));
    snn->n_neuron_input_synapses       = (size_t*)malloc(N * sizeof(size_t));
    snn->neuron_input_synapses_offset  = (size_t*)malloc(N * sizeof(size_t));

    // allocate memory for synapses
    snn->w                             = (float*)malloc(S * sizeof(float));
    snn->init_w                        = (float*)malloc(S * sizeof(float));
    snn->dw                            = (float*)malloc(S * sizeof(float));
    snn->delay                         = (int*)malloc(S * sizeof(int));
    snn->lr                            = (int*)malloc(S * sizeof(int));
    snn->pre_neuron_index              = (size_t*)malloc(S * sizeof(size_t));
    snn->post_neuron_index             = (size_t*)malloc(S * sizeof(size_t));
    snn->pre_fired                     = (char*)malloc(S * sizeof(char));
    snn->pre_trace                     = (float*)malloc(S * sizeof(float));

    
    // allocate memory for spike matrix spk matrix
    snn->spk_matrix                    = (char*)malloc((N + iN) * mD * sizeof(char)); // time steps is temporal

    // return allocated structure
    return snn;
}

void deallocate_snn_str(GPU_SNN_t *snn){

    // deallocate neuron model dependant parameters
    switch(snn->neuron_type){
        case 1:
            snn->deallocate_neurons(snn);
            break;
        default:
            printf(" > Default neuron model not implemented!\n");
            exit(1);
            break;
    }

    // deallocate paramteres shared by all neuron models
    if(snn->post_fired) free(snn->post_fired);
    if(snn->post_trace) free(snn->post_trace);
    if(snn->n_neuron_input_synapses) free(snn->n_neuron_input_synapses);
    if(snn->neuron_input_synapses_offset) free(snn->neuron_input_synapses_offset);

    // deallocate synapses parameters
    if(snn->w) free(snn->w);
    if(snn->init_w) free(snn->init_w);
    if(snn->dw) free(snn->dw);
    if(snn->delay) free(snn->delay);
    if(snn->lr) free(snn->lr);
    if(snn->pre_neuron_index) free(snn->pre_neuron_index);
    if(snn->post_neuron_index) free(snn->post_neuron_index);
    if(snn->pre_fired) free(snn->pre_fired);
    if(snn->pre_trace) free(snn->pre_trace);

    // deallocate spike matrix spk matrix
    if(snn->spk_matrix) free(snn->spk_matrix);

    // deallocate snn pointer
    if(snn) free(snn);
}

/// @brief Function to set the pointer functions depending on the neuron type
/// @param snn SNN network structure
/// @param conf Structure containing the configuration data of the simulation
void set_network_neuron_ptr(GPU_SNN_t *snn, simulation_configuration_t *conf){

    switch(conf->neuron_type){
        case 1: // LIF
            set_LIF_ptrs(snn);
            break;
        default: // 

            printf(" Default neuron type not defined\n");
            exit(1);
            break;      
    }
}

void initialize_neurons_CPU(GPU_SNN_t *snn, topology_t *topology, simulation_configuration_t *conf){

    size_t i;
    
    // initialize neuron model dependant variables
    snn->init_neurons(snn, topology, conf);

    // initialize variables shared by all neuron models
    for(i = 0; i<snn->n_neurons; i++){

        snn->n_neuron_input_synapses[i] = 0;
        snn->neuron_input_synapses_offset[i] = 0;
        snn->post_fired[i] = 0;
        snn->post_trace[i] = 0.0;
    }
}

void initialize_synapses_CPU(GPU_SNN_t *snn, topology_t *topology, simulation_configuration_t *conf){

    size_t i;

    for(i = 0; i<snn->n_synapses; i++){

        // read synapse parameters from lists
        snn->w[i]      = topology->synapses.w[i];
        snn->delay[i]  = topology->synapses.delay[i];
        snn->init_w[i] = snn->w[i];
        snn->dw[i]     = 0.0;
        snn->lr[i]     = topology->synapses.lr[i]; // [TODO]: this will call the function in an array of function pointers

        // initialized later
        snn->pre_neuron_index[i]  = 0;
        snn->post_neuron_index[i] = 0;
        
        // control variables
        snn->pre_fired[i] = 0;
        snn->pre_trace[i] = 0.0;
    }
}

void connect_network_input_criteria(GPU_SNN_t *snn, topology_t *topology, simulation_configuration_t *conf){

    // count the number of input synapses and the offsets for each neuron
    count_neurons_input_synapses_and_conpute_offsets(snn, topology);

    // set pre and post neurons indexes for each synapse
    connect_synapses_to_network(snn, topology);
}

size_t get_max_delay(GPU_SNN_t *snn){

    size_t i = 0;

    size_t max_delay = 0;
    for(i=0; i<snn->n_synapses; i++){

        if(snn->delay[i] > max_delay){
            max_delay = snn->delay[i];
        }
    }

    // set max delay
    snn->LT = max_delay;

    return max_delay;
}

void init_batch_snn(GPU_SNN_t *snn, simulation_configuration_t *conf){

    size_t i, j, t, p, B, LT, N, iN, S;

    // temporal variables to store values
    float *tmp_w, *tmp_dw, *tmp_pre_trace, *tmp_post_trace;
    char *tmp_prefired, *tmp_pstfired, *tmp_spk_matrix;

    // [TODO]: improve
    size_t padding = 32; // 32 bytes (256 bits) of padding for copied variables (avoid segfault in vectorized code)

    B = conf->batch_size;
    LT = snn->LT;
    N = snn->n_neurons;
    iN = snn->n_input_neurons;
    S = snn->n_synapses;

    // copy neuron models variables batch size times if they will be written
    switch (conf->neuron_type){
        case 1:
            init_batch_LIF(snn, conf);
            break;
        default:
            init_batch_LIF(snn, conf);
            break;
    }

    tmp_pstfired    = snn->post_fired;
    tmp_post_trace  = snn->post_trace;

    tmp_w           = snn->w;
    tmp_dw          = snn->dw;
    tmp_pre_trace   = snn->pre_trace;
    tmp_prefired    = snn->pre_fired;

    tmp_spk_matrix  = snn->spk_matrix;

    // allocate memory for copied variables
    snn->post_fired = (char*)malloc(N * B * sizeof(char) + padding);
    snn->post_trace = (float*)malloc(N * B * sizeof(float) + padding);
    
    snn->w          = (float*)malloc(S * B * sizeof(float) + padding);
    snn->dw         = (float*)malloc(S * B * sizeof(float) + padding);
    snn->pre_trace  = (float*)malloc(S * B * sizeof(float) + padding);
    snn->pre_fired  = (char*)malloc(S * B * sizeof(char) + padding); 
    
    // reallocate spk matrix
    snn->spk_matrix = (char*)malloc((iN + N) * B * LT * sizeof(char) + padding);


    // initialize neurons using old values
    for(i = 0; i<N; i++){

        for(j = 0; j<B; j++){
            snn->post_fired[i * B + j] = tmp_pstfired[i];
            snn->post_trace[i * B + j] = tmp_post_trace[i];
        }
    }

    // initialize synapses using old values
    for(i = 0; i<S; i++){

        for(j = 0; j<B; j++){

            // copy values
            snn->w[i * B + j] = tmp_w[i];
            snn->dw[i * B + j] = tmp_dw[i];
            snn->pre_trace[i * B + j] = tmp_pre_trace[i];
            snn->pre_fired[i * B + j] = tmp_prefired[i];
        }
    }

    // initialize spk matrix to 0
    for(t = 0; t<LT; t++){
        
        for(i = 0; i<iN + N; i++){
            
            for(j = 0; j<B; j++){

                snn->spk_matrix[(iN + N) * B * t + B * i + j] = 0;
            }
        }
    }


    // initialize padding positions to 0
    for(p = N * B ; p<(padding / sizeof(float)) + N * B; p++){
                      
        snn->post_trace[p] = 0.0;               
    }

    for(p = N * B ; p<(padding / sizeof(char)) + N * B; p++){

        snn->post_fired[p] = 0;     
    }
    for(p = S * B ; p<(padding / sizeof(float)) + S * B; p++){

        snn->w[p] = 0.0;                    
        snn->dw[p] = 0.0;  
        snn->pre_trace[p] = 0.0;                
    }
    for(p = S * B ; p<(padding / sizeof(char)) + S * B; p++){

        snn->pre_fired[p] = 0;
    }

    for(p = (iN + N) * B * LT ; p<(padding / sizeof(char)) + (iN + N) * B * LT; p++){

        snn->spk_matrix[p] = 0;
    }

    // deallocate original arrays
    if(tmp_pstfired) free(tmp_pstfired);
    if(tmp_post_trace) free(tmp_post_trace);

    if(tmp_w) free(tmp_w);
    if(tmp_dw) free(tmp_dw);
    if(tmp_pre_trace) free(tmp_pre_trace);
    if(tmp_prefired) free(tmp_prefired);
    
    if(tmp_spk_matrix) free(tmp_spk_matrix);
}


/* [GET SIZE] */ // [TODO]: revise

double get_snn_size(GPU_SNN_t *snn){
    
    size_t N = snn->n_neurons;
    size_t iN = snn->n_input_neurons;
    size_t S = snn->n_synapses;
    size_t LT = snn->LT;

    // n_input neurons or synapses?
    return (
        (sizeof(GPU_SNN_t) +
        sizeof(size_t) * 10 + // scalars
        
        // neurons 
        sizeof(float) * 5 * N + // neurons floats
        sizeof(int) * 3 * N + // neurons integers
        sizeof(char) * 1 * N + // neurons chars
        sizeof(size_t) * 2 * N + // neurons size_t
        
        // synapses
        sizeof(float) * 4 * S + // synapses floats
        sizeof(int) * 2 * S + // synapses integers
        sizeof(char) * 1 * S + // synapses chars
        sizeof(size_t) * 2 * S + // synapses size_t
        
        // spk matrix
        sizeof(char) * (iN + N) * LT) / 8.0
    );
}

double get_snn_cpy_size(GPU_SNN_t *snn){

    size_t N = snn->n_neurons;
    size_t iN = snn->n_input_neurons;
    size_t S = snn->n_synapses;
    size_t LT = snn->LT;

    return (

        // neurons
        (sizeof(float) * 3 * N + // copied float arrays
        sizeof(int) * 1 * N + // copied int arrays
        sizeof(char) * 1 * N + // copied char arrays

        // synapses
        sizeof(float) * 3 * S + // copied float arrays 
        sizeof(char) * 1 * S + // copied char arrays
        
        // spk matrix
        sizeof(char) * LT * (iN + N)) / 8.0
    );
}

/* [PRINT] */

void print_network(GPU_SNN_t *snn){

    printf(" > Function for printing a network not implemented yet!\n");
}

void print_networks(GPU_SNN_t *snn, simulation_configuration_t *conf){

    printf(" > Function for printing networks not implemented yet!\n");
}

/* [STORE] */

void store_network_in_file(GPU_SNN_t *snn, simulation_configuration_t *conf){
 
    printf(" > Function for storing a network not implemented yet!\n");
}