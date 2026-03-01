#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>

#include "toml_c/toml.h"

#include "config/config_loader.h"
#include "networks/snn.h"
#include "utils.h"
#include "neuron_models/neuron_models.h"


/* [PRIVATE] */

/// @brief Structure containing the initial values of the network parameters
typedef struct network_construction_lists_t {
    
    // general network information
    size_t n_neurons, n_input_neurons, n_output_neurons, n_synapses, max_delay;

    // connectivity
    int *neuron_excitatory, *training_zones; // change to uint in the future
    int *delay_list;
    double *weight_list;
    int **synaptic_connections;
    
    // LIF neuron
    double *v_list, *v_thres_list, *v_rest_list, *R_list; // parameters for all neurons [n_neurons]
    int *r_time_list; // refractory times for all neurons [n_neurons] 

} network_construction_lists_t;


/// @brief Load the information of the neuron from the input files to a intermediate structure.
/// @param conf Configuration structure with information about the network
/// @return Structure with the arrays storing the data of the neurons, synapses and connectivity to initialize the network
network_construction_lists_t* load_network_information_in_lists_from_file(simulation_configuration_t *conf){

    FILE *f = NULL, *f_neurons, *f_synapses;
    char errbuf[100];
    size_t i, j;

    // struct to store arrays
    network_construction_lists_t *lists = (network_construction_lists_t*)calloc(1, sizeof(network_construction_lists_t));

    // define table and parameters variables
    toml_table_t *tbl, *tbl_general, *tbl_neurons, *tbl_synapses;

    toml_array_t *behaviour_lst, *v_thres_lst, *v_rest_lst, *t_refract_lst, *res_lst, *latency_lst, 
                    *weights_lst, *training_zones_lst, *connection_lst_lst, *connection_lst;

    toml_value_t n_neurons, n_input_neurons, n_output_neurons, n_synapses, v_thres, v_rest, t_refract, res, latency, 
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
    
    // load information in lists structure
    lists->n_neurons = N;
    lists->n_input_neurons = iN;
    lists->n_output_neurons = oN;
    lists->n_synapses = S;

    // check that all the information has been loaded correctly
    if(!(n_neurons.ok && n_input_neurons.ok && n_output_neurons.ok && n_synapses.ok)){
        printf("The number of neurons, input neurons, output neurons, synapses, input synapses and output synapses must be provided in the network file!");
        fflush(stdout);
        exit(1);
    }


    // read files with the information of the neurons and synapses
    open_file(&f_neurons, conf->network_neurons_file); // TOML file
    open_file(&f_synapses, conf->network_synapses_file); // TOML file


    // allocate memory for arrays related to neurons data
    lists->v_thres_list = (double *)malloc(N * sizeof(double)); // thresholds
    lists->v_rest_list = (double *)malloc(N * sizeof(double)); // resting potentials
    lists->r_time_list = (int *)malloc(N * sizeof(int)); // refractary times
    lists->R_list = (double *)malloc(N * sizeof(double)); // resistances


    /* load [Neurons] section */
    // load whether different parameters are included in the input file
    v_thres = toml_table_int(tbl_neurons, "v_thres");
    v_rest = toml_table_int(tbl_neurons, "v_rest");
    t_refract = toml_table_int(tbl_neurons, "t_refract");
    res = toml_table_int(tbl_neurons, "resistance");

    float tmp = 0.0;
    int tmp_int = 0;

    // load thresholds
    if(v_thres.ok && v_thres.u.i == 1){
        for(i=0; i<N; i++)
            fscanf(f_neurons, "%lf", &((lists->v_thres_list)[i]));
    }
    else if(v_thres.ok && v_thres.u.i == 2){
        fscanf(f_neurons, "%lf", &(tmp));
        for(i=0; i<N; i++)
            lists->v_thres_list[i] = tmp;
    }
    else{
        for(i=0; i<N; i++)
            lists->v_thres_list[i] = 1.0;
    }

    // load resting potentials
    if(v_rest.ok && v_rest.u.i == 1){
        for(i=0; i<N; i++)
            fscanf(f_neurons, "%lf", &((lists->v_rest_list)[i]));
    }    
    else if(v_rest.ok && v_rest.u.i == 2){
        fscanf(f_neurons, "%lf", &(tmp));
        for(i=0; i<N; i++)
            lists->v_rest_list[i] = tmp;
    }
    else{
        for(i=0; i<N; i++)
            lists->v_rest_list[i] = 0.0;
    }

    // load refractory periods
    if(t_refract.ok && t_refract.u.i == 1){
        for(i=0; i<N; i++)
            fscanf(f_neurons, "%d", &((lists->r_time_list)[i]));
    }    
    else if(t_refract.ok && t_refract.u.i == 2){
        fscanf(f_neurons, "%d", &(tmp_int));
        for(i=0; i<N; i++)
            lists->r_time_list[i] = (int)(tmp_int);
    }
    else{
        for(i=0; i<N; i++)
            lists->r_time_list[i] = 1;
    }

    // load resistances
    if(res.ok && res.u.i == 1){
        for(i=0; i<N; i++)
            fscanf(f_neurons, "%lf", &((lists->R_list)[i]));
    }    
    else if(res.ok && res.u.i == 2){
        fscanf(f_neurons, "%lf", &(tmp));
        for(i=0; i<N; i++)
            lists->R_list[i] = tmp;
    }
    else{
        for(i=0; i<N; i++)
            lists->R_list[i] = 1.0;
    }

    fclose(f_neurons);
    

    /* Synapses section */

    // allocate memory for arrays related to synapses data
    lists->weight_list = (double *)malloc(S * sizeof(double)); // weights
    lists->delay_list = (int *)malloc(S * sizeof(int)); // delays
    lists->training_zones = (int *)malloc(S * sizeof(int)); // learnin rules

    // load whether different parameters are included in the input file
    latency = toml_table_int(tbl_synapses, "delay");
    weight = toml_table_int(tbl_synapses, "weight");
    training_zone = toml_table_int(tbl_synapses, "training_zone");


    // load delays and compute the maximum value
    if(latency.ok && latency.u.i == 1){
    
        lists->max_delay = 0;
        for(i=0; i<S; i++){
            fscanf(f_synapses, "%d", &(lists->delay_list[i]));

            // store maximum delay value
            if(lists->delay_list[i] > lists->max_delay)
                lists->max_delay = lists->delay_list[i];
        }
        lists->max_delay += 1; // store max_delay
    }
    else if(latency.ok && latency.u.i == 2){
        fscanf(f_synapses, "%d", &(tmp_int));
        for(i = 0; i<S; i++){
            lists->delay_list[i] = tmp_int;
        }
        lists->max_delay = tmp_int + 1;
    }
    else{
        for(i=0; i<S; i++){
            lists->delay_list[i] = 1;
        }
        lists->max_delay = 2;
    }

    // load weights
    if(weight.ok && weight.u.i == 1){
        for(i=0; i<S; i++){
            fscanf(f_synapses, "%lf", &(lists->weight_list[i]));
        }
    }
    else if(weight.ok && weight.u.i == 2){
        fscanf(f_synapses, "%lf", &(tmp));
        for(i = 0; i<S; i++)
            lists->weight_list[i] = tmp;
    }
    else{
        for(i = 0; i<S; i++)
            lists->weight_list[i] = 0.25;
    }

    // load training zones
    if(training_zone.ok && training_zone.u.i == 1){
        for(i=0; i<S; i++){
            fscanf(f_synapses, "%d", &(lists->training_zones[i]));
        }
    }
    else if(training_zone.ok && training_zone.u.i == 2){
        fscanf(f_synapses, "%d", &(tmp_int));
        for(i = 0; i<S; i++)
            lists->training_zones[i] = tmp_int;
    }
    else{
        for(i = 0; i<S; i++)
            lists->training_zones[i] = 0;
    }

    printf(" >> Synapses loaded from file!\n");
    fflush(stdout);


    /* Connectivity */

    // allocate memory for connectibity
    lists->synaptic_connections = (int **)malloc(N * sizeof(int *)); // +1, since the input layer is stored too

    int number_connections;
    for(i=0; i<N; i++){ // network input synapses are loaded first and output synapses last
        fscanf(f_synapses, "%d", &number_connections);

        // alloc memory
        (lists->synaptic_connections)[i] = (int*)malloc((number_connections * 2 + 1) * sizeof(int)); // for each connection the neuron id and the number of synapses must be stored
        (lists->synaptic_connections)[i][0] = number_connections;

        for(j = 0; j<number_connections; j++){
            fscanf(f_synapses, "%d", &((lists->synaptic_connections)[i][j * 2 + 1])); // number of synapses connected to that neuron
            fscanf(f_synapses, "%d", &((lists->synaptic_connections)[i][j * 2 + 2])); // number of synapses connected to that neuron
        }
    }
    fclose(f_synapses);

    printf(" >> Synapses section loaded\n");
    fflush(stdout);

    return lists;
}

/// @brief Deallocate intermediate structure that stores data of the neurons and synapses
/// @param data Structure of arrays that stores the values of the neurons and synapses
void deallocate_arrays_str(network_construction_lists_t *data){

    // deallocate internal arrays
    if(data->v_list)
        free(data->v_list);
    if(data->v_thres_list)
        free(data->v_thres_list);
    if(data->v_rest_list)
        free(data->v_rest_list);
    if(data->R_list)
        free(data->R_list);
    if(data->r_time_list)
        free(data->r_time_list);

    if(data->weight_list)
        free(data->weight_list);
    if(data->training_zones)
        free(data->training_zones);
    if(data->delay_list)
        free(data->delay_list);
    
    for(size_t i = 0; i<data->n_neurons; i++){
        if(data->synaptic_connections[i])
            free(data->synaptic_connections[i]);
    }
    if(data->synaptic_connections)
        free(data->synaptic_connections);

    // deallocate struct
    if(data)
        free(data);
}


/* [Neuron models allocations, deallocations and initializations] */


/// @brief Function to set the pointer functions depending on the neuron type
/// @param snn SNN network structure
/// @param conf Structure containing the configuration data of the simulation
void set_network_neuron_ptr(GPU_SNN_t *snn, simulation_configuration_t *conf){

    switch(conf->neuron_type){
        case 1: // LIF
            snn->reinit_neurons = &reinitialize_LIF_neurons_batch; // LIF neurons reinitializations
            snn->neuron_dynamics = &compute_LIF_V_batch; // LIF neuron dynamics computation
            break;
        default: // LIF
            snn->reinit_neurons = &reinitialize_LIF_neurons_batch; // LIF neurons reinitializations
            snn->neuron_dynamics = &compute_LIF_V_batch; // LIF neuron dynamics computation
            break;      
    }
}

/// @brief Allocate memory for  LIF neuron arrays
/// @param snn SNN structure to initialize neurons in
/// @param N Number of neurons in the network
/// @param conf structure with configuration information
void allocate_memory_for_LIF_neurons(GPU_SNN_t *snn, size_t N, simulation_configuration_t *conf){

    // allocate memory for neurons and synapses
    snn->v               = (float*)malloc(N * conf->batch_size * sizeof(float));
    snn->v_thresh        = (float*)malloc(N * conf->batch_size * sizeof(float));
    snn->v_rest          = (float*)malloc(N * conf->batch_size * sizeof(float));
    snn->arrI            = (float*)malloc(N * conf->batch_size * sizeof(float));
    snn->r_period        = (int*)malloc(N * conf->batch_size * sizeof(int));
    snn->r_period_remain = (int*)malloc(N * conf->batch_size * sizeof(int));
    snn->res             = (int*)malloc(N * conf->batch_size * sizeof(int));
}

/// @brief Initialize LIF neuron arrays
/// @param snn SNN structure to initialize neurons in
/// @param data structure with the values for initializing the neurons
/// @param conf structure with configuration information
/// @return SNN structure
void initialize_LIF_neurons(GPU_SNN_t *snn, network_construction_lists_t *data, simulation_configuration_t *conf){

    size_t i;
    
    // initialize neuron parameters from array
    for(i = 0; i<snn->n_neurons; i++){

        snn->v_thresh[i] = data->v_thres_list[i];
        snn->v_rest[i] = data->v_rest_list[i]; // this or the next one?
        snn->v[i] = snn->v_rest[i]; 
        snn->res[i] = data->R_list[i];
        snn->r_period[i] = data->r_time_list[i];
        snn->r_period_remain[i] = -1;
        snn->arrI[i] = 0.0;
    }
}

/// @brief Initialize synaptic arrays
/// @param snn SNN structure to initialize synapses in
/// @param data structure with the values for initializing the synapses
/// @param conf structure with configuration information
/// @return SNN structure
void initialize_neurons_CPU(GPU_SNN_t *snn, network_construction_lists_t *data, simulation_configuration_t *conf){

    size_t i;
    
    // init variables that depend on the neuron type
    switch(conf->neuron_type){
        case 1:
            initialize_LIF_neurons(snn, data, conf);
        default:
            initialize_LIF_neurons(snn, data, conf);
    }
    
    // initialize shared variables by all neurons
    for(i = 0; i<snn->n_neurons; i++){

        snn->n_neuron_input_synapses[i] = 0;
        snn->neuron_input_synapses_offset[i] = 0;
        snn->post_fired[i] = 0;
        snn->post_trace[i] = 0.0;
    }
}


/* [PUBLIC] */
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

GPU_SNN_t* allocate_memory_for_SNN(size_t N, size_t iN, size_t S, size_t mD, simulation_configuration_t *conf){

    GPU_SNN_t *snn = (GPU_SNN_t *)malloc(sizeof(GPU_SNN_t));

    // allocate memory for neurons and synapses
    switch (conf->neuron_type){
        case 1:
            allocate_memory_for_LIF_neurons(snn, N, conf);
        default:
            allocate_memory_for_LIF_neurons(snn, N, conf);
    }
    
    // general neuron properties
    snn->post_fired                    = (char*)malloc(N * conf->batch_size * sizeof(char));
    snn->post_trace                    = (float*)malloc(N * conf->batch_size * sizeof(float));
    snn->n_neuron_input_synapses       = (size_t*)malloc(N * conf->batch_size * sizeof(size_t));
    snn->neuron_input_synapses_offset  = (size_t*)malloc(N * conf->batch_size * sizeof(size_t));

    // allocate memory for synapses
    snn->w                             = (float*)malloc(S * conf->batch_size * sizeof(float));
    snn->init_w                        = (float*)malloc(S * conf->batch_size * sizeof(float));
    snn->dw                            = (float*)malloc(S * conf->batch_size * sizeof(float));
    snn->delay                         = (int*)malloc(S * conf->batch_size * sizeof(int));
    snn->lr                            = (int*)malloc(S * conf->batch_size * sizeof(int));
    snn->pre_neuron_index              = (size_t*)malloc(S * conf->batch_size * sizeof(size_t));
    snn->post_neuron_index             = (size_t*)malloc(S * conf->batch_size * sizeof(size_t));
    snn->pre_fired                     = (char*)malloc(S * conf->batch_size * sizeof(char));
    snn->pre_trace                     = (float*)malloc(S * conf->batch_size * sizeof(float));

    
    // allocate memory for spike matrix spk matrix
    snn->spk_matrix                    = (char*)malloc((N + iN) * mD * conf->batch_size * sizeof(char)); // time steps is temporal

    // return allocated structure
    return snn;
}

void deallocate_snn_str(GPU_SNN_t *snn){

    // deallocate internal arrays
    if(snn->v)
        free(snn->v);
    if(snn->v_thresh)
        free(snn->v_thresh);
    if(snn->v_rest)
        free(snn->v_rest);
    if(snn->r_period)
        free(snn->r_period);
    if(snn->r_period_remain)
        free(snn->r_period_remain);
    if(snn->res)
        free(snn->res);
    if(snn->post_fired)
        free(snn->post_fired);
    if(snn->post_trace)
        free(snn->post_trace);
    if(snn->arrI)
        free(snn->arrI);
    if(snn->n_neuron_input_synapses)
        free(snn->n_neuron_input_synapses);
    if(snn->neuron_input_synapses_offset)
        free(snn->neuron_input_synapses_offset);
    
    if(snn->w)
        free(snn->w);
    if(snn->init_w)
        free(snn->init_w);
    if(snn->dw)
        free(snn->dw);
    if(snn->delay)
        free(snn->delay);
    if(snn->lr)
        free(snn->lr);
    if(snn->pre_neuron_index)
        free(snn->pre_neuron_index);
    if(snn->post_neuron_index)
        free(snn->post_neuron_index);
    if(snn->pre_fired)
        free(snn->pre_fired);
    if(snn->pre_trace)
        free(snn->pre_trace);

    if(snn->spk_matrix)
        free(snn->spk_matrix);

    // deallocate struct
    if(snn)
        free(snn);
}

/// @brief Initialize synaptic arrays
/// @param snn SNN structure to initialize synapses in
/// @param data structure with the values for initializing the synapses
/// @param conf structure with configuration information
void initialize_synapses_CPU(GPU_SNN_t *snn, network_construction_lists_t *data, simulation_configuration_t *conf){

    size_t i;

    for(i = 0; i<snn->n_synapses; i++){

        // read synapse parameters from lists
        snn->w[i]      = data->weight_list[i];
        snn->delay[i]  = data->delay_list[i];
        snn->init_w[i] = snn->w[i];
        snn->dw[i]     = 0.0;
        snn->lr[i]     = data->training_zones[i]; // [TODO]: this will call the function in an array of function pointers

        // initialized later
        snn->pre_neuron_index[i]  = 0;
        snn->post_neuron_index[i] = 0;
        
        // control variables
        snn->pre_fired[i] = 0;
        snn->pre_trace[i] = 0.0;
    }
}

/// @brief Compute the input synapses offsets for each neuron
/// @param snn SNN structure
/// @param data Structure containing helper values
void count_neurons_input_synapses_and_conpute_offsets(GPU_SNN_t *snn, network_construction_lists_t *data){
    
    size_t i, j, off = 0;

    // loop over neurons and count number of input synapses
    for(i = 0; i<snn->n_neurons; i++){

        // loop over input synapses of the neuron
        for(j = 0; j<(size_t)data->synaptic_connections[i][0]; j++){

            // set the number of input synapses
            // [i * 2 + 1] is the index of the neuron
            // [i * 2 + 2] the number of synapses from that neuron to the current 
            snn->n_neuron_input_synapses[i] += (size_t)data->synaptic_connections[i][j*2+2];
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
/// @param data Structure containing helper values
void connect_synapses_to_network(GPU_SNN_t *snn, network_construction_lists_t *data){

    size_t i, j, l, next_syn = 0;

    // loop over the neurons
    for(i = 0; i<snn->n_neurons; i++){

        // loop over the input synapses of the neuron
        for(j = 0; j<(size_t)data->synaptic_connections[i][0]; j++){

            // set the input and output neurons of the synapse
            for(l = 0; l<(size_t)data->synaptic_connections[i][j * 2 + 2]; l++){
                
                snn->post_neuron_index[next_syn] = i + snn->n_input_neurons; // first [n_input_neurons] neurons are virtual input
                snn->pre_neuron_index[next_syn] = (size_t)data->synaptic_connections[i][j * 2 + 1];
                next_syn ++;
            }
        }
    }
}

/// @brief Connect the network following an input criteria, where, for each neuron, its input synapses are stored
/// @param snn SNN structure
/// @param data Structure containing helper values
/// @param conf structure containing the configuration of the simulation
void connect_network_input_criteria(GPU_SNN_t *snn, network_construction_lists_t *data, simulation_configuration_t *conf){

    // count the number of input synapses and the offsets for each neuron
    count_neurons_input_synapses_and_conpute_offsets(snn, data);

    // set pre and post neurons indexes for each synapse
    connect_synapses_to_network(snn, data);
}


/* [PUBLIC FUNCTIONS] */

GPU_SNN_t* initialize_network_cpu(simulation_configuration_t *conf){

    // load network information into intermedaite arrays
    network_construction_lists_t *data; 


    // store information of network, neurons and synapses in an intermediate structure
    if(conf->load_network == 0){ // load from file
        data = load_network_information_in_lists_from_file(conf);
    }
    else if(conf->load_network == 1){ // generate
        
        // TODO
    }
    else if(conf->load_network == 2){ // other?
        
        // TODO
    }

    // allocate memory for the SNN structure arrays
    GPU_SNN_t* snn = allocate_memory_for_SNN(data->n_neurons, data->n_input_neurons, data->n_synapses, data->max_delay, conf);
    printf(" Memory allocated");
    fflush(stdout);

    // load information from intermediate structure to the SNN

    // initialize general info
    snn->n_neurons = data->n_neurons;
    snn->n_input_neurons = data->n_input_neurons;
    snn->n_output_neurons = data->n_output_neurons;
    snn->n_synapses = data->n_synapses;
    snn->LT = data->max_delay;

    // set function pointers for managing neurons reinitialization and dynamics simualtion
    set_network_neuron_ptr(snn, conf);

    // initialize the network neurons

    initialize_neurons_CPU(snn, data, conf);

    // initialize synapses
    initialize_synapses_CPU(snn, data, conf);

    // connect the synapses and the neurons
    connect_network_input_criteria(snn, data, conf);

    // deallocate intermeadite arrays structure used for initializing the neurons and synapses
    deallocate_arrays_str(data);

    // return the initialized SNN structure
    return snn;
}

void cpy_snn(GPU_SNN_t *snn, simulation_configuration_t *conf){

    size_t i, j, t, p, B, LT, N, iN, S;

    // temporal variables to store values
    float *tmp_v, *tmp_arrI, *tmp_w, *tmp_dw, *tmp_pre_trace, *tmp_post_trace;
    int *tmp_rperiod_remain;
    char *tmp_pstfired, *tmp_prefired, *tmp_spk_matrix;

    // [TODO]: improve
    size_t padding = 32; // 32 bytes (256 bits) of padding for copied variables (avoid segfault in vectorized code)

    B = conf->batch_size;
    LT = snn->LT;
    N = snn->n_neurons;
    iN = snn->n_input_neurons;
    S = snn->n_synapses;

    // store all values for copying later
    tmp_v                         = snn->v;
    tmp_rperiod_remain            = snn->r_period_remain;
    tmp_pstfired                  = snn->post_fired;
    tmp_post_trace                = snn->post_trace;
    tmp_arrI                      = snn->arrI;

    tmp_w                         = snn->w;
    tmp_dw                        = snn->dw;
    tmp_pre_trace                 = snn->pre_trace;
    tmp_prefired                  = snn->pre_fired;

    tmp_spk_matrix                = snn->spk_matrix;

    // allocate memory for copied variables
    snn->v                            = (float*)malloc(N * B * sizeof(float) + padding);
    snn->arrI                         = (float*)malloc(N * B * sizeof(float) + padding);
    snn->r_period_remain              = (int*)malloc(N * B * sizeof(int) + padding);
    snn->post_fired                   = (char*)malloc(N * B * sizeof(char) + padding);
    snn->post_trace                   = (float*)malloc(N * B * sizeof(float) + padding);

    snn->w                            = (float*)malloc(S * B * sizeof(float) + padding);
    snn->dw                           = (float*)malloc(S * B * sizeof(float) + padding);
    snn->pre_trace                    = (float*)malloc(S * B * sizeof(float) + padding);
    snn->pre_fired                    = (char*)malloc(S * B * sizeof(char) + padding); 
    
    // reallocate spk matrix
    snn->spk_matrix                   = (char*)malloc((iN + N) * B * LT * sizeof(char) + padding);


    // initialize neurons using old values
    for(i = 0; i<N; i++){

        for(j = 0; j<B; j++){

            // copy values
            snn->v[i * B + j] = tmp_v[i];
            snn->r_period_remain[i * B + j] = tmp_rperiod_remain[i];
            snn->post_fired[i * B + j] = tmp_pstfired[i];
            snn->post_trace[i * B + j] = tmp_post_trace[i];
            snn->arrI[i * B + j] = tmp_arrI[i];
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
    for(p = N * B ; p<padding / sizeof(float); p++){

        snn->v[p] = 0.0;                       
        snn->arrI[p] = 0.0;                         
        snn->post_trace[p] = 0.0;               
    }
    for(p = N * B ; p<padding / sizeof(int); p++){

        snn->r_period_remain[p] = 0;              
    }
    for(p = N * B ; p<padding / sizeof(char); p++){

        snn->post_fired[p] = 0;     
    }
    for(p = S * B ; p<padding / sizeof(float); p++){

        snn->w[p] = 0.0;                    
        snn->dw[p] = 0.0;  
        snn->pre_trace[p] = 0.0;                
    }
    for(p = S * B ; p<padding / sizeof(char); p++){

        snn->pre_fired[p] = 0;
    }

    for(p = (iN + N) * B * LT ; p<padding / sizeof(char); p++){

        snn->spk_matrix[p] = 0;
    }

    // deallocate original arrays
    if(tmp_v)
        free(tmp_v);
    if(tmp_rperiod_remain)
        free(tmp_rperiod_remain);
    if(tmp_pstfired)
        free(tmp_pstfired);
    if(tmp_post_trace)
        free(tmp_post_trace);
    if(tmp_arrI)
        free(tmp_arrI);
    if(tmp_w)
        free(tmp_w);
    if(tmp_dw)
        free(tmp_dw);
    if(tmp_pre_trace)
        free(tmp_pre_trace);
    if(tmp_prefired)
        free(tmp_prefired);
    if(tmp_spk_matrix)
        free(tmp_spk_matrix);
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