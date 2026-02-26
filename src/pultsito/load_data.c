#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h> 
#include <ctype.h>
#include <toml_c/toml-c.h>
#include <time.h>

#include "snn_library.h"
#include "load_data.h"


/* [PRIVATE] */

int open_file(FILE **f, const char *file_name){
    
    *f = fopen(file_name, "r");
    if (*f == NULL){
        printf("Error opening the file %s\n", file_name);
        return 1;
    }    
    //printf("File openned!\n");
    return 0;
}

int open_file_w(FILE **f, const char *file_name){
    
    *f = fopen(file_name, "a"); // append mode, no overwriting
    if (*f == NULL){
        printf(" > Error opening the file %s\n", file_name);
        return 1;
    }    
    return 0;
}



// TODO: revise

void open_results_files(simulation_configuration_t *conf){

    // open files for storing results

    if(conf->store_generated_spikes == 1){
        conf->f_gs = fopen(conf->generated_spikes_file, "w");
    }

    if(conf->store_n_spikes == 1){        
        conf->f_ns = fopen(conf->n_spikes_per_neuron_file, "w");
    }

    if(conf->store_execution_times == 1){
        conf->f_et = fopen(conf->execution_times_file, "w");
    }

    if(conf->store_network == 1){
        conf->f_sn = fopen(conf->store_network_file, "w");
    }
}

// TODO: revise

void close_results_files(simulation_configuration_t *conf){

    // open files for storing results

    if(conf->store_generated_spikes == 1){
        fclose(conf->f_gs);
    }

    if(conf->store_n_spikes == 1){        
        fclose(conf->f_ns);
    }

    if(conf->store_execution_times == 1){
        fclose(conf->f_et);
    }

    if(conf->store_network == 1){
        fclose(conf->f_sn);
    }
}


simulation_configuration_t* load_general_section_from_toml(simulation_configuration_t *conf, toml_table_t *tbl){

    // [general] section
    toml_value_t neuron_type, n_process, cuda, multigpu, learn, batch_size, thrN, load_network, load_dataset;

    /* read [general] section, and set default values in case data is not provided */
    neuron_type   =  toml_table_int(tbl, "neuron_type"); // neuron type (0: combined, 1: LIF)
    n_process     =  toml_table_int(tbl, "n_process"); // number of CPU threads (min 1)
    cuda          =  toml_table_int(tbl, "cuda"); // whether simulation should be performed using GPU (0: all, N: number of devices)
    multigpu      =  toml_table_int(tbl, "multi_gpu"); // number of GPUs
    learn         =  toml_table_int(tbl, "learn"); // inference or training
    batch_size    =  toml_table_int(tbl, "batch_size");
    thrN          =  toml_table_int(tbl, "thrN");
    load_network  =  toml_table_int(tbl, "load_network");
    load_dataset  =  toml_table_int(tbl, "load_dataset");


    if(!neuron_type.ok)
        neuron_type.u.i = 0; // LIF neuron
    if(!n_process.ok)
        n_process.u.i = 1; // serial execution
    if(!cuda.ok)
        cuda.u.i = 0; // no cuda
    if(!multigpu.ok)
        multigpu.u.i = 1; // no multigpu
    if(!learn.ok)
        learn.u.i = 0; // inference
    if(!batch_size.ok)
        batch_size.u.i = 1; // only one sample in the batch
    if(!thrN.ok)
        thrN.u.i = 1;
    if(!load_network.ok)
        load_network.u.i = 1; // load
    if(!load_dataset.ok)
        load_dataset.u.i = 1; // load


    // load information in configuration struct
    conf->neuron_type = neuron_type.u.i;
    conf->n_process = (size_t)n_process.u.i;
    conf->cuda = cuda.u.i;
    conf->multi_gpu = multigpu.u.i;
    conf->learn = learn.u.i;
    conf->batch_size = (size_t)batch_size.u.i;
    conf->thrN = (size_t)thrN.u.i;
    conf->load_network = load_network.u.i;
    conf->load_dataset = load_dataset.u.i;

    return conf;
}

simulation_configuration_t* load_simulation_section_from_toml(simulation_configuration_t *conf, toml_table_t *tbl){

    // [simulation] section
    toml_value_t time_steps, max_spikes, max_input_spikes;

    /* read [simulation] section */
    time_steps = toml_table_int(tbl, "time_steps"); // simulation time steps
    max_input_spikes = toml_table_int(tbl, "max_input_spikes"); // length for input and output neurons

    if(!time_steps.ok){

        printf(" >> It is necessary to provide simulation time steps!\n");
        fflush(stdout);
        exit(1);
    }
    if(!max_input_spikes.ok)
        max_input_spikes.u.i = time_steps.u.i;

    // load information in configuration struct
    conf->time_steps = (size_t)time_steps.u.i;
    conf->max_input_spikes = (size_t)max_input_spikes.u.i;

    return conf;
}

simulation_configuration_t* load_dataset_section_from_toml(simulation_configuration_t *conf, toml_table_t *tbl){

    // [dataset] section
    toml_value_t train_provided, train_set, train_labels, n_train, 
                test_provided, test_set, test_labels, n_test,
                dataset_name, n_classes, epochs, input_size, shuffle_samples;
    
    
    /* read [dataset] section */
    train_provided  =  toml_table_int(tbl, "train_provided");
    train_set       =  toml_table_string(tbl, "train_set");
    train_labels    =  toml_table_string(tbl, "train_labels");
    n_train         =  toml_table_int(tbl, "n_train");

    test_provided   =  toml_table_int(tbl, "test_provided");
    test_set        =  toml_table_string(tbl, "test_set");
    test_labels     =  toml_table_string(tbl, "test_labels");
    n_test          =  toml_table_int(tbl, "n_test");

    dataset_name    =  toml_table_string(tbl, "dataset_name");
    n_classes       =  toml_table_int(tbl, "num_classes");
    input_size      =  toml_table_int(tbl, "input_size");

    shuffle_samples =  toml_table_int(tbl, "shuffle_samples");
    epochs          =  toml_table_int(tbl, "epochs");

    // check that all is correctly loaded
    if(!train_set.ok){
        
        printf(" >> File to load the train / simulation dataset must be provided!\n");
        fflush(stdout);
        exit(1);
    }
    if(!n_train.ok){
        
        printf(" >> The number of samples in the set must be provided!\n");
        fflush(stdout);
        exit(1);    
    }
    if(!train_labels.ok){
        
        printf(" >> WARNING: labels not provided for the training set!\n");
        fflush(stdout);
    }

    if(!test_set.ok){
        
        printf(" >> WARNING: test set not provided!\n");
        fflush(stdout);
    }
    else{
        if(!n_test.ok){
        
            printf(" >> The number of samples in the test set must be provided!\n");
            fflush(stdout);
            exit(1);    
        }
        if(!test_labels.ok){
            
            printf(" >> WARNING: labels not provided for test set!\n");
            fflush(stdout);
        }
    }

    if(!dataset_name.ok){
        dataset_name.u.s = "Dataset";
    }
    if(!n_classes.ok){ // [TODO]: count?
        
        n_classes.u.i = 1;
    }
    if(!epochs.ok){
        
        epochs.u.i = 1;
    }
    if(!input_size.ok){
        printf(" >> The input size is necessary!\n");
        fflush(stdout);
        exit(1);
    }
    if(!shuffle_samples.ok){
        shuffle_samples.u.i = 0; // do not shuffle
    }

    // store dataset information
    conf->train_provided = train_provided.u.i;
    conf->train_set = train_set.u.s;
    conf->train_labels = train_labels.u.s;
    conf->n_train = (size_t)n_train.u.i;

    conf->test_provided = test_provided.u.i;
    conf->test_set = test_set.u.s;
    conf->test_labels = test_labels.u.s;
    conf->n_test = (size_t)n_test.u.i;

    conf->dataset_name = dataset_name.u.s;
    conf->n_classes = (size_t)n_classes.u.i;
    conf->epochs = (size_t)epochs.u.i;
    conf->input_size = (size_t)input_size.u.i;

    conf->shuffle_samples = shuffle_samples.u.i;

    return conf;
}

// TODO: revise
simulation_configuration_t* load_output_section_from_toml(simulation_configuration_t *conf, toml_table_t *tbl){

    // [output] section
    toml_value_t generated_spikes_file, execution_times_file, n_spikes_per_neuron_file, store_network, 
                store_generated_spikes, store_n_spikes, store_execution_times, store_network_file;

    /* read [outputt] section */
    store_generated_spikes   =  toml_table_int(tbl, "store_generated_spikes");
    generated_spikes_file    =  toml_table_string(tbl, "generated_spikes");
    store_n_spikes           =  toml_table_int(tbl, "store_n_spikes");
    n_spikes_per_neuron_file =  toml_table_string(tbl, "spikes_per_neuron");
    store_execution_times    =  toml_table_int(tbl, "store_execution_times");
    execution_times_file     =  toml_table_string(tbl, "execution_times");
    store_network            =  toml_table_int(tbl, "store_network"); // whether to store or not the file
    store_network_file       =  toml_table_string(tbl, "store_network_file"); // file to store the network

    // check that everything is loaded
    if(!store_generated_spikes.ok){
        store_generated_spikes.u.i = 0;
    }
    if(store_generated_spikes.u.i == 1 && !generated_spikes_file.ok){
        
        printf("A file to store generated spikes must be provided!\n");
        fflush(stdout);
        exit(1);
    }
    if(!store_execution_times.ok){
        store_execution_times.u.i = 0;
    }
    if(store_execution_times.u.i == 1 && !execution_times_file.ok){
        
        printf("A file to store execution times must be provided!\n");
        fflush(stdout);
        exit(1);
    }
    if(!store_n_spikes.ok){
        store_n_spikes.u.i = 0;
    }
    if(store_n_spikes.u.i == 1 && !n_spikes_per_neuron_file.ok){
        
        printf("A file to store the number of spikes generated per neuron must be provided!\n");
        fflush(stdout);
        exit(1);
    }
    if(!store_network.ok){
        
        if(conf->learn == 0)
            store_network.u.i = 0; // do not store since it is inference
        else
            store_network.u.i = 1; // network is trained, so store
    }
    if(store_network.u.i == 1 && !store_network_file.ok){

        printf("The file name to store the final network must be provided!\n");
        fflush(stdout);
        exit(1);
    }

    conf->store_generated_spikes = store_generated_spikes.u.i;
    conf->generated_spikes_file = generated_spikes_file.u.s;

    conf->store_execution_times = store_execution_times.u.i;
    conf->execution_times_file = execution_times_file.u.s;

    conf->store_n_spikes = store_n_spikes.u.i;
    conf->n_spikes_per_neuron_file = n_spikes_per_neuron_file.u.s;

    conf->store_network = store_network.u.i;
    if(conf->store_network == 1)
        conf->store_network_file = store_network_file.u.s;

    return conf;
}

simulation_configuration_t* load_network_section_from_toml(simulation_configuration_t *conf, toml_table_t *tbl){

    // [network] section
    toml_value_t network_file, network_neurons_file, network_synapses_file, behaviours, delays, weights,
                training_zones, thresh, v_rest, t_refract, R;

    
    /* read [network] section */
    network_file           =  toml_table_string(tbl, "network_file");
    network_neurons_file   =  toml_table_string(tbl, "network_neurons_file");
    network_synapses_file  =  toml_table_string(tbl, "network_synapses_file");

    // check if all is in configuration file
    if(!network_file.ok)
    {
        printf(" > The file to load the network from must be provided!\n");
        fflush(stdout);
        exit(1);
    }

    if(network_file.ok){
        conf->network_file = network_file.u.s;
    }
    if(network_neurons_file.ok){
        conf->network_neurons_file = network_neurons_file.u.s;
    }
    if(network_synapses_file.ok){
        conf->network_synapses_file = network_synapses_file.u.s;
    }

    return conf;
}


/* [PUBLIC] */

simulation_configuration_t* load_configuration_params_from_toml(const char *file_name){
    
    FILE *f = NULL;
    char errbuf[1000];
    int l_file_names = 300;
    simulation_configuration_t *conf = (simulation_configuration_t*)malloc(sizeof(simulation_configuration_t));

    // define tables and variables to store data from configuration file (toml format)
    toml_table_t *tbl, *tbl_general, *tbl_simulation, *tbl_dataset, *tbl_output, *tbl_network;
    
    // open configuration file and convert to TOML structure
    open_file(&f, file_name);
    tbl = toml_parse_file(f, errbuf, l_file_names);
    fclose(f);


    /* get sections from file */
    tbl_general = toml_table_table(tbl, "general");
    tbl_simulation = toml_table_table(tbl, "simulation");
    tbl_dataset = toml_table_table(tbl, "dataset");
    tbl_output = toml_table_table(tbl, "output");
    tbl_network = toml_table_table(tbl, "network");

    // load sections
    conf = load_general_section_from_toml(conf, tbl_general);
    conf = load_simulation_section_from_toml(conf, tbl_simulation);
    conf = load_dataset_section_from_toml(conf, tbl_dataset);
    conf = load_output_section_from_toml(conf, tbl_output);
    conf = load_network_section_from_toml(conf, tbl_network);
    
    // return configuration struct
    return conf;
}

network_construction_lists_t* load_network_information_in_lists(simulation_configuration_t *conf){

    FILE *f = NULL, *f_neurons, *f_synapses;
    char errbuf[100];
    size_t i, j;

    // struct to store arrays
    network_construction_lists_t *lists = (network_construction_lists_t*)malloc(sizeof(network_construction_lists_t));

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

    // allocate memory for arrays related to synapses data
    lists->weight_list = (double *)malloc(S * sizeof(double)); // weights
    lists->delay_list = (int *)malloc(S * sizeof(int)); // delays
    lists->training_zones = (int *)malloc(S * sizeof(int)); // learnin rules

    // allocate memory for connectibity
    lists->synaptic_connections = (int **)malloc(N * sizeof(int *)); // +1, since the input layer is stored too


    /* load [Neurons] section */
    // load whether different parameters are included in the input file
    v_thres = toml_table_int(tbl_neurons, "v_thres");
    v_rest = toml_table_int(tbl_neurons, "v_rest");
    t_refract = toml_table_int(tbl_neurons, "t_refract");
    res = toml_table_int(tbl_neurons, "resistance");

    if(v_thres.ok && v_thres.u.i == 1){
        for(i=0; i<N; i++)
            fscanf(f_neurons, "%lf", &((lists->v_thres_list)[i]));
    }
    else{
        for(i=0; i<N; i++)
            lists->v_thres_list[i] = 1.0;
    }

    if(v_rest.ok && v_rest.u.i == 1){
        for(i=0; i<N; i++)
            fscanf(f_neurons, "%lf", &((lists->v_rest_list)[i]));
    }
    else{
        for(i=0; i<N; i++)
            lists->v_rest_list[i] = 0.0;
    }

    if(t_refract.ok && t_refract.u.i == 1){
        for(i=0; i<N; i++)
            fscanf(f_neurons, "%d", &((lists->r_time_list)[i]));
    }
    else{
        for(i=0; i<N; i++)
            lists->r_time_list[i] = 1;
    }

    if(res.ok && res.u.i == 1){
        for(i=0; i<N; i++)
            fscanf(f_neurons, "%lf", &((lists->R_list)[i]));
    }
    else{
        for(i=0; i<N; i++)
            lists->R_list[i] = 1.0;
    }

    fclose(f_neurons);
    

    /* Synapses section */
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
    else{
        for(i=0; i<S; i++){
            lists->delay_list[i] = 1;
        }
        lists->max_delay = 1;
    }

    if(weight.ok && weight.u.i == 1){
        for(i=0; i<S; i++){
            fscanf(f_synapses, "%lf", &(lists->weight_list[i]));
        }
    }
    else{
        lists->weight_list[i] = 0.25;
    }

    if(training_zone.ok && training_zone.u.i == 1){
        for(i=0; i<S; i++){
            fscanf(f_synapses, "%d", &(lists->training_zones[i]));
        }
    }
    else{
        lists->training_zones[i] = 0;
    }

    printf(" >> Synapses loaded from file!\n");
    fflush(stdout);


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

GPU_dataset_t* load_dataset_from_file_cpu(const char *file_name, const char *labels_file_name, size_t n_samples, simulation_configuration_t *conf){

    size_t i, j, l;
    FILE *f = NULL;

    open_file(&f, file_name);

    // allocate memory for dataset
    GPU_dataset_t *dataset = (GPU_dataset_t*)malloc(sizeof(GPU_dataset_t));

    // load general dataset information from configuration struct
    //dataset->type = 0;
    dataset->type = 1; // frequencies
    
    dataset->n_classes = conf->n_classes;
    dataset->n_samples = n_samples;
    dataset->n_features = conf->input_size;
    dataset->n_spikes = 0; // counted during loading

    // allocate memory for arrays
    dataset->n_spikes_per_feature = (size_t*)malloc(dataset->n_samples * dataset->n_features * sizeof(size_t)); 
    dataset->sample_offset = (size_t*)malloc(dataset->n_samples * sizeof(size_t));
    dataset->feature_offset = (size_t*)malloc(dataset->n_samples * dataset->n_features * sizeof(size_t));
    

    // count number of spikes in the dataset
    size_t **tmp_spikes = (size_t **)malloc(dataset->n_samples * dataset->n_features * sizeof(size_t*));
    size_t offset = 0;
    size_t n_spikes;
    for(i = 0; i<dataset->n_samples; i++){
        
        // set sample offset
        dataset->sample_offset[i] = offset;

        // loop over features
        for(j = 0; j<dataset->n_features; j++){

            // set feature offset
            dataset->feature_offset[i * dataset->n_features + j] = offset;

            // scan number of spikes of the feature
            fscanf(f, "%zu", &(n_spikes));

            // update offset
            offset += n_spikes;
            dataset->n_spikes += n_spikes;

            // store number of spikes
            dataset->n_spikes_per_feature[i * dataset->n_features + j] = n_spikes;

            // store the spike times of the feature in tmp_spikes
            tmp_spikes[i * dataset->n_features + j] = (size_t*)malloc(n_spikes * sizeof(size_t));
            for(l = 0; l < n_spikes; l++){
                fscanf(f, "%zu", &(tmp_spikes[i * dataset->n_features + j][l]));
            }
        }
    }


    // copy spikes to the dataset struct, compute frequencies and store first spike time
    dataset->spikes = (size_t*)malloc(dataset->n_spikes * sizeof(size_t));
    dataset->freq = (size_t*)malloc(dataset->n_samples * dataset->n_features * sizeof(size_t));
    dataset->first_spk = (size_t*)malloc(dataset->n_samples * dataset->n_features * sizeof(size_t));

    size_t next_spike = 0;
    for(i = 0; i<dataset->n_samples; i++){

        // loop over features
        for(j = 0; j<dataset->n_features; j++){


            n_spikes = dataset->n_spikes_per_feature[i * dataset->n_features + j];

            for(l = 0; l < n_spikes; l++){

                dataset->spikes[next_spike] = tmp_spikes[i * dataset->n_features + j][l];
                next_spike ++;
            }

            // store first spike and frequency
            if(n_spikes > 0){
                dataset->freq[i * dataset->n_features + j] = conf->max_input_spikes / n_spikes; // spikes each freq time steps
                dataset->first_spk[i * dataset->n_features + j] = dataset->spikes[dataset->feature_offset[i * dataset->n_features + j]];
            }
            else{
                dataset->freq[i * dataset->n_features + j] = 0; // spikes each freq time steps
                dataset->first_spk[i * dataset->n_features + j] = 0;
            }

            // deallocate memory
            free(tmp_spikes[i * dataset->n_features + j]);
        }
    }

    // free temporaly allocated memory
    free(tmp_spikes);

    // return dataset
    return dataset;
}