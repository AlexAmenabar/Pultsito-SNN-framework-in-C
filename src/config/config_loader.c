#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h> 
#include <ctype.h>
#include <time.h>

#include "toml_c/toml.h"

#include "config/config_loader.h"
#include "utils.h"

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
        neuron_type.u.i = 1; // LIF neuron
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
        load_network.u.i = 0; // load
    if(!load_dataset.ok)
        load_dataset.u.i = 0; // load


    // load information in configuration struct
    conf->neuron_type  = neuron_type.u.i;
    conf->n_process    = (size_t)n_process.u.i;
    conf->cuda         = cuda.u.i;
    conf->multi_gpu    = multigpu.u.i;
    conf->learn        = learn.u.i;
    conf->batch_size   = (size_t)batch_size.u.i;
    conf->thrN         = (size_t)thrN.u.i;
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
    conf->time_steps       = (size_t)time_steps.u.i;
    conf->max_input_spikes = (size_t)max_input_spikes.u.i;

    return conf;
}

simulation_configuration_t* load_dataset_section_from_toml(simulation_configuration_t *conf, toml_table_t *tbl){

    // [dataset] section
    toml_value_t dataset, labels, n_samples, input_size, labels_provided,
                n_classes, epochs, shuffle_samples;
    
    
    /* read [dataset] section */
    labels_provided =  toml_table_int(tbl, "labels_provided");
    dataset         =  toml_table_string(tbl, "dataset");
    labels          =  toml_table_string(tbl, "labels");
    n_samples       =  toml_table_int(tbl, "n_samples");

    n_classes       =  toml_table_int(tbl, "num_classes");
    input_size      =  toml_table_int(tbl, "input_size");

    shuffle_samples =  toml_table_int(tbl, "shuffle_samples");
    epochs          =  toml_table_int(tbl, "epochs");

    // check that all is correctly loaded
    if(conf->load_dataset == 0){
        if(!dataset.ok){
            
            printf(" >> File to load dataset from must be provided when load_dataset is 0! Exiting.\n");
            fflush(stdout);
            exit(1);
        }
        if(!labels_provided.ok){
        
            printf(" >> WARNING: labels not provided!\n");
            labels_provided.u.i = 0;
            fflush(stdout);
        }
        if(!labels.ok && labels_provided.u.i == 1){
            
            printf(" >> Dataset labels must be provided when load_datase and labels_provided are 0! Exiting.\n");
            fflush(stdout);
            exit(1);   
        }
    }
    // [TODO]
    else if(conf->load_dataset == 1){

        printf(" > load_dataset = 1 not implemented yet. Exiting.\n");
    }
    else if(conf->load_dataset == 2){

        printf(" > load_dataset = 2 not implemented yet. Exiting.\n");
    }

    if(!n_samples.ok){
        
        printf(" >> The number of samples in the dataset must be provided when load_dataset is 1! Exiting\n");
        fflush(stdout);
        exit(1);    
    }
    
    if(!n_classes.ok && conf->load_dataset == 0 && labels_provided.u.i == 1){
        
        printf(" >> WARNING: number of classes in the dataset not provided. The application will try to count them automatically.");
        n_classes.u.i = 0;
    }
    if(!epochs.ok){
        
        printf(" >> WARNING: number of epochs to simulate not provided. Setting 1 as default value.");
        epochs.u.i = 1;
    }

    if(!input_size.ok && conf->load_dataset == 0){
        printf(" >> Number of features in the dataset not provided! Exiting.\n");
        fflush(stdout);
        exit(1);
    }
    if(!shuffle_samples.ok){
        shuffle_samples.u.i = 0; // do not shuffle
    }

    // store dataset information
    conf->labels_provided = labels_provided.u.i;
    conf->dataset         = dataset.u.s;
    conf->labels          = labels.u.s;
    conf->n_samples       = (size_t)n_samples.u.i;

    conf->n_classes       = (size_t)n_classes.u.i;
    conf->input_size      = (size_t)input_size.u.i;

    conf->shuffle_samples = shuffle_samples.u.i;
    conf->epochs          = (size_t)epochs.u.i;

    return conf;
}

// [TODO]: revise
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

    conf->store_generated_spikes   = store_generated_spikes.u.i;
    conf->generated_spikes_file    = generated_spikes_file.u.s;

    conf->store_execution_times    = store_execution_times.u.i;
    conf->execution_times_file     = execution_times_file.u.s;

    conf->store_n_spikes           = store_n_spikes.u.i;
    conf->n_spikes_per_neuron_file = n_spikes_per_neuron_file.u.s;

    conf->store_network            = store_network.u.i;
    if(conf->store_network == 1)
        conf->store_network_file   = store_network_file.u.s;

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
    if(!network_file.ok && conf->load_network == 0)
    {
        printf(" > The file to load the network from must be provided when load_network is 1! Exiting.\n");
        fflush(stdout);
        exit(1);
    }

    if(conf->load_network == 0){
        
        if(network_file.ok){
            conf->network_file = network_file.u.s;
        }
        else{
            printf(" > network_file for loading the general network properties not provided! Exiting.\n");
            exit(1);
        }
        if(network_neurons_file.ok){
            conf->network_neurons_file = network_neurons_file.u.s;
        }
        else{
            
            printf(" > network_file_neurons for loading the neurons properties not provided! Exiting.\n");
            exit(1);
        }
        if(network_synapses_file.ok){
            conf->network_synapses_file = network_synapses_file.u.s;
        }
        else{
            
            printf(" > network_file_synapses for loading the synapses properties not provided! Exiting.\n");
            exit(1);
        }
    }
    // [TODO]
    else if(conf->load_network == 1){

        printf(" > load_network = 1 not implemented yet. Exiting.\n");
    }
    // [TODO]
    else if(conf->load_network == 2){

        printf(" > load_network = 2 not implemented yet. Exiting.\n");
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