#include "snn_library.h"
#include "load_data.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h> 
#include <ctype.h>
#include <toml_c/toml-c.h>

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


// Function to load a network into the simulation from a network definition file
void load_network_information(const char *file_name, spiking_nn_t *snn, network_construction_lists_t *lists, simulation_configuration_t *conf) {
    
    FILE *f = NULL, *f_neurons, *f_synapses;
    char errbuf[100], *original_file_name, *tmp_file_name, *copied_file_name;
    int i, j, l;

    // define table and parameters variables
    toml_table_t *tbl, *tbl_general, *tbl_neurons, *tbl_synapses;

    toml_array_t *behaviour_lst, *v_thres_lst, *v_rest_lst, *t_refract_lst, *res_lst, *input_synapses_lst, 
                *output_synapses_lst, *latency_lst, *weights_lst, *training_zones_lst,
                *connection_lst_lst, *connection_lst;

    toml_value_t n_neurons, n_input_neurons, n_output_neurons, n_synapses, n_input_synapses, n_output_synapses,
                behaviour, v_thres, v_rest, t_refract, res, latency, training_zone, weight,
                n_connections, neuron_id, n_synapses_to_neuron, network_is_separated;


    // open TOML file
    open_file(&f, file_name); // TOML file

    // read TOML file
    tbl = toml_parse_file(f, errbuf, 100);
    
    // close the file as the information has been readed
    fclose(f);
    
    // get sections from file
    tbl_general = toml_table_table(tbl, "general");
    tbl_neurons = toml_table_table(tbl, "neurons");
    tbl_synapses = toml_table_table(tbl, "synapsis");
    printf(" >> General information loaded\n");
    fflush(stdout);
    
    /* General section */

    n_neurons = toml_table_int(tbl_general, "neurons");
    n_input_neurons = toml_table_int(tbl_general, "input_neurons");
    n_output_neurons = toml_table_int(tbl_general, "output_neurons");
    n_synapses = toml_table_int(tbl_general, "synapsis");
    n_input_synapses = toml_table_int(tbl_general, "input_synapsis");
    n_output_synapses = toml_table_int(tbl_general, "output_synapsis");
    network_is_separated = toml_table_int(tbl_general, "network_is_separated"); // it indicates if the network is separated in several files
    printf(" >> General section loaded\n");
    fflush(stdout);

    // NETWORK IS SEPARATED IS TEMPORAL UNTIL SOMETHING BETTER IS FOUND
    // if network_is_separated == 1 network must be loaded from more than one file
    if(network_is_separated.ok && network_is_separated.u.i == 1){

        // read files
        open_file(&f_neurons, conf->network_neurons_file); // TOML file
        open_file(&f_synapses, conf->network_synapses_file); // TOML file
    }


    // check that all the information has been loaded correctly
    if(!(n_neurons.ok && n_input_neurons.ok && n_output_neurons.ok && n_synapses.ok && n_input_synapses.ok && n_output_synapses.ok)){
        printf("The number of neurons, input neurons, output neurons, synapses, input synapses and output synapses must be provided in the network file!");
        exit(1);
    }

    // if correctly readed, copy to the snn structure
    snn->n_neurons = n_neurons.u.i;
    snn->n_input = n_input_neurons.u.i;
    snn->n_output = n_output_neurons.u.i;
    snn->n_synapses = n_synapses.u.i;
    snn->n_input_synapses = n_input_synapses.u.i;
    snn->n_output_synapses = n_output_synapses.u.i;


#ifdef DEBUG
    printf("n_neurons = %d\n", snn->n_neurons);
    printf("n_input_neurons = %d\n", snn->n_input);
    printf("n_output_neurons = %d\n", snn->n_output);
    printf("n_synapses = %d\n", snn->n_synapses);
    printf("n_input_synapses = %d\n", snn->n_input_synapses);
    printf("n_output_synapses = %d\n", snn->n_output_synapses);
#endif


    // reserve memory for lists related to neurons
    lists->neuron_excitatory = (int *)malloc(snn->n_neurons * sizeof(int));
    lists->r_time_list = (int *)malloc(snn->n_neurons * sizeof(int));
    lists->v_thres_list = (double *)malloc(snn->n_neurons * sizeof(double));
    lists->v_rest_list = (double *)malloc(snn->n_neurons * sizeof(double));
    lists->R_list = (double *)malloc(snn->n_neurons * sizeof(double));

    lists->weight_list = (double *)malloc(snn->n_synapses * sizeof(double));
    lists->delay_list = (int *)malloc(snn->n_synapses * sizeof(int));
    lists->training_zones = (int *)malloc(snn->n_synapses * sizeof(int));

    lists->synaptic_connections = (int **)malloc((snn->n_neurons + 1) * sizeof(int *)); // + 2, one row input layer and the other the output layet



    /* Neurons section */

    behaviour = toml_table_int(tbl_neurons, "behaviour");
    v_thres = toml_table_double(tbl_neurons, "v_thres");
    v_rest = toml_table_double(tbl_neurons, "v_rest");
    t_refract = toml_table_int(tbl_neurons, "t_refract");
    res = toml_table_double(tbl_neurons, "resistance");
    printf(" >> Neurons section loaded\n");
    fflush(stdout);


    // if network information is not separated into more than one files
    if(!network_is_separated.ok || network_is_separated.ok && network_is_separated.u.i != 1){
        
        // read information from TOML file
        behaviour_lst = toml_table_array(tbl_neurons, "behaviour_list");
        v_thres_lst = toml_table_array(tbl_neurons, "v_thres_list");
        v_rest_lst = toml_table_array(tbl_neurons, "v_rest_list");
        t_refract_lst = toml_table_array(tbl_neurons, "t_refract_list");
        res_lst = toml_table_array(tbl_neurons, "res_list");
        input_synapses_lst = toml_table_array(tbl_neurons, "input_synapsis");
        output_synapses_lst = toml_table_array(tbl_neurons, "output_synapsis");

        // load information into snn structure
        for(i=0; i<snn->n_neurons; i++){
            // load data into value variables
            behaviour = toml_array_int(behaviour_lst, i);
            t_refract = toml_array_int(t_refract_lst, i);
            v_thres = toml_array_double(v_thres_lst, i);
            v_rest = toml_array_double(v_rest_lst, i);
            res = toml_array_double(res_lst, i);

            // analyze if this parameters are not provided when it is supposed that they are provided
            // if data was provided load it directly on the network
            if(!behaviour.ok && conf->behaviours_provided == 1 || conf->behaviours_provided == 0){
                printf("Following configuration file, behaviours for neurons must be provided, setting 1 (excitatory)\n");
                behaviour.u.i = 1;
            }
            if(!v_thres.ok && conf->thresholds_provided == 1 || conf->thresholds_provided == 0){
                printf("Following configuration file, thresholds for neurons must be provided, setting 150\n");
                v_thres.u.d = 150;
            }
            if(!v_rest.ok && conf->v_rests_provided == 1 || conf->v_rests_provided == 0){
                printf("Following configuration file, resting potentials for neurons must be provided, setting 50\n");
                v_rest.u.d = 50;
            }
            if(!res.ok && conf->res_provided == 1 || conf->res_provided == 0){
                printf("Following configuration file, resistances for neurons not proveided, setting 1\n");
                res.u.d = 1;
            }
            if(!t_refract.ok && conf->refract_times_provided == 1 || conf->refract_times_provided == 0){
                printf("Following configuration file, refractary times for neurons must be provided, setting 3\n");
                t_refract.u.i = 3;
            }

            // load data into lists structure
            lists->neuron_excitatory[i] = behaviour.u.i;
            lists->v_thres_list[i] = v_thres.u.d;
            lists->v_rest_list[i] = v_rest.u.d;
            lists->r_time_list[i] = t_refract.u.i;
            lists->R_list[i] = res.u.d; 
        }
    }
    // if it is separated, read the information from the other file
    else{
        // TODO: I am not handling the case in which some parameter is not provided
        for(i=0; i<snn->n_neurons; i++){
            fscanf(f_neurons, "%d", &((lists->neuron_excitatory)[i]));
        }
        for(i=0; i<snn->n_neurons; i++){
            fscanf(f_neurons, "%lf", &((lists->v_thres_list)[i]));
        }
        for(i=0; i<snn->n_neurons; i++){
            fscanf(f_neurons, "%lf", &((lists->v_rest_list)[i]));
        }
        for(i=0; i<snn->n_neurons; i++){
            fscanf(f_neurons, "%d", &((lists->r_time_list)[i]));
        }

        // TODO: not correct
        for(i = 0; i<snn->n_neurons; i++)
            lists->R_list[i] = 1;

        fclose(f_neurons);
    }
    

    /* Synapses section */

    // if network information is not separated into more than one files
    if(!network_is_separated.ok || network_is_separated.ok && network_is_separated.u.i != 1){
        latency_lst = toml_table_array(tbl_synapses, "latency_list"); // load latencies
        weights_lst = toml_table_array(tbl_synapses, "weights"); // load weights
        training_zones_lst = toml_table_array(tbl_synapses, "training_zones_list"); // load training zones
        connection_lst_lst = toml_table_array(tbl_synapses, "connections"); // load connections
        

        // load information into snn structure
        for(i=0; i<snn->n_synapses; i++){
            // load data into value variables
            latency = toml_array_int(latency_lst, i);
            weight = toml_array_double(weights_lst, i);
            training_zone = toml_array_int(training_zones_lst, i);

            // analyze if this parameters are not provided when it is supposed that they are provided
            // if data was provided load it directly on the network
            if(!latency.ok && conf->delays_provided == 1 || conf->delays_provided == 0){
                printf("Following configuration file, latencies for synapses must be provided, setting 1\n");
                latency.u.i = 1;
            }
            if(!weight.ok && conf->weights_provided == 1 || conf->weights_provided == 0){
                printf("Following configuration file, weights for synapses must be provided, setting 100\n");
                weight.u.d = 100;
            }
            if(!training_zone.ok && conf->training_zones_provided == 1 || conf->training_zones_provided == 0){
                printf("Following configuration file, training zones for neurons must be provided, setting 0 (normal STDP)\n");
                training_zone.u.i = 0;
            }

            // load data into lists structure
            (lists->weight_list)[i] = (double)weight.u.d;
            (lists->delay_list)[i] = (int)latency.u.i;
            (lists->training_zones)[i] = (int)training_zone.u.i;
        }
    }
    // if it is separated, read the information from the other file
    else{
        for(i=0; i<snn->n_synapses; i++){
            fscanf(f_synapses, "%d", &((lists->delay_list)[i]));
        }
        for(i=0; i<snn->n_synapses; i++){
            fscanf(f_synapses, "%lf", &((lists->weight_list)[i]));
        }
        for(i=0; i<snn->n_synapses; i++){
            fscanf(f_synapses, "%d", &((lists->training_zones)[i]));
        }
    }


    // load connectivity information
    if(!network_is_separated.ok || network_is_separated.ok && network_is_separated.u.i != 1){
        
        // load connections (list of lists)
        for(i=0; i<snn->n_neurons+1; i++){

            connection_lst = toml_array_array(connection_lst_lst, i);
            n_connections = toml_array_int(connection_lst, 0);
            
            // check that data has been correctly loaded
            if(!n_connections.ok){
                printf(" > Connection list is incorrect. Exiting.\n");
                exit(1);
            }

            //
            (lists->synaptic_connections)[i] = malloc((n_connections.u.i * 2 + 1) * sizeof(int)); // for each connection the neuron id and the number of synapses must be stored
            (lists->synaptic_connections)[i][0] = n_connections.u.i;
        
            for(int j = 0; j<n_connections.u.i; j++){
                neuron_id = toml_array_int(connection_lst, j * 2 + 1);
                n_synapses_to_neuron = toml_array_int(connection_lst, j * 2 + 2);

                // check that the information have been correctly loaded
                if(!(neuron_id.ok && n_synapses_to_neuron.ok)){
                    printf("Connection list data is incorrect. Exiting\n");
                    exit(1);
                }

                (lists->synaptic_connections)[i][j * 2 + 1] = neuron_id.u.i;
                (lists->synaptic_connections)[i][j * 2 + 2] = n_synapses_to_neuron.u.i;
            } 
        }
    }
    // if network information is separated into more than one file
    else{
        int number_connections;
        for(i=0; i<(snn->n_neurons + 1); i++){ // network input synapses are loaded first and output synapses last
            fscanf(f_synapses, "%d", &number_connections);

            // alloc memory
            (lists->synaptic_connections)[i] = malloc((number_connections * 2 + 1) * sizeof(int)); // for each connection the neuron id and the number of synapses must be stored
            (lists->synaptic_connections)[i][0] = number_connections;

            for(j = 0; j<number_connections; j++){
                fscanf(f_synapses, "%d", &((lists->synaptic_connections)[i][j * 2 + 1])); // number of synapses connected to that neuron
                fscanf(f_synapses, "%d", &((lists->synaptic_connections)[i][j * 2 + 2])); // number of synapses connected to that neuron
            }
        }
        fclose(f_synapses);
    }
    printf(" >> Synapses section loaded\n");
    fflush(stdout);
    
    // free memory
    /*toml_free(tbl);
    toml_free(tbl_general);
    toml_free(tbl_neurons);
    toml_free(tbl_synapses);*/
}



// I THINK THIS SHOULDN'T BE HERE // NOW WITH NEURON IT IS NOT CORRECT FOR GENERALIZING
/// @brief Function to load data into the SNN structure 
/// @param file_name File name to load spikes from
/// @param snn SNN structure to load spikes in
void load_input_spike_trains_on_snn(const char *file_name, spiking_nn_t *snn){
    
    FILE *f = NULL;
    int i, j, n_spikes;

    // open file
    open_file(&f, file_name);


    // the first synapses of the network are the input ones
    for(i = 0; i<snn->n_input; i++){

        // read number of spikes
        fscanf(f, "%d", &n_spikes);

        // load spikes for i neuron
        for(j=0; j<n_spikes; j++){
            fscanf(f, "%d", &(snn->input_lif_neurons[i].spike_times_arr[j]));
        }

        // refresh spikes index for synapse // TODO: this should be introduced as stream?
        //snn->input_lif_neurons[i].spike_times_arr += n_spikes;
        snn->input_lif_neurons[i].last_spike = n_spikes;
    }

    // close file
    fclose(f);
}


/// @brief [Deprecated] Function to read a file name
/// @param f 
/// @param max_length 
/// @param file_name 
/// @return 
int read_file_name(FILE *f, int max_length, char *file_name){
    char ch;
    int length = 0;

    // skip any leading whitespace
    while ((ch = fgetc(f)) != EOF && isspace(ch));

    // read word and count number of characters
    while (ch != EOF && !isspace(ch) && length < max_length) {
        file_name[length] = ch;
        ch = fgetc(f);
        length++;
    }

    printf("Length of first word: %d\n", length);
    printf("Word: %s\n", file_name);
    if(length >= max_length) return 1; // error

    return 0;
}

/// @brief [Deprecated] Function to load configuration parameters
/// @param file_name 
/// @param conf 
/// @return 
int load_configuration_params(const char *file_name, simulation_configuration_t *conf){
    FILE *f = NULL;

    int err = 0;

    // allocate memory for file names
    conf->network_file = malloc(100 * sizeof(char));
    conf->spike_times_file = malloc(100 * sizeof(char));
    conf->times_file = malloc(100 * sizeof(char));
    conf->n_spikes_file = malloc(100 * sizeof(char));
    conf->final_network_file = malloc(100 * sizeof(char));

    // read parameters
    open_file(&f, file_name);

    fscanf(f, "%d", &(conf->simulation_type));
    fscanf(f, "%d", &(conf->neuron_type));
    err = read_file_name(f, 100, conf->network_file);
    if(err == 1) 
        return err; //err


    fscanf(f, "%d", &(conf->n_process));
    fscanf(f, "%d", &(conf->store));
    if(conf->store == 1){
        err = read_file_name(f, 100, conf->final_network_file);
        if(err == 1)
            return err; //err
    }

    // output files
    err = read_file_name(f, 100, conf->spike_times_file);
    if(err == 1) 
        return err; //err

    err = read_file_name(f, 100, conf->times_file);
    if(err == 1) 
        return err; //err

    err = read_file_name(f, 100, conf->n_spikes_file);
    if(err == 1) 
        return err; //err    
    
    // IFDEF SIMULATION...
    conf->input_spikes_file = malloc(100 * sizeof(char));
    err = read_file_name(f, 100, conf->input_spikes_file);
    if(err == 1) 
        return err; //err    
    printf("Hasta aquío bien\n");

    fscanf(f, "%d", &(conf->time_steps));
    printf("Vaya por dios \n");
    // ELSE...
}

// Function to load the simulation configuration file
int load_configuration_params_from_toml(const char *file_name, simulation_configuration_t *conf){
    
    FILE *f = NULL;
    char errbuf[1000];
    int i, l, l_file_names = 300;

    // define table and variables to store configuration parameters
    toml_table_t *tbl, *tbl_general, *tbl_simulation, *tbl_samples, *tbl_output, *tbl_network;
    
    toml_value_t execution_type, neuron_type, execution_obj, n_process, cuda, learn, 
                time_steps, input_file,
                dataset, dataset_name, num_classes, epochs, n_samples,
                generated_spikes, execution_times, spikes_per_neuron, store_network, store_network_file,
                network, network_neurons, network_synapses, behaviours, delays, weights, training_zones, 
                thresh, v_rest, t_refract, R;


    // open configuration 
    open_file(&f, file_name); // TOML file

    // read TOML field in the file
    tbl = toml_parse_file(f, errbuf, l_file_names);

    // close file
    fclose(f);


    /* get sections from file */
    tbl_general = toml_table_table(tbl, "general");
    tbl_simulation = toml_table_table(tbl, "simulation");
    tbl_samples = toml_table_table(tbl, "samples");
    tbl_output = toml_table_table(tbl, "output");
    tbl_network = toml_table_table(tbl, "network");

    // read subsections in general section
    execution_type = toml_table_int(tbl_general, "execution_type");
    neuron_type = toml_table_int(tbl_general, "neuron_type");
    execution_obj = toml_table_int(tbl_general, "execution_obj");
    n_process = toml_table_int(tbl_general, "n_process");
    cuda = toml_table_int(tbl_general, "cuda");
    learn = toml_table_int(tbl_general, "learn");

    // if something is missing in configuration file, set default values
    if(!execution_type.ok)
        execution_type.u.i = 0; // clock-based

    if(!neuron_type.ok)
        neuron_type.u.i = 0; // LIF neuron

    if(!execution_obj.ok)
        execution_obj.u.i = 0; // biological simulation

    if(!n_process.ok)
        n_process.u.i = 1; // serial execution

    if(!cuda.ok)
        cuda.u.i = 0; // no cuda

    if(!learn.ok)
        learn.u.i = 1; // not learn


    // load information into configuration structure
    conf->simulation_type = execution_type.u.i;
    conf->neuron_type = neuron_type.u.i;
    conf->simulation_obj = execution_obj.u.i;
    conf->n_process = n_process.u.i;
    conf->cuda = cuda.u.i;
    conf->learn = learn.u.i;
    



    /* simulation section in TOML file */ 
    // load data from TOML file
    time_steps = toml_table_int(tbl_simulation, "time_steps");
    input_file = toml_table_string(tbl_simulation, "input_file");
    n_samples = toml_table_int(tbl_simulation, "n_samples");
    epochs = toml_table_int(tbl_simulation, "epochs");
    printf(" > Simulation section loaded\n");
    fflush(stdout);

    // set default values if something is missing
    if(!time_steps.ok)
    {
        printf(" >> Amount of time-steps haven't been provided, setting default value: 1000\n");
        time_steps.u.i = 1000;
    }

    if(!input_file.ok){
        printf(" >> A input file must be provided for input spikes!\n");
        exit(1);
    }

    if(!n_samples.ok)
    {
        printf(" >> Setting number of samples in 1\n");
        n_samples.u.i = 1;
    }

    if(!epochs.ok){
        printf(" >> Number of epochs not provided! Setting 1.\n");
        epochs.u.i = 1;
    }
    
    // load information into config structure
    conf->time_steps = time_steps.u.i;
    l = input_file.u.sl; // get string length
    conf->input_spikes_file = malloc(l * sizeof(char)); // allocate memory for string
    conf->input_spikes_file = input_file.u.s; // store file name
    conf->n_samples = n_samples.u.i;
    conf->epochs = epochs.u.i;

    printf(" > Simulation section copied\n");
    fflush(stdout);


    /* output section in TOML file */
    generated_spikes = toml_table_string(tbl_output, "generated_spikes");
    execution_times = toml_table_string(tbl_output, "execution_times");
    spikes_per_neuron = toml_table_string(tbl_output, "spikes_per_neuron");
    store_network = toml_table_int(tbl_output, "store_network"); // whether to store or not the file
    store_network_file = toml_table_string(tbl_output, "store_network_file"); // file to store the network

    // check that everything is loaded
    if(!generated_spikes.ok){
        printf("A file to store generated spikes must be provided!\n");
        exit(1);
    }
    if(!execution_times.ok){
        printf("A file to store execution times must be provided!\n");
        exit(1);
    }
    if(!spikes_per_neuron.ok){
        printf("A file to store the number of spikes generated per neuron must be provided!\n");
        exit(1);
    }

    if(!store_network.ok){
        store_network.u.i = 0; // default option
    }
    else if(store_network.u.i == 1 && !store_network_file.ok){
        printf("The file name to store the final network must be provided!\n");
        exit(1);
    }

    printf(" > Output section loaded\n");
    fflush(stdout);


    // allocate memory for strings and store them into lists structure
    l = generated_spikes.u.sl;
    conf->spike_times_file = malloc(l * sizeof(char));
    conf->spike_times_file = generated_spikes.u.s;

    l = execution_times.u.sl;
    conf->times_file = malloc(l * sizeof(char));
    conf->times_file = execution_times.u.s;

    l = spikes_per_neuron.u.sl;
    conf->n_spikes_file = malloc(l * sizeof(char));
    conf->n_spikes_file = spikes_per_neuron.u.s;

    conf->store = store_network.u.i;
    if(conf->store == 1){
        l = store_network_file.u.sl;
        conf->final_network_file = malloc(l * sizeof(char));
        conf->final_network_file = store_network_file.u.s;
    }

    printf(" > Output section copied\n");
    fflush(stdout);


    /* network section */
    // parameters that are provided for the network
    network = toml_table_string(tbl_network, "network_file");
    network_neurons = toml_table_string(tbl_network, "network_neurons_file");
    network_synapses = toml_table_string(tbl_network, "network_synapses_file");
    behaviours = toml_table_int(tbl_network, "behaviours");
    delays = toml_table_int(tbl_network, "delays");
    weights = toml_table_int(tbl_network, "weights");
    training_zones = toml_table_int(tbl_network, "training_zones");
    thresh = toml_table_int(tbl_network, "thresh");
    v_rest = toml_table_int(tbl_network, "v_rests");
    t_refract = toml_table_int(tbl_network, "t_refract");
    R = toml_table_int(tbl_network, "resistance");
    // TODO: more parameters should be added in the future

    printf(" > Network section loaded\n");
    fflush(stdout);

    // if something is missing in configuration file, set default value
    if(!network.ok)
    {
        printf("The file to load the network must be provided\n");
        exit(1);
    }

    // if there is not in the file, by default set as not proveded
    if(!behaviours.ok){
        behaviours.u.i = 0;
    }
    
    if(!delays.ok){
        delays.u.i = 0;
    }

    if(!weights.ok){
        weights.u.i = 0;
    }    
    
    if(!training_zones.ok){
        training_zones.u.i = 0;
    }    
    
    if(!thresh.ok){
        thresh.u.i = 0;
    }    

    if(!v_rest.ok){
        v_rest.u.i = 0;
    }    

    if(!R.ok){
        R.u.i = 0; 
    }

    if(!t_refract.ok){
        t_refract.u.i = 0;
    }    
    
    // copy information into config structure
    l = network.u.sl;
    conf->network_file = malloc(l * sizeof(char));
    conf->network_neurons_file = malloc(l * sizeof(char));
    conf->network_synapses_file = malloc(l * sizeof(char));
    conf->network_file = network.u.s;
    if(network_neurons.ok)
        conf->network_neurons_file = network_neurons.u.s;
    if(network_synapses.ok)
        conf->network_synapses_file = network_synapses.u.s;
    
    conf->behaviours_provided = behaviours.u.i;
    conf->delays_provided = delays.u.i;
    conf->weights_provided = weights.u.i;
    conf->training_zones_provided = training_zones.u.i;
    conf->thresholds_provided = thresh.u.i;
    conf->v_rests_provided = v_rest.u.i;
    conf->refract_times_provided = t_refract.u.i;
    conf->res_provided = R.u.i;


    // free memory
    /*toml_free(tbl);
    toml_free(tbl_general);
    toml_free(tbl_simulation);
    toml_free(tbl_samples);
    toml_free(tbl_output);
    toml_free(tbl_network);*/
}



// TODO: this functions should be adapted for several samples

// Functions to store results and data into files
void store_results(simulation_results_t *results, simulation_configuration_t *conf, spiking_nn_t *snn){

    // store
    store_generated_spikes(results, conf, snn);
    //store_network();
    store_number_of_spikes(results, conf,snn);
    store_times(results, conf, snn);

}

void store_generated_spikes(simulation_results_t *results, simulation_configuration_t *conf, spiking_nn_t *snn){
    
    int i,j;

    // TODO: In this moment this function only stores the first sample results
    simulation_results_per_sample_t *results_per_sample = &(results->results_per_sample[0]);
    FILE *f;

    // file to store generated spikes
    f = fopen(conf->spike_times_file, "w");
    if(f == NULL){
        printf("Error opening the file %s. \n", conf->spike_times_file);
        exit(1);
    }

    // write generated spikes // TODO: for each sample (a file?)
    for (i = 0; i<snn->n_neurons; i++)
    {
        for(j = 0; j<conf->time_steps; j++)
            fprintf(f, "%c", results_per_sample->generated_spikes[i][j]);
        
        fprintf(f, "\n");
    }

    // close file
    fclose(f);
}


void store_network(simulation_results_t *results, simulation_configuration_t *conf, spiking_nn_t *snn){
    
    int i,j;

    // TODO: In this moment this function only stores the first sample results
    simulation_results_per_sample_t *results_per_sample = &(results->results_per_sample[0]);
    FILE *f;

    // file to store generated spikes
    f = fopen(conf->spike_times_file, "w");
    if(f == NULL){
        printf("Error opening the file %s. \n", conf->spike_times_file);
        exit(1);
    }

    // store network

    // close file
    fclose(f);
}


void store_number_of_spikes(simulation_results_t *results, simulation_configuration_t *conf, spiking_nn_t *snn){
    
    int i,j;

    // TODO: In this moment this function only stores the first sample results
    simulation_results_per_sample_t *results_per_sample = &(results->results_per_sample[0]);
    FILE *f;

    // file to store generated spikes
    f = fopen(conf->n_spikes_file, "w");
    if(f == NULL){
        printf("Error opening the file %s. \n", conf->n_spikes_file);
        exit(1);
    }

    // store number of spikes
    for(i = 0; i<snn->n_neurons; i++)
        fprintf(f, "%d ", results->results_per_sample[0].n_spikes_per_neuron[i]);

    // close file
    fclose(f);
}


void store_times(simulation_results_t *results, simulation_configuration_t *conf, spiking_nn_t *snn){
    
    int i,j;

    // TODO: In this moment this function only stores the first sample results
    simulation_results_per_sample_t *results_per_sample = &(results->results_per_sample[0]);
    FILE *f;

    // file to store generated spikes
    f = fopen(conf->times_file, "w");
    if(f == NULL){
        printf("Error opening the file %s. \n", conf->spike_times_file);
        exit(1);
    }

    // store number of spikes
    fprintf(f, "%lf ", results_per_sample->elapsed_time);
    fprintf(f, "%lf ", results_per_sample->elapsed_time_neurons);
    fprintf(f, "%lf ", results_per_sample->elapsed_time_neurons_input);
    fprintf(f, "%lf ", results_per_sample->elapsed_time_neurons_output);
    fprintf(f, "%lf ", results_per_sample->elapsed_time_synapses);
    fprintf(f, "%lf ", results_per_sample->elapsed_time_synapses_input);
    fprintf(f, "%lf ", results_per_sample->elapsed_time_synapses_output);
    fprintf(f, "%lf \n", results_per_sample->elapsed_time_learning);


    // close file
    fclose(f);
}