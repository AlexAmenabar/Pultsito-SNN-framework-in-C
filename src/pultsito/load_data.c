#include "snn_library.h"
#include "load_data.h"
#include<unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h> 
#include <ctype.h>
#include <toml_c/toml-c.h>
#include <time.h>



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

simulation_configuration_t* load_configuration_params_from_toml(const char *file_name){
    
    FILE *f = NULL;
    char errbuf[1000];
    int l_file_names = 300;
    simulation_configuration_t *conf = (simulation_configuration_t*)malloc(sizeof(simulation_configuration_t));

    // define tables and variables to store data from configuration file (toml format)
    toml_table_t *tbl, *tbl_general, *tbl_simulation, *tbl_dataset, *tbl_output, *tbl_network;
    
    // [general] section
    toml_value_t execution_type, neuron_type, execution_obj, n_process, n_inner_process, cuda, multigpu, 
                learn, encode, batch_size;
    
    // [simulation] section
    toml_value_t time_steps, max_spikes, max_input_spikes;

    // [dataset] section
    toml_value_t train_provided, train_set, train_labels, n_train, 
                test_provided, test_set, test_labels, n_test,
                dataset_name, n_classes, epochs, input_size, shuffle_samples;
    
    // [output] section
    toml_value_t generated_spikes_file, execution_times_file, n_spikes_per_neuron_file, store_network, 
                store_generated_spikes, store_n_spikes, store_execution_times, store_network_file;

    // [network] section
    toml_value_t network_file, network_neurons_file, network_synapses_file, behaviours, delays, weights,
                training_zones, thresh, v_rest, t_refract, R;


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


    /* read [general] section, and set default values in case data is not provided */
    execution_type = toml_table_int(tbl_general, "execution_type"); // clock / event-driven
    neuron_type = toml_table_int(tbl_general, "neuron_type"); // neuron type
    execution_obj = toml_table_int(tbl_general, "execution_obj"); // execution obj. (ML,...)
    n_process = toml_table_int(tbl_general, "n_process"); // number of CPU processes
    n_inner_process = toml_table_int(tbl_general, "n_inner_process");
    cuda = toml_table_int(tbl_general, "cuda"); // simulated on cuda
    multigpu = toml_table_int(tbl_general, "multi_gpu"); // simulated on cuda
    learn = toml_table_int(tbl_general, "learn"); // inference or training
    encode = toml_table_int(tbl_general, "encode"); // encode input or not
    batch_size = toml_table_int(tbl_general, "batch_size");
    
    if(!execution_type.ok)
        execution_type.u.i = 0; // clock-based
    if(!neuron_type.ok)
        neuron_type.u.i = 0; // LIF neuron
    if(!execution_obj.ok)
        execution_obj.u.i = 0; // biological simulation
    if(!n_process.ok)
        n_process.u.i = 1; // serial execution
    if(!n_inner_process.ok)
        n_inner_process.u.i = 1;
    if(!cuda.ok)
        cuda.u.i = 0; // no cuda
    if(!multigpu.ok)
        multigpu.u.i = 0; // no multigpu
    if(!learn.ok)
        learn.u.i = 0; // inference
    if(!encode.ok)
        encode.u.i = 0; // input already encoded
    if(!batch_size.ok)
        batch_size.u.i = 1; // only one sample in the batch

    // load information in configuration struct
    conf->simulation_type = execution_type.u.i;
    conf->neuron_type = neuron_type.u.i;
    conf->simulation_obj = execution_obj.u.i;
    conf->n_process = (size_t)n_process.u.i;
    conf->n_inner_process = (size_t)n_inner_process.u.i;
    conf->cuda = cuda.u.i;
    conf->multi_gpu = multigpu.u.i;
    conf->learn = learn.u.i;
    conf->encode = encode.u.i;
    conf->batch_size = (size_t)batch_size.u.i;



    /* read [simulation] section */
    time_steps = toml_table_int(tbl_simulation, "time_steps"); // simulation time steps
    max_spikes = toml_table_int(tbl_simulation, "max_spikes"); // max spikes for middle spike arrays
    max_input_spikes = toml_table_int(tbl_simulation, "max_input_spikes"); // length for input and output neurons

    if(!time_steps.ok){

        printf(" >> It is necessary to provide simulation time steps!\n");
        fflush(stdout);
        exit(1);
    }
    if(!max_spikes.ok)
        max_spikes.u.i = time_steps.u.i * 10; 
    if(!max_input_spikes.ok)
        max_input_spikes.u.i = time_steps.u.i * 10;

    // load information in configuration struct
    conf->time_steps = (size_t)time_steps.u.i;
    conf->max_spikes = (size_t)max_spikes.u.i; // REVISE THIS AND THE NEXT
    conf->max_input_spikes = (size_t)max_input_spikes.u.i;


    /* read [dataset] section */
    train_provided = toml_table_int(tbl_dataset, "train_provided");
    train_set = toml_table_string(tbl_dataset, "train_set");
    train_labels = toml_table_string(tbl_dataset, "train_labels");
    n_train = toml_table_int(tbl_dataset, "n_train");

    test_provided = toml_table_int(tbl_dataset, "test_provided");
    test_set = toml_table_string(tbl_dataset, "test_set");
    test_labels = toml_table_string(tbl_dataset, "test_labels");
    n_test = toml_table_int(tbl_dataset, "n_test");

    dataset_name = toml_table_string(tbl_dataset, "dataset_name");
    n_classes = toml_table_int(tbl_dataset, "num_classes");
    epochs = toml_table_int(tbl_dataset, "epochs");
    input_size = toml_table_int(tbl_dataset, "input_size");

    shuffle_samples = toml_table_int(tbl_dataset, "shuffle_samples");

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
    // test dataset is only necessary on ML tasks, not in biological simulations
    if(!test_provided.ok && conf->simulation_obj == 1){

        if(!test_set.ok){
            
            printf(" >> File to load the test dataset must be provided for ML!\n");
            fflush(stdout);
            exit(1);
        }
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
    if(!n_classes.ok && conf->simulation_obj == 1){
        
        printf(" >> In ML simulations, the number of classes in the dataset is necessary!\n");
        fflush(stdout);
        exit(1);
    }
    if(!epochs.ok && conf->simulation_obj == 1){
        
        printf(" >> In ML simulations, the number of epochs to be simulated is necessary!\n");
        fflush(stdout);
        exit(1);
    }
    // if there is no learn it makes no sense to simulate more than one epochs
    if(conf->simulation_obj == 0 && conf->learn == 0){
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


    /* read [outputt] section */
    store_generated_spikes = toml_table_int(tbl_output, "store_generated_spikes");
    generated_spikes_file = toml_table_string(tbl_output, "generated_spikes");
    store_n_spikes = toml_table_int(tbl_output, "store_n_spikes");
    n_spikes_per_neuron_file = toml_table_string(tbl_output, "spikes_per_neuron");
    store_execution_times = toml_table_int(tbl_output, "store_execution_times");
    execution_times_file = toml_table_string(tbl_output, "execution_times");
    store_network = toml_table_int(tbl_output, "store_network"); // whether to store or not the file
    store_network_file = toml_table_string(tbl_output, "store_network_file"); // file to store the network

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


    /* read [network] section */ // TODO: not provide network, but properties and generated automatically??
    network_file = toml_table_string(tbl_network, "network_file");
    network_neurons_file = toml_table_string(tbl_network, "network_neurons_file");
    network_synapses_file = toml_table_string(tbl_network, "network_synapses_file");
    behaviours = toml_table_int(tbl_network, "behaviours");
    delays = toml_table_int(tbl_network, "delays");
    weights = toml_table_int(tbl_network, "weights");
    training_zones = toml_table_int(tbl_network, "training_zones");
    thresh = toml_table_int(tbl_network, "thresh");
    v_rest = toml_table_int(tbl_network, "v_rests");
    t_refract = toml_table_int(tbl_network, "t_refract");
    R = toml_table_int(tbl_network, "R");
    // TODO: more parameters should be added in the future, such as STDP parameters...

    // check if all is in configuration file
    if(!network_file.ok)
    {
        printf(" > The file to load the network from must be provided!\n");
        fflush(stdout);
        exit(1);
    }

    // if there is not in the file, by default set as not provided
    if(!behaviours.ok)
        behaviours.u.i = 0;
    if(!delays.ok)
        delays.u.i = 0;
    if(!weights.ok)
        weights.u.i = 0;
    if(!training_zones.ok)
        training_zones.u.i = 0;
    if(!thresh.ok)
        thresh.u.i = 0;
    if(!v_rest.ok)
        v_rest.u.i = 0;
    if(!R.ok)
        R.u.i = 0;
    if(!t_refract.ok)
        t_refract.u.i = 0;



    conf->store_generated_spikes = store_generated_spikes.u.i;
    conf->generated_spikes_file = generated_spikes_file.u.s;

    conf->store_execution_times = store_execution_times.u.i;
    conf->execution_times_file = execution_times_file.u.s;

    conf->store_n_spikes = store_n_spikes.u.i;
    conf->n_spikes_per_neuron_file = n_spikes_per_neuron_file.u.s;

    conf->store_network = store_network.u.i;
    if(conf->store_network == 1)
        conf->store_network_file = store_network_file.u.s;

    conf->network_file = network_file.u.s;
    if(network_neurons_file.ok)
        conf->network_neurons_file = network_neurons_file.u.s;
    if(network_synapses_file.ok)
        conf->network_synapses_file = network_synapses_file.u.s;
    
    conf->behaviours_provided = behaviours.u.i;
    conf->thresholds_provided = thresh.u.i;
    conf->v_rests_provided = v_rest.u.i;
    conf->refract_times_provided = t_refract.u.i;
    conf->R_provided = R.u.i;

    conf->delays_provided = delays.u.i;
    conf->weights_provided = weights.u.i;
    conf->training_zones_provided = training_zones.u.i;

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

    toml_value_t n_neurons, n_input_neurons, n_output_neurons, n_synapses, n_input_synapses, n_output_synapses,
                    behaviour, v_thres, v_rest, t_refract, res, latency, training_zone, weight,
                    n_connections, neuron_id, n_synapses_to_neuron, network_is_separated;


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
    n_input_synapses = toml_table_int(tbl_general, "input_synapsis");
    n_output_synapses = toml_table_int(tbl_general, "output_synapsis");
    network_is_separated = toml_table_int(tbl_general, "network_is_separated"); // it indicates if the network is separated in several files

    size_t N, iN, oN, S, iS, oS;
    N = (size_t)n_neurons.u.i;
    iN = (size_t)n_input_neurons.u.i;
    oN = (size_t)n_output_neurons.u.i;
    S = (size_t)n_synapses.u.i;
    iS = (size_t)n_input_synapses.u.i;
    oS = (size_t)n_output_synapses.u.i;
    
    // load information in lists structure
    lists->n_neurons = N;
    lists->n_input_neurons = iN;
    lists->n_output_neurons = oN;
    lists->n_synapses = S;
    lists->n_input_synapses = iS;
    lists->n_output_synapses = oS;

    // if the network is separated in several files, load the rest
    if(network_is_separated.ok && network_is_separated.u.i == 1){

        // read files
        open_file(&f_neurons, conf->network_neurons_file); // TOML file
        open_file(&f_synapses, conf->network_synapses_file); // TOML file
    }

    // check that all the information has been loaded correctly
    if(!(n_neurons.ok && n_input_neurons.ok && n_output_neurons.ok && n_synapses.ok && n_input_synapses.ok && n_output_synapses.ok)){
        printf("The number of neurons, input neurons, output neurons, synapses, input synapses and output synapses must be provided in the network file!");
        fflush(stdout);
        exit(1);
    }

    // allocate memory for arrays related to neurons data
    lists->neuron_excitatory = (int *)malloc(N * sizeof(int)); // I think this can be removed
    lists->r_time_list = (int *)malloc(N * sizeof(int)); // refractary times
    lists->v_thres_list = (double *)malloc(N * sizeof(double)); // thresholds
    lists->v_rest_list = (double *)malloc(N * sizeof(double)); // resting potentials
    lists->R_list = (double *)malloc(N * sizeof(double)); // resistances

    // allocate memory for arrays related to synapses data
    lists->weight_list = (double *)malloc(S * sizeof(double)); // weights
    lists->delay_list = (int *)malloc(S * sizeof(int)); // delays
    lists->training_zones = (int *)malloc(S * sizeof(int)); // learnin rules

    // allocate memory for connectibity
    lists->synaptic_connections = (int **)malloc(N * sizeof(int *)); // +1, since the input layer is stored too


    /* load Neurons section */
    behaviour = toml_table_int(tbl_neurons, "behaviour");
    v_thres = toml_table_double(tbl_neurons, "v_thres");
    v_rest = toml_table_double(tbl_neurons, "v_rest");
    t_refract = toml_table_int(tbl_neurons, "t_refract");
    res = toml_table_double(tbl_neurons, "resistance");
    printf(" >> Neurons section loaded\n");
    fflush(stdout);


    // if network information is not separated into more than one files
    if(!network_is_separated.ok || (network_is_separated.ok && network_is_separated.u.i != 1)){
        
        // read information from TOML file
        behaviour_lst = toml_table_array(tbl_neurons, "behaviour_list");
        v_thres_lst = toml_table_array(tbl_neurons, "v_thres_list");
        v_rest_lst = toml_table_array(tbl_neurons, "v_rest_list");
        t_refract_lst = toml_table_array(tbl_neurons, "t_refract_list");
        //res_lst = toml_table_array(tbl_neurons, "res_list");

        // load information into snn structure
        for(i=0; i<N; i++){
            // load data into value variables
            behaviour = toml_array_int(behaviour_lst, i);
            t_refract = toml_array_int(t_refract_lst, i);
            v_thres = toml_array_double(v_thres_lst, i);
            v_rest = toml_array_double(v_rest_lst, i);
            res = toml_array_double(res_lst, i);

            // analyze if this parameters are not provided when it is supposed that they are provided
            // if data was provided load it directly on the network
            if((!behaviour.ok && conf->behaviours_provided == 1) || conf->behaviours_provided == 0){
                printf("Following configuration file, behaviours for neurons must be provided, setting 1 (excitatory)\n");
                behaviour.u.i = 1;
            }
            if((!v_thres.ok && conf->thresholds_provided == 1) || conf->thresholds_provided == 0){
                printf("Following configuration file, thresholds for neurons must be provided, setting 150\n");
                v_thres.u.d = 150;
            }
            if((!v_rest.ok && conf->v_rests_provided == 1) || conf->v_rests_provided == 0){
                printf("Following configuration file, resting potentials for neurons must be provided, setting 50\n");
                v_rest.u.d = 50;
            }
            //if((!res.ok && conf->R_provided == 1) || conf->R_provided == 0){
            //    printf("Following configuration file, resistances for neurons not proveided, setting 1\n");
            //    res.u.d = 1;
            //}
            if((!t_refract.ok && conf->refract_times_provided == 1) || conf->refract_times_provided == 0){
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
        for(i=0; i<N; i++){
            fscanf(f_neurons, "%d", &((lists->neuron_excitatory)[i]));
        }
        for(i=0; i<N; i++){
            fscanf(f_neurons, "%lf", &((lists->v_thres_list)[i]));
        }
        for(i=0; i<N; i++){
            fscanf(f_neurons, "%lf", &((lists->v_rest_list)[i]));
        }
        for(i=0; i<N; i++){
            fscanf(f_neurons, "%d", &((lists->r_time_list)[i]));
        }

        // TODO: not correct
        for(i = 0; i<N; i++)
            lists->R_list[i] = 1;

        fclose(f_neurons);
    }
    

    /* Synapses section */

    lists->max_delay = 0;

    // if network information is not separated into more than one files
    if(!network_is_separated.ok || (network_is_separated.ok && network_is_separated.u.i != 1)){

               
        latency_lst = toml_table_array(tbl_synapses, "latency_list"); // load latencies
        weights_lst = toml_table_array(tbl_synapses, "weights"); // load weights
        training_zones_lst = toml_table_array(tbl_synapses, "training_zones_list"); // load training zones
        //connection_lst_lst = toml_table_array(tbl_synapses, "connections"); // load connections
        

        // load information into snn structure
        for(i=0; i<S; i++){
            // load data into value variables
            latency = toml_array_int(latency_lst, i);
            weight = toml_array_double(weights_lst, i);
            training_zone = toml_array_int(training_zones_lst, i);

            // analyze if this parameters are not provided when it is supposed that they are provided
            // if data was provided load it directly on the network
            if((!latency.ok && conf->delays_provided == 1) || conf->delays_provided == 0){
                printf("Following configuration file, latencies for synapses must be provided, setting 1\n");
                latency.u.i = 1;
            }
            if((!weight.ok && conf->weights_provided == 1) || conf->weights_provided == 0){
                printf("Following configuration file, weights for synapses must be provided, setting 100\n");
                weight.u.d = 100;
            }
            if((!training_zone.ok && conf->training_zones_provided == 1) || conf->training_zones_provided == 0){
                printf("Following configuration file, training zones for neurons must be provided, setting 0 (normal STDP)\n");
                training_zone.u.i = 0;
            }

            // load data into lists structure
            lists->weight_list[i] = weight.u.d;
            lists->delay_list[i] = latency.u.i;
            lists->training_zones[i] = training_zone.u.i;
        }
    }
    // if it is separated, read the information from the other file
    else{

        for(i=0; i<S; i++){
            fscanf(f_synapses, "%d", &(lists->delay_list[i]));

            // store maximum delay value
            if(lists->delay_list[i] > lists->max_delay)
                lists->max_delay = lists->delay_list[i];
        }
        lists->max_delay += 1; // store max_delay

        for(i=0; i<S; i++){
            fscanf(f_synapses, "%lf", &(lists->weight_list[i]));
        }
        for(i=0; i<S; i++){
            fscanf(f_synapses, "%d", &(lists->training_zones[i]));
        }
    }
    
    printf(" >> Synapses loaded from file!\n");
    fflush(stdout);


    // load connectivity information
    if(!network_is_separated.ok || (network_is_separated.ok && network_is_separated.u.i != 1)){
        
        // load connections (list of lists)
        for(i=0; i<N+1; i++){

            connection_lst = toml_array_array(connection_lst_lst, i);
            n_connections = toml_array_int(connection_lst, 0);
            
            // check that data has been correctly loaded
            if(!n_connections.ok){
                printf(" > Connection list is incorrect. Exiting.\n");
                fflush(stdout);
                exit(1);
            }

            //
            (lists->synaptic_connections)[i] = malloc((n_connections.u.i * 2 + 1) * sizeof(int)); // for each connection the neuron id and the number of synapses must be stored
            (lists->synaptic_connections)[i][0] = n_connections.u.i;
        
            for(j = 0; j<(size_t)n_connections.u.i; j++){
                neuron_id = toml_array_int(connection_lst, j * 2 + 1);
                n_synapses_to_neuron = toml_array_int(connection_lst, j * 2 + 2);

                // check that the information have been correctly loaded
                if(!(neuron_id.ok && n_synapses_to_neuron.ok)){
                    printf("Connection list data is incorrect. Exiting\n");
                    fflush(stdout);
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
    }
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
                dataset->freq[i * dataset->n_features + j] = conf->max_spikes / n_spikes; // spikes each freq time steps
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


// [DEPRECATED]
// Function to load a network into the simulation from a network definition file
void load_network_information(const char *file_name, spiking_nn_t *snn, network_construction_lists_t *lists, simulation_configuration_t *conf) {
    
    FILE *f = NULL, *f_neurons, *f_synapses;
    char errbuf[100];
    int i, j;

    // define table and parameters variables
    toml_table_t *tbl, *tbl_general, *tbl_neurons, *tbl_synapses;

    toml_array_t *behaviour_lst, *v_thres_lst, *v_rest_lst, *t_refract_lst, *res_lst, *latency_lst, 
                    *weights_lst, *training_zones_lst, *connection_lst_lst, *connection_lst;

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

    // if the network is separated in several files, load the rest
    if(network_is_separated.ok && network_is_separated.u.i == 1){

        // read files
        open_file(&f_neurons, conf->network_neurons_file); // TOML file
        open_file(&f_synapses, conf->network_synapses_file); // TOML file
    }


    // check that all the information has been loaded correctly
    if(!(n_neurons.ok && n_input_neurons.ok && n_output_neurons.ok && n_synapses.ok && n_input_synapses.ok && n_output_synapses.ok)){
        printf("The number of neurons, input neurons, output neurons, synapses, input synapses and output synapses must be provided in the network file!");
        fflush(stdout);
        exit(1);
    }

    // if correctly readed, copy the information to the snn structure
    snn->n_neurons = n_neurons.u.i;
    snn->n_input = n_input_neurons.u.i;
    snn->n_output = n_output_neurons.u.i;
    snn->n_synapses = n_synapses.u.i;
    snn->n_input_synapses = n_input_synapses.u.i;
    snn->n_output_synapses = n_output_synapses.u.i;

    lists->n_neurons = n_neurons.u.i;
    lists->n_input_neurons = n_input_neurons.u.i;
    lists->n_output_neurons = n_output_neurons.u.i;
    lists->n_synapses = n_synapses.u.i;
    lists->n_input_synapses = n_input_synapses.u.i;
    lists->n_output_synapses = n_output_synapses.u.i;


    // allocate memory for arrays related to neurons data
    if(conf->behaviours_provided == 1)
        lists->neuron_excitatory = (int *)malloc(snn->n_neurons * sizeof(int)); // I think this can be removed
    
    if(conf->refract_times_provided == 1)
        lists->r_time_list = (int *)malloc(snn->n_neurons * sizeof(int)); // refractary times
    
    if(conf->thresholds_provided == 1)
        lists->v_thres_list = (double *)malloc(snn->n_neurons * sizeof(double)); // thresholds
    
    if(conf->v_rests_provided == 1)
        lists->v_rest_list = (double *)malloc(snn->n_neurons * sizeof(double)); // resting potentials
    
    if(conf->R_provided == 1)
        lists->R_list = (double *)malloc(snn->n_neurons * sizeof(double)); // resistances

    // allocate memory for arrays related to synapses data
    lists->weight_list = (double *)malloc(snn->n_synapses * sizeof(double)); // weights
    lists->delay_list = (int *)malloc(snn->n_synapses * sizeof(int)); // delays
    lists->training_zones = (int *)malloc(snn->n_synapses * sizeof(int)); // learnin rules

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
    if(!network_is_separated.ok || (network_is_separated.ok && network_is_separated.u.i != 1)){
        
        // read information from TOML file
        behaviour_lst = toml_table_array(tbl_neurons, "behaviour_list");
        v_thres_lst = toml_table_array(tbl_neurons, "v_thres_list");
        v_rest_lst = toml_table_array(tbl_neurons, "v_rest_list");
        t_refract_lst = toml_table_array(tbl_neurons, "t_refract_list");
        //res_lst = toml_table_array(tbl_neurons, "res_list");

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
            if((!behaviour.ok && conf->behaviours_provided == 1) || conf->behaviours_provided == 0){
                printf("Following configuration file, behaviours for neurons must be provided, setting 1 (excitatory)\n");
                behaviour.u.i = 1;
            }
            if((!v_thres.ok && conf->thresholds_provided == 1) || conf->thresholds_provided == 0){
                printf("Following configuration file, thresholds for neurons must be provided, setting 150\n");
                v_thres.u.d = 150;
            }
            if((!v_rest.ok && conf->v_rests_provided == 1) || conf->v_rests_provided == 0){
                printf("Following configuration file, resting potentials for neurons must be provided, setting 50\n");
                v_rest.u.d = 50;
            }
            //if((!res.ok && conf->R_provided == 1) || conf->R_provided == 0){
            //    printf("Following configuration file, resistances for neurons not proveided, setting 1\n");
            //    res.u.d = 1;
            //}
            if((!t_refract.ok && conf->refract_times_provided == 1) || conf->refract_times_provided == 0){
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

    struct timespec start, end; // to measure simulation complete time
    clock_gettime(CLOCK_MONOTONIC, &start);

    // if network information is not separated into more than one files
    if(!network_is_separated.ok || (network_is_separated.ok && network_is_separated.u.i != 1)){

               
        latency_lst = toml_table_array(tbl_synapses, "latency_list"); // load latencies
        weights_lst = toml_table_array(tbl_synapses, "weights"); // load weights
        training_zones_lst = toml_table_array(tbl_synapses, "training_zones_list"); // load training zones
        //connection_lst_lst = toml_table_array(tbl_synapses, "connections"); // load connections
        

        // load information into snn structure
        for(i=0; i<snn->n_synapses; i++){
            // load data into value variables
            latency = toml_array_int(latency_lst, i);
            weight = toml_array_double(weights_lst, i);
            training_zone = toml_array_int(training_zones_lst, i);

            // analyze if this parameters are not provided when it is supposed that they are provided
            // if data was provided load it directly on the network
            if((!latency.ok && conf->delays_provided == 1) || conf->delays_provided == 0){
                printf("Following configuration file, latencies for synapses must be provided, setting 1\n");
                latency.u.i = 1;
            }
            if((!weight.ok && conf->weights_provided == 1) || conf->weights_provided == 0){
                printf("Following configuration file, weights for synapses must be provided, setting 100\n");
                weight.u.d = 100;
            }
            if((!training_zone.ok && conf->training_zones_provided == 1) || conf->training_zones_provided == 0){
                printf("Following configuration file, training zones for neurons must be provided, setting 0 (normal STDP)\n");
                training_zone.u.i = 0;
            }

            // load data into lists structure
            lists->weight_list[i] = weight.u.d;
            lists->delay_list[i] = latency.u.i;
            lists->training_zones[i] = training_zone.u.i;
        }
    }
    // if it is separated, read the information from the other file
    else{

        for(i=0; i<snn->n_synapses; i++){
            fscanf(f_synapses, "%d", &(lists->delay_list[i]));
        }
        for(i=0; i<snn->n_synapses; i++){
            fscanf(f_synapses, "%lf", &(lists->weight_list[i]));
        }
        for(i=0; i<snn->n_synapses; i++){
            fscanf(f_synapses, "%d", &(lists->training_zones[i]));
        }
    }
    clock_gettime(CLOCK_MONOTONIC, &end);

    double elpt = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;
    printf(" Elapsed time reading synapses data: %lf\n", elpt);
    
    printf(" >> Synapses loaded from file!\n");
    fflush(stdout);


    // load connectivity information
    if(!network_is_separated.ok || (network_is_separated.ok && network_is_separated.u.i != 1)){
        
        // load connections (list of lists)
        for(i=0; i<snn->n_neurons+1; i++){

            connection_lst = toml_array_array(connection_lst_lst, i);
            n_connections = toml_array_int(connection_lst, 0);
            
            // check that data has been correctly loaded
            if(!n_connections.ok){
                printf(" > Connection list is incorrect. Exiting.\n");
                fflush(stdout);
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
                    fflush(stdout);
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


// [DEPRECATED]
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


// [DEPRECATED]
void load_dataset_from_file(input_data_t *dataset, const char *file_name, const char *labels_file_name, int n_samples, simulation_configuration_t *conf){

    int i, j, l, n_spikes;
    
    // copy dataset information from conf to dataset
    dataset->image_size = conf->input_size;
    dataset->n_classes = conf->n_classes;
    dataset->n_samples = n_samples;

    // load dataset from file
    FILE *f = NULL, *f_labels = NULL;
    open_file(&f, file_name);

    // if ML, load labels
    if(conf->simulation_obj == 1)
        open_file(&f_labels, labels_file_name);

    // allocate memory for samples and labels
    dataset->samples = (sample_t *)malloc(dataset->n_samples * sizeof(sample_t));

    // allocate memory for labels
    if(conf->simulation_obj == 1)
        dataset->labels = (int *)malloc(dataset->n_samples * sizeof(int));

    // load samples
    for(i = 0; i<dataset->n_samples; i++){

        // load label
        if(conf->simulation_obj == 1)
            fscanf(f_labels, "%d", &(dataset->labels[i]));
        
        // load each element of the sample
        dataset->samples[i].st = (spike_train_t *)malloc(dataset->image_size * sizeof(spike_train_t));
        
        for(j = 0; j<dataset->image_size; j++){

            // read number of spikes
            fscanf(f, "%d", &n_spikes);
            dataset->samples[i].st[j].n_spikes = n_spikes;
            dataset->samples[i].st[j].stimes = (int *)malloc(n_spikes * sizeof(int));

            // load spike train for each element of the sample
            for(l = 0; l<n_spikes; l++){
                fscanf(f, "%d", &(dataset->samples[i].st[j].stimes[l]));
            }
        }
    }
}


// TODO: all the functions below must be revised!!!

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

// Functions to store results and data into files
void store_results(simulation_results_t *results, simulation_configuration_t *conf, spiking_nn_t *snn, input_data_t *dataset, int *batches, int batch){

    // store only if it is the last epoch
    if(conf->store_generated_spikes == 1){
        store_generated_spikes(results, conf, snn, dataset, batches, batch);
    }

    if(conf->store_n_spikes == 1){        
        store_number_of_spikes(results, conf,snn, dataset, batches, batch);
    }

    if(conf->store_execution_times == 1){
        store_times(results, conf, snn, dataset);
    }

    if(conf->store_network == 1){
        //store_network();
    }
}

void store_epoch_results(simulation_results_t *results, simulation_configuration_t *conf, spiking_nn_t *snn, input_data_t *dataset){

}

void store_generated_spikes(simulation_results_t *results, simulation_configuration_t *conf, spiking_nn_t *snn, input_data_t *dataset, int *batches, int batch){
    
    int i, j, l, s;
    simulation_results_per_sample_t *results_per_sample;
    FILE *f = conf->f_gs;

    /*for(i = 0; i<dataset->n_samples; i++){
        
        results_per_sample = &(results->results_per_sample[i]);
        for(j = 0; j<snn->n_neurons; j++){
            
            for(l = 0; l<conf->time_steps; l++){

                fprintf(f, "%c", results_per_sample->generated_spikes[j][l]);
            }
            fprintf(f, "\n");
        }
        fprintf(f, "\n");
    }*/

    for(i = 0; i<conf->batch_size; i++){
        
        s = batches[batch * conf->batch_size + i];

        if(s != -1){
            results_per_sample = &(results->results_per_sample[s]);        
            
            for(j = 0; j<snn->n_neurons; j++){
                
                for(l = 0; l<conf->time_steps; l++){

                    fprintf(f, "%c", results_per_sample->generated_spikes[j][l]);
                }
                fprintf(f, "\n");
            }
        }
        fprintf(f, "\n");
    }
}


void store_network(simulation_results_t *results, simulation_configuration_t *conf, spiking_nn_t *snn){
    
    // TODO: In this moment this function only stores the first sample results
    simulation_results_per_sample_t *results_per_sample = &(results->results_per_sample[0]);
    FILE *f = conf->f_sn;

    // store network

}


void store_number_of_spikes(simulation_results_t *results, simulation_configuration_t *conf, spiking_nn_t *snn, input_data_t *dataset, int *batches, int batch){
    
    int i, j, s;
    simulation_results_per_sample_t *results_per_sample;
    FILE *f = conf->f_ns;

    
    // store all samples number of spikes
    /*for(i = 0; i<dataset->n_samples; i++){

        results_per_sample = &(results->results_per_sample[i]);

        // store number of spikes
        for(j = 0; j<snn->n_neurons; j++){
            fprintf(f, "%d ", results_per_sample->n_spikes_per_neuron[j]);
        }

        fprintf(f, "\n");
    }
    fprintf(f, "\n");*/

    for(i = 0; i<conf->batch_size; i++){

        s = batches[batch * conf->batch_size + i];
        if(s != -1){
            results_per_sample = &(results->results_per_sample[s]);

            // store number of spikes
            for(j = 0; j<snn->n_neurons; j++){
                fprintf(f, "%d ", results_per_sample->n_spikes_per_neuron[j]);
            }

            fprintf(f, "\n");
        }
    }
    fprintf(f, "\n");
}


void store_times(simulation_results_t *results, simulation_configuration_t *conf, spiking_nn_t *snn, input_data_t *dataset){
    
    int i;
    simulation_results_per_sample_t *results_per_sample;
    FILE *f = conf->f_et;


    // times for the epoch
    fprintf(f, "%lf ", results->elapsed_time_epoch); // store total elapsed time
    fprintf(f, "%lf ", results->elapsed_time_neurons_input); // store total elapsed time
    fprintf(f, "%lf ", results->elapsed_time_neurons_output); // store total elapsed time
    fprintf(f, "%lf ", results->elapsed_time_learning); // store total elapsed time
    fprintf(f, "%lf ", results->elapsed_time_re_neurons); // store total elapsed time
    fprintf(f, "%lf ", results->elapsed_time_re_synapses); // store total elapsed time
    fprintf(f, "%lf\n", results->elapsed_time_load_sample); // store total elapsed time
    fprintf(f, "\n");


    /*for(i=0; i<dataset->n_samples; i++){

        results_per_sample = &(results->results_per_sample[i]);
        fprintf(f, "%lf ", results_per_sample->elapsed_time_neurons_input);
        fprintf(f, "%lf ", results_per_sample->elapsed_time_neurons_output);
        fprintf(f, "%lf ", results_per_sample->elapsed_time_learning);
        fprintf(f, "%lf ", results_per_sample->elapsed_time_re_neurons);
        fprintf(f, "%lf ", results_per_sample->elapsed_time_re_synapses);
        fprintf(f, "%lf\n", results_per_sample->elapsed_time_load_sample);
    }
    fprintf(f, "\n");*/
}