#include "snn_library.h"

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>

#include "training_rules/stdp.h"
#include "neuron_models/lif_neuron.h"
#include "encoders/image_encoders.h"

void initialize_network_neurons(spiking_nn_t *snn, network_construction_lists_t *lists){
    
    int i, j;
    int *neurons_input_synapses, *neurons_output_synapses;
    
    // initialize neurons
    if(snn->neuron_type == 0){ // lif neurons
        snn->lif_neurons = (lif_neuron_t *)malloc(snn->n_neurons * sizeof(lif_neuron_t)); // reserve memory for lif neurons
    }
    // else {}


    // alloc memory to store the number of input and output synapses for each neuron
    neurons_input_synapses = calloc(snn->n_neurons, sizeof(int)); 
    neurons_output_synapses = calloc(snn->n_neurons, sizeof(int));
    
    
    // count neuron input synapses
    
    // input layer (stored in the first row of lists)
    for(i=0; i<lists->synaptic_connections[0][0]; i++){
        
        // [i * 2 + 1] is the index of the neuron, [i * 2 + 2] the number of input synapses 
        neurons_input_synapses[lists->synaptic_connections[0][i*2+1]] += lists->synaptic_connections[0][i*2+2];
    }

    // rest of layers
    for(i=1; i<snn->n_neurons+1; i++){
        
        for(j=0; j<lists->synaptic_connections[i][0]; j++){
            
            // check if it is an output synapse
            if(lists->synaptic_connections[i][j*2+1] != -1){
                neurons_input_synapses[lists->synaptic_connections[i][j*2+1]] += lists->synaptic_connections[i][j*2+2];
                neurons_output_synapses[i-1] += lists->synaptic_connections[i][j*2+2];
            }
            // is only an output connection
            else{
                neurons_output_synapses[i-1] += lists->synaptic_connections[i][j*2+2];
            }
        }
    }


    // initialize neurons
    for(i=0; i<snn->n_neurons; i++){

        if(snn->neuron_type == 0){ // lif neurons
            initialize_lif_neuron(snn, i, lists, neurons_input_synapses[i], neurons_output_synapses[i]);
        }
    }

    // free memory
    free(neurons_input_synapses);
    free(neurons_output_synapses);
}


void initialize_synapse(synapse_t *synapse, network_construction_lists_t *lists, spiking_nn_t *snn, int synapse_id){
    
    int i;

    // read synapse parameters from lists
    synapse->w = lists->weight_list[synapse_id];
    synapse->delay = lists->delay_list[synapse_id];

    if(synapse_id < snn->n_input_synapses)
        synapse->delay = 0;

    // TODO: REVISE REVISE REVISE REVISE 
    // TODO: IN GENERAL, THE WAY I MANAGE THE INPUT MUST BE REVISED
    //if(synapse_id < snn->n_input_synapses || synapse_id >= snn->n_synapses - snn->n_output_synapses)
    //    synapse->max_spikes = INPUT_MAX_SPIKES;
    //else    
    //    synapse->max_spikes = MAX_SPIKES;

    // initialize control parameters
    synapse->t_last_post_spike = -1;
    synapse->t_last_pre_spike = -1;
    
    synapse->post_neuron_computed = -1;
    synapse->pre_neuron_computed = -1;
    
    //synapse->last_spike = 0; 
    //synapse->next_spike = 0;

    // initialize array of spikes for the synapse
    //synapse->l_spike_times = (int *)malloc(synapse->max_spikes * sizeof(int));
    //for(i = 0; i<synapse->max_spikes; i++)
    //    synapse->l_spike_times[i] = -1; // no spikes yet


    // set training rule
    switch (lists->training_zones[synapse_id]) // get synapse training zone from list
    {
    case 0:
        synapse->learning_rule = &add_stdp;//(void (*)())&add_stdp;
        break;
    case 1:
        synapse->learning_rule = &mult_stdp;//(void (*)())&mult_stdp;
        break;
    case 2:
        synapse->learning_rule = &anti_stdp;//(void (*)())&anti_stdp;
        break;
    //case 3:
    //    synapse->learning_rule = &triplet_stdp;//(void (*)())&triplet_stdp;
    //    break;*/
    default:
        synapse->learning_rule = &add_stdp;//(void (*)())&add_stdp;
        break;
    }
}



void re_initialize_synapses(spiking_nn_t *snn){
    
    // reinitialize all synapses
    for(int i = 0; i<snn->n_synapses; i++)
        re_initialize_synapse(&(snn->synapses[i]));
}


void re_initialize_synapse(synapse_t *synapse){
    
    // reinitialize synapse control variables
    synapse->t_last_post_spike = -1;
    synapse->t_last_pre_spike = -1;
    
    synapse->post_neuron_computed = -1;
    synapse->pre_neuron_computed = -1;
    
    //synapse->last_spike = 0; 
    //synapse->next_spike = 0;
}


void initialize_network_synapses(spiking_nn_t *snn, int n_synapses, network_construction_lists_t *lists){
    
    int i;

    // allocate memory for synapses
    snn->synapses = (synapse_t *)malloc(n_synapses * sizeof(synapse_t));

    // initialize synapses
    for(i = 0; i<n_synapses; i++)
        initialize_synapse(&(snn->synapses[i]), lists, snn, i);
}


void add_input_synapse_to_neuron(spiking_nn_t *snn, int neuron_index, int synapse_index){
    
    // LIF neuron
    if(snn->neuron_type == 0)
        add_input_synapse_to_lif_neuron(&snn->lif_neurons[neuron_index], &(snn->synapses[synapse_index]), synapse_index);
    // else{}
}


void add_output_synapse_to_neuron(spiking_nn_t *snn, int neuron_index, int synapse_index){
    
    // LIF neuron
    if(snn->neuron_type == 0)
        add_output_synapse_to_lif_neuron(&snn->lif_neurons[neuron_index], &(snn->synapses[synapse_index]), synapse_index);
    // else{}
}


void connect_neurons_and_synapses(spiking_nn_t *snn, int **synaptic_connections){
    
    int n_neurons, n_input, n_output, i, j, l, i_synapse;
    
    // copy data for legibility
    n_neurons = snn->n_neurons;
    n_input = snn->n_input;
    n_output = snn->n_output;
    i_synapse = 0; // index of the actual synapse


    // add network input synapses (the firt neurons and synapses)
    for(i = 0; i<(synaptic_connections[0][0]); i++){ 
        
        // 
        for(j=0; j<(synaptic_connections[0][i*2+2]); j++){
            add_input_synapse_to_neuron(snn, synaptic_connections[0][i*2+1], i_synapse);
            i_synapse++;
        }
    }
    
    // connect neurons with output and input synapses
    for(i = 0; i<n_neurons; i++){ 
        for(j = 0; j<(synaptic_connections[i+1][0]); j++){
            if(synaptic_connections[i+1][j*2+1] != -1){ // network output synapse
                for(l=0; l<synaptic_connections[i+1][j*2+2]; l++){
                    add_input_synapse_to_neuron(snn, synaptic_connections[i+1][j*2+1], i_synapse);
                    add_output_synapse_to_neuron(snn, i, i_synapse); // neuron identifier is i, but in file neuron i is located on i+1 row
                    i_synapse++;
                }
            }
        }
    }

    // connect output neurons to output synapses
    for(i = 0; i<n_neurons; i++){
        // loop over neurons lists
        for(j = 0; j<(synaptic_connections[i+1][0]); j++){ // i+1 since first row is input layer
            // check if it is an output connection
            if(synaptic_connections[i+1][j*2+1] == -1){ 
                for(l=0; l<synaptic_connections[i+1][j*2+2]; l++){
                    add_output_synapse_to_neuron(snn, i, i_synapse);   
                    i_synapse++;
                }
            }    
        }
    }
}


/* Function to initialize pointers to functions that are used during the simulation */
void initialize_network_function_pointers(spiking_nn_t *snn){
    snn->neuron_type = 0;
    switch (snn->neuron_type)
    { 
        case 0: // LIF neuron
            snn->complete_step = &lif_neuron_step;
            snn->neuron_initializer = &initialize_lif_neuron;
            snn->neuron_re_initializer = &re_initialize_lif_neuron;
            snn->input_step = &lif_neuron_compute_input_synapses;
            snn->output_step = &lif_neuron_compute_output_synapses;

            // allocate memory for neurons
            snn->lif_neurons = (lif_neuron_t *)calloc(snn->n_neurons, sizeof(lif_neuron_t));

            break;
        default:
            snn->complete_step = &lif_neuron_step;
            snn->neuron_initializer = &initialize_lif_neuron;
            snn->neuron_re_initializer = &re_initialize_lif_neuron;
            snn->input_step = &lif_neuron_compute_input_synapses;
            snn->output_step = &lif_neuron_compute_output_synapses;

            // allocate memory for neurons
            snn->lif_neurons = (lif_neuron_t *)calloc(snn->n_neurons, sizeof(lif_neuron_t));
            break;
    }
}


void initialize_network(spiking_nn_t *snn, simulation_configuration_t *conf, network_construction_lists_t *lists){ //void *neuron_initializer()){
    
    int i, j;

    // initialize general information
    snn->neuron_type = conf->neuron_type;

    // initialize step function depending on neuron type
    switch (snn->neuron_type)
    {
        case 0: // LIF - Leaky-Integrate-and-Fire
            snn->complete_step = &lif_neuron_step;
            snn->neuron_initializer = &initialize_lif_neuron;
            snn->neuron_re_initializer = &re_initialize_lif_neuron;
            snn->input_step = &lif_neuron_compute_input_synapses;
	        snn->output_step = &lif_neuron_compute_output_synapses;
	    break;
        case 1: // HH
            break;
        default:
            break;
    }

    // initialize neurons of the network
    initialize_network_neurons(snn, lists);
    printf(" >> Neurons initialized!\n");
    fflush(stdout);

    // initialize synapses
    initialize_network_synapses(snn, snn->n_synapses, lists);
    printf(" >> Synapses initialized!\n");
    fflush(stdout);

    // connect neurons and synapses
    connect_neurons_and_synapses(snn, lists->synaptic_connections);
    printf(" >> Connected!\n");
    fflush(stdout);



    // free and reallocate all spikes arrays for neurons
    for(i = 0; i<snn->n_neurons; i++){
    
        int highest_delay = 0;

        // find for each neuron the highest postsynaptic delay value
        for(j = 0; j<snn->lif_neurons[i].n_output_synapse; j++){

            int d = snn->synapses[snn->lif_neurons[i].output_synapse_indexes[j]].delay;
            if(d > highest_delay)
                highest_delay = d;
        }
        snn->lif_neurons[i].highest_post_delay = highest_delay;

        
        int L = highest_delay + 1;
        L = L*10;
        // reallocate spike array
        //if(snn->lif_neurons[i].is_input_neuron != 1){
            
        free(snn->lif_neurons[i].spike_times_arr);
        snn->lif_neurons[i].spike_times_arr = (int *)malloc(L * sizeof(int));
        snn->lif_neurons[i].max_spikes = L;

        for(j = 0; j<L; j++)
            snn->lif_neurons[i].spike_times_arr[j] = -1;
        //}
    }

    // initialize input spikes
    snn->input_lif_neurons = (lif_neuron_t *)malloc(snn->n_input * sizeof(lif_neuron_t));
    for(i = 0; i<snn->n_input; i++){
        
        lif_neuron_t *neuron = &(snn->input_lif_neurons[i]);

        neuron->spike_times_arr = (int *)calloc(INPUT_MAX_SPIKES, sizeof(int));
        neuron->max_spikes = INPUT_MAX_SPIKES;

        // connect input layer neurons with the rest 
        neuron->n_output_synapse = 1;
        neuron->output_synapse_indexes = (int *)malloc(neuron->n_output_synapse * sizeof(int));
        neuron->output_synapse_indexes[0] = i;
        snn->synapses[i].pre_synaptic_lif_neuron = neuron;

        neuron->r_time = 0;
        neuron->r_time_rest = 0;
        neuron->last_spike = 0;
        neuron->t_last_spike = 0;
    }
}


// This function reorders the list of synapses following the input criterion
void reorder_synapse_list(spiking_nn_t *snn){
    
    int *new_order_indexes, *map_indexes;
    int i, j, l, next_pos = 0, n_synap;
    lif_neuron_t *lif_neuron, *pre_neuron;
    synapse_t *new_list, *synapse;    

    // allocate memory
    new_order_indexes = (int *)malloc(snn->n_synapses * sizeof(int));
    map_indexes = (int *)malloc(snn->n_synapses * sizeof(int));
    new_list = (synapse_t *)malloc(snn->n_synapses * sizeof(synapse_t));
    

    // input synapses are not moved
    for(i=0; i<snn->n_input_synapses; i++){
        new_order_indexes[next_pos] = i;
        map_indexes[i] = next_pos; // store mapping
        next_pos ++;
    }

    // order by NEURONS input synapses
    for(i=0; i<snn->n_neurons; i++){
        //printf("Computing neuron %d\n", i);
        if(snn->neuron_type == 0){ // LIF neurons
            //printf("Neuron %d\n", i);
            lif_neuron = &(snn->lif_neurons[i]);
            n_synap = snn->lif_neurons[i].n_input_synapse;
            for(j=0; j<n_synap; j++){
                if(lif_neuron->input_synapse_indexes[j] >= snn->n_input_synapses){
                    new_order_indexes[next_pos] = lif_neuron->input_synapse_indexes[j]; // store synapse index in the list
                    map_indexes[lif_neuron->input_synapse_indexes[j]] = next_pos; // store mapping
                    next_pos ++;
                }
            }
        }
    }


    // add NETWORK output synapses
    for(i = snn->n_synapses - snn->n_output_synapses; i< snn->n_synapses; i++){
        new_order_indexes[next_pos] = i;
        map_indexes[i] = next_pos; // store mapping
        next_pos ++;
    }
    
    
    // copy the new order
    for(i = 0; i<snn->n_synapses; i++){
        new_list[i] = snn->synapses[new_order_indexes[i]];
    }

    // copy again but in snn structure
    for(i = 0; i<snn->n_synapses; i++){
        snn->synapses[i] = new_list[i];
    }    

    for(i=0; i<snn->n_neurons; i++){
        lif_neuron = &(snn->lif_neurons[i]);
        n_synap = snn->lif_neurons[i].n_input_synapse;

        for(j=0; j<n_synap; j++){
            lif_neuron->input_synapse_indexes[j] = map_indexes[lif_neuron->input_synapse_indexes[j]];
        }

        n_synap = snn->lif_neurons[i].n_output_synapse;
        for(j=0; j<n_synap; j++){
            lif_neuron->output_synapse_indexes[j] = map_indexes[lif_neuron->output_synapse_indexes[j]];
        }
    }

#ifdef DEBUG
    printf(" > Index lists by neurons\n");
    for(i=0; i<snn->n_neurons;i++){
        printf(" >> Printing neuron %d lists\n", i);
        lif_neuron = &(snn->lif_neurons[i]);

        printf("        -");
        for(j=0; j<lif_neuron->n_input_synapse; j++){
            printf("%d ,", lif_neuron->input_synapse_indexes[j]);
        }
        printf("\n");

        printf("        -");
        for(j=0; j<lif_neuron->n_output_synapse; j++){
            printf("%d ,", lif_neuron->output_synapse_indexes[j]);
        }
        printf("\n");
    }
#endif

    free(new_list);
    free(new_order_indexes);
}


void initialize_results_struct(simulation_results_t *results, results_configuration_t *conf){
    
    // initialize struct to store simulation configuration and outputs
    results->results_per_sample = (simulation_results_per_sample_t *)malloc(conf->n_samples * sizeof(simulation_results_per_sample_t));
    
    for(int i = 0; i<conf->n_samples; i++)
        initialize_sample_results_struct(&(results->results_per_sample[i]), conf);
}

void initialize_sample_results_struct(simulation_results_per_sample_t *results_per_sample, results_configuration_t *conf){
    
    int i;

    // initialize results struct
    results_per_sample->elapsed_time = 0;
    results_per_sample->elapsed_time_neurons = 0;
    results_per_sample->elapsed_time_neurons_input = 0;
    results_per_sample->elapsed_time_neurons_output = 0;  
    results_per_sample->elapsed_time_synapses = 0;  
    results_per_sample->elapsed_time_synapses_input = 0;  
    results_per_sample->elapsed_time_synapses_output = 0;  
    results_per_sample->elapsed_time_learning = 0;  
    
    results_per_sample->generated_spikes = (unsigned char **)calloc(conf->n_neurons, sizeof(unsigned char *));
    results_per_sample->n_spikes_per_neuron = (int *)calloc(conf->n_neurons, sizeof(int));

    // initialize generated spikes
    for (i = 0; i<conf->n_neurons; i++){
        results_per_sample->generated_spikes[i] = (unsigned char *)calloc((conf->time_steps), sizeof(unsigned char));
        results_per_sample->n_spikes_per_neuron[i] = 0;
    }
}

void reinitialize_results_struct(simulation_results_t *results, results_configuration_t *conf){
    
    int i;

    // initialize struct to store simulation configuration and outputs    
    for(i = 0; i<conf->n_samples; i++)
        reinitialize_sample_results_struct(&(results->results_per_sample[i]), conf);
}


void reinitialize_sample_results_struct(simulation_results_per_sample_t *results_per_sample, results_configuration_t *conf){
    
    int i, j;
    
    results_per_sample->elapsed_time = 0;
    results_per_sample->elapsed_time_neurons = 0;
    results_per_sample->elapsed_time_neurons_input = 0;
    results_per_sample->elapsed_time_neurons_output = 0;  
    results_per_sample->elapsed_time_synapses = 0;  
    results_per_sample->elapsed_time_synapses_input = 0;  
    results_per_sample->elapsed_time_synapses_output = 0;  
    results_per_sample->elapsed_time_learning = 0;    
    
    for (i = 0; i<conf->n_neurons; i++){
        for(j = 0; j<conf->time_steps; j++)
            results_per_sample->generated_spikes[i][j] = 0;
        results_per_sample->n_spikes_per_neuron[i] = 0;
    }
}

void free_results_struct_memory(simulation_results_t *results, results_configuration_t *conf){
    for(int i = 0; i<conf->n_samples; i++)
        free_sample_results_struct_memory(&(results->results_per_sample[i]), conf);
    free(results->results_per_sample);
}

void free_sample_results_struct_memory(simulation_results_per_sample_t *results_per_sample, results_configuration_t *conf){
    for (int i = 0; i<conf->n_neurons; i++)
        free(results_per_sample->generated_spikes[i]);

    free(results_per_sample->generated_spikes);
    free(results_per_sample->n_spikes_per_neuron);
}


void free_lists_memory(network_construction_lists_t *lists, spiking_nn_t *snn){
    int i;

    for(i=0; i<(snn->n_neurons+1); i++){
        free(lists->synaptic_connections[i]);
    }
    free(lists->synaptic_connections);
    free(lists->weight_list);
    free(lists->delay_list);
    free(lists->training_zones);
}


// TODO: this should be refactorized
// Synapse weight should depend on thresholds...
void initialize_network_weights(spiking_nn_t *snn){
    int i; 
    for(i=0; snn->n_synapses; i++){
        snn->synapses[i].w = (double)rand() / RAND_MAX;

        // TODO: this is temporal, should be changed in the near future, using some intelligent weight initialization or a distribution
        if(snn->synapses[i].w < 10)
            snn->synapses[i].w = 10;
        if(snn->synapses[i].w > 200)
            snn->synapses[i].w = 200;

        // inhibitory / excitatory
    }
}


void print_lif_neuron_information(lif_neuron_t *lif_neuron){
    int i;
    
    printf("      > Excitatory (1) / inhibitory (0) = %d; ", lif_neuron->excitatory);
    printf("Neuron membrane potential = %f; ", lif_neuron->v);
    printf("Neuron membrane potential threshold = %f; ", lif_neuron->v_tresh);
    printf("Neuron resting membrane potential = %f; ", lif_neuron->v_rest);
    printf("Neuron resistance = %f; ", lif_neuron->r);
    printf("Neuron refractory time = %d; ", lif_neuron->r_time);
    printf("Neuron is input = %d; ", lif_neuron->is_input_neuron);
    printf("Neuron is output = %d; ", lif_neuron->is_output_neuron);
    printf("Neuron number of input synapses = %d; ", lif_neuron->n_input_synapse);
    printf("Neuron number of output synapses = %d\n", lif_neuron->n_output_synapse);
    printf("      > Neuron input synapse indexes = ");
    for(i=0; i<lif_neuron->n_input_synapse; i++){
        printf("%d ", lif_neuron->input_synapse_indexes[i]);
    }
    printf("\n");

    printf("      > Neuron output synapse indexes = ");
    for(i=0; i<lif_neuron->n_output_synapse; i++){
        printf("%d ", lif_neuron->output_synapse_indexes[i]);
    }
    printf("\n");
    printf("\n");
}

void print_lif_neurons_information(spiking_nn_t *snn){
    int i;
    
    for(i = 0; i<snn->n_neurons; i++){
        print_lif_neuron_information(&(snn->lif_neurons[i]));
    }
}

void print_neurons_information(spiking_nn_t *snn){
    switch(snn->neuron_type)
    {
        case 0:
            print_lif_neurons_information(snn);
            break;
        default:
            print_lif_neurons_information(snn);
            break;
    }
}


void print_input_synapse_information(synapse_t *synapse){
    printf("        > Weight = %f; ", synapse->w);
    printf("Delay = %d\n", synapse->delay);
}

void print_synapse_information(synapse_t *synapse){
    printf("        > Weight = %f; ", synapse->w);
    printf("Delay = %d\n", synapse->delay);
}

void print_output_synapse_information(synapse_t *synapse){
    printf("        > Weight = %f; ", synapse->w);
    printf("Delay = %d\n", synapse->delay);
}

void print_synapses_information(spiking_nn_t *snn){
    int i;

    // print input synapses
    printf("      > Printing input synapses:\n");
    for(i=0; i<snn->n_input_synapses; i++){
        print_input_synapse_information(&(snn->synapses[i]));
    }
    printf("\n");

    // print synapses (not input or output)
    printf("      > Printing middle synapses:\n");
    for(i = snn->n_input_synapses; i<snn->n_synapses - snn->n_output_synapses; i++){
        print_synapse_information(&(snn->synapses[i]));
    }
    printf("\n");

    // print output synapses
    printf("       > Printing output synapses:\n");
    for(i=snn->n_synapses - snn->n_output_synapses; i<snn->n_synapses; i++){
        print_output_synapse_information(&(snn->synapses[i]));
    }
    printf("\n");
}

void print_network_information(spiking_nn_t *snn){
    printf("    > Printing general information:\n");
    printf("      > Type of the neurons = %d\n", snn->neuron_type);
    printf("      > Number of neurons = %d\n", snn->n_neurons);
    printf("      > Number of input neurons = %d\n", snn->n_input);
    printf("      > Number of output neurons = %d\n", snn->n_output);
    printf("      > Number of synapses = %d\n", snn->n_synapses);
    printf("      > Number of input synapses = %d\n", snn->n_input_synapses);
    printf("      > Number of output synapses = %d\n\n", snn->n_output_synapses);

    printf("    > Printing neurons information:\n");
    print_neurons_information(snn);
    printf("    > All neurons printed!\n\n");

    //printf("    > Printing synapses information:\n");
    //print_synapses_information(snn);
    //printf("    > All synapses printed!\n\n");

    
}

// Functions to free resources

void free_snn_struct_memory(spiking_nn_t *snn){
    // free neurons
    //printf(" >>>>>>>> Deallocating neurons memory\n");
    //fflush( stdout );

    if(snn->neuron_type == 0){
        free_lif_neurons(snn);
        //printf(" >>>>>>> Checkpoint\n");
        fflush( stdout );
        free(snn->lif_neurons);
    }

    // free synapses
    //printf(" >>>>>>>> Deallocating synapses memory\n");
    free_synapses(snn);
    free(snn->synapses);

    //printf(" >>>>>>> Deallocating snn pointer\n");
    free(snn);
    //printf(" >>>>>>> Deallocated\n");
    snn = NULL;
}

void free_lif_neurons(spiking_nn_t *snn){
    //printf(" >>>>>>>>> n_neurons = %d\n", snn->n_neurons);
    fflush(stdout);
    for(int i=0; i<snn->n_neurons; i++){
        free(snn->lif_neurons[i].input_synapse_indexes);
        free(snn->lif_neurons[i].output_synapse_indexes);
    }
}

void free_synapses(spiking_nn_t *snn){
    //for(int i = 0; i<snn->n_synapses; i++){
    //    free(snn->synapses[i].l_spike_times);
    //}
    printf(" Not impelemented yet");
}