#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>

#include "snn_library.h"
#include "training_rules/stdp.h"
#include "neuron_models/lif_neuron.h"


void initialize_network(spiking_nn_t *snn, simulation_configuration_t *conf, network_construction_lists_t *data){ 
    
    int i;

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
            snn->load_sample = &load_sample_in_LIF;
	    break;
        case 1: // HH
            break;
        default:
            break;
    }

    // initialize neurons of the network
    initialize_neurons(snn, conf, data);
    printf(" >> Neurons initialized!\n");
    fflush(stdout);

    // initialize synapses
    initialize_synapses(snn, snn->n_synapses, conf, data);
    printf(" >> Synapses initialized!\n");
    fflush(stdout);

    // connect neurons and synapses
    connect_neurons_and_synapses(snn, data->synaptic_connections);
    printf(" >> Connected!\n");
    fflush(stdout);



    // free and reallocate all spikes arrays for neurons (TODO: rethink about L and generalize)
    /*for(i = 0; i<snn->n_neurons; i++){
        
        free(snn->lif_neurons[i].spike_times_arr);
        snn->lif_neurons[i].spike_times_arr = (int *)malloc(100 * sizeof(int));
        snn->lif_neurons[i].max_spikes = 100;

        for(j = 0; j<snn->lif_neurons[i].max_spikes; j++)
            snn->lif_neurons[i].spike_times_arr[j] = -1;
    }*/

    /*for(i = 0; i<snn->n_neurons; i++){
    
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
    }*/


    // initialize input spikes and virtual input neurons // TODO: refactorize
    snn->input_lif_neurons = (lif_neuron_t *)malloc(snn->n_input * sizeof(lif_neuron_t));
    for(i = 0; i<snn->n_input; i++){
        
        lif_neuron_t *neuron = &(snn->input_lif_neurons[i]);

        neuron->max_spikes = conf->max_input_spikes;
        neuron->spike_times_arr = (int *)malloc(conf->max_input_spikes * sizeof(int));
        for(int j = 0; j<neuron->max_spikes; j++){
            neuron->spike_times_arr[j] = -1;
        }

        // connect input layer neurons with the rest 
        neuron->n_input_synapse = 0;

        neuron->n_output_synapse = 1;
        neuron->output_synapse_indexes = (int *)malloc(neuron->n_output_synapse * sizeof(int));
        neuron->output_synapse_indexes[0] = i;
        snn->synapses[i].pre_synaptic_lif_neuron = neuron;
        snn->synapses[i].pre_neuron_index = i; // this is in the array of input neurons

        neuron->r_time = 0;
        neuron->r_time_rest = 0;
        neuron->last_spike = 0;

        // TODO: TEMP: Not necessary in input neurons????
        //neuron->n_last_spikes = 3;
        //neuron->t_last_spikes = (int*)malloc(neuron->n_last_spikes * sizeof(int));
        //for(int j = 0; j<neuron->n_last_spikes; j++){
        //    neuron->t_last_spikes[j] = -1;
        //}
        //neuron->next_last_spike = 0;   
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
            snn->load_sample = &load_sample_in_LIF;

            // allocate memory for neurons
            snn->lif_neurons = (lif_neuron_t *)calloc(snn->n_neurons, sizeof(lif_neuron_t));

            break;
        default:
            snn->complete_step = &lif_neuron_step;
            snn->neuron_initializer = &initialize_lif_neuron;
            snn->neuron_re_initializer = &re_initialize_lif_neuron;
            snn->input_step = &lif_neuron_compute_input_synapses;
            snn->output_step = &lif_neuron_compute_output_synapses;
            snn->load_sample = &load_sample_in_LIF;

            // allocate memory for neurons
            snn->lif_neurons = (lif_neuron_t *)calloc(snn->n_neurons, sizeof(lif_neuron_t));
            break;
    }
}


void initialize_neurons(spiking_nn_t *snn, simulation_configuration_t  *conf, network_construction_lists_t *data){
    
    int i, j;
    int *neurons_input_synapses, *neurons_output_synapses;
    
    // initialize neurons
    if(snn->neuron_type == 0){ // lif neurons
        snn->lif_neurons = (lif_neuron_t *)malloc(snn->n_neurons * sizeof(lif_neuron_t)); // reserve memory for lif neurons
    }
    // else {}

    // alloc memory to store the number of input and output synapses for each neuron
    neurons_input_synapses = (int *)calloc(snn->n_neurons, sizeof(int)); 
    neurons_output_synapses = (int *)calloc(snn->n_neurons, sizeof(int));
    
    
    // count neuron input synapses
    
    // input layer (stored in the first row of lists)
    for(i=0; i<data->synaptic_connections[0][0]; i++){
        
        // [i * 2 + 1] is the index of the neuron, [i * 2 + 2] the number of input synapses 
        neurons_input_synapses[data->synaptic_connections[0][i*2+1]] += data->synaptic_connections[0][i*2+2];
    }

    // rest of layers
    for(i=1; i<snn->n_neurons+1; i++){
        
        for(j=0; j<data->synaptic_connections[i][0]; j++){
            
            // check if it is an output synapse
            if(data->synaptic_connections[i][j*2+1] != -1){
                neurons_input_synapses[data->synaptic_connections[i][j*2+1]] += data->synaptic_connections[i][j*2+2];
                neurons_output_synapses[i-1] += data->synaptic_connections[i][j*2+2];
            }
            // is only an output connection
            else{
                neurons_output_synapses[i-1] += data->synaptic_connections[i][j*2+2];
            }
        }
    }


    // initialize neurons
    for(i=0; i<snn->n_neurons; i++){

        if(snn->neuron_type == 0){ // lif neurons
            initialize_lif_neuron(snn, i, data, neurons_input_synapses[i], neurons_output_synapses[i], conf->max_spikes);
        }
    }

    // free memory
    free(neurons_input_synapses);
    free(neurons_output_synapses);
}


void initialize_synapses(spiking_nn_t *snn, int n_synapses, simulation_configuration_t *conf, network_construction_lists_t *data){
    
    int i;

    // allocate memory for synapses
    snn->synapses = (synapse_t *)calloc(n_synapses, sizeof(synapse_t));

    // initialize synapses
    for(i = 0; i<n_synapses; i++)
        initialize_synapse(&(snn->synapses[i]), data, snn, i);
}

void initialize_synapse(synapse_t *synapse, network_construction_lists_t *data, spiking_nn_t *snn, int synapse_id){
    
    // read synapse parameters from lists
    synapse->w = data->weight_list[synapse_id];
    synapse->dw = 0; // initialize difference in weights
    synapse->init_w = synapse->w;
    synapse->delay = data->delay_list[synapse_id];
    synapse->stdp_steps = 3;
    synapse->next_pre_spike = 0;
    synapse->next_post_spike = 0;

    // input synapses do not have delay
    if(synapse_id < snn->n_input_synapses) 
        synapse->delay = 0;

    // initialize control parameters
    //synapse->t_last_post_spike = -1;
    //synapse->t_last_pre_spike = -1;
        
    // set training rule
    synapse->lr = data->training_zones[synapse_id];
    switch (data->training_zones[synapse_id]) // get synapse training zone from list
    {
        case 0:
            synapse->learning_rule = &addSTDP;//(void (*)())&add_stdp;
            break;
        case 1:
            synapse->learning_rule = &mltSTDP;//(void (*)())&mult_stdp;
            break;
        case 2:
            synapse->learning_rule = &antiSTDP;//(void (*)())&anti_stdp;
            break;
        //case 3:
        //    synapse->learning_rule = &triplet_stdp;//(void (*)())&triplet_stdp;
        //    break;*/
        default:
            synapse->learning_rule = &mltSTDP;//(void (*)())&add_stdp;
            break;
    }
}


// I think there is nothing to reinitialize actually
void re_initialize_synapses(spiking_nn_t *snn){
    
    // reinitialize all synapses
    for(int i = 0; i<snn->n_synapses; i++)
        re_initialize_synapse(&(snn->synapses[i]));
}


void re_initialize_synapse(synapse_t *synapse){
    
    // reinitialize weight
    synapse->w = synapse->init_w;
    synapse->next_pre_spike = 0;
    synapse->next_post_spike = 0;
}

void update_weights(spiking_nn_t *snn){

    int i;

    for(i=0; i<snn->n_synapses; i++)
        update_weight(&(snn->synapses[i]));
}

void update_weight(synapse_t *synapse){

    int sign = synapse->w > 0;
    synapse->w += synapse->dw;

    // I have to check if this penalises a lot the performance
    if((synapse->w > 0) != sign){
        if(synapse->w > 0) synapse->w = -0.01;
        else synapse->w = 0.01; 
    }
}


void connect_neurons_and_synapses(spiking_nn_t *snn, int **synaptic_connections){
    
    int n_neurons, i, j, l, i_synapse;
    
    // copy data for legibility
    n_neurons = snn->n_neurons;
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


void add_input_synapse_to_neuron(spiking_nn_t *snn, int neuron_index, int synapse_index){
    
    // LIF neuron
    if(snn->neuron_type == 0)
        add_input_synapse_to_lif_neuron(snn, neuron_index, synapse_index);
    // else{}
}


void add_output_synapse_to_neuron(spiking_nn_t *snn, int neuron_index, int synapse_index){
    
    // LIF neuron
    if(snn->neuron_type == 0)
        add_output_synapse_to_lif_neuron(snn, neuron_index, synapse_index);
    // else{}
}


void cp_network(spiking_nn_t *cp_snn, spiking_nn_t *or_snn, simulation_configuration_t *conf){

    // cp general information
    cp_snn->neuron_type = or_snn->neuron_type;
    
    cp_snn->n_neurons = or_snn->n_neurons;
    cp_snn->n_input = or_snn->n_input;
    cp_snn->n_output = or_snn->n_output;

    cp_snn->n_synapses = or_snn->n_synapses;
    cp_snn->n_input_synapses = or_snn->n_input_synapses;
    cp_snn->n_output_synapses = or_snn->n_output_synapses;

    // function pointers
    cp_snn->neuron_initializer = or_snn->neuron_initializer;
    cp_snn->neuron_re_initializer = or_snn->neuron_re_initializer;
    cp_snn->complete_step = or_snn->complete_step;
    cp_snn->input_step = or_snn->input_step;
    cp_snn->output_step = or_snn->output_step;
    cp_snn->load_sample = or_snn->load_sample;


    // cp neurons
    cp_neurons(cp_snn, or_snn);

    // cp input neurons
    cp_input_neurons(cp_snn, or_snn);

    // cp synapses
    cp_synapses(cp_snn, or_snn);
}

void cp_neurons(spiking_nn_t *cp_snn, spiking_nn_t *or_snn){

    switch(or_snn->neuron_type){
        case 0:
            // allocate memory
            cp_snn->lif_neurons = (lif_neuron_t *)malloc(or_snn->n_neurons * sizeof(lif_neuron_t));
            // copy
            cp_lif_neurons(cp_snn->lif_neurons, or_snn->lif_neurons, or_snn->n_neurons);
        break;
        default:
            // allocate memory
            cp_snn->lif_neurons = (lif_neuron_t *)malloc(or_snn->n_neurons * sizeof(lif_neuron_t));
            // copy
            cp_lif_neurons(cp_snn->lif_neurons, or_snn->lif_neurons, or_snn->n_neurons);
        break;
    }
    
}

void cp_input_neurons(spiking_nn_t *cp_snn, spiking_nn_t *or_snn){

    switch(or_snn->neuron_type){
        case 0:
            // allocate memory
            cp_snn->input_lif_neurons = (lif_neuron_t *)malloc(or_snn->n_input * sizeof(lif_neuron_t));
            // copy
            cp_lif_neurons(cp_snn->input_lif_neurons, or_snn->input_lif_neurons, or_snn->n_input);
        break;
        default:
            // allocate memory
            cp_snn->input_lif_neurons = (lif_neuron_t *)malloc(or_snn->n_input * sizeof(lif_neuron_t));
            // copy
            cp_lif_neurons(cp_snn->input_lif_neurons, or_snn->input_lif_neurons, or_snn->n_input);        
        break;
    }
    
}

void cp_synapses(spiking_nn_t *cp_snn, spiking_nn_t *or_snn){

    int i;
    synapse_t *cp_synapse, *or_synapse;

    // allocate memory for synapses
    cp_snn->synapses = (synapse_t *)malloc(cp_snn->n_synapses * sizeof(synapse_t));
    
    for(i = 0; i<or_snn->n_synapses; i++){

        or_synapse = &(or_snn->synapses[i]);
        cp_synapse = &(cp_snn->synapses[i]);

        cp_synapse->w = or_synapse->w;
        cp_synapse->init_w = or_synapse->init_w;
        cp_synapse->dw = or_synapse->dw; // should be 0
        cp_synapse->delay = or_synapse->delay;
        cp_synapse->lr = or_synapse->lr;
        cp_synapse->learning_rule = or_synapse->learning_rule; // copy function pointer
        cp_synapse->stdp_steps = or_synapse->stdp_steps;
        cp_synapse->next_pre_spike = 0;
        cp_synapse->next_post_spike = 0;

        cp_synapse->pre_neuron_index = or_synapse->pre_neuron_index;
        cp_synapse->post_neuron_index = or_synapse->post_neuron_index;

        // reference pre and post synaptic neurons // This does not work with input neurons...
        switch(or_snn->neuron_type){
            case 0:
                if(i < or_snn->n_input)
                    cp_synapse->pre_synaptic_lif_neuron = &(cp_snn->input_lif_neurons[cp_synapse->pre_neuron_index]);
                else
                    cp_synapse->pre_synaptic_lif_neuron = &(cp_snn->lif_neurons[cp_synapse->pre_neuron_index]);
                
                if(i < or_snn->n_synapses - or_snn->n_output_synapses)
                    cp_synapse->post_synaptic_lif_neuron = &(cp_snn->lif_neurons[cp_synapse->post_neuron_index]);
            break;
            default:
                if(i < or_snn->n_input)
                    cp_synapse->pre_synaptic_lif_neuron = &(cp_snn->input_lif_neurons[cp_synapse->pre_neuron_index]);
                else
                    cp_synapse->pre_synaptic_lif_neuron = &(cp_snn->lif_neurons[cp_synapse->pre_neuron_index]);
                cp_synapse->post_synaptic_lif_neuron = &(cp_snn->lif_neurons[cp_synapse->post_neuron_index]);
            break;
        }
    }
}


// This function reorders the list of synapses following the input criterion
void reorder_synapse_list(spiking_nn_t *snn){
    
    int *new_order_indexes, *map_indexes;
    int i, j, next_pos = 0, n_synap;
    lif_neuron_t *lif_neuron;
    synapse_t *new_list;    

    // allocate memory
    new_order_indexes = (int *)malloc(snn->n_synapses * sizeof(int));
    map_indexes = (int *)malloc(snn->n_synapses * sizeof(int));
    new_list = (synapse_t *)malloc(snn->n_synapses * sizeof(synapse_t));
    

    // input synapses are not moved
    /*for(i=0; i<snn->n_input_synapses; i++){
        new_order_indexes[next_pos] = i;
        map_indexes[i] = next_pos; // store mapping
        next_pos ++;
    }*/

    // order by NEURONS input synapses
    for(i=0; i<snn->n_neurons; i++){
        //printf("Computing neuron %d\n", i);
        if(snn->neuron_type == 0){ // LIF neurons
            //printf("Neuron %d\n", i);
            lif_neuron = &(snn->lif_neurons[i]);
            n_synap = snn->lif_neurons[i].n_input_synapse;
            for(j=0; j<n_synap; j++){
                //if(lif_neuron->input_synapse_indexes[j] >= snn->n_input_synapses){
                    new_order_indexes[next_pos] = lif_neuron->input_synapse_indexes[j]; // store synapse index in the list
                    map_indexes[lif_neuron->input_synapse_indexes[j]] = next_pos; // store mapping
                    next_pos ++;
                //}
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



// DEPRECATED
/* Mapping functions */ // TODO: Change the L to a configuration parameter
GPU_SNN_t* SNN_CPU2GPU_mapping(spiking_nn_t *snn, simulation_configuration_t *conf){

    int i, j, offset = 0;
    GPU_SNN_t *gpu_snn;

    // allocate memory for GPU SNN structure
    gpu_snn = (GPU_SNN_t*)malloc(sizeof(GPU_SNN_t));

    gpu_snn->n_neurons = snn->n_neurons;
    gpu_snn->n_input_neurons = snn->n_input;
    gpu_snn->n_output_neurons = snn->n_output;
    gpu_snn->n_synapses = snn->n_synapses;
    gpu_snn->n_input_synapses = snn->n_input_synapses;
    gpu_snn->n_output_synapses = snn->n_output_synapses;
    gpu_snn->max_spikes = conf->max_input_spikes;

    // neurons - allocate and cpy
    gpu_snn->v = (float *)malloc((snn->n_neurons) * sizeof(float));
    gpu_snn->v_thresh = (float *)malloc((snn->n_neurons) * sizeof(float));
    gpu_snn->v_rest = (float *)malloc((snn->n_neurons) * sizeof(float));
    gpu_snn->r_period = (int *)malloc((snn->n_neurons) * sizeof(int));
    gpu_snn->r_period_remain = (int *)calloc(snn->n_neurons, sizeof(int)); // 0
    gpu_snn->res = (int *)malloc((snn->n_neurons) * sizeof(int));
    gpu_snn->pre_fired = (int *)calloc((snn->n_synapses), sizeof(int));
    gpu_snn->post_fired = (int *)calloc((snn->n_neurons), sizeof(int));
    gpu_snn->arrI = (float *)malloc(snn->n_neurons * sizeof(float));
    gpu_snn->inR = (int *)calloc(snn->n_neurons, sizeof(int));

    //gpu_snn->n_last_spikes = 3; // TEMPORAL // TODO
    //gpu_snn->t_last_spikes = (int*)calloc(snn->n_neurons * gpu_snn->n_last_spikes, sizeof(int)); // 0
    //gpu_snn->next_last_spike = (int*)calloc(snn->n_neurons, sizeof(int)); // '
    gpu_snn->n_neuron_input_synapses = (size_t*)malloc(snn->n_neurons * sizeof(size_t));
    gpu_snn->neuron_input_synapses_offset = (size_t*)malloc(snn->n_neurons * sizeof(size_t));


    // TODO: generalize??
    for(i = 0; i<snn->n_neurons; i++){
        
        // copy information from CPU
        gpu_snn->v[i] = (float)(snn->lif_neurons[i].v);
        gpu_snn->v_thresh[i] = (float)(snn->lif_neurons[i].v_tresh);
        gpu_snn->v_rest[i] = (float)(snn->lif_neurons[i].v_rest);
        gpu_snn->r_period[i] = snn->lif_neurons[i].r_time;
        gpu_snn->r_period_remain[i] = 0;
        gpu_snn->res[i] = snn->lif_neurons[i].r;

        // initialize t_last_spikes
        //for(j = 0; j<gpu_snn->n_last_spikes; j++)
        //    gpu_snn->t_last_spikes[i * gpu_snn->n_last_spikes + j] = -1; // initialize to -1, no spikes

        // input synapses and offsets
        gpu_snn->n_neuron_input_synapses[i] = snn->lif_neurons[i].n_input_synapse;
        gpu_snn->neuron_input_synapses_offset[i] = offset;

        // update offset for the next iteration
        offset += gpu_snn->n_neuron_input_synapses[i];    
    }

    // synapses - allocate and cpy
    gpu_snn->w = (float*)malloc(snn->n_synapses * sizeof(float));
    gpu_snn->init_w = (float*)malloc(snn->n_synapses * sizeof(float));
    gpu_snn->dw = (float*)calloc(snn->n_synapses, sizeof(float)); // 0
    gpu_snn->delay = (int*)malloc(snn->n_synapses * sizeof(int));
    gpu_snn->lr = (int*)malloc(snn->n_synapses * sizeof(int));
    gpu_snn->pre_neuron_index = (size_t*)malloc(snn->n_synapses * sizeof(size_t));
    gpu_snn->post_neuron_index = (size_t*)malloc(snn->n_synapses * sizeof(size_t));
    //gpu_snn->next_pre_spike = (int*)malloc(snn->n_synapses * sizeof(int));
    //gpu_snn->next_post_spike = (int*)malloc(snn->n_synapses * sizeof(int));
    gpu_snn->pre_trace = (float*)calloc(snn->n_synapses, sizeof(float));
    gpu_snn->post_trace = (float*)calloc(snn->n_synapses, sizeof(float));

    for(i = 0; i<snn->n_synapses; i++){

        gpu_snn->w[i] = (float)(snn->synapses[i].w);
        gpu_snn->init_w[i] = (float)(snn->synapses[i].init_w);
        gpu_snn->dw[i] = 0.0;
        gpu_snn->delay[i] = snn->synapses[i].delay;
        gpu_snn->lr[i] = snn->synapses[i].lr;
        gpu_snn->pre_neuron_index[i] = (size_t)snn->synapses[i].pre_neuron_index;
        gpu_snn->post_neuron_index[i] = (size_t)snn->synapses[i].post_neuron_index;
        //gpu_snn->next_pre_spike[i] = 0;
        //gpu_snn->next_post_spike[i] = 0;

        // if it is an input synapse, then it's pre neuron is a virtual one. Map indexes
        if(gpu_snn->delay[i] > 0) // to detect input synapses
            gpu_snn->pre_neuron_index[i] += snn->n_input; // n_input_synapses??
        gpu_snn->post_neuron_index[i] += snn->n_input;
    }

    // spike matrix
    gpu_snn->spk_matrix = (char*)malloc((snn->n_neurons + snn->n_input) * conf->max_input_spikes * sizeof(char));
    for(i = 0; i<snn->n_neurons + snn->n_input; i++)
        for(j = 0; j<conf->max_input_spikes; j++)
            gpu_snn->spk_matrix[i * conf->max_input_spikes + j] = -1; // no spikes yet


    return gpu_snn;
}


// DEPRECATED
GPU_dataset_t* dataset_CPU2GPU_mapping(input_data_t *dataset, simulation_configuration_t *conf){

    size_t i, j, l;
    size_t next_spike, next_feature;
    size_t n_samples, n_features;

    GPU_dataset_t *gpu_dataset;
    gpu_dataset = (GPU_dataset_t *)malloc(sizeof(GPU_dataset_t));
    
    // copy first to a temporal struct in the CPU
    gpu_dataset->type = dataset->type;
    gpu_dataset->n_classes = dataset->n_classes;
    gpu_dataset->n_samples = dataset->n_samples;
    gpu_dataset->n_features = dataset->image_size; // TODO: generalize
    gpu_dataset->n_spikes_per_feature = (size_t*)malloc(dataset->n_samples * dataset->image_size * sizeof(size_t)); 

    gpu_dataset->sample_offset = (size_t*)malloc(gpu_dataset->n_samples * sizeof(size_t));
    gpu_dataset->feature_offset = (size_t*)malloc(gpu_dataset->n_samples * gpu_dataset->n_features * sizeof(size_t));

    gpu_dataset->n_spikes = 0;


    // loop over samples to count the total number of spikes in the dataset
    for(i = 0; i<dataset->n_samples; i++){
        
        // loop over the features of each sample
        for(j = 0; j<dataset->image_size; j++){ // TODO: Generalize
        
            // count the number of spikes
            gpu_dataset->n_spikes += (size_t)(dataset->samples[i].st[j].n_spikes);
        }
    }

    // allocate memory to store all the spikes in a 1D array
    gpu_dataset->spikes = (size_t*)malloc(gpu_dataset->n_spikes * sizeof(size_t));


    // compute offsets for samples and features
    next_spike = 0;
    next_feature = 0;
    n_samples = gpu_dataset->n_samples;
    n_features = gpu_dataset->n_features;

    for(i = 0; i<n_samples; i++){

        gpu_dataset->sample_offset[i] = next_spike;
        //next_feature = 0; // local for sample

        // loop over features
        for(j = 0; j<n_features; j++){

            gpu_dataset->feature_offset[i * n_features + j] = next_feature;

            // store the number of spikes of the feature
            gpu_dataset->n_spikes_per_feature[i * n_features + j] = dataset->samples[i].st[j].n_spikes;
            //printf(" > Sample %d, feature %d, n. spikes = %d\n", i, j, dataset->samples[i].st[j].n_spikes);

            // loop over spikes
            for(l = 0; l<dataset->samples[i].st[j].n_spikes; l++){

                gpu_dataset->spikes[next_spike] = dataset->samples[i].st[j].stimes[l];
                next_spike ++;
                next_feature ++;
            }
        }
    }

    return gpu_dataset;
}

void free_gpu_snn_in_CPU(GPU_SNN_t *gpu_snn){

    // deallocate internal arrays
    free(gpu_snn->v);
    free(gpu_snn->v_thresh);
    free(gpu_snn->v_rest);
    free(gpu_snn->r_period);
    free(gpu_snn->r_period_remain);
    free(gpu_snn->res);
    //free(gpu_snn->t_last_spikes);
    //free(gpu_snn->next_last_spike);
    free(gpu_snn->n_neuron_input_synapses);
    free(gpu_snn->neuron_input_synapses_offset);
    
    free(gpu_snn->w);
    free(gpu_snn->init_w);
    free(gpu_snn->dw);
    free(gpu_snn->delay);
    free(gpu_snn->lr);
    free(gpu_snn->pre_neuron_index);
    free(gpu_snn->post_neuron_index);
    //free(gpu_snn->next_pre_spike);
    //free(gpu_snn->next_post_spike);

    free(gpu_snn->spk_matrix);

    // deallocate struct
    free(gpu_snn);
}

void free_gpu_dataset_in_CPU(GPU_dataset_t *gpu_dataset){
    
    // deallocate internal arrays
    free(gpu_dataset->spikes);
    free(gpu_dataset->n_spikes_per_feature);
    free(gpu_dataset->sample_offset);
    free(gpu_dataset->feature_offset);

    // deallocate struct
    free(gpu_dataset);
}


void free_gpu_results_in_CPU(GPU_results_t *gpu_results){

    free(gpu_results->n_spks);
    free(gpu_results->gnt_spks);

    free(gpu_results);
}


/////


void initialize_results_struct(simulation_results_t *results, simulation_configuration_t *conf, int n_samples, int n_neurons){
    
    int i;

    if(conf->store_execution_times == 1){
        results->elapsed_time_epoch = 0;
        results->elapsed_time = 0;
    }

    // initialize struct to store simulation configuration and outputs
    results->results_per_sample = (simulation_results_per_sample_t *)malloc(n_samples * sizeof(simulation_results_per_sample_t));
        
    for(i = 0; i<n_samples; i++)
        initialize_sample_results_struct(&(results->results_per_sample[i]), conf, n_samples, n_neurons);
}

void initialize_sample_results_struct(simulation_results_per_sample_t *results_per_sample, simulation_configuration_t *conf, int n_samples, int n_neurons){
    
    int i;

    // initialize simulation time variables
    if(conf->store_execution_times == 1){
        results_per_sample->elapsed_time_neurons = 0;
        results_per_sample->elapsed_time_neurons_input = 0;
        results_per_sample->elapsed_time_neurons_output = 0;  
        results_per_sample->elapsed_time_synapses = 0;  
        results_per_sample->elapsed_time_synapses_input = 0;  
        results_per_sample->elapsed_time_synapses_output = 0;  
        results_per_sample->elapsed_time_learning = 0;  
        results_per_sample->elapsed_time_sample = 0;
        results_per_sample->elapsed_time_load_sample = 0;  
        results_per_sample->elapsed_time_re_neurons = 0;
        results_per_sample->elapsed_time_re_synapses = 0;
    }


    //TODO

    // spikes info
    if(conf->store_generated_spikes == 1)
        results_per_sample->generated_spikes = (unsigned char **)calloc(n_neurons, sizeof(unsigned char *));
    
    if(conf->store_n_spikes == 1)
        results_per_sample->n_spikes_per_neuron = (int *)calloc(n_neurons, sizeof(int));

    // initialize generated spikes
    for (i = 0; i<n_neurons; i++){
        
        if(conf->store_generated_spikes == 1){
            results_per_sample->generated_spikes[i] = (unsigned char *)calloc((conf->time_steps), sizeof(unsigned char));
            results_per_sample->generated_spikes[i][0] = ' '; // first time step is not simulated
        }
        
        if(conf->store_n_spikes == 1)
            results_per_sample->n_spikes_per_neuron[i] = 0;
    }
}

void reinitialize_results_struct(simulation_results_t *results, simulation_configuration_t *conf, int n_samples, int n_neurons){
    
    int i;

    // initialize struct to store simulation configuration and outputs    
    if(conf->store_execution_times == 1){
        results->elapsed_time_epoch = 0;//(double *)calloc(conf->n_epochs, sizeof(double));
        results->elapsed_time = 0;
    }

    for(i = 0; i<n_samples; i++)
        reinitialize_sample_results_struct(&(results->results_per_sample[i]), conf, n_samples, n_neurons);
}


void reinitialize_sample_results_struct(simulation_results_per_sample_t *results_per_sample, simulation_configuration_t *conf, int n_samples, int n_neurons){
    
    int i, j;
    
    if(conf->store_execution_times == 1){
        results_per_sample->elapsed_time_neurons = 0;
        results_per_sample->elapsed_time_neurons_input = 0;
        results_per_sample->elapsed_time_neurons_output = 0;  
        results_per_sample->elapsed_time_synapses = 0;  
        results_per_sample->elapsed_time_synapses_input = 0;  
        results_per_sample->elapsed_time_synapses_output = 0;  
        results_per_sample->elapsed_time_learning = 0;    
        results_per_sample->elapsed_time_sample = 0;
        results_per_sample->elapsed_time_load_sample = 0;  
        results_per_sample->elapsed_time_re_neurons = 0;
        results_per_sample->elapsed_time_re_synapses = 0;    
    }
   
    for (i = 0; i<n_neurons; i++){
        for(j = 0; j<conf->time_steps; j++){
            if(conf->store_generated_spikes == 1){
                results_per_sample->generated_spikes[i][j] = 0;
            }
        }
        if(conf->store_n_spikes == 1){
            results_per_sample->n_spikes_per_neuron[i] = 0;
        }
    }
}



void free_results_struct_memory(simulation_results_t *results, simulation_configuration_t *conf, int n_samples, int n_neurons){
    for(int i = 0; i<n_samples; i++)
        free_sample_results_struct_memory(&(results->results_per_sample[i]), conf, n_samples, n_neurons);
    free(results->results_per_sample);
}

void free_sample_results_struct_memory(simulation_results_per_sample_t *results_per_sample, simulation_configuration_t *conf, int n_samples, int n_neurons){
    for (int i = 0; i<n_neurons; i++)
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




////// NEW FUNCTIONS FOR THE NEW DATA STRUCTS

/// @brief Allocate memory for a SNN structure
/// @param snn SNN structure
/// @param conf configuration structure with helpful information
void allocate_memory_for_SNN(GPU_SNN_t *snn, simulation_configuration_t *conf){

    // allocate memory for neurons and synapses

    snn->v                             = (float*)malloc(snn->n_neurons * conf->batch_size * sizeof(float));
    snn->v_thresh                      = (float*)malloc(snn->n_neurons * conf->batch_size * sizeof(float));
    snn->v_rest                        = (float*)malloc(snn->n_neurons * conf->batch_size * sizeof(float));
    snn->arrI                          = (float*)malloc(snn->n_neurons * conf->batch_size * sizeof(float));
    snn->r_period                      = (int*)malloc(snn->n_neurons * conf->batch_size * sizeof(int));
    snn->r_period_remain               = (int*)malloc(snn->n_neurons * conf->batch_size * sizeof(int));
    snn->res                           = (int*)malloc(snn->n_neurons * conf->batch_size * sizeof(int));
    snn->post_fired                    = (char*)malloc(snn->n_neurons * conf->batch_size * sizeof(char));
    snn->post_trace                    = (float*)malloc(snn->n_neurons * conf->batch_size * sizeof(float));
    snn->inR                           = (int*)malloc(snn->n_neurons * conf->batch_size * sizeof(int));
    snn->n_neuron_input_synapses       = (size_t*)malloc(snn->n_neurons * conf->batch_size * sizeof(size_t));
    snn->neuron_input_synapses_offset  = (size_t*)malloc(snn->n_neurons * conf->batch_size * sizeof(size_t));
    snn->n_neuron_output_synapses      = (size_t*)malloc(snn->n_neurons * conf->batch_size * sizeof(size_t));
    snn->neuron_output_synapses_offset = (size_t*)malloc(snn->n_neurons * conf->batch_size * sizeof(size_t));

    // allocate memory for neurons and synapses

    snn->w                             = (float*)malloc(snn->n_synapses * conf->batch_size * sizeof(float));
    snn->init_w                        = (float*)malloc(snn->n_synapses * conf->batch_size * sizeof(float));
    snn->dw                            = (float*)malloc(snn->n_synapses * conf->batch_size * sizeof(float));
    snn->delay                         = (int*)malloc(snn->n_synapses * conf->batch_size * sizeof(int));
    snn->lr                            = (int*)malloc(snn->n_synapses * conf->batch_size * sizeof(int));
    snn->pre_fired                     = (char*)malloc(snn->n_synapses * conf->batch_size * sizeof(char));
    snn->pre_trace                     = (float*)malloc(snn->n_synapses * conf->batch_size * sizeof(float));
    snn->pre_neuron_index              = (size_t*)malloc(snn->n_synapses * conf->batch_size * sizeof(size_t));
    snn->post_neuron_index             = (size_t*)malloc(snn->n_synapses * conf->batch_size * sizeof(size_t));
    

    // spk matrix
    snn->spk_matrix                    = (char*)malloc((snn->n_neurons + snn->n_input_neurons) * snn->LT * conf->batch_size * sizeof(char)); // time steps is temporal
}

void cpy_snn(GPU_SNN_t *snn, simulation_configuration_t *conf){

    size_t i, j, t, B, LT, N, iN, S;

    // temporal variables to store values
    float *tmp_v, *tmp_arrI, *tmp_w, *tmp_dw, *tmp_pre_trace, *tmp_post_trace;
    int *tmp_rperiod_remain;
    char *tmp_pstfired, *tmp_prefired;

    size_t padding = 32; // 32 bytes (256 bits) of padding for copied variables // TODO: improve

    //char *tmp_spk_mtrx;

    //float *tmp_thresh, *tmp_rest, *tmp_init_w;
    //int *tmp_rperiod, *tmp_res, *tmp_inR, *tmp_delay, *tmp_lr;
    //size_t *tmp_in_off, *tmp_pre_neuron_index, *tmp_post_neuron_index, *tmp_n_insyn;

    B = conf->batch_size;
    LT = snn->LT;
    N = snn->n_neurons;
    iN = snn->n_input_neurons;
    S = snn->n_synapses;

    // allocate neural properties
    tmp_v                         = snn->v;
    //tmp_thresh                    = snn->v_thresh;
    //tmp_rest                      = snn->v_rest;
    //tmp_rperiod                   = snn->r_period;
    tmp_rperiod_remain            = snn->r_period_remain;
    //tmp_res                       = snn->res;
    tmp_pstfired                  = snn->post_fired;
    tmp_post_trace                = snn->post_trace;
    tmp_arrI                      = snn->arrI;
    //tmp_inR                       = snn->inR;
    //tmp_n_insyn                   = snn->n_neuron_input_synapses;
    //tmp_in_off                    = snn->neuron_input_synapses_offset;

    tmp_w                         = snn->w;
    //tmp_init_w                    = snn->init_w;
    tmp_dw                        = snn->dw;
    tmp_pre_trace                 = snn->pre_trace;
    tmp_prefired                  = snn->pre_fired;
    //tmp_delay                     = snn->delay;
    //tmp_lr                        = snn->lr;
    //tmp_pre_neuron_index          = snn->pre_neuron_index;
    //tmp_post_neuron_index         = snn->post_neuron_index;

    //tmp_spk_mtrx                  = snn->spk_matrix;


    /* Memory allocation */
    snn->v                            = (float*)malloc(N * B * sizeof(float) + padding);
    snn->arrI                         = (float*)malloc(N * B * sizeof(float) + padding);
    snn->r_period_remain              = (int*)malloc(N * B * sizeof(int) + padding);
    snn->post_fired                   = (char*)malloc(N * B * sizeof(char) + padding);
    snn->post_trace                   = (float*)malloc(N * B * sizeof(float) + padding);


    snn->w                            = (float*)malloc(S * B * sizeof(float) + padding);
    snn->dw                           = (float*)malloc(S * B * sizeof(float) + padding);
    snn->pre_trace                    = (float*)malloc(S * B * sizeof(float) + padding);
    snn->pre_fired                    = (char*)malloc(S * B * sizeof(char) + padding); 
    

    //snn->v_thresh                     = (float*)malloc(N * B * sizeof(float));
    //snn->v_rest                       = (float*)malloc(N * B * sizeof(float));
    //snn->r_period                     = (int*)malloc(N * B * sizeof(int));
    //snn->res                          = (int*)malloc(N * B * sizeof(int));
    //snn->inR                          = (int*)malloc(N * B * sizeof(int));
    //snn->n_neuron_input_synapses      = (size_t*)malloc(N * B * sizeof(size_t));
    //snn->neuron_input_synapses_offset = (size_t*)malloc(N * B * sizeof(size_t));

    //snn->init_w                       = (float*)malloc(S * B * sizeof(float));
    //snn->delay                        = (int*)malloc(S * B * sizeof(int));
    //snn->lr                           = (int*)malloc(S * B * sizeof(int));
    //snn->pre_neuron_index             = (size_t*)malloc(S * B * sizeof(size_t));
    //snn->post_neuron_index            = (size_t*)malloc(S * B * sizeof(size_t));

    // reallocate spk matrix
    snn->spk_matrix                   = (char*)malloc((iN + N) * B * LT * sizeof(char) + padding);


    /* Initialization */
    for(i = 0; i<N; i++){

        for(j = 0; j<B; j++){

            // copy values
            snn->v[i * B + j] = tmp_v[i];
            //snn->v_thresh[i * B + j] = tmp_thresh[i];
            //snn->v_rest[i * B + j] = tmp_rest[i];
            //snn->r_period[i * B + j] = tmp_rperiod[i];
            snn->r_period_remain[i * B + j] = tmp_rperiod_remain[i];
            //snn->res[i * B + j] = tmp_res[i];
            snn->post_fired[i * B + j] = tmp_pstfired[i];
            snn->post_trace[i * B + j] = tmp_post_trace[i];
            //snn->inR[i * B + j] = tmp_inR[i];
            snn->arrI[i * B + j] = tmp_arrI[i];
            //snn->n_neuron_input_synapses[i * B + j] = tmp_n_insyn[i];
            //snn->neuron_input_synapses_offset[i * B + j] = tmp_in_off[i];
        }
    }

    for(i = 0; i<S; i++){

        for(j = 0; j<B; j++){

            // copy values
            snn->w[i * B + j] = tmp_w[i];
            //snn->init_w[i * B + j] = tmp_init_w[i];
            snn->dw[i * B + j] = tmp_dw[i];
            snn->pre_trace[i * B + j] = tmp_pre_trace[i];
            snn->pre_fired[i * B + j] = tmp_prefired[i];
            //snn->delay[i * B + j] = tmp_delay[i];
            //snn->lr[i * B + j] = tmp_lr[i];
            //snn->pre_neuron_index[i * B + j] = tmp_pre_neuron_index[i];
            //snn->post_neuron_index[i * B + j] = tmp_post_neuron_index[i];
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

    // initialize padding for floats
    for(size_t p = N * B ; p<padding / sizeof(float); p++){

        snn->v[p] = 0.0;                       
        snn->arrI[p] = 0.0;                         
        snn->post_trace[p] = 0.0;               
    }
    for(size_t p = N * B ; p<padding / sizeof(int); p++){

        snn->r_period_remain[p] = 0;              
    }
    for(size_t p = N * B ; p<padding / sizeof(char); p++){

        snn->post_fired[p] = 0;     
    }
    for(size_t p = S * B ; p<padding / sizeof(float); p++){

        snn->w[p] = 0.0;                    
        snn->dw[p] = 0.0;  
        snn->pre_trace[p] = 0.0;                
    }
    for(size_t p = S * B ; p<padding / sizeof(char); p++){

        snn->pre_fired[p] = 0;
    }

    for(size_t p = (iN + N) * B * LT ; p<padding / sizeof(char); p++){

        snn->spk_matrix[p] = 0;
    }
    // TODO
    // free original arrays memory
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
    snn->max_delay = max_delay;

    return max_delay;
}


// function to estimate structure size
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

/*double get_neuron_size(){
    
}*/

// TODO
double get_neuron_size(){
    
    double size = 0.0;
    // get memory related to one neuron
    

    return size;
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

double get_dataset_size(GPU_dataset_t *dataset){
    
    size_t nS = dataset->n_samples;
    size_t nF = dataset->n_features;
    size_t nSpks = dataset->n_spikes;

    return (
        // general and scalars
        (sizeof(GPU_dataset_t) +
        sizeof(int) + 
        sizeof(size_t) * 4 + 
        // arrays
        sizeof(size_t) * (nS + nS * nF * 4 + nSpks)) / 8.0
    );
}

double get_results_size(size_t N, size_t nS, size_t T){

    return (
        (sizeof(GPU_results_t) + 
        sizeof(int) * N * nS + // number of spike generated by each neuron
        sizeof(double) * 7 // variables to store execution times
        ) / 8.0
    );
}

double get_tmp_batch_size(size_t iN, size_t nS, size_t T){

    return (
        (sizeof(tmp_batch_cpu_t) + 
        sizeof(size_t) * iN * nS * T 
        ) / 8.0
    );
}



// this should allocate the arrays of LIF, instead of previously : generalization for new neuron models
void initialize_lif_neuron_CPU(GPU_SNN_t *snn, network_construction_lists_t *data, simulation_configuration_t *conf){

    size_t i;
    
    // initialize neuron parameters from array
    for(i = 0; i<snn->n_neurons; i++){

        snn->v_thresh[i] = data->v_thres_list[i];
        snn->v_rest[i] = data->v_rest_list[i]; // this or the next one?
        snn->v[i] = snn->v_rest[i]; 
        snn->res[i] = data->R_list[i];
        snn->r_period[i] = data->r_time_list[i];
        snn->r_period_remain[i] = -1;
        snn->n_neuron_input_synapses[i] = 0;
        snn->neuron_input_synapses_offset[i] = 0;

        // control variables
        snn->post_fired[i] = 0;
        snn->post_trace[i] = 0.0;
        snn->arrI[i] = 0.0;
        snn->inR[i] = 0;
    }
}

// TODO: simplify
void initialize_synapses_CPU(GPU_SNN_t *snn, network_construction_lists_t *data, simulation_configuration_t *conf){

    size_t i;

    for(i = 0; i<snn->n_synapses; i++){

        // read synapse parameters from lists
        snn->w[i] = data->weight_list[i];
        snn->init_w[i] = snn->w[i];
        snn->dw[i] = 0.0; 
        snn->delay[i] = data->delay_list[i];
        snn->lr[i] = 0; // temporal: this will call the function in an array of function pointers

        // initialized later
        snn->pre_neuron_index[i] = 0;
        snn->post_neuron_index[i] = 0;
        
        // control variables
        snn->pre_fired[i] = 0;
        snn->pre_trace[i] = 0.0;
    }
}


void count_neurons_input_synapses_and_conpute_offsets(GPU_SNN_t *snn, network_construction_lists_t *data){
    
    size_t i, j, off = 0;

    // loop over neurons
    for(i = 0; i<snn->n_neurons; i++){

        for(j = 0; j<(size_t)data->synaptic_connections[i][0]; j++){

            // [i * 2 + 1] is the index of the neuron, [i * 2 + 2] the number of input synapses 
            snn->n_neuron_input_synapses[i] += (size_t)data->synaptic_connections[i][j*2+2];
        }
    }

    // loop over neurons
    for(i = 0; i<snn->n_neurons; i++){

        snn->neuron_input_synapses_offset[i] = off;
        off += snn->n_neuron_input_synapses[i];
    }
}

void connect_synapses_to_network(GPU_SNN_t *snn, network_construction_lists_t *data){

    size_t i, j, l, next_syn = 0;

    for(i = 0; i<snn->n_neurons; i++){

        for(j = 0; j<(size_t)data->synaptic_connections[i][0]; j++){

            for(l = 0; l<(size_t)data->synaptic_connections[i][j * 2 + 2]; l++){
                

                snn->post_neuron_index[next_syn] = i + snn->n_input_neurons; // first [n_input_neurons] neurons are virtual input
                snn->pre_neuron_index[next_syn] = (size_t)data->synaptic_connections[i][j * 2 + 1];
                next_syn ++;
            }
        }
    }
}


void connect_network_input_criteria(GPU_SNN_t *snn, network_construction_lists_t *data, simulation_configuration_t *conf){

    // count input synapses for each neuron
    count_neurons_input_synapses_and_conpute_offsets(snn, data);

    // sets pre and post neurons indexes for each synapse
    connect_synapses_to_network(snn, data);
}

void reallocate_spk_matrix(GPU_SNN_t *snn, size_t N, size_t B, size_t T){

    free(snn->spk_matrix);
    snn->spk_matrix = (char*)calloc(N * B * T, sizeof(char));
}

GPU_SNN_t* initialize_network_cpu(simulation_configuration_t *conf, network_construction_lists_t *data){

    // allocate memory for SNN structure
    GPU_SNN_t *snn = (GPU_SNN_t*)malloc(sizeof(GPU_SNN_t));

    // initialize general info
    snn->n_neurons = data->n_neurons;
    snn->n_input_neurons = data->n_input_neurons;
    snn->n_output_neurons = data->n_output_neurons;

    snn->n_synapses = data->n_synapses;
    snn->n_input_synapses = data->n_input_synapses;
    snn->n_output_synapses = data->n_output_synapses;

    snn->max_spikes = conf->max_spikes; // DEPRECATED
    snn->max_delay = data->max_delay; // max_delay == LT
    snn->LT = data->max_delay;

    // allocate memory for the SNN structure arrays
    allocate_memory_for_SNN(snn, conf); // TODO: generalize?

    // initialize neurons
    initialize_lif_neuron_CPU(snn, data, conf);
    
    // initialize synapses
    initialize_synapses_CPU(snn, data, conf);

    // connect network
    connect_network_input_criteria(snn, data, conf);

    // return SNN structure
    return snn;
}

void deallocate_network_cpu(GPU_SNN_t *snn){

    // deallocate internal arrays
    free(snn->v);                            
    free(snn->v_thresh);                      
    free(snn->v_rest);                     
    free(snn->arrI);                      
    free(snn->r_period);                     
    free(snn->r_period_remain);              
    free(snn->res);             
    free(snn->post_fired);                   
    free(snn->inR);                  
    free(snn->n_neuron_input_synapses);      
    free(snn->neuron_input_synapses_offset); 
    free(snn->n_neuron_output_synapses);      
    free(snn->neuron_output_synapses_offset); 

    free(snn->w);
    free(snn->init_w);
    free(snn->dw);
    free(snn->pre_trace);
    free(snn->post_trace);
    free(snn->delay);
    free(snn->lr);
    free(snn->pre_neuron_index);
    free(snn->post_neuron_index);
    free(snn->pre_fired);
    
    free(snn->spk_matrix);

    // deallocate general structure
    free(snn);
}

void print_network(GPU_SNN_t *snn){

    size_t i;

    printf(" > Printing network:\n\n");

    printf(" > > Network general information:\n");
    printf(" > >>  N: %zu\n", snn->n_neurons);
    printf(" > >> iN: %zu\n", snn->n_input_neurons);
    printf(" > >> oN: %zu\n", snn->n_output_neurons);
    printf(" > >>  S: %zu\n", snn->n_synapses);
    printf(" > >> iS: %zu\n", snn->n_input_synapses);
    printf(" > >> oS: %zu\n", snn->n_output_synapses);
    printf(" > >> mS: %zu\n", snn->max_spikes);
    //printf(" > >> nN: %zu\n\n", snn->n_networks);

    printf(" > > Neurons data:\n");

    printf(" > >> V: [");
    for(i = 0; i<snn->n_neurons-1; i++){

        printf("%f, ", snn->v[i]);
    }
    printf("%f]\n", snn->v[i]);

    printf(" > >> Thresholds: [");
    for(i = 0; i<snn->n_neurons-1; i++){

        printf("%f, ", snn->v_thresh[i]);
    }
    printf("%f]\n", snn->v_thresh[i]);

    printf(" > >> Rests: [");
    for(i = 0; i<snn->n_neurons-1; i++){

        printf("%f, ", snn->v_rest[i]);
    }
    printf("%f]\n", snn->v_rest[i]);

    printf(" > >> Refract: [");
    for(i = 0; i<snn->n_neurons-1; i++){

        printf("%d (%d), ", snn->r_period[i], snn->r_period_remain[i]);
    }
    printf("%d (%d)]\n", snn->r_period[i], snn->r_period_remain[i]);

    printf(" > >> R: [");
    for(i = 0; i<snn->n_neurons-1; i++){

        printf("%d, ", snn->res[i]);
    }
    printf("%d]\n", snn->res[i]);

    printf(" > >> Post fired: [");
    for(i = 0; i<snn->n_neurons-1; i++){

        printf("%d, ", snn->post_fired[i]);
    }
    printf("%d]\n", snn->post_fired[i]);

    printf(" > >> I: [");
    for(i = 0; i<snn->n_neurons-1; i++){

        printf("%f, ", snn->arrI[i]);
    }
    printf("%f]\n", snn->arrI[i]);

    printf(" > >> N input synapses: [");
    for(i = 0; i<snn->n_neurons-1; i++){

        printf("%zu, ", snn->n_neuron_input_synapses[i]);
    }
    printf("%zu]\n", snn->n_neuron_input_synapses[i]);

    printf(" > >> Input synapses off: [");
    for(i = 0; i<snn->n_neurons-1; i++){

        printf("%zu, ", snn->neuron_input_synapses_offset[i]);
    }
    printf("%zu]\n", snn->neuron_input_synapses_offset[i]);



    printf(" > > Synapses data:\n");

    printf(" > >> W: [");
    for(i = 0; i<snn->n_synapses-1; i++){

        printf("%f, ", snn->w[i]);
    }
    printf("%f]\n", snn->w[i]);

    printf(" > >> Init W: [");
    for(i = 0; i<snn->n_synapses-1; i++){

        printf("%f, ", snn->init_w[i]);
    }
    printf("%f]\n", snn->init_w[i]);

    printf(" > >> dW: [");
    for(i = 0; i<snn->n_synapses-1; i++){

        printf("%f, ", snn->dw[i]);
    }
    printf("%f]\n", snn->dw[i]);

    printf(" > >> Delay: [");
    for(i = 0; i<snn->n_synapses-1; i++){

        printf("%d, ", snn->delay[i]);
    }
    printf("%d]\n", snn->delay[i]);

    printf(" > >> Lr: [");
    for(i = 0; i<snn->n_synapses-1; i++){

        printf("%d, ", snn->lr[i]);
    }
    printf("%d]\n", snn->lr[i]);

    printf(" > >> Pre neuron: [");
    for(i = 0; i<snn->n_synapses-1; i++){

        printf("%zu, ", snn->pre_neuron_index[i]);
    }
    printf("%zu]\n", snn->pre_neuron_index[i]);

    printf(" > >> Post neuron: [");
    for(i = 0; i<snn->n_synapses-1; i++){

        printf("%zu, ", snn->post_neuron_index[i]);
    }
    printf("%zu]\n", snn->post_neuron_index[i]);

    printf(" > >> Pre fired: [");
    for(i = 0; i<snn->n_synapses-1; i++){

        printf("%d, ", snn->pre_fired[i]);
    }
    printf("%d]\n", snn->pre_fired[i]);

    printf(" > >> Pre trace: [");
    for(i = 0; i<snn->n_synapses-1; i++){

        printf("%f, ", snn->pre_trace[i]);
    }
    printf("%f]\n", snn->pre_trace[i]);

    printf(" > >> Post trace: [");
    for(i = 0; i<snn->n_synapses-1; i++){

        printf("%f, ", snn->post_trace[i]);
    }
    printf("%f]\n", snn->post_trace[i]);

    fflush(stdout);
}

void print_networks(GPU_SNN_t *snn, simulation_configuration_t *conf){

    size_t i, b, B;
    B = conf->batch_size;

    printf(" > Printing network:\n\n");

    printf(" > > Network general information:\n");
    printf(" > >>  N: %zu\n", snn->n_neurons);
    printf(" > >> iN: %zu\n", snn->n_input_neurons);
    printf(" > >> oN: %zu\n", snn->n_output_neurons);
    printf(" > >>  S: %zu\n", snn->n_synapses);
    printf(" > >> iS: %zu\n", snn->n_input_synapses);
    printf(" > >> oS: %zu\n", snn->n_output_synapses);
    printf(" > >> mS: %zu\n", snn->max_spikes);
    //printf(" > >> nN: %zu\n\n", snn->n_networks);

    printf(" > > Neurons data:\n");

    printf(" > >> V: [");
    for(i = 0; i<snn->n_neurons-1; i++){
        
        printf("[");
        for(b = 0; b<B-1; b++)
            printf("%f, ", snn->v[i * B + b]);

        printf("%f]", snn->v[i * B + B - 1]);
    }
    printf("%f]\n", snn->v[i * B + B-1]);

    printf(" > >> Thresholds: [");
    for(i = 0; i<snn->n_neurons-1; i++){

        printf("[");
        for(b = 0; b<B-1; b++)
            printf("%f, ", snn->v_thresh[i * B + b]);

        printf("%f]", snn->v_thresh[i * B + B - 1]);    
    }
    printf("%f]\n", snn->v_thresh[i * B + B - 1]);

    printf(" > >> Rests: [");
    for(i = 0; i<snn->n_neurons-1; i++){

        printf("[");
        for(b = 0; b<B-1; b++)
            printf("%f, ", snn->v_rest[i * B + b]);

        printf("%f]", snn->v_rest[i * B + B - 1]);      
    }
    printf("%f]\n", snn->v_rest[i * B + B - 1]);

    printf(" > >> Refract: [");
    for(i = 0; i<snn->n_neurons-1; i++){

        printf("[");
        for(b = 0; b<B-1; b++)
            printf("%d (%d), ", snn->r_period[i * B + b], snn->r_period_remain[i * B + b]);

        printf("%d (%d)]", snn->r_period[i * B + B-1], snn->r_period_remain[i * B + B-1]);      
    }
    printf("%d (%d)]\n", snn->r_period[i * B + B-1], snn->r_period_remain[i * B + B-1]);

    printf(" > >> R: [");
    for(i = 0; i<snn->n_neurons-1; i++){

        printf("[");
        for(b = 0; b<B-1; b++)
            printf("%d, ", snn->res[i * B + b]);

        printf("%d]", snn->res[i * B + B-1]);      
    }
    printf("%d]\n", snn->res[i * B + B-1]);

    printf(" > >> Post fired: [");
    for(i = 0; i<snn->n_neurons-1; i++){

        printf("[");
        for(b = 0; b<B-1; b++)
            printf("%d, ", snn->post_fired[i * B + b]);

        printf("%d]", snn->post_fired[i * B + B-1]);      
    }
    printf("%d]\n", snn->post_fired[i * B + B-1]);

    printf(" > >> I: [");
    for(i = 0; i<snn->n_neurons-1; i++){

        printf("[");
        for(b = 0; b<B-1; b++)
            printf("%f, ", snn->arrI[i * B + b]);

        printf("%f]", snn->arrI[i * B + B-1]);      
    }
    printf("%f]\n", snn->arrI[i * B + B-1]);

    printf(" > >> N input synapses: [");
    for(i = 0; i<snn->n_neurons-1; i++){

        printf("[");
        for(b = 0; b<B-1; b++)
            printf("%zu, ", snn->n_neuron_input_synapses[i * B + b]);

        printf("%zu]", snn->n_neuron_input_synapses[i * B + B-1]);      
    }
    printf("%zu]\n", snn->n_neuron_input_synapses[i * B + B-1]);

    printf(" > >> Input synapses off: [");
    for(i = 0; i<snn->n_neurons-1; i++){

        printf("[");
        for(b = 0; b<B-1; b++)
            printf("%zu, ", snn->neuron_input_synapses_offset[i * B + b]);

        printf("%zu]", snn->neuron_input_synapses_offset[i * B + B-1]);      
    }
    printf("%zu]\n", snn->neuron_input_synapses_offset[i * B + B-1]);



    printf(" > > Synapses data:\n");

    printf(" > >> W: [");
    for(i = 0; i<snn->n_synapses-1; i++){

        printf("[");
        for(b = 0; b<B-1; b++)
            printf("%f, ", snn->w[i * B + b]);

        printf("%f]", snn->w[i * B + B-1]);      
    }
    printf("%f]\n", snn->w[i * B + B-1]);

    printf(" > >> Init W: [");
    for(i = 0; i<snn->n_synapses-1; i++){

        printf("[");
        for(b = 0; b<B-1; b++)
            printf("%f, ", snn->init_w[i * B + b]);

        printf("%f]", snn->init_w[i * B + B-1]);      
    }
    printf("%f]\n", snn->init_w[i * B + B-1]);

    printf(" > >> dW: [");
    for(i = 0; i<snn->n_synapses-1; i++){

        printf("[");
        for(b = 0; b<B-1; b++)
            printf("%f, ", snn->dw[i * B + b]);

        printf("%f]", snn->dw[i * B + B-1]);      
    }
    printf("%f]\n", snn->dw[i * B + B-1]);

    printf(" > >> Delay: [");
    for(i = 0; i<snn->n_synapses-1; i++){

        printf("[");
        for(b = 0; b<B-1; b++)
            printf("%d, ", snn->delay[i * B + b]);

        printf("%d]", snn->delay[i * B + B-1]);      
    }
    printf("%d]\n", snn->delay[i * B + B-1]);

    printf(" > >> Lr: [");
    for(i = 0; i<snn->n_synapses-1; i++){

        printf("[");
        for(b = 0; b<B-1; b++)
            printf("%d, ", snn->lr[i * B + b]);

        printf("%d]", snn->lr[i * B + B-1]);      
    }
    printf("%d]\n", snn->lr[i * B + B-1]);

    printf(" > >> Pre neuron: [");
    for(i = 0; i<snn->n_synapses-1; i++){

        printf("[");
        for(b = 0; b<B-1; b++)
            printf("%zu, ", snn->pre_neuron_index[i * B + b]);

        printf("%zu]", snn->pre_neuron_index[i * B + B-1]);      
    }
    printf("[");
    for(b = 0; b<B-1; b++)
        printf("%zu, ", snn->pre_neuron_index[i * B + b]);
    printf("%zu]\n", snn->pre_neuron_index[i * B + B-1]);

    printf(" > >> Post neuron: [");
    for(i = 0; i<snn->n_synapses-1; i++){

        printf("[");
        for(b = 0; b<B-1; b++)
            printf("%zu, ", snn->post_neuron_index[i * B + b]);

        printf("%zu]", snn->post_neuron_index[i * B + B-1]);      
    }
    printf("%zu]\n", snn->post_neuron_index[i * B + B-1]);

    printf(" > >> Pre fired: [");
    for(i = 0; i<snn->n_synapses-1; i++){

        printf("[");
        for(b = 0; b<B-1; b++)
            printf("%d, ", snn->pre_fired[i * B + b]);

        printf("%d]", snn->pre_fired[i * B + B-1]);      
    }
    printf("%d]\n", snn->pre_fired[i * B + B-1]);

    printf(" > >> Pre trace: [");
    for(i = 0; i<snn->n_synapses-1; i++){

        printf("[");
        for(b = 0; b<B-1; b++)
            printf("%f, ", snn->pre_trace[i * B + b]);

        printf("%f]", snn->pre_trace[i * B + B-1]);      
    }
    printf("%f]\n", snn->pre_trace[i * B + B-1]);

    printf(" > >> Post trace: [");
    for(i = 0; i<snn->n_synapses-1; i++){

        printf("[");
        for(b = 0; b<B-1; b++)
            printf("%f, ", snn->post_trace[i * B + b]);

        printf("%f]", snn->post_trace[i * B + B-1]);      
    }
    printf("%f]\n", snn->post_trace[i * B + B-1]);

    fflush(stdout);
}



void print_dataset(GPU_dataset_t *dataset){

    size_t i, j, l, next = 0;

    printf(" > Printing dataset: \n");
    for(i = 0; i<dataset->n_samples; i++){

        printf(" > > Sample %zu (offset = %zu)\n", i, dataset->sample_offset[i]);

        // loop over features
        for(j = 0; j<dataset->n_features; j++){

            printf(" > >> Feature %zu (offset = %zu): [", j, dataset->feature_offset[i * dataset->n_features + j]);

            // print spikes in feature
            for(l = 0; l<dataset->n_spikes_per_feature[i * dataset->n_features + j]-1; l++){

                printf("%zu, ", dataset->spikes[next]);
                next++;
            }
            printf("%zu]\n", dataset->spikes[next]);
            next++;
        }
    }
}

// initialize struct to store batch results
GPU_results_t* initialize_batch_results_cpu(simulation_configuration_t *conf, size_t N, size_t batch_size, size_t T){

    GPU_results_t *results = (GPU_results_t*)malloc(sizeof(GPU_results_t));

    // number of spikes generated by each neuron processing n_samples
    results->n_spks = (int*)calloc(N * batch_size, sizeof(int));

    // initialize execution times
    results->t = 0.0; // total execution time
    results->t_in = 0.0; // time processing input spikes
    results->t_v = 0.0; // time processing neurons dynamics
    results->t_out = 0.0; // time processing output spikes
    results->t_learn = 0.0; // time learning
    results->t_reinit = 0.0; // network reinitialization
    results->t_load = 0.0; // loading sample or batch in network

    return results;
}


void reinitialize_batch_results_cpu(GPU_results_t *results, simulation_configuration_t *conf, size_t N, size_t batch_size, size_t T){

    // number of spikes generated by each neuron processing n_samples
    for(size_t i = 0; i<N; i++){
        
        for(size_t j = 0; j<batch_size; j++){

            results->n_spks[i * batch_size + j] = 0;
        }
    }

    // initialize execution times
    results->t = 0.0; // total execution time
    results->t_in = 0.0; // time processing input spikes
    results->t_v = 0.0; // time processing neurons dynamics
    results->t_out = 0.0; // time processing output spikes
    results->t_learn = 0.0; // time learning
    results->t_reinit = 0.0; // network reinitialization
    results->t_load = 0.0; // loading sample or batch in network
}

void deallocate_results_struct(GPU_results_t *results, simulation_configuration_t *conf){

    // deallocate internal structs
    free(results->n_spks);

    // deallocate general struct
    free(results);
}

void acc_batch_execution_times(GPU_results_t *results, GPU_results_t *batch_results){

    results->t        += batch_results->t; // total execution time
    results->t_in     += batch_results->t_in; // time processing input spikes
    results->t_v      += batch_results->t_v; // time processing neurons dynamics
    results->t_out    += batch_results->t_out; // time processing output spikes
    results->t_learn  += batch_results->t_learn; // time learning
    results->t_reinit += batch_results->t_reinit; // network reinitialization
    results->t_load   += batch_results->t_load; // loading sample or batch in network
}




// CUDA
void compute_input_synapses_subsets(cuda_info_t *cuda_info, GPU_SNN_t *snn, simulation_configuration_t *conf){

    size_t i, j;

    // get available bytes in shared memory per each block
    double shared_memory = cuda_info->shared_memory_mem[0];

    // store max synapses to be computed per thread
    size_t max_synapses_per_thr = 100; // temporal


    // define arrays to store offsets and numbers of offsets
    size_t *n_subsets_neuron; // number of subsets of input synapses for each neuron
    size_t *off_subsets_neuron; // first subset of the neuron (if neuron[i-1] has 3 subsets, the neuron[i] subset is 3)

    size_t *n_input_synapses_subset; // number of input synapses for each subset
    size_t *off_input_synapses_subset; // number of synapses processed before this subset

    size_t n_subsets = 0; // total number of subsets (the sum of the number of subsets of each neuron)
    size_t r_subset = 0; // helper variable

    // allocate memory
    off_subsets_neuron = (size_t*)calloc(snn->n_neurons, sizeof(size_t));
    n_subsets_neuron = (size_t*)calloc(snn->n_neurons, sizeof(size_t));


    // count total number of subsets first, and set number of subsets per neuron and the offset of the neuron
    for(i = 0; i<snn->n_neurons; i++){

        // subset offset is the number of previous offsets
        off_subsets_neuron[i] = n_subsets; 

        // compute number of subsets of the neuron and the number of synapses in the last subset
        n_subsets_neuron[i] = snn->n_neuron_input_synapses[i] / max_synapses_per_thr;

        // if there are remaining subsets, add one more subset
        if((snn->n_neuron_input_synapses[i] % max_synapses_per_thr) > 0)
            n_subsets_neuron[i] += 1;

        // compute total number of subsets
        n_subsets += n_subsets_neuron[i];
    }


    // allocate memory synapses information of each offset
    off_input_synapses_subset = (size_t*)calloc(n_subsets, sizeof(size_t)); // offset of the first synapse on each subset
    n_input_synapses_subset = (size_t*)calloc(n_subsets, sizeof(size_t)); // number of synapses on each subset
    
    // control variables
    size_t next_subset = 0;
    size_t n_synapses = 0; // total number of synapses processed to be used as offset

    // loop over input synapses of the neuron to form the subsets
    for(i = 0; i<snn->n_neurons; i++){

        // all subsets until the last has the maximum number of synapses
        for(j = 0; j<n_subsets_neuron[i]-1; j++){

            off_input_synapses_subset[next_subset] = n_synapses; // set the offset of the subset in the array of synapses
            n_input_synapses_subset[next_subset] = max_synapses_per_thr; // subset has the maximum number of synapses in it

            // update control variables
            n_synapses += max_synapses_per_thr; // processed synapses
            next_subset ++; // synapse counter
        }

        // set the last subset
        off_input_synapses_subset[next_subset] = n_synapses;

        // if r_subset > 0, the last subset of the neuron has less synapses than the max 
        r_subset = snn->n_neuron_input_synapses[i] % max_synapses_per_thr; 
        if(r_subset == 0){
            
            n_input_synapses_subset[next_subset] = max_synapses_per_thr;
            n_synapses += max_synapses_per_thr;
        }
        else{
            
            n_input_synapses_subset[next_subset] = r_subset;
            n_synapses += r_subset;
        }

        next_subset ++;
    }
}


void configure_cuda_simulation(cuda_info_t *cuda_info, GPU_SNN_t *snn, GPU_dataset_t *dataset, simulation_configuration_t *conf){

    double available_memory, total_mem_snn;

    // store information in cuda info struct
    cuda_info->n_samples = dataset->n_samples;
    cuda_info->time_steps = conf->time_steps;
    cuda_info->multi_gpu_allowed = conf->multi_gpu;
    cuda_info->batch_size = conf->batch_size;

    // get memory occupation data for GPU by data structures
    cuda_info->network_size = get_snn_size(snn);
    cuda_info->network_cpy_size = get_snn_cpy_size(snn);
    cuda_info->dataset_size = get_dataset_size(dataset); 
    cuda_info->results_size = get_results_size(snn->n_neurons, dataset->n_samples, conf->time_steps);


    // compute available memory for copies after copying network, datasets, results, and reserving some memory for cuda 
    available_memory = cuda_info->gpu_free_mem[0];
    available_memory -= cuda_info->dataset_size;
    available_memory -= cuda_info->results_size;
    available_memory -= cuda_info->network_size;
    available_memory -= 500 * 1024 * 1024; // reserve 500MB for cuda // too much?

    // usable memory
    printf(" > Available memory for SNN copies in GPU = %.2fMB\n", available_memory / 1024.0 / 1024.0);


    // check if all batch copies can be stored in the device
    total_mem_snn = cuda_info->network_size + cuda_info->network_cpy_size * (float)conf->batch_size;
    cuda_info->dev_batch_size = conf->batch_size;

    if(total_mem_snn < available_memory){

        cuda_info->n_networks_per_dev = conf->batch_size;
        cuda_info->dev_batch_size = conf->batch_size;
    }
    else{

        // TODO
    }

    // TODO: multi-gpu / single-GPU differenciation


    // simple configuration: one thread per each neuron / synapse in batch
    cuda_info->n_thr_per_blk_neurons_x = 512;
    cuda_info->n_thr_per_blk_neurons_y = 1;
    cuda_info->n_thr_per_blk_neurons_z = 1;
    cuda_info->n_blk_neurons_x = ((snn->n_neurons * conf->batch_size) / cuda_info->n_thr_per_blk_neurons_x) + 1;
    cuda_info->n_blk_neurons_y = 1;
    cuda_info->n_blk_neurons_z = 1;

    cuda_info->n_threads_per_blk_synapses_x = 512;
    cuda_info->n_threads_per_blk_synapses_y = 1;
    cuda_info->n_threads_per_blk_synapses_z = 1;
    cuda_info->n_blk_synapses_x = (snn->n_synapses * conf->batch_size) / cuda_info->n_threads_per_blk_synapses_x + 1;
    cuda_info->n_blk_synapses_y = 1;
    cuda_info->n_blk_synapses_z = 1;

    cuda_info->n_thr_per_blk_in_neurons_x = 512;
    cuda_info->n_thr_per_blk_in_neurons_y = 1;
    cuda_info->n_thr_per_blk_in_neurons_z = 1;
    cuda_info->n_blk_in_neurons_x = snn->n_input_neurons * conf->batch_size / cuda_info->n_thr_per_blk_in_neurons_x + 1;
    cuda_info->n_blk_in_neurons_y = 1;
    cuda_info->n_blk_in_neurons_z = 1;

    cuda_info->n_thr_per_blk_all_neurons_x = 512;
    cuda_info->n_thr_per_blk_all_neurons_y = 1;
    cuda_info->n_thr_per_blk_all_neurons_z = 1;
    cuda_info->n_blk_all_neurons_x = (snn->n_input_neurons + snn->n_neurons) * conf->batch_size / cuda_info->n_thr_per_blk_all_neurons_x + 1;
    cuda_info->n_blk_all_neurons_y = 1;
    cuda_info->n_blk_all_neurons_z = 1;
    
    cuda_info->n_thr_per_blk_uw_x = 512;
    cuda_info->n_thr_per_blk_uw_y = 1;
    cuda_info->n_thr_per_blk_uw_z = 1;
    cuda_info->n_blk_uw_x = snn->n_synapses / cuda_info->n_thr_per_blk_uw_x + 1;
    cuda_info->n_blk_uw_y = 1;
    cuda_info->n_blk_uw_z = 1;


    // {N * batch_size * thrN} threads, max per block = 512
    size_t neurons_is = snn->n_neurons * thrN;

    cuda_info->n_thr_per_blk_is_x = 512;
    cuda_info->n_thr_per_blk_is_y = 1;
    cuda_info->n_thr_per_blk_is_z = 1;
    cuda_info->n_blk_is_x = (neurons_is * conf->batch_size) / cuda_info->n_thr_per_blk_is_x + 1;
    cuda_info->n_blk_is_y = 1;
    cuda_info->n_blk_is_z = 1;


    // find the maximum number of possible copies

    // memory of each kernel for thread (bits)
    
    /*
    int max_copies = (int)(cuda_info->gpu_usable_mem / cuda_info->network_cpy_size);
    int n_cps = max_copies;
    int valid = 0;

    // cost of each kernel (in Bytes)
    double reinit_kernel = 64 / 8;
    double load_sample_kernel = 672 / 8;
    double reinit_neurons_kernel = 256 / 8;
    double lif_in = 832 / 8;
    double lif_out = 480 / 8;

    int n_rk_threads, n_lsk_threads, n_neur_threads;

    // compute the number of copies that can be stored for the network
    while(valid == 0){

        n_rk_threads = ((((gpu_snn_in_cpu->n_neurons + gpu_snn_in_cpu->n_input_neurons) * gpu_snn_in_cpu->max_spikes * n_cps) / 1024) + 1) * 1024;
        n_lsk_threads = ((gpu_snn_in_cpu->n_input_neurons * gpu_snn_in_cpu->max_spikes * n_cps / 1024) + 1) * 1024;
        n_neur_threads = (((gpu_snn_in_cpu->n_neurons * n_cps) / 516) + 1) * 516;
        
        if((n_rk_threads * reinit_kernel) < cuda_info->gpu_usable_mem && (n_lsk_threads * load_sample_kernel) < cuda_info->gpu_usable_mem &&
            (n_neur_threads * lif_in) < cuda_info->gpu_usable_mem){
                valid = 1;
        }
        else{
            n_cps --;
        }
    }
    
*/
    /* compute the number of networks and samples per batch per device */
/*
    // multi gpu not allowed
    if(cuda_info->multi_gpu_allowed == 0){

        if(conf->batch_size < n_cps){
            
            // the maximum number of copies we need is the size of the batch
            cuda_info->n_networks = conf->batch_size;

            // revise this, only works for unigpu
            cuda_info->nDevices = 1;
            cuda_info->n_batch_per_dev = cuda_info->batch_size;
            cuda_info->n_networks_per_dev = cuda_info->n_networks;
        }

        else{

            // number of networks should be multiple of batch_size, so compute
            cuda_info->n_networks = conf->batch_size / 2;
            
            // loop until the number of copies is bigger than the maximum
            while(cuda_info->n_networks > n_cps){
                cuda_info->n_networks /= 2;
            }

            // store
            cuda_info->nDevices = 1;
            cuda_info->n_batch_per_dev = cuda_info->batch_size; // only one device
            cuda_info->n_networks_per_dev = cuda_info->n_networks;
        }
    }
    // number of devices to use specified
    else if(cuda_info->multi_gpu_allowed != 0){

        // check if there are enough devices available
        if(cuda_info->multi_gpu_allowed > 0){

            if(cuda_info->multi_gpu_allowed < cuda_info->nDevices){
                cuda_info->nDevices = cuda_info->multi_gpu_allowed;
            }
            else if(cuda_info->nDevices < cuda_info->multi_gpu_allowed){
                printf(" > WARNING: Only %d devices available!\n", cuda_info->nDevices);
            }
        }


        // compute number of batches and networks per dev
        cuda_info->n_batch_per_dev = cuda_info->batch_size / cuda_info->nDevices; // number of samples per batch computed in each device

        // if the number of samples per batch is bigger than the maximum number of network copies per dev, compute number of cpies
        cuda_info->n_networks = cuda_info->n_batch_per_dev;
        
        if(cuda_info->n_batch_per_dev > n_cps){

            cuda_info->n_networks = cuda_info->n_batch_per_dev / 2;

            while(cuda_info->n_networks > n_cps){
                cuda_info->n_networks /= 2;
            }
        }

        cuda_info->n_networks_per_dev = cuda_info->n_networks;
    }

    // set n networks for the snn
    gpu_snn_in_cpu->n_networks = cuda_info->n_networks_per_dev;


    printf("\n\n > Printing cuda simulation info:\n");
    printf(" >> Batch size = %d\n", cuda_info->batch_size);
    printf(" >> Number of networks = %d\n", cuda_info->n_networks);
    printf(" >> Number of devices = %d\n", cuda_info->nDevices);
    printf(" >> Number of samples in batch simulated by each device = %d\n", cuda_info->n_batch_per_dev);
    printf(" >> Number of networks in each device = %d\n\n", cuda_info->n_networks_per_dev);

    // compute kernel configuration for each execution part
    // The execution is done by a 3D grid, where each dimension of the grid
    // contains the blocks and thredas to simulate a cpy


    // REVISE THIS FOR MULTIGPU
    cuda_info->n_threads_per_blk_rsm_x = 516;
    cuda_info->n_threads_per_blk_rsm_y = 1;
    cuda_info->n_threads_per_blk_rsm_z = 1;
    cuda_info->n_blk_rsm_x =
            ((unsigned int)(gpu_snn_in_cpu->n_neurons + gpu_snn_in_cpu->n_input_neurons) *
            (unsigned int)gpu_snn_in_cpu->max_spikes *
            (unsigned int)gpu_snn_in_cpu->n_networks) /
            (unsigned int)cuda_info->n_threads_per_blk_rsm_x + 1;
    cuda_info->n_blk_rsm_y = 1;
    cuda_info->n_blk_rsm_z = 1;

    cuda_info->n_threads_per_blk_ls_x = 516;
    cuda_info->n_threads_per_blk_ls_y = 1;
    cuda_info->n_threads_per_blk_ls_z = 1;
    cuda_info->n_blk_ls_x = ((unsigned int)gpu_snn_in_cpu->n_input_neurons * (unsigned int)gpu_snn_in_cpu->n_networks) / (unsigned int)cuda_info->n_threads_per_blk_ls_x + 1;
    cuda_info->n_blk_ls_y = 1;
    cuda_info->n_blk_ls_z = 1;

    cuda_info->n_threads_per_blk_nrs_x = 516;
    cuda_info->n_threads_per_blk_nrs_y = 1;
    cuda_info->n_threads_per_blk_nrs_z = 1;
    cuda_info->n_blk_nrs_x = ((unsigned int)gpu_snn_in_cpu->n_neurons * (unsigned int)gpu_snn_in_cpu->n_networks) / (unsigned int)cuda_info->n_threads_per_blk_nrs_x + 1;
    cuda_info->n_blk_nrs_y = 1;
    cuda_info->n_blk_nrs_z = 1;

    cuda_info->n_threads_per_blk_synapses_x = 516;
    cuda_info->n_threads_per_blk_synapses_y = 1;
    cuda_info->n_threads_per_blk_synapses_z = 1;
    cuda_info->n_blk_synapses_x = ((unsigned int)gpu_snn_in_cpu->n_synapses * (unsigned int)gpu_snn_in_cpu->n_networks) / (unsigned int)cuda_info->n_threads_per_blk_synapses_x + 1;
    cuda_info->n_blk_synapses_y = 1;
    cuda_info->n_blk_synapses_z = 1;

    cuda_info->n_threads_per_blk_uw_x = 516;
    cuda_info->n_threads_per_blk_uw_y = 1;
    cuda_info->n_threads_per_blk_uw_z = 1;
    cuda_info->n_blk_uw_x = ((unsigned int)gpu_snn_in_cpu->n_synapses) / (unsigned int)cuda_info->n_threads_per_blk_synapses_x + 1;
    cuda_info->n_blk_uw_y = 1;
    cuda_info->n_blk_uw_z = 1;
    */
}