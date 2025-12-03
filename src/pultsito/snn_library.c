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
    snn->Dw = (double *)calloc(snn->n_synapses, sizeof(double));
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
        neuron->n_last_spikes = 3;
        neuron->t_last_spikes = (int*)malloc(neuron->n_last_spikes * sizeof(int));
        for(int j = 0; j<neuron->n_last_spikes; j++){
            neuron->t_last_spikes[j] = -1;
        }
        neuron->next_last_spike = 0;   
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

    // weights change
    cp_snn->Dw = (double*)calloc(cp_snn->n_synapses, sizeof(double));

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
    gpu_snn->n_last_spikes = 3; // TEMPORAL // TODO
    gpu_snn->t_last_spikes = (int*)calloc(snn->n_neurons * gpu_snn->n_last_spikes, sizeof(int)); // 0
    gpu_snn->next_last_spike = (int*)calloc(snn->n_neurons, sizeof(int)); // '
    gpu_snn->n_neuron_input_synapses = (int*)malloc(snn->n_neurons * sizeof(int));
    gpu_snn->neuron_input_synapses_offset = (int*)malloc(snn->n_neurons * sizeof(int));

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
        for(j = 0; j<gpu_snn->n_last_spikes; j++)
            gpu_snn->t_last_spikes[i * gpu_snn->n_last_spikes + j] = -1; // initialize to -1, no spikes

        // input synapses and offsets
        gpu_snn->n_neuron_input_synapses[i] = snn->lif_neurons[i].n_input_synapse;
        gpu_snn->neuron_input_synapses_offset[i] = offset;

        // update offset for the next iteration
        offset += gpu_snn->n_neuron_input_synapses[i];    
    }

    // synapses - allocate and cpy
    gpu_snn->w = (float*)malloc(snn->n_synapses * sizeof(float));
    gpu_snn->dw = (float*)calloc(snn->n_synapses, sizeof(float)); // 0
    gpu_snn->delay = (int*)malloc(snn->n_synapses * sizeof(int));
    gpu_snn->lr = (int*)malloc(snn->n_synapses * sizeof(int));
    gpu_snn->pre_neuron_index = (int*)malloc(snn->n_synapses * sizeof(int));
    gpu_snn->post_neuron_index = (int*)malloc(snn->n_synapses * sizeof(int));
    gpu_snn->next_pre_spike = (int*)malloc(snn->n_synapses * sizeof(int));
    gpu_snn->next_post_spike = (int*)malloc(snn->n_synapses * sizeof(int));

    for(i = 0; i<snn->n_synapses; i++){

        gpu_snn->w[i] = (float)(snn->synapses[i].w);
        gpu_snn->dw[i] = 0.0;
        gpu_snn->delay[i] = snn->synapses[i].delay;
        gpu_snn->lr[i] = snn->synapses[i].lr;
        gpu_snn->pre_neuron_index[i] = snn->synapses[i].pre_neuron_index;
        gpu_snn->post_neuron_index[i] = snn->synapses[i].post_neuron_index;
        gpu_snn->next_pre_spike[i] = 0;
        gpu_snn->next_post_spike[i] = 0;

        // if it is an input synapse, then it's pre neuron is a virtual one. Map indexes
        if(gpu_snn->delay[i] > 0) // to detect input synapses
            gpu_snn->pre_neuron_index[i] += snn->n_input; // n_input_synapses??
        gpu_snn->post_neuron_index[i] += snn->n_input;
    }

    // spike matrix
    gpu_snn->spk_matrix = (int*)malloc((snn->n_neurons + snn->n_input) * conf->max_input_spikes * sizeof(int));
    for(i = 0; i<snn->n_neurons + snn->n_input; i++)
        for(j = 0; j<conf->max_input_spikes; j++)
            gpu_snn->spk_matrix[i * conf->max_input_spikes + j] = -1; // no spikes yet


    return gpu_snn;
}

GPU_dataset_t* dataset_CPU2GPU_mapping(input_data_t *dataset, simulation_configuration_t *conf){

    int i, j, l;
    size_t next_spike, next_feature;
    int n_samples, n_features;

    GPU_dataset_t *gpu_dataset;
    gpu_dataset = (GPU_dataset_t *)malloc(sizeof(GPU_dataset_t));
    
    // copy first to a temporal struct in the CPU
    gpu_dataset->type = dataset->type;
    gpu_dataset->n_classes = dataset->n_classes;
    gpu_dataset->n_samples = dataset->n_samples;
    gpu_dataset->n_features = dataset->image_size; // TODO: generalize
    gpu_dataset->n_spikes_per_feature = (int*)malloc(dataset->n_samples * dataset->image_size * sizeof(int)); 

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
    gpu_dataset->spikes = (int*)malloc(gpu_dataset->n_spikes * sizeof(int));


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


// function to estimate structure size
double get_gpu_snn_size(GPU_SNN_t *gpu_snn){
    
    // n_input neurons or synapses?
    return (
        (sizeof(GPU_SNN_t) +
        sizeof(int) * gpu_snn->n_neurons * 6 + sizeof(float) * gpu_snn->n_neurons * 3 + 
        sizeof(int) * (gpu_snn->n_synapses - gpu_snn->n_output_synapses) +
        sizeof(float) * gpu_snn->n_synapses * 2 + sizeof(int) * gpu_snn->n_synapses * 4 + 
        sizeof(int) * (gpu_snn->n_neurons + gpu_snn->n_input_neurons) * gpu_snn->max_spikes) / 8.0
    );
}

double get_gpu_snn_cpy_size(GPU_SNN_t *gpu_snn){

    // n_input neurons or synapses?
    // v, r_period_remain, t_last_spike, w, dw, spk_matrix
    return (
        (sizeof(int) * gpu_snn->n_neurons * 2 + sizeof(float) * gpu_snn->n_neurons + 
        sizeof(float) * gpu_snn->n_synapses * 2 + 
        sizeof(int) * (gpu_snn->n_neurons + gpu_snn->n_input_neurons) * gpu_snn->max_spikes) / 8.0
    );
}

double get_gpu_dataset_size(GPU_dataset_t *gpu_dataset){
    
    return (
        (sizeof(GPU_dataset_t) +
        sizeof(int) * gpu_dataset->n_spikes +
        sizeof(int) * gpu_dataset->n_samples + 
        sizeof(int) * gpu_dataset->n_samples * gpu_dataset->n_features) / 8.0
    );
}

double get_gpu_results_size(int n_neurons, int n_samples, int L){

    return (
        (sizeof(GPU_results_t) + 
        sizeof(int) * n_neurons * n_samples +
        sizeof(char) * n_neurons * n_samples * L) / 8.0
    );
}


void configure_cuda_simulation(cuda_info_t *cuda_info, GPU_SNN_t *gpu_snn_in_cpu, GPU_dataset_t *gpu_dataset_in_cpu, simulation_configuration_t *conf){

    // store information in cuda info struct
    cuda_info->n_samples = gpu_dataset_in_cpu->n_samples;
    cuda_info->time_steps = conf->time_steps;
    cuda_info->batch_size = conf->batch_size;
    cuda_info->multi_gpu_allowed = conf->multi_gpu;

    // get memory occupation data for GPU by data structures
    cuda_info->network_size = get_gpu_snn_size(gpu_snn_in_cpu);
    cuda_info->network_cpy_size = get_gpu_snn_cpy_size(gpu_snn_in_cpu);
    cuda_info->dataset_size = get_gpu_dataset_size(gpu_dataset_in_cpu); 
    cuda_info->results_size = get_gpu_results_size(gpu_snn_in_cpu->n_neurons, gpu_dataset_in_cpu->n_samples, conf->max_input_spikes);

    printf(" > Free memory in GPU = %.2f\n", cuda_info->gpu_usable_mem / 1024.0 / 1024.0);

    printf(" >> SNN size = %.2f MB (cpy size = %.2f MB)\n >> Dataset size = %.2f MB\n >> Results struct size = %.2f MB\n",
            cuda_info->network_size / 1024.0 / 1024.0, cuda_info->network_cpy_size / 1024.0 / 1024.0,
            cuda_info->dataset_size / 1024.0 / 1024.0, cuda_info->results_size / 1024.0 / 1024.0);
    fflush(stdout);


    // get free memory
    cuda_info->gpu_usable_mem -= cuda_info->dataset_size;
    cuda_info->gpu_usable_mem -= cuda_info->results_size;
    cuda_info->gpu_usable_mem -= cuda_info->network_size;

    // store memory for cuda (500 MB)
    cuda_info->gpu_usable_mem -= (500 * 1024 * 1024); // B

    // usable memory
    printf(" > Free memory in GPU = %.2f\n", cuda_info->gpu_usable_mem / 1024.0 / 1024.0);


    // find the maximum number of possible copies

    // memory of each kernel for thread (bits)
    
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
    
    // compute the number of networks
    if(conf->batch_size < n_cps) {
        cuda_info->n_networks = conf->batch_size;

        // revise this, only works for unigpu
        cuda_info->nDevices = 1;
        cuda_info->n_batch_per_dev = cuda_info->batch_size;
        cuda_info->n_networks_per_dev = cuda_info->n_networks;
    }
    else { // conf->batch_size > n_cps
        // batch_size should be multiple of batch_size
        cuda_info->n_networks = conf->batch_size / 2;
        
        while(cuda_info->n_networks > n_cps){
            cuda_info->n_networks /= 2;
        }



        // multi gpu
        if(cuda_info->multi_gpu_allowed == 1){

            int n_required_devices = cuda_info->batch_size / cuda_info->n_networks;
        
            // set the number of devices that will be used
            if(n_required_devices < cuda_info->nDevices){
                cuda_info->nDevices = n_required_devices;
            }

            // compute the number of networks per GPU
            cuda_info->n_batch_per_dev = cuda_info->batch_size / cuda_info->nDevices;

            // compute the number of networks per each device (ideally, the number of samples to be computed)
            cuda_info->n_networks_per_dev = cuda_info->n_batch_per_dev;

            // if the number of network per device is bigger than the maximum, change
            if(cuda_info->n_networks_per_dev > cuda_info->n_networks)
                cuda_info->n_networks_per_dev = cuda_info->n_networks;
        }
        // balance work in all available gpus
        else if(cuda_info->multi_gpu_allowed == 2){

            // compute how much samples are simulated on each device
            cuda_info->n_batch_per_dev = cuda_info->batch_size / cuda_info->nDevices; 

            // try to use one network per each sample
            cuda_info->n_networks_per_dev = cuda_info->n_batch_per_dev;
            // if a devoce do not have enough memory to store the ideal number of copies, reduce the number
            if(cuda_info->n_networks_per_dev > cuda_info->n_networks){
                
                int valid = 0;

                while(valid == 0){
                    
                    cuda_info->n_networks_per_dev /= 2;

                    if(cuda_info->n_networks_per_dev <= cuda_info->n_networks)
                        valid = 1;
                }
            }

        }

        // if multigpu not allowed, then copy all information for unigpu
        else{

            cuda_info->nDevices = 1;
            cuda_info->n_batch_per_dev = cuda_info->batch_size;
            cuda_info->n_networks_per_dev = cuda_info->n_networks;

        }
    }

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
    cuda_info->n_blk_rsm_x = ((gpu_snn_in_cpu->n_neurons + gpu_snn_in_cpu->n_input_neurons) * gpu_snn_in_cpu->max_spikes * gpu_snn_in_cpu->n_networks) / cuda_info->n_threads_per_blk_rsm_x + 1;
    cuda_info->n_blk_rsm_y = 1;
    cuda_info->n_blk_rsm_z = 1;

    cuda_info->n_threads_per_blk_ls_x = 516;
    cuda_info->n_threads_per_blk_ls_y = 1;
    cuda_info->n_threads_per_blk_ls_z = 1;
    cuda_info->n_blk_ls_x = (gpu_snn_in_cpu->n_input_neurons * gpu_snn_in_cpu->n_networks) / cuda_info->n_threads_per_blk_ls_x + 1;
    cuda_info->n_blk_ls_y = 1;
    cuda_info->n_blk_ls_z = 1;

    cuda_info->n_threads_per_blk_nrs_x = 516;
    cuda_info->n_threads_per_blk_nrs_y = 1;
    cuda_info->n_threads_per_blk_nrs_z = 1;
    cuda_info->n_blk_nrs_x = (gpu_snn_in_cpu->n_neurons * gpu_snn_in_cpu->n_networks) / cuda_info->n_threads_per_blk_nrs_x + 1;
    cuda_info->n_blk_nrs_y = 1;
    cuda_info->n_blk_nrs_z = 1;

    cuda_info->n_threads_per_blk_synapses_x = 516;
    cuda_info->n_threads_per_blk_synapses_y = 1;
    cuda_info->n_threads_per_blk_synapses_z = 1;
    cuda_info->n_blk_synapses_x = (gpu_snn_in_cpu->n_synapses * gpu_snn_in_cpu->n_networks) / cuda_info->n_threads_per_blk_synapses_x + 1;
    cuda_info->n_blk_synapses_y = 1;
    cuda_info->n_blk_synapses_z = 1;
}



void free_gpu_snn_in_CPU(GPU_SNN_t *gpu_snn){

    // deallocate internal arrays
    free(gpu_snn->v);
    free(gpu_snn->v_thresh);
    free(gpu_snn->v_rest);
    free(gpu_snn->r_period);
    free(gpu_snn->r_period_remain);
    free(gpu_snn->res);
    free(gpu_snn->t_last_spikes);
    free(gpu_snn->next_last_spike);
    free(gpu_snn->n_neuron_input_synapses);
    free(gpu_snn->neuron_input_synapses_offset);
    
    free(gpu_snn->w);
    free(gpu_snn->init_w);
    free(gpu_snn->dw);
    free(gpu_snn->delay);
    free(gpu_snn->lr);
    free(gpu_snn->pre_neuron_index);
    free(gpu_snn->post_neuron_index);

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

    free(gpu_results->nspk);
    free(gpu_results->gs);

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