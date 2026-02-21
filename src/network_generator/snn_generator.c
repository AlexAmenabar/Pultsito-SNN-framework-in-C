//#include "network_generator/network_generator.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <malloc.h>

#include <toml_c/toml-c.h>



/*
Configuration file format:

[general]
- Layered: int (0/1)
- n_input: int
- n_neurons: int
- n_output_neurons: int
- max_pair_neurons_connections: int

[layered]
- n_layers: int
- n_neurons_per_layer: [int]

- fully_connected: int (0/1)

[no_fully_connected]
- TODO: connectivity?

[non_layered]
- n_synapses: int
- randomization: int

[neurons]
- min_v: float (min v)
- max_v: float (max v)

- min_v_thresh: float (min v)
- max_v_thresh: float (max v)

- min_v_rest: float (min v)
- max_v_rest: float (max v)

- min_R: float (min v)
- max_R: float (max v)

- min_rft_per: float (min v)
- max_rft_per: float (max v)

[synapses]
- min_w: float (min v)
- max_w: float (max v)

- min_lr: int
- max_lr: int

- min_delay: int
- max_delay: int

[IO]
- output: char*
- output_neurons: char*
- output_synapses: char*
- output_is_separated: int
- criteria: (0, 1, 2)
*/

typedef struct {

    // general
    int layered;
    size_t n_input;
    size_t n_neurons;
    size_t n_output_neurons;
    size_t max_pair_neurons_connections;

    // layered
    size_t n_layers;
    size_t *n_neurons_per_layer;
    int fully_connected;

    // non-layered
    size_t n_synapses;
    float randomization;

    // neurons
    float v_min, v_max; 
    float v_thresh_min, v_thresh_max;
    float v_rest_min, v_rest_max;
    int R_min, R_max;
    int rft_per_min, rft_per_max;
    
    // synapses
    float w_min, w_max;
    int delay_min, delay_max;
    int lr_min, lr_max;

    // IO
    char *output_file;
    char *output_file_neurons;
    char *output_file_synapses;

    char *output_file_out;
    char *output_file_neurons_out;
    char *output_file_synapses_out;

    int output_is_separated;
    int criteria;

} generator_conf_t;


// struct to store the values of the neurons and synapses
typedef struct {

    // neurons properties
    float *v;
    float *v_thresh;
    float *v_rest;
    int *R;
    int *rft_per;


} lif_neurons_t;

typedef struct {

    float *w;
    int *lr;
    int *delay;

} synapses_t;

// struct to store topologies
typedef struct {

    size_t n_neurons, n_output_neurons, n_input;
    size_t n_synapses;

    size_t n_layers;
    size_t *n_neurons_per_layer;

    lif_neurons_t neurons;
    synapses_t synapses;

    // topology
    size_t **input_neurons_per_neuron;
    size_t **output_neurons_per_neuron;

} topology_t;


size_t rand_lim(int limit) {
/* return a random number between 0 and limit inclusive.
 */

    int divisor = RAND_MAX/(limit+1);
    int retval;

    do { 
        retval = rand() / divisor;
    } while (retval > limit);

    return (size_t)retval;
}

int rand_int_range(int min, int max){

    // Taking current time as seed
    unsigned int seed = time(0);

    // Generate a random number in the range [min, max]
    int random_number = rand() % (max - min + 1) + min;

    return random_number;
}

float rand_float_range(float min, float max){

    // Taking current time as seed
    unsigned int seed = time(0);

    // Generate a random number in the range [min, max]
    float random_number = (float)rand() / RAND_MAX * max;  
    
    return random_number;
}

void selectionSort(size_t *arr) {
    
    size_t n = arr[0] * 2;

    for (size_t i = 1; i < n; i+=2) {
        
        size_t min = i; // i index minimum
        for (size_t j = i + 2; j < n; j+=2) {
            if (arr[j] < arr[min])
                min = j;
        }

        if (min != i) {
            size_t temp = arr[min];
            size_t temp2 = arr[min + 1];

            arr[min] = arr[i];
            arr[min + 1] = arr[i + 1];

            arr[i] = temp;
            arr[i + 1] = temp2;
        }
    }
}


/* IO functions */

generator_conf_t* load_general_section(generator_conf_t *conf, toml_table_t *tbl){
    
    toml_value_t layered, n_input, n_neurons, n_output_neurons, max_pair_neurons_connections;

    // load section
    layered = toml_table_int(tbl, "layered");
    n_input = toml_table_int(tbl, "n_input");
    n_neurons = toml_table_int(tbl, "n_neurons");
    n_output_neurons = toml_table_int(tbl, "n_output_neurons");
    max_pair_neurons_connections = toml_table_int(tbl, "max_pair_neurons_connections");

    // check provided values and set default values
    if(!layered.ok){
        
        layered.u.i = 0; // non-layered by default
    }
    if(!n_input.ok){
        
        printf(" > Number of input spike trains must be provided! Exiting.\n");
        exit(1);
    }
    if(!n_neurons.ok){
        
        printf(" > Number of neurons must be provided! Exiting.\n");
        exit(1);
    }
    if(!n_output_neurons.ok){

        printf(" > WARNING: number of output neurons not provided. Setting 0.");
        n_output_neurons.u.i = 0;
    }
    if(!max_pair_neurons_connections.ok){

        printf(" > WARNING: number of max pair neurons connections. Setting 1.");
        max_pair_neurons_connections.u.i = 1;
    }   

    // move to configuration struct
    conf->layered = layered.u.i;
    conf->n_input = n_input.u.i;
    conf->n_neurons = n_neurons.u.i;
    conf->n_output_neurons = n_output_neurons.u.i;
    conf->max_pair_neurons_connections = max_pair_neurons_connections.u.i;

    return conf;
}

generator_conf_t* load_layered_section(generator_conf_t *conf, toml_table_t *tbl){
    
    toml_value_t n_layers, fully_connected, n_neurons_in_layer;
    toml_array_t *n_neurons_per_layer;

    // load section
    n_layers = toml_table_int(tbl, "n_layers");
    fully_connected = toml_table_int(tbl, "fully_connected");
    n_neurons_per_layer = toml_table_array(tbl, "n_neurons_per_layer"); // load latencies


    // check provided values and set default values
    if(!n_layers.ok){
                
        printf(" > Number of layers must be provided! Exiting.\n");
        exit(1);    
    }
    if(!fully_connected.ok){
        
        printf(" > WARNING: Setting fully connected.\n");
    }

    // copy values to conf file
    conf->n_layers = n_layers.u.i;
    conf->fully_connected = fully_connected.u.i;

    // allocate memory for n neurons per each layer
    conf->n_neurons_per_layer = (size_t*)malloc(conf->n_layers * sizeof(size_t));

    // load number of neurons per layer
    for(size_t i = 0; i<n_layers.u.i; i++){

        n_neurons_in_layer = toml_array_int(n_neurons_per_layer, i);

        if(!n_neurons_in_layer.ok){

            printf(" > Number of neurons in layer %zu not provided. Exiting.", i);
            exit(1);
        }

        // copy to conf struct
        conf->n_neurons_per_layer[i] = n_neurons_in_layer.u.i;
    }

    return conf;
}


generator_conf_t* load_no_fully_connected_section(generator_conf_t *conf, toml_table_t *tbl){

    // TODO
    toml_value_t connectivity;

    return conf;
}

generator_conf_t* load_no_layered_section(generator_conf_t *conf, toml_table_t *tbl){

    toml_value_t n_synapses, randomization;

    // load section
    n_synapses = toml_table_int(tbl, "n_synapses");
    randomization = toml_table_double(tbl, "randomization");

    // check provided values and set default values
    if(!n_synapses.ok){
                
        printf(" > Number of synapses must be provided! Exiting.\n");
        exit(1);    
    }
    if(!randomization.ok){
        
        randomization.u.d = 0.0; 
    }

    // copy values to conf file
    conf->n_synapses = n_synapses.u.i;
    conf->randomization = (float)randomization.u.d;

    return conf;
}

generator_conf_t* load_neurons_section(generator_conf_t *conf, toml_table_t *tbl){

    toml_value_t v_min, v_max, v_thresh_min, v_thresh_max, v_rest_min, v_rest_max, 
                 R_min, R_max, rft_per_min, rft_per_max;
    
    // load section
    v_max = toml_table_double(tbl, "v_max");
    v_min = toml_table_double(tbl, "v_min");
    v_thresh_max = toml_table_double(tbl, "v_thresh_max");
    v_thresh_min = toml_table_double(tbl, "v_thresh_min");
    v_rest_max = toml_table_double(tbl, "v_rest_max");
    v_rest_min = toml_table_double(tbl, "v_rest_min");
    R_max = toml_table_int(tbl, "R_max");
    R_min = toml_table_int(tbl, "R_min");
    rft_per_max = toml_table_int(tbl, "rft_per_max");
    rft_per_min = toml_table_int(tbl, "rft_per_min");
    
    // check provided values and set default values
    if(!v_max.ok){
                
        printf(" > v max not provided! Exiting.\n");
        exit(1);    
    }
    if(!v_min.ok){
            
        printf(" > v min not provided! Exiting.\n");
        exit(1);    
    }
    if(!v_thresh_max.ok){
                
        printf(" > v_thresh max not provided! Exiting.\n");
        exit(1);    
    }
    if(!v_thresh_min.ok){
                
        printf(" > v_thresh min not provided! Exiting.\n");
        exit(1);    
    }
    if(!v_rest_max.ok){
                
        printf(" > v_rest max not provided! Exiting.\n");
        exit(1);    
    }
    if(!v_rest_min.ok){
                
        printf(" > v_rest min not provided! Exiting.\n");
        exit(1);    
    }
    if(!R_max.ok){
                
        printf(" > R max not provided! Exiting.\n");
        exit(1);    
    }
    if(!R_min.ok){
                
        printf(" > R min not provided! Exiting.\n");
        exit(1);    
    }
    if(!rft_per_max.ok){
                
        printf(" > Refractary period max not provided! Exiting.\n");
        exit(1);    
    }
    if(!rft_per_min.ok){
                
        printf(" > Refractary period min not provided! Exiting.\n");
        exit(1);    
    }

    // copy values to conf file
    conf->v_max = (float)v_max.u.d;
    conf->v_min = (float)v_min.u.d;
    conf->v_thresh_max = (float)v_thresh_max.u.d;
    conf->v_thresh_min = (float)v_thresh_min.u.d;
    conf->v_rest_max = (float)v_rest_max.u.d;
    conf->v_rest_min = (float)v_rest_min.u.d;
    conf->R_max = R_max.u.i;
    conf->R_min = R_min.u.i;
    conf->rft_per_max = rft_per_max.u.i;
    conf->rft_per_min = rft_per_min.u.i;
    
    return conf;
}

generator_conf_t* load_synapses_section(generator_conf_t *conf, toml_table_t *tbl){

    toml_value_t w_min, w_max, delay_min, delay_max, lr_min, lr_max;
    
    // load section
    w_max = toml_table_double(tbl, "w_max");
    w_min = toml_table_double(tbl, "w_min");
    delay_max = toml_table_int(tbl, "delay_max");
    delay_min = toml_table_int(tbl, "delay_min");
    lr_max = toml_table_int(tbl, "lr_max");
    lr_min = toml_table_int(tbl, "lr_min");
    
    // check provided values and set default values
    if(!w_max.ok){
                
        printf(" > w max not provided! Exiting.\n");
        exit(1);    
    }
    if(!w_min.ok){
            
        printf(" > w min not provided! Exiting.\n");
        exit(1);    
    }
    if(!delay_max.ok){
                
        printf(" > delay max not provided! Setting 1.\n");
        delay_max.u.i = 1;
    }
    if(!delay_min.ok){
                
        printf(" > delay min not provided! Setting 1.\n");
        delay_min.u.i = 1;
    }
    if(!lr_max.ok){
                
        printf(" > lr max not provided! Setting 0.\n");
        lr_min.u.i = 0;
    }
    if(!lr_min.ok){
                
        printf(" > lr min not provided! Setting 0.\n");
        lr_min.u.i = 0;    
    }

    // copy values to conf file
    conf->w_max = (float)w_max.u.d;
    conf->w_min = (float)w_min.u.d;
    conf->delay_max = delay_max.u.i;
    conf->delay_min = delay_min.u.i;
    conf->lr_max = lr_max.u.i;
    conf->lr_min = lr_min.u.i;

    return conf;
}

generator_conf_t* load_IO_section(generator_conf_t *conf, toml_table_t *tbl){

    toml_value_t output_file, output_file_neurons, output_file_synapses, output_file_out, output_file_neurons_out, output_file_synapses_out,
                output_is_separated, criteria;
    
    // load section
    output_file = toml_table_string(tbl, "output_file");
    output_file_neurons = toml_table_string(tbl, "output_file_neurons");
    output_file_synapses = toml_table_string(tbl, "output_file_synapses");

    output_file_out = toml_table_string(tbl, "output_file_out");
    output_file_neurons_out = toml_table_string(tbl, "output_file_neurons_out");
    output_file_synapses_out = toml_table_string(tbl, "output_file_synapses_out");

    output_is_separated = toml_table_int(tbl, "output_is_separated");
    criteria = toml_table_int(tbl, "criteria");

    // check provided values and set default values
    if(!output_file.ok){
                
        printf(" > Output file must be provided! Exiting.\n");
        exit(1);    
    }
    if(!output_is_separated.ok){
                
        printf(" > WARNING: Setting output not separated.\n");
        output_is_separated.u.i = 0;
    }

    // copy values to conf file
    conf->output_file = output_file.u.s;
    conf->output_file_out = output_file_out.u.s;

    conf->output_is_separated = output_is_separated.u.i;

    if(conf->output_is_separated == 1 && !output_file_neurons.ok){
                
        printf(" > Output file for neurons must be provided! Exiting.\n");
        exit(1);    
    }
    if(conf->output_is_separated == 1 && !output_file_synapses.ok){
                
        printf(" > Output file for synapses must be provided! Exiting.\n");
        exit(1);    
    }

    if(conf->output_is_separated == 1 && !output_file_neurons_out.ok){
                
        printf(" > Output file for neurons must be provided! Exiting.\n");
        exit(1);    
    }
    if(conf->output_is_separated == 1 && !output_file_synapses_out.ok){
                
        printf(" > Output file for synapses must be provided! Exiting.\n");
        exit(1);    
    }

    if(!criteria.ok){
             
        criteria.u.i=0;
        printf(" > Criteria not provided. Default value 0.\n");
    }


    conf->output_file_neurons = output_file_neurons.u.s;
    conf->output_file_synapses = output_file_synapses.u.s;

    conf->output_file_neurons_out = output_file_neurons_out.u.s;
    conf->output_file_synapses_out = output_file_synapses_out.u.s;

    conf->criteria = criteria.u.i;

    return conf;
}

generator_conf_t* read_configuration_file(char* conf_file){

    generator_conf_t *conf = (generator_conf_t*)malloc(sizeof(generator_conf_t));

    // open configuration file
    FILE *f = fopen(conf_file, "r");;
    if (f == NULL){
        perror("Error opening the file!\n");
        exit(1);
    }    

    // define TOML library variables
    toml_table_t *tbl, *tbl_general, *tbl_layered, *tbl_no_fully_connected, *tbl_no_layered, *tbl_neurons, *tbl_synapses, *tbl_IO;
    

    toml_value_t configuration_file, output_file;


    // load toml file
    char errbuf[100];
    tbl = toml_parse_file(f, errbuf, 100);
    // close file
    fclose(f);

    // get sections from toml file
    tbl_general = toml_table_table(tbl, "general");
    tbl_layered = toml_table_table(tbl, "layered");
    tbl_no_fully_connected = toml_table_table(tbl, "no_fully_connected");
    tbl_no_layered = toml_table_table(tbl, "no_layered");
    tbl_neurons = toml_table_table(tbl, "neurons");
    tbl_synapses = toml_table_table(tbl, "synapses");
    tbl_IO = toml_table_table(tbl, "IO");

    // load general section
    conf = load_general_section(conf, tbl_general);

    if(conf->layered == 1){
        conf = load_layered_section(conf, tbl_layered);
        conf = load_no_fully_connected_section(conf, tbl_no_fully_connected);
    }
    else{
        conf = load_no_layered_section(conf, tbl_no_layered);
    }
    conf = load_neurons_section(conf, tbl_neurons);
    conf = load_synapses_section(conf, tbl_synapses);
    conf = load_IO_section(conf, tbl_IO);

    // return configuration struct
    return conf;
}


/* Function for generating networks */

topology_t generate_layered_topology(generator_conf_t *conf){

    size_t n_neurons, n_input, n_output_neurons, n_synapses, n_layers;
    size_t *n_neurons_per_layer;
    size_t **input_neurons_per_neuron;

    n_neurons = conf->n_neurons;
    n_input = conf->n_input;
    n_output_neurons = conf->n_output_neurons;
    n_layers = conf->n_layers;
    n_neurons_per_layer = conf->n_neurons_per_layer;
    n_synapses = 0;


    // TODO
    // n neurons is the sum of all neuron in all layers, and n_unput
    // should be the neurons in the first layer? The last is incorrect
    //n_input = n_neurons_per_layer[0];

    input_neurons_per_neuron = (size_t **)malloc(n_neurons * sizeof(size_t*));

    // process layer by layer: input neurons of layer[i] neurons are layer[i-1]
    
    // input layer: TEMPORAL: 1:1
    size_t next_neuron = 0;
    /*for(size_t i = 0; i<n_neurons_per_layer[0]; i++){

        input_neurons_per_neuron[next_neuron] = (size_t*)malloc((1 * 2 + 1) * sizeof(size_t));
        input_neurons_per_neuron[next_neuron][0] = 1;
        input_neurons_per_neuron[next_neuron][1] = i;
        input_neurons_per_neuron[next_neuron][2] = 1;

        n_synapses ++;
        next_neuron ++;
    }*/
    
    // all input neurons to all input spike trains
    for(size_t i = 0; i<n_neurons_per_layer[0]; i++){

        input_neurons_per_neuron[next_neuron] = (size_t*)malloc((n_input * 2 + 1) * sizeof(size_t));
        input_neurons_per_neuron[next_neuron][0] = n_input; // each neuron connected to all input spike trains

        for(size_t j = 0; j<n_input; j++){
        
            input_neurons_per_neuron[next_neuron][j * 2 + 1] = j;
            input_neurons_per_neuron[next_neuron][j * 2 + 2] = 1;
            n_synapses ++;
        }

        next_neuron ++;
    }

    for(size_t i = 1; i<n_layers; i++){

        for(size_t j = 0; j<n_neurons_per_layer[i]; j++){

            input_neurons_per_neuron[next_neuron] = (size_t*)malloc((n_neurons_per_layer[i-1] * 2 + 1) * sizeof(size_t));
            input_neurons_per_neuron[next_neuron][0] = n_neurons_per_layer[i-1];
        
            for(size_t l = 0; l<n_neurons_per_layer[i-1]; l++){
                input_neurons_per_neuron[next_neuron][l * 2 + 1] = next_neuron - j - n_neurons_per_layer[i-1] + l + n_input;
                input_neurons_per_neuron[next_neuron][l * 2 + 2] = 1;
            }

            n_synapses += n_neurons_per_layer[i-1];
            next_neuron ++;
        }
    }

    // initialize and return topology
    topology_t topology;
    topology.n_neurons = n_neurons;
    topology.n_input = n_input;
    topology.n_output_neurons = n_output_neurons;
    topology.n_synapses = n_synapses;

    topology.input_neurons_per_neuron = input_neurons_per_neuron;

    topology.n_layers = n_layers;
    topology.n_neurons_per_layer = n_neurons_per_layer;     

    return topology;
}

topology_t generate_non_layered_topology(generator_conf_t *conf){

    size_t n_neurons, n_input, n_output_neurons, n_synapses;
    size_t n_mean_input_connections, *n_input_connections_per_neuron, **input_neurons_per_neuron;

    float randomization;

    size_t i, j;

    n_neurons = conf->n_neurons;
    n_input = conf->n_input;
    n_output_neurons = conf->n_output_neurons;
    n_synapses = conf->n_synapses;
    randomization = conf->randomization;

    // compute mean number of connections per each neuron
    n_mean_input_connections = n_synapses / n_neurons;
    input_neurons_per_neuron = (size_t**)malloc(n_neurons * sizeof(size_t*));
    n_input_connections_per_neuron = (size_t*)malloc(n_neurons * sizeof(size_t));
    
    // set number of connections per each neuron using randomization
    for(i=0; i<n_neurons; i++){

        n_input_connections_per_neuron[i] = n_mean_input_connections;

        // randomize
        size_t n = (size_t)(rand() % (size_t)((float)n_neurons * randomization)) / 100;
        size_t dec = (size_t)(rand() % 2);

        if(dec == 0){
            
            if(n > n_input_connections_per_neuron[i])
                n = n_input_connections_per_neuron[i];

            n_input_connections_per_neuron[i] -= n;
        }
        else if(dec == 1){
            
            n_input_connections_per_neuron[i] += n;
        }
    }

    // manage connection with input spike trains (TEMPORAL, 1:1 connections)
    for(i = 0; i<n_input; i++){

        n_input_connections_per_neuron[i] += 1;
    }

    // generate topology
    size_t generated_n_synapses = 0; 
    for(i = 0; i<n_neurons; i++){

        // allocate memory to store input neurons
        input_neurons_per_neuron[i] = (size_t*)malloc((n_input_connections_per_neuron[i] * 2 + 1) * sizeof(size_t));
        input_neurons_per_neuron[i][0] = 0;

        size_t remaining_input_connections = n_input_connections_per_neuron[i]; 
        
        j = 0;

        // TEMPORAL: check if neuron is connected to an input neuron (at least one)
        if(i < n_input){

            // add one connection to input neuron i
            input_neurons_per_neuron[i][0] ++;
            input_neurons_per_neuron[i][1] = i;
            input_neurons_per_neuron[i][2] = 1;
            generated_n_synapses += 1;

            // update remaining input connections
            remaining_input_connections --;
            j++;
        }
        
        while(remaining_input_connections > 0){

            // select a neuron
            size_t in_neuron = rand_lim(n_neurons + n_input - 1); // [0, n_neurons + n_input]
            
            size_t n = rand_lim(conf->max_pair_neurons_connections); // number of connections to that neuron
            if(n == 0) 
                n = 1;

            // check if the neuron has already been selected previously
            int valid = 1;
            for(size_t l = 0; l<j; l++){

                if(input_neurons_per_neuron[i][l * 2 +1] == in_neuron){

                    // already in the array, update number of connections
                    valid = 0;
                    input_neurons_per_neuron[i][l*2+2] += n;
                }
            }

            // neuron not encountered in the array, so add it
            if(valid == 1){

                input_neurons_per_neuron[i][j * 2 + 1] = in_neuron;
                input_neurons_per_neuron[i][j * 2 + 2] = n;
                
                // new input neuron added
                input_neurons_per_neuron[i][0] ++;

                // update j
                j++;
            }   

            // update number of synapses generated in the network
            generated_n_synapses += n;

            // reduce remaining connection to be done
            remaining_input_connections -= n;
        }

        // sort array
        selectionSort(input_neurons_per_neuron[i]);
    }

    // free memory
    free(n_input_connections_per_neuron);

    
    // initialize and return topology
    topology_t topology;
    topology.n_neurons = n_neurons;
    topology.n_input = n_input;
    topology.n_output_neurons = n_output_neurons;
    topology.n_synapses = generated_n_synapses;

    topology.input_neurons_per_neuron = input_neurons_per_neuron;

    topology.n_layers = 1;
    topology.n_neurons_per_layer = 0; // no used 


    // return topology
    return topology;
}

topology_t generate_topology(generator_conf_t *conf){

    topology_t topology;

    if(conf->layered == 1){
        topology = generate_layered_topology(conf);
    }
    else{
        topology = generate_non_layered_topology(conf);
    }

    // update number of synapses in configuration
    conf->n_synapses = topology.n_synapses;

    return topology;
}


lif_neurons_t initialize_neurons(generator_conf_t *conf){

    lif_neurons_t neurons;

    // allocate memory for neurons properties
    neurons.v        = (float*)malloc(conf->n_neurons * sizeof(float));
    neurons.v_thresh = (float*)malloc(conf->n_neurons * sizeof(float));
    neurons.v_rest   = (float*)malloc(conf->n_neurons * sizeof(float));
    neurons.R        = (int*)malloc(conf->n_neurons * sizeof(int));
    neurons.rft_per  = (int*)malloc(conf->n_neurons * sizeof(int));

    // initialize values
    for(size_t i = 0; i<conf->n_neurons; i++){

        neurons.v[i] = rand_float_range(conf->v_min, conf->v_max);
        neurons.v_thresh[i] = rand_float_range(conf->v_thresh_min, conf->v_thresh_max);
        neurons.v_rest[i] = rand_float_range(conf->v_rest_min, conf->v_rest_max);
        neurons.R[i] = rand_int_range(conf->R_min, conf->R_max);
        neurons.rft_per[i] = rand_int_range(conf->rft_per_min, conf->rft_per_max);
    }
    
    // return neurons
    return neurons;
}

synapses_t initialize_synapses(generator_conf_t *conf){

    synapses_t synapses;

    printf(" n_synapses = %zu\n", conf->n_synapses);
    // allocate memory for synapses properties
    synapses.w     = (float*)malloc(conf->n_synapses * sizeof(float));
    synapses.delay = (int*)malloc(conf->n_synapses * sizeof(int));
    synapses.lr    = (int*)malloc(conf->n_synapses * sizeof(int));

    // initialize values
    for(size_t i = 0; i<conf->n_synapses; i++){

        synapses.w[i] = rand_float_range(conf->w_min, conf->w_max);
        synapses.delay[i] = rand_int_range(conf->delay_min, conf->delay_max);
        synapses.lr[i] = rand_int_range(conf->lr_min, conf->lr_max);
    }

    // return synapses
    return synapses;
}



void map_network_to_output_criteria(topology_t *topology, generator_conf_t *conf){

    size_t i, j, k;


    // store old values
    float *tmp_w = topology->synapses.w;
    int *tmp_lr = topology->synapses.lr;
    int *tmp_delay = topology->synapses.delay;

    // allocate memory for reordered arrays
    topology->synapses.w = (float*)malloc(topology->n_synapses * sizeof(float));
    topology->synapses.lr = (int*)malloc(topology->n_synapses * sizeof(int));
    topology->synapses.delay = (int*)malloc(topology->n_synapses * sizeof(int));

    // count
    size_t n_total_neurons = topology->n_neurons + topology->n_input; 

    size_t *n_output_neurons_per_neuron = (size_t*)CALLOC(n_total_neurons, sizeof(size_t));


    for(i = 0; i<topology->n_neurons; i++){

        for(j = 0; j<topology->input_neurons_per_neuron[i][0]; j++){

            size_t in_neuron = topology->input_neurons_per_neuron[i][j * 2 + 1];
            n_output_neurons_per_neuron[in_neuron] += 1;//topology->input_neurons_per_neuron[i][j * 2 + 2];
        }
    }

    // allocate
    topology->output_neurons_per_neuron = (size_t**)CALLOC(n_total_neurons, sizeof(size_t*));
    for(i = 0; i<n_total_neurons; i++){

        topology->output_neurons_per_neuron[i] = (size_t*)malloc((n_output_neurons_per_neuron[i] * 2 + 1) * sizeof(size_t));
        topology->output_neurons_per_neuron[i][0] = n_output_neurons_per_neuron[i];
    }

    // map connections
    size_t *next = (size_t*)CALLOC(n_total_neurons, sizeof(size_t));
    size_t *n_output_synapses = (size_t*)CALLOC(n_total_neurons, sizeof(size_t));

    for(i = 0; i<topology->n_neurons; i++){

        size_t out = i;

        for(j = 0; j<topology->input_neurons_per_neuron[i][0]; j++){

            size_t in = topology->input_neurons_per_neuron[i][j * 2 + 1];
            size_t n = topology->input_neurons_per_neuron[i][j * 2 + 2];
            topology->output_neurons_per_neuron[in][next[in] * 2 + 1] = out + topology->n_input;
            topology->output_neurons_per_neuron[in][next[in] * 2 + 2] = n;
            next[in] ++;

            n_output_synapses[in] += n;
        }
    }

    // compute offsets
    size_t *synapses_offset = (size_t*)CALLOC(n_total_neurons, sizeof(size_t));
    size_t tmp_offset = 0;
    for(i = 0; i<n_total_neurons; i++){

        synapses_offset[i] = tmp_offset;
        tmp_offset += n_output_synapses[i];
    }

    // map synapses properties
    size_t *neuron_synapse_offset = (size_t*)CALLOC(n_total_neurons, sizeof(size_t));
    size_t synapse_index = 0;
    for(i=0; i<topology->n_neurons; i++){

        size_t out = i;

        for(j = 0; j<topology->input_neurons_per_neuron[i][0]; j++){

            size_t in = topology->input_neurons_per_neuron[i][j * 2 + 1];
            size_t n = topology->input_neurons_per_neuron[i][j * 2 + 2];

            for(k = 0; k<n; k++){

                topology->synapses.w[synapses_offset[i] + neuron_synapse_offset[i]] = tmp_w[synapse_index];
                topology->synapses.lr[synapses_offset[i] + neuron_synapse_offset[i]] = tmp_lr[synapse_index];
                topology->synapses.delay[synapses_offset[i] + neuron_synapse_offset[i]] = tmp_delay[synapse_index];

                neuron_synapse_offset[i] += n;
                synapse_index ++;
            }
        }
    }
}



void store_network(topology_t *topology, generator_conf_t *conf, int criteria){

    // open files to store the network
    FILE *f, *f_out;

    f = fopen(conf->output_file, "w"); // write
    if (f == NULL){
        perror("Error opening the file\n");
        exit(1);
    }    


    FILE *f_neurons, *f_neurons_out;
    FILE *f_synapses, *f_synapses_out;

    if(conf->output_is_separated == 1){
    
        f_neurons = fopen(conf->output_file_neurons, "w"); // write
        if (f_neurons == NULL){
            perror("Error opening the file\n");
            exit(1);
        }    

        if(criteria == 0){
        
            f_synapses = fopen(conf->output_file_synapses, "w"); // write
            if (f_synapses == NULL){
                perror("Error opening the file\n");
                exit(1);
            }    
        }

        if(criteria == 1){
        
            f_synapses_out = fopen(conf->output_file_synapses_out, "w"); // write
            if (f_synapses_out == NULL){
                perror("Error opening the file\n");
                exit(1);
            }    
        }

    }
    
    // write general information
    printf(" > Storing general information...\n");

    fprintf(f, "[general]\n");
    fprintf(f, "    neurons = %zu\n", topology->n_neurons);
    fprintf(f, "    input_neurons = %zu\n", topology->n_input);
    fprintf(f, "    output_neurons = %zu\n", topology->n_output_neurons);
    fprintf(f, "    synapsis = %zu\n", topology->n_synapses);
    fprintf(f, "    input_synapsis = %zu\n", 0); // DEPRECATED
    fprintf(f, "    output_synapsis = %zu\n", 0);
    fprintf(f, "    network_is_separated = %d\n\n", conf->output_is_separated);

    printf(" > General information stored!\n");

    // write neurons information
    printf(" > Storing neurons information...\n");

    fprintf(f, "[neurons]\n");
    fprintf(f, "    behaviour = 1\n"); // TODO: temporal, var deprecated
    fprintf(f, "\n");

    // if is separated == 1, store neurons information as plain text
    //if(output_is_separated == 1){
    size_t i;
    for(i=0; i<topology->n_neurons; i++){ // excitatory
        fprintf(f_neurons, "%d ", 1);
    }
    fprintf(f_neurons, "\n");

    for(i=0; i<topology->n_neurons; i++){
        fprintf(f_neurons, "%lf ", topology->neurons.v_thresh[i]);
    }
    fprintf(f_neurons, "\n");

    for(i=0; i<topology->n_neurons; i++){
        fprintf(f_neurons, "%lf ", topology->neurons.v_rest[i]);
    }
    fprintf(f_neurons, "\n");

    for(i=0; i<topology->n_neurons; i++){
        fprintf(f_neurons, "%d ", topology->neurons.rft_per[i]);
    }
    fprintf(f_neurons, "\n");


    printf(" > Neurons information stored!\n");

    // write synapses information
    printf(" > Storing synapses information...\n");
    
    // if is separated == 1, store neurons information in plain text
    //if(output_is_separated == 1){
    for(i=0; i<topology->n_synapses; i++){
        
        if(criteria == 0)
            fprintf(f_synapses, "%d ", topology->synapses.delay[i]);

        if(criteria == 1)
            fprintf(f_synapses_out, "%d ", topology->synapses.delay[i]);
    }
    fprintf(f_synapses, "\n");

    for(i=0; i<topology->n_synapses; i++){

        if(criteria == 0)
            fprintf(f_synapses, "%lf ", topology->synapses.w[i]);
        
        if(criteria == 1)
            fprintf(f_synapses_out, "%lf ", topology->synapses.w[i]);
    }
    fprintf(f_synapses, "\n");

    for(i=0; i<topology->n_synapses; i++){
        if(criteria == 0)
            fprintf(f_synapses, "%d ", topology->synapses.lr[i]);

        if(criteria == 1)
            fprintf(f_synapses_out, "%d ", topology->synapses.lr[i]);
    }
    fprintf(f_synapses, "\n\n");

    // store connectivity


    printf(" Storing connectivity\n");
    if(criteria == 0){
        
        for(size_t i = 0; i<topology->n_neurons; i++){

            fprintf(f_synapses, "%zu ", topology->input_neurons_per_neuron[i][0]);
            
            for(size_t j = 0; j<topology->input_neurons_per_neuron[i][0]; j++){

                fprintf(f_synapses, "%zu ", topology->input_neurons_per_neuron[i][j * 2 + 1]);
                fprintf(f_synapses, "%zu ", topology->input_neurons_per_neuron[i][j * 2 + 2]);
            }

            fprintf(f_synapses,"\n");

        }
    }
    if(criteria == 1){
        
        for(size_t i = 0; i<topology->n_neurons + topology->n_input; i++){

            fprintf(f_synapses_out, "%zu ", topology->output_neurons_per_neuron[i][0]);
            
            for(size_t j = 0; j<topology->output_neurons_per_neuron[i][0]; j++){

                fprintf(f_synapses_out, "%zu ", topology->output_neurons_per_neuron[i][j * 2 + 1]);
                fprintf(f_synapses_out, "%zu ", topology->output_neurons_per_neuron[i][j * 2 + 2]);
            }

            fprintf(f_synapses_out,"\n");

        }
    }
    fclose(f);
    fclose(f_neurons);
    fclose(f_synapses);
}




int main(int argc, char *argv[]) {
    
    // randomize exection to create different networks each time
    srand(time(NULL));

    // load configuration file
    generator_conf_t *conf = read_configuration_file(argv[1]);

    // initialize topology
    topology_t topology = generate_topology(conf);

    // initialize neurons and synapses
    topology.neurons = initialize_neurons(conf);
    topology.synapses = initialize_synapses(conf);

    // store generated network
    store_network(&topology, conf, 0);
/*
    // map to output criteria
    map_network_to_output_criteria(&topology, conf);
    
    // store generated network
    store_network(&topology, conf, 1);*/
}