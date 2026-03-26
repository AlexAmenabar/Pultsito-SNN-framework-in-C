#ifndef SNN_GENERATOR_H
#define SNN_GENERATOR_H

#include <stddef.h>

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

[neurons] (LIF for now)
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
- store_in_file: int (0: no, 1: yes)
- output: char*
- output_neurons: char*
- output_synapses: char*
- output_is_separated: int
- criteria: (0, 1, 2)
*/


#ifdef __cplusplus
extern "C" {
#endif

/// @brief Structure to store all the data relative to the generation configuration
typedef struct generator_conf_t {

    // general
    int neuron_type; // 1: LIF
    int layered; // 0: non-layered, 1: layered
    size_t n_input; // number of input neurons
    size_t n_neurons; // number of neurons
    size_t n_output_neurons; // number of output neurons
    size_t max_pair_neurons_connections; // maximum number of synapses connecting a pair of neurons

    // layered
    size_t n_layers; // number of layers
    size_t *n_neurons_per_layer; // number of neurons per each layer
    int fully_connected; // whether network is fully-connected

    // non-layered
    size_t n_synapses; // number of synapses in the network
    float randomization; // randomization

    // neurons: maximum and minimum values
    float v_thresh_min, v_thresh_max;
    float v_rest_min, v_rest_max;
    float R_min, R_max;
    int rft_per_min, rft_per_max;
    
    // synapses: maximum and minimum values
    float w_min, w_max;
    int delay_min, delay_max;
    int lr_min, lr_max;

    // IO: path for storing the networks
    int store_in_file;

    char *output_file;
    char *output_file_neurons;
    char *output_file_synapses;

    char *output_file_out;
    char *output_file_neurons_out;
    char *output_file_synapses_out;

    int output_is_separated; // whether network should be stored in several files
    int criteria; // storage criteria

} generator_conf_t;

/// @brief Structure to store all the data relative to the neurons
typedef struct neurons_t {

    // neurons properties
    float *v_thresh;
    float *v_rest;
    float *R;
    int *rft_per;


} neurons_t;

/// @brief Structure to store all the data relative to the synapses
typedef struct synapses_t {

    float *w;
    int *lr;
    int *delay;

} synapses_t;

/// @brief Structure to store all the data of the network, including the topology and neuronal and synaptic values
typedef struct topology_t {

    size_t neuron_type, n_neurons, n_output_neurons, n_input;
    size_t n_synapses;

    size_t n_layers;
    size_t *n_neurons_per_layer;

    neurons_t neurons;
    synapses_t synapses;

    // topology
    size_t **input_neurons_per_neuron;
    size_t **output_neurons_per_neuron;

} topology_t;


/// @brief Function to load the configuration file information in the memory
/// @param conf_file path of the configuration file
/// @return Structure with the configuration data
generator_conf_t* read_configuration_file(char* conf_file);

/// @brief Generate the network topology following the instructions in the configuration structure
/// @param conf Structure that describes how to create the topology
/// @return Structure that describes the topology
topology_t generate_topology(generator_conf_t *conf);

/// @brief Function to initialize the network neurons following the instructions in the configuration structure
/// @param conf Structure with the neurons configuration
/// @return structure with the neurons data
neurons_t initialize_neurons(generator_conf_t *conf);

/// @brief Function to initialize the network synapses following the instructions in the configuration structure
/// @param conf Structure with the synapses configuration
/// @return structure with the synapses data
synapses_t initialize_synapses(generator_conf_t *conf);

/// @brief Function to store the generated network in a file
/// @param topology Structure that contains all the information about the generated network
/// @param conf Structure with the configuration data that indicates where to store the network
/// @param criteria Input or output criteria that defines the order of the synapses relative neurons
void store_network(topology_t *topology, generator_conf_t *conf, int criteria);

/// @brief Deallocate memory allocated by the topology structure
/// @param topology Structure to deallocate
void deallocate_topology_str(topology_t *topology);

#ifdef __cplusplus
}
#endif

#endif