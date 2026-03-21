#ifndef PRIV_SNN_H
#define PRIV_SNN_H

#ifdef __cplusplus
extern "C" {
#endif

typedef struct GPU_SNN_t GPU_SNN_t; // forward declaration
typedef struct simulation_configuration_t simulation_configuration_t; // forward declaration
typedef struct topology_t topology_t; // forward declaration
typedef struct toml_table_t toml_table_t; // forward declaration


/// @brief Function to allocate memory for a SNN structure
/// @param snn SNN network structure
/// @param conf Structure containing the configuration data of the simulation
GPU_SNN_t* allocate_memory_for_SNN(size_t N, size_t iN, size_t S, size_t mD, simulation_configuration_t *conf);

/// @brief Load the information of the neuron from the input files to a intermediate structure.
/// @param conf Configuration structure with information about the network
/// @return Structure with the arrays storing the data of the neurons, synapses and connectivity to initialize the network
topology_t* load_network_information_in_topology_from_file(simulation_configuration_t *conf);

/// @brief Function to set the pointer functions depending on the neuron type
/// @param snn SNN network structure
/// @param conf Structure containing the configuration data of the simulation
void set_network_neuron_ptr(GPU_SNN_t *snn, simulation_configuration_t *conf);

/// @brief Initialize synaptic arrays
/// @param snn SNN structure to initialize synapses in
/// @param topology structure with the values for initializing the neurons
/// @param conf structure with configuration information
/// @return SNN structure
void initialize_neurons_CPU(GPU_SNN_t *snn, topology_t *topology, simulation_configuration_t *conf);

/// @brief Initialize synaptic arrays
/// @param snn SNN structure to initialize synapses in
/// @param topology structure with the values for initializing the synapses
/// @param conf structure with configuration information
void initialize_synapses_CPU(GPU_SNN_t *snn, topology_t *topology, simulation_configuration_t *conf);

/// @brief Connect the network following an input criteria, where, for each neuron, its input synapses are stored
/// @param snn SNN structure
/// @param data Structure containing helper values
/// @param conf structure containing the configuration of the simulation
void connect_network_input_criteria(GPU_SNN_t *snn, topology_t *topology, simulation_configuration_t *conf);

#ifdef __cplusplus
}
#endif

#endif