#ifndef PRIV_NEURON_MODELS_H
#define PRIV_NEURON_MODELS_H


// Forwared declarations
typedef struct GPU_SNN_t GPU_SNN_t;
typedef struct simulation_configuration_t simulation_configuration_t; 
typedef struct topology_t topology_t;
typedef struct toml_table_t toml_table_t;
typedef struct neurons_t neurons_t;

/* General functions (used by several neuron models) */

/// @brief Function for computing the input currents received by the neurons in time step t
/// @param snn structure that stores the network to be simulated
/// @param conf stores the configuration data for guiding the simulation
/// @param t time step index in the spike matrix
/// @param gt time step of the simulation
void compute_input_current_batch(GPU_SNN_t *snn, simulation_configuration_t *conf, size_t t, size_t gt);

/// @brief Function to simulate neurons firing process
/// @param snn stores the network to be simulated
/// @param conf stores the configuration data for guiding the simulation
/// @param results structure to store the neccessary results
/// @param t local time step to simulate
/// @param gt time step of the simulation
void process_neuron_firing_batch(GPU_SNN_t *snn, simulation_configuration_t *conf, GPU_results_t *results, size_t t, size_t gt);


/* LIF neuron model functions */

/// @brief Function to load the LIF model data from the toml table to the topology structure
/// @param topology Structure to store data
/// @param conf Structure that contains configuration data
/// @param tbl Structure that contains the data read from the toml file
/// @return Structure with the data related to the LIF neurons
neurons_t load_LIF_neurons_from_file(topology_t *topology, simulation_configuration_t *conf, toml_table_t *tbl);

/// @brief Set function pointers for the LIF neuron model simulation
/// @param snn structure
void set_LIF_ptrs(GPU_SNN_t *snn);

/// @brief Allocate memory for  LIF neuron arrays
/// @param snn SNN structure to initialize neurons in
/// @param N Number of neurons in the network
/// @param conf structure with configuration information
void allocate_memory_for_LIF_neurons(GPU_SNN_t *snn, size_t N, simulation_configuration_t *conf);

/// @brief Deallocate memory of LIF neuron arrays
/// @param snn SNN structure to initialize neurons in
void deallocate_memory_for_LIF_neurons(GPU_SNN_t *snn);

/// @brief Initialize LIF neuron arrays
/// @param snn SNN structure to initialize neurons in
/// @param data structure with the values for initializing the neurons
/// @param conf structure with configuration information
/// @return SNN structure
void initialize_LIF_neurons(GPU_SNN_t *snn, topology_t *topology, simulation_configuration_t *conf);

/// @brief Copy each LIF neuron batch size times for enabling batch simulation
/// @param snn SNN structure to expand
/// @param conf Structure with the configuration structure
void init_batch_LIF(GPU_SNN_t *snn, simulation_configuration_t *conf);

/// @brief Reinitialize LIF neurons in network
/// @param snn stores the network to be simulated
/// @param conf stores the configuration data for guiding the simulation
void reinitialize_LIF_neurons_batch(GPU_SNN_t *snn, simulation_configuration_t *conf);

/// @brief Function to compute LIF neuron dynamics
/// @param snn stores the network to be simulated
/// @param conf stores the configuration data for guiding the simulation
/// @param t position to get spikes in the spike matrix
/// @param gt time step of the simulation
void compute_LIF_V_batch(GPU_SNN_t *snn, simulation_configuration_t *conf, size_t t, size_t gt);

#endif