#ifndef NEURON_MODELS_H
#define NEURON_MODELS_H

typedef struct GPU_SNN_t GPU_SNN_t;
typedef struct simulation_configuration_t simulation_configuration_t; 
typedef struct topology_t topology_t;

/* LIF neuron */

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

/// @brief Copy neurons of the SNN for batch simulation
/// @param snn SNN structure to expand
/// @param conf Structure with the configuration structure
void cpy_LIF(GPU_SNN_t *snn, simulation_configuration_t *conf);

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


/* Hodwking-Huxley */



/* ... */

#endif