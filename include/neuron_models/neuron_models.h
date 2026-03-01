#ifndef NEURON_MODELS_H
#define NEURON_MODELS_H

typedef struct GPU_SNN_t GPU_SNN_t;
typedef struct simulation_configuration_t simulation_configuration_t; 

/* LIF neuron */

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