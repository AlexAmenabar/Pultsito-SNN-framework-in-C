#ifndef STDP_H
#define STDP_H

/// Forwarded declarations
typedef struct GPU_SNN_t GPU_SNN_t;
typedef struct simulation_configuration_t simulation_configuration_t;

/// @brief Process trace-based STDP for a simulated batch time step
/// @param snn stores the network to be simulated
/// @param conf stores the configuration data for guiding the simulation
void process_trace_based_STDP_batch(GPU_SNN_t *snn, simulation_configuration_t *conf);

#endif