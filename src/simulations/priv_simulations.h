#ifndef PRIV_SIMULATIONS_H
#define PRIV_SIMULATIONS_H

/// Forwarded declarations
typedef struct GPU_SNN_t GPU_SNN_t;
typedef struct GPU_dataset_t GPU_dataset_t;
typedef struct simulation_configuration_t simulation_configuration_t;
typedef struct GPU_results_t GPU_results_t;
typedef struct tmp_batch_cpu_t tmp_batch_cpu_t;


/* Intermediate structure */

/// @brief Function to initialize the intermediate structure for storing the bitmap of the samples in the simulated batch
/// @param snn stores the network to be simulated
/// @param dataset stores the encoded dataset to simulate
/// @param conf stores the configuration data for guiding the simulation
/// @param results structure to store the neccessary results
/// @param bidx indicates the index of the batch to be simulated in the dataset
tmp_batch_cpu_t* initialize_batch_matrix(GPU_SNN_t *snn, GPU_dataset_t *dataset, simulation_configuration_t *conf, size_t bidx);

/// @brief Function to deallocate batch matrix struct
/// @param batchmatrix Structure to deallocate
void deallocate_batch_matrix(tmp_batch_cpu_t* batch_matrix);


/* Simulation */

/// @brief Function to load a time step from the intermediate structure to the spike matrix of the SNN
/// @param snn stores the network to be simulated
/// @param dataset stores the encoded dataset to simulate
/// @param conf stores the configuration data for guiding the simulation
/// @param results structure to store the neccessary results
/// @param bidx indicates the index of the batch to be simulated in the dataset
/// @param t time step to load
void load_batch_time_step_in_SNN_batch(GPU_SNN_t *snn, GPU_dataset_t *dataset, simulation_configuration_t *conf, tmp_batch_cpu_t *batch_data, size_t bidx, size_t t);

/// @brief Function to reinitialize the neurons in the network
/// @param snn Structure that stores the network
/// @param conf Structure that stores the configuration 
void reinitialize_neurons_batch(GPU_SNN_t *snn, simulation_configuration_t *conf);

/// @brief Function to reinitialize the synapses in the network
/// @param snn stores the network to be simulated
/// @param conf stores the configuration data for guiding the simulation
void reinitialize_synapses_batch(GPU_SNN_t *snn, simulation_configuration_t *conf);

/// @brief Function to reinitialize the spike matrix (set 0 values)
/// @param snn stores the network to be simulated
/// @param conf stores the configuration data for guiding the simulation
void reinitialize_spk_matrix_batch(GPU_SNN_t *snn, simulation_configuration_t *conf);


/* Print */

void print_spk_matrix_2D(GPU_SNN_t *snn);

void print_spk_matrix_3D(GPU_SNN_t *snn, size_t B);

void print_weights_batch(GPU_SNN_t *snn, size_t B, size_t b);

void print_weights_3D(GPU_SNN_t *snn, size_t B);

void print_dw_3D(GPU_SNN_t *snn, size_t B);


#endif