#ifndef SIMULATIONS_H
#define SIMULATIONS_H

/// Forwarded declarations
typedef struct GPU_SNN_t GPU_SNN_t;
typedef struct GPU_dataset_t GPU_dataset_t;
typedef struct simulation_configuration_t simulation_configuration_t;
typedef struct GPU_results_t GPU_results_t;


/// @brief structure to store the samples in a batch as a bitmap
typedef struct tmp_batch_cpu_t {

    char *spikes;

} tmp_batch_cpu_t;


/* General simulation functions */

/// @brief Function for simulating a dataset in batches.
/// @param snn stores the network to be simulated
/// @param dataset stores the encoded dataset to simulate
/// @param conf stores the configuration data for guiding the simulation
/// @param results is used to store the necessary results during and after the simulation
void simulate_batches(GPU_SNN_t *snn, GPU_dataset_t *dataset, simulation_configuration_t *conf, GPU_results_t *results);

/// @brief Function for simulating a batch
/// @param snn stores the network to be simulated
/// @param dataset stores the encoded dataset to simulate
/// @param conf stores the configuration data for guiding the simulation
/// @param results is used to store the results of the batch simulation
/// @param bidx indicates the index of the batch to be simulated in the dataset
/// @param print_data indicates whether data should be printed during simulation (execution times...)
void simulate_batch_CPU(GPU_SNN_t *snn, GPU_dataset_t *dataset, simulation_configuration_t *conf, GPU_results_t *results, size_t bidx, int print_data);

/// @brief Function to update network weights after simulating a batch
/// @param snn stores the network to be simulated
/// @param conf stores the configuration data for guiding the simulation
void update_weights_cpu(GPU_SNN_t *snn, simulation_configuration_t *conf);

/// @brief Function to get the size of the intermediate structure that stores the bitmap of the batch spikes
/// @param N number of neurons
/// @param nS number of synapses
/// @param T simulated time steps
/// @return size in bytes
double get_tmp_batch_size(size_t iN, size_t nS, size_t T);



#endif