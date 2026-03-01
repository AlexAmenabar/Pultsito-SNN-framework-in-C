#ifndef SIMULATIONS_H
#define SIMULATIONS_H

/// Forwarded declarations
typedef struct GPU_SNN_t GPU_SNN_t;
typedef struct GPU_dataset_t GPU_dataset_t;
typedef struct simulation_configuration_t simulation_configuration_t;
typedef struct GPU_results_t GPU_results_t;


/// @brief struct to store the samples in a batch
typedef struct {

    char *spikes;

} tmp_batch_cpu_t;

/* General simulation functions */

/// @brief Function for simulating a dataset in batches
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
/// @param print_data indicates whether data should be printed during simulation
void simulate_batch_CPU(GPU_SNN_t *snn, GPU_dataset_t *dataset, simulation_configuration_t *conf, GPU_results_t *results, size_t bidx, int print_data);


/* Function to compute input currents and neurons firing */

/// @brief Function for computing the input currents
/// @param snn stores the network to be simulated
/// @param conf stores the configuration data for guiding the simulation
/// @param t local time step to simulate
/// @param gt time step of the simulation
void compute_input_current_batch(GPU_SNN_t *snn, simulation_configuration_t *conf, size_t t, size_t gt);

/// @brief Function to simulate neurons firing process
/// @param snn stores the network to be simulated
/// @param conf stores the configuration data for guiding the simulation
/// @param results structure to store the neccessary results
/// @param t local time step to simulate
/// @param gt time step of the simulation
void process_neuron_firing_batch(GPU_SNN_t *snn, simulation_configuration_t *conf, GPU_results_t *results, size_t t, size_t gt);


/* Functions for managing dataset  */

/// @brief Function to initialize the intermediate structure for storing the data of the samples in the simulated batch
/// @param snn stores the network to be simulated
/// @param dataset stores the encoded dataset to simulate
/// @param conf stores the configuration data for guiding the simulation
/// @param results structure to store the neccessary results
/// @param bidx indicates the index of the batch to be simulated in the dataset
tmp_batch_cpu_t* initialize_batch_matrix(GPU_SNN_t *snn, GPU_dataset_t *dataset, simulation_configuration_t *conf, size_t bidx);

/// @brief Function to get the size of the intermediate structure that stores the bitmap of spikes
/// @param N number of neurons
/// @param nS number of synapses
/// @param T simulated time steps
/// @return size in bytes
double get_tmp_batch_size(size_t iN, size_t nS, size_t T);

/// @brief Function to load a time step from the intermediate structure to the spike matrix of the SNN
/// @param snn stores the network to be simulated
/// @param dataset stores the encoded dataset to simulate
/// @param conf stores the configuration data for guiding the simulation
/// @param results structure to store the neccessary results
/// @param bidx indicates the index of the batch to be simulated in the dataset
/// @param t time step to load
void load_batch_time_step_in_SNN_batch(GPU_SNN_t *snn, GPU_dataset_t *dataset, simulation_configuration_t *conf, tmp_batch_cpu_t *batch_data, size_t bidx, size_t t);


/* Functions to reinitialize the network(s) */

/// @brief Function to reinitialize the synapses in the network
/// @param snn stores the network to be simulated
/// @param conf stores the configuration data for guiding the simulation
void reinitialize_synapses_batch(GPU_SNN_t *snn, simulation_configuration_t *conf);

/// @brief Function to reinitialize the spike matrix (set 0 values)
/// @param snn stores the network to be simulated
/// @param conf stores the configuration data for guiding the simulation
void reinitialize_spk_matrix_batch(GPU_SNN_t *snn, simulation_configuration_t *conf);


/* Weights updating */

/// @brief Function to update network weights after simulating a batch
/// @param snn stores the network to be simulated
/// @param conf stores the configuration data for guiding the simulation
void update_weights_cpu(GPU_SNN_t *snn, simulation_configuration_t *conf);

#endif