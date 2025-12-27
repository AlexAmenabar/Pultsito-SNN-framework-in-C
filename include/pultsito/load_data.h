#ifndef LOAD_DATA_H
#define LOAD_DATA_H

#include <stdio.h>
#include "snn_library.h"


/**
* Functions to manage de Input/Output: load configuration files, networks, store results...
*/

/* General */

/// @brief Function to open a file
/// @param f FILE pointer to stream file 
/// @param file_name File name
/// @return Execution code
int open_file(FILE **f, const char *file_name);

/// @brief Function to open a file in write mode without overwriting current contents
/// @param f File to write in
/// @param file_name Name of the file to be opened
/// @return Execution code
int open_file_w(FILE **f, const char *file_name);


/* Input*/

/// @brief Read configuration parameters from configuration TOML file
/// @param file_name Name of the file to load configuration information from
/// @return Struct containing all the configuration data
simulation_configuration_t* load_configuration_params_from_toml(const char *file_name);

/// @brief Load network information from toml file to a struct of arrays
/// @param conf struct containing configuration data
/// @return Struct containing the arrays to initialize the network
network_construction_lists_t* load_network_information_in_lists(simulation_configuration_t *conf);

/// @brief Load dataset
/// @param file_name file name to load data from
/// @param labels_file_name file containing the labels
/// @param n_samples number of samples in the dataset
GPU_dataset_t* load_dataset_from_file_cpu(const char *file_name, const char *labels_file_name, size_t n_samples, simulation_configuration_t *conf);


/* Output */




/* Deprecated */

/// [DEPRECATED]
/// @brief Load network data from file (number of neurons, number of synapses, connections, weights of synapses...)
/// @param file_name File to load information from
/// @param snn Spiking nueral network structure to load information in
/// @param lists Lists used to create the SNN
/// @param conf Structure that contains the simulation configuration
void load_network_information(const char *file_name, spiking_nn_t *snn, network_construction_lists_t *lists, simulation_configuration_t *conf);

// [DEPRECATED]: usable for old SNN structure
/// @brief Load input spike trains inside the SNN structure
/// @param file_name File to load spike trains from
/// @param snn Spiking Neural Network structure to store input spike trains on input synapses
void load_input_spike_trains_on_snn(const char *file_name, spiking_nn_t *snn);

/// [DEPRECATED]
/// @brief Read configuration parameters from configuration file
/// @param file_name Configuration file name
/// @param conf struct to store configuration parameters
int load_configuration_params(const char *file_name, simulation_configuration_t *conf);

/// [DEPREACTED]: used for the old dataset struct
/// @brief Load dataset from file in struct
/// @param file_name 
/// @param dataset 
void load_dataset_from_file(input_data_t *dataset, const char *file_name, const char *labels_file_name, int n_samples, simulation_configuration_t *conf);

/// [DREPRECATED]
void load_network_information(const char *file_name, spiking_nn_t *snn, network_construction_lists_t *lists, simulation_configuration_t *conf);


void open_results_files(simulation_configuration_t *conf);
void close_results_files(simulation_configuration_t *conf);

/// @brief Function that stores the obtained results in files
/// @param results Structure that contains the results of the simulation
/// @param conf Structure to read where to store the results
/// @param snn Spiking neural network structure
void store_results(simulation_results_t *results, simulation_configuration_t *conf, spiking_nn_t *snn, input_data_t *dataset, int *batches, int batch);
void store_generated_spikes(simulation_results_t *results, simulation_configuration_t *conf, spiking_nn_t *snn, input_data_t *dataset, int *batches, int batch);
void store_network(simulation_results_t *results, simulation_configuration_t *conf, spiking_nn_t *snn);
void store_number_of_spikes(simulation_results_t *results, simulation_configuration_t *conf, spiking_nn_t *snn, input_data_t *dataset, int *batches, int batch);
void store_times(simulation_results_t *results, simulation_configuration_t *conf, spiking_nn_t *snn, input_data_t *dataset);

#endif