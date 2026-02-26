#ifndef LOAD_DATA_H
#define LOAD_DATA_H

#include "snn_library.h"

/* [Input management]*/
// Functions to load configuration files, networks, datasets...

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


/* [Output management] */
// Function to store results, networks...


#endif