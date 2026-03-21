#ifndef DATASETS_H
#define DATASETS_H

/// TODO: chunks management
/// TODO: include metadata of the dataset in the dataset file (number of samples...), currently they are indicated in the simulation configuration file

#ifdef __cplusplus
extern "C" {
#endif

/// Forwarded declarations
typedef struct simulation_configuration_t simulation_configuration_t;

/// @brief Structure to store the dataset
typedef struct GPU_dataset_t {

    // [Dataset information]
    int type; // dataset type
    size_t n_classes; // number of classes in the dataset
    size_t n_samples; // number of samples in the dataset
    size_t n_features; // number of features of the samples

    // [Spike times storage]
    size_t *n_spikes_per_feature; // [n_samples * n_features]: number of spikes per each feature 
    size_t *sample_offset; // [n_samples]: offset that indicates where each sample starts in "spikes" array   
    size_t *feature_offset; // [n_samples * n_features]: offset that indicates where each feature starts in "spikes" array
    
    size_t *spikes; // entire dataset (spike times)
    size_t n_spikes; // total number of spikes in the dataset

    // [Spike train frequecy coding]
    size_t *freq; // spike trains described by frequencies

    // [First spike coding]
    size_t *first_spk; // first spike time

} GPU_dataset_t;


/// @brief Load dataset from file
/// @param file_name file name to load data from
/// @param labels_file_name file containing the labels
/// @param n_samples number of samples in the dataset
GPU_dataset_t* load_dataset_from_file_cpu(const char *file_name, const char *labels_file_name, size_t n_samples, simulation_configuration_t *conf);

/// @brief Function to get the size of the dataset in memory
/// @param dataset dataset structure
/// @return size in bytes
double get_dataset_size(GPU_dataset_t *dataset);

/// @brief Funtiong for printing the dataset
/// @param dataset struct storing the dataset
void print_dataset(GPU_dataset_t *dataset);

#ifdef __cplusplus
}
#endif


#endif