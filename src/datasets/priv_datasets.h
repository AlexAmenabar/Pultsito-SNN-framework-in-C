#ifndef PRIV_DATASETS_H
#define PRIV_DATASETS_H

/// Forwarded declarations
typedef struct GPU_dataset_t GPU_dataset_t;
typedef struct simulation_configuration_t simulation_configuration_t;

/// @brief Allcoate memory for storing a dataset. If an input parameter is 0, the arrays related to it are not allocated
/// @param n_samples Number of samples in the dataset
/// @param n_features Number of features on each sample
/// @param n_classes Number of classes in the dataset
/// @param n_spikes Number of spikes in the dataset, including all samples and features on the samples
/// @return dataset struct
GPU_dataset_t* allocate_dataset_str(size_t n_samples, size_t n_features, size_t n_classes, size_t n_spikes);

/// @brief Deallocate struct for storing the dataset
/// @param dataset dataset structure to be deallocated
void deallocate_dataset_str(GPU_dataset_t *dataset);

#endif