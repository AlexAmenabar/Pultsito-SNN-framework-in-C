#ifndef METRICS_H
#define METRICS_H

#include "snn_library.h"


// shuffle indexes
void shuffle_sample_indexes(int *indexes, int n);

// creates batches from the dataset
int* create_batches(int *n_batches, int *extra_samples, input_data_t *dataset, int n_samples, int batch_size);


float accuracy_by_output_neurons(simulation_results_t *results, input_data_t *dataset, int *batches, int batch_size, int batch, int n_out);


#endif