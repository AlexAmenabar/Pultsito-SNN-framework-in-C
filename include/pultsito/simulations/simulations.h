/// FUNCTIONS WITH SIMULATION TYPES
#include "snn_library.h"
#include "load_data.h"
#include "helpers.h"
#include "training_rules/stdp.h"

#include "neuron_models/lif_neuron.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

/// @brief 
/// @param snn 
/// @param conf 
/// @param results 
/// @param dataset 
void train_network(spiking_nn_t *snn, simulation_configuration_t *conf, simulation_results_t *results, input_data_t *dataset);

/// @brief 
/// @param snn 
/// @param conf 
/// @param results 
/// @param dataset 
void test_network(spiking_nn_t *snn, simulation_configuration_t *conf, simulation_results_t *results, input_data_t *dataset);

/// @brief 
/// @param snn 
/// @param conf 
/// @param results 
/// @param dataset 
void simulate_samples(spiking_nn_t *snn, simulation_configuration_t *conf, simulation_results_t *results, input_data_t *dataset);

/// @brief Function to simulate an SNN where input is stream data // TODO
void stream_simulation(spiking_nn_t *snn, simulation_configuration_t *conf, simulation_results_t *results);









/// Function to loop over all input synapses of a neuron and calculate the input current
float compute_neuron_input_current(GPU_SNN_t *snn, size_t t, size_t i);
void compute_input_current(GPU_SNN_t *snn, simulation_configuration_t *conf, size_t t);
void compute_LIF_V(GPU_SNN_t *snn, simulation_configuration_t *conf, size_t t);
void process_neuron_firing(GPU_SNN_t *snn, simulation_configuration_t *conf, GPU_results_t *results, size_t t);


///
void simulate_sample_CPU(GPU_SNN_t *snn, GPU_dataset_t *dataset, simulation_configuration_t *conf, GPU_results_t *results, size_t sidx);
void simulate_batch_CPU(GPU_SNN_t *snn, GPU_dataset_t *dataset, simulation_configuration_t *conf, GPU_results_t *results, size_t bidx);
void update_weights_cpu(GPU_SNN_t *snn, size_t batch_size);
void simulate_batches(GPU_SNN_t *snn, GPU_dataset_t *dataset, simulation_configuration_t *conf, GPU_results_t *results);
