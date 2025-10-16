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


/// @brief Function to simulate an SNN where all the input is already introduced to the network
/// @param snn SNN structure
/// @param conf Structure that contains information about the simulation configuration
/// @param results Structure to store the results of the simulation
void simulate(spiking_nn_t *snn, simulation_configuration_t *conf, simulation_results_t *results);

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