#ifndef SNN_H
#define SNN_H

typedef struct GPU_SNN_t GPU_SNN_t; // forward declaration
typedef struct simulation_configuration_t simulation_configuration_t; // forward declaration
typedef struct topology_t topology_t; // forward declaration
typedef struct GPU_results_t GPU_results_t; // forward declaration


/// @brief SNN structure
typedef struct GPU_SNN_t {

    // char is used to reduce memory usage
    size_t neuron_type;
    size_t n_neurons;
    size_t n_input_neurons;
    size_t n_output_neurons;
    size_t n_synapses;
    size_t n_networks; // [TODO]: n_networks < batch_size in case batch can not be stored in the GPU
    size_t LT; // matrix L dimension // [DEPRECATED] ???

    // neccessary parameters for all neurons
    float *v; // [n_neurons]: membrane potential of the neuron
    float *v_thresh; // [n_neurons]: threshold of the neuron
    float *v_rest; // [n_neurons]: resting potential of the neuron
    int *r_period; // [n_neurons]: refractary period of the neuron
    int *r_period_remain; // [n_neurons]: remaining refractory period of the neuron
    int *res; // [n_neurons]: resistance of the neuron
    char *post_fired; // [n_neurons]: describes if neuron fired on time t
    float *post_trace; // [n_neurons]: postsynaptic trace
    float *arrI; // [n_neurons]: input current in the time step
    
    size_t *n_neuron_input_synapses; // [n_neurons]: number of input synapses for the neuron
    size_t *neuron_input_synapses_offset; // [n_neurons]: index of the first input synapse for each neuron
    
    // function pointer to generalize neurons management
    void (*allocate_neurons)(GPU_SNN_t *snn, size_t N, simulation_configuration_t *conf);
    void (*init_neurons)(GPU_SNN_t *snn, topology_t *topology, simulation_configuration_t *conf);
    void (*cpy_neurons)(GPU_SNN_t *snn, simulation_configuration_t *conf);
    void (*reinit_neurons)(GPU_SNN_t *snn, simulation_configuration_t *conf); 
    void (*compute_input_current)(GPU_SNN_t *snn, simulation_configuration_t *conf, size_t t, size_t gt);
    void (*compute_dynamics)(GPU_SNN_t *snn, simulation_configuration_t *conf, size_t t, size_t gt);
    void (*compute_firing)(GPU_SNN_t *snn, simulation_configuration_t *conf, GPU_results_t *results, size_t t, size_t gt);
    void (*deallocate_neurons)(GPU_SNN_t *snn);
    
    // synapses
    float *w; // [n_synapses]: weight of the synapse
    float *init_w; // [n_synapses]: initial weight of the synapse
    float *dw; // [n_synapses]: weight difference of the synapse
    float *acc_dw; // [n_synapses]
    int *delay; // [n_synapses]: delay or latency of the synapse
    int *lr; // [n_synapses]: learning rule of the synapse
    size_t *pre_neuron_index; // [n_synapses]: index of the presynaptic neuron
    size_t *post_neuron_index; // [n_synapses]: index of the postsynaptic neuron
    char *pre_fired; // [n_synapses]: describes if the presynaptic neuron fired
    float *pre_trace; // [n_synapses]: presynaptic trace

    char *spk_matrix; // [(n_input_synapses + n_neurons) * t_len * batch_size]

} GPU_SNN_t;


/// @brief Function to initialize the SNN structure. Load from file or generate following the instructions in the configuration file
/// @param conf configuration file with information about the network
/// @return SNN structure initialized
GPU_SNN_t* initialize_network_cpu(simulation_configuration_t *conf);

/// @brief Function to allocate memory for a SNN structure
/// @param snn SNN network structure
void deallocate_snn_str(GPU_SNN_t *snn);

/// @brief Function to copy the required SNN parameters for enabling processing in batches
/// @param snn SNN structure
/// @param conf str with the configuration information
void init_batch_snn(GPU_SNN_t *snn, simulation_configuration_t *conf);

/// @brief Function to get the maximum delay value in the network
/// @param snn SNN network structure
/// @return Maximum delay value
size_t get_max_delay(GPU_SNN_t *snn);

/// @brief Function to get the number of bytes occuped by the SNN strcuture
/// @param snn SNN structure
/// @return size in bytes
double get_snn_size(GPU_SNN_t *snn);

/// @brief Function to get the size of each snn copy
/// @param snn SNN structure
/// @return size in bytes
double get_snn_cpy_size(GPU_SNN_t *snn);

/// @brief Funtiong for printing the SNN information
/// @param snn structure
void print_network(GPU_SNN_t *snn);

/// @brief Funtiong for printing the SNN and its batch_size copies
/// @param snn structure
/// @param conf configuration of the simulation
void print_networks(GPU_SNN_t *snn, simulation_configuration_t *conf);

/// @brief Store the network in a file
/// @param snn Structure that contains the network to be stored
/// @param conf Structure with configuration information
void store_network_in_file(GPU_SNN_t *snn, simulation_configuration_t *conf);


#endif