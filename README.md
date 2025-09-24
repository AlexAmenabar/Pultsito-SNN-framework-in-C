# Pultsito

Framework to build, simulate and train SNNs. This framework offers features to simulate SNNs of several biological plausible degrees. This networks can be used either run biological simulations and machine learning tasks.


## Features

The framework offers the following features:

### Neuron models

Actually, only the LIF neuron model is implemented. However, implementation allows to easily integrate more neuron models.

In addition, neurons can be excitatory or inhibitory.

#### Leaky-Integrate-and-Fire

The equations that simulate the LIF neurons are the following:


Additionally, the following parameters are included in the structures to simulate this models:
- Membrane potential.
- Membrane potential threshold.
- Resting potential.
- Neuron resistance.
- Refractory period.

### Synapses

Synapses include the following properties:
- Weight.
- Latency / Delay
- Learning rule.

Including the learning rule as a synaptic parameter allows to incorporate several learning rules in the same network. However, not all learning rules are compatible.

### Learning rules

Actually, 3 STDP variants are implemented:
- Additive-STDP
- Multiplicative-STDP
- Anti-STDP


### Simulation schema

SNNs can be simulated either by clock-based and event-driven approaches. Actually, only clock-based ones are implemented.


### Parallelization strategy

Simulations are accelerated in CPU using OpenMP. However, a GPU implementation is in progress.


### More tools: network generator, input spike trains generator



## Configuration files

Simulations are configured using several configuration files. This file especifies 


### Simulation configuration files

The files to set the configure the simulations has a .toml format. It is organized with the following fields:


| concept     | parameter             | description                                         | values range | type            |
| ------------|-----------------------|-----------------------------------------------------|--------------|-----------------|
| general     | execution_type        | clock-based / event-driven                          | 0 / 1        | int             |
|             | neuron_type           | neuron type for simulations                         | 0            | int             |  
|             | execution_obj         | biological simulations / machine learning           | 0 / 1        | int             |
|             | n_process             | number of parallel threads                          | n > 0        | int             |
|             | learn                 | learning done or not                                | 0 / 1        | int             |
|-------------|-----------------------|-----------------------------------------------------|--------------|-----------------|
| simulation  | time_steps            | time steps for the simulation                       | n > 0        | int             |
|             | input_file            | file path containing input spikes                   | *            | char[]          |
|-------------|-----------------------|-----------------------------------------------------|--------------|-----------------|
| samples     | dataset               | file path containing the dataset to simulate        | *            | char[]          |
|             | dataset_name          | dataset name                                        | *            | char[]          |
|             | num_classes           | number of classes in the dataset                    | n >= 0       | int             |
|             | epochs                | number of epochs to train the network               | n >= 0       | int             |
|-------------|-----------------------|-----------------------------------------------------|--------------|-----------------|
| output      | generated_spikes      | file path to store the spikes generated             | *            | char[]          |
|             | final_weights         | file path to store the weights after train          | *            | char[]          |
|             | execution_times       | file path to store the simulation times             | *            | char[]          |
|             | spikes_per_neuron     | file path to store the n spikes per neuron          | *            | char[]          |
|             | store_network         | not store final network / store                     | 0 / 1        | int             |
|             | store_network_file    | file path to store final network                    | *            | char[]          |
|-------------|-----------------------|-----------------------------------------------------|--------------|-----------------|
| network     | network_file          | file path to load the input network                 | *            | char[]          |
|             | network_neurons_file  | file path to load neurons data                      | *            | char[]          |
|             | network_synapses_file | file path to load synapses data                     | *            | char[]          |
|             | behaviours            | whether neurons behaviours are in network_file      | 0 / 1        | int             |
|             | delays                | whether synapses delays are in network_file         | 0 / 1        | int             |
|             | training_zones        | whehter learning rules are in network_file          | 0 / 1        | int             |
|             | thresh                | whether thresholds are in network_file              | 0 / 1        | int             |
|             | t_refract             | whether refractory periods are in network_file      | 0 / 1        | int             |
|-------------|-----------------------|-----------------------------------------------------|--------------|-----------------|

The input SNN file format is especified in the next section of this documentation.



### SNN file

The files containing the information of the networks is a .toml file. It is organized with the following fields:


| concept  | parameter           | description                                  | values range | type            |
| ---------|---------------------|----------------------------------------------|--------------|-----------------|
| general  | neurons             | number of neurons  [int]                     | n > 0        | int             |
|          | input_neurons       | number of input neurons                      | n > 0        | int             |  
|          | output neurons      | number of output neurons                     | n >= 0       | int             |
|          | synapsis            | number of synapses                           | n > 0        | int             |
|          | input synapsis      | number of input synapses                     | n > 0        | int             |
|          | output_synapsis     | number of output synapses                    | n >= 0       | int             |
|----------|---------------------|----------------------------------------------|--------------|-----------------|
| neurons  | behaviour           | excitatory / inhibitory                      | 1 / 0        | int             |
|          | behaviour_list      | array of behaviours for all neurons          | 1 / 0        | int[neurons]    |
|          | v_thres             | unique threshold for all neurons             | n            | double          |
|          | v_thres_list        | array of thresholds for all neurons          | n            | double[neurons] |
|          | v_rest              | resting potential for all neurons            | n > 0        | double          |
|          | v_rest_list         | resting potentials for each neuron           | n > 0        | double[neurons] |
|          | t_refract           | refratory period for all neurons             | n > 0        | int             |
|          | t_reafract_list     | array of refractory periods                  | n > 0        | int[neurons]    |
|          | input_synapsis      | n. input synapses for each neuron            | n > 0        | int[neurons]    | 
|          | output_synapsis     | n. output synapses for each neuron           | n > 0        | int[neurons]    |
|----------|---------------------|----------------------------------------------|--------------|-----------------|
| synapsis | latency             | latency for all synapses                     | n > 0        | int             |
|          | latency_list        | latencies for each synapse                   | n > 0        | int[synapsis]   |
|          | weights             | weights for each synapse                     | n            | double[synapsis]|
|          | training_zones      | training_zones                               | n > 0        | int             |
|          | training_zones_list | training_zones_list                          | n > 0        | int[synapsis]   |
|          | connections         | indexes of motifs each motif connected to    | n > 0        | int[neurons[]]  |
|----------|---------------------|----------------------------------------------|--------------|-----------------|


The last field (connections) is an array of arrays. The length of the first dimension of arrays is the number of neurons, while each neuron has the length of the number of synapses of the neuron.



## Usage

The source code can be compiled using the compile.sh file with the command ./compile.sh. However, a Makefile will be available soon.

To run the code, it is necessary to provide a configuration file. Some examples of configuration and network files are available on the examples directory.



## Implementation details

This section provides details about the implementation in case it is desired or necessary to expand or modify the features of the library.

### TODO