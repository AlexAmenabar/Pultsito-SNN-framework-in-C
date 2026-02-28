# ARCEUS

Framework to build, simulate and train Spiking Neural Networks (SNNs). This framework offers features to simulate SNNs of several biological plausible degrees, especially focused on Machine Learning oriented simulations, providing features to efficiently simulate datasets. 


## Features

The framework offers the following features:

### Neuron models

Actually, only the LIF neuron model [ref] is implemented. However, implementation allows to easily integrate new neuron models without changing the core code.


#### Leaky-Integrate-and-Fire (LIF)

The equations that simulate the LIF neurons are the following:



Additionally, the following parameters are included in the structures to simulate this models:
- Membrane potential (V).
- Membrane potential threshold (V_{thresh}).
- Resting potential (V_{rest}).
- Neuron resistance (R).
- Refractory period (rftp).

### Synapses

Synaptic connections (connections between neurons) include the following properties:
- Weight.
- Delay
- Learning rule.

The learning rule allows to combine different learning paradigms in the same network.

### Learning rules

Currently only the trace-based STDP learning rule is included [ref]


### Simulation schema

SNNs can be simulated either following clock-based and event-driven approaches. This framework incorporates the clock-based approach, where in all time steps all neurons and synapses are processed. This approach allows to efficiently simulate samples in batches, which is very common in Machine Learning simulations.


### CPU simulations

CPU simulations are accelerated combining both parallelization and vectorization. The former is achieved using OpenMP, while the latter takes advantage of AVX512.

### GPU simulations

The simulation loop is also implemented in the GPU, allowing both single-GPU and multi-GPU simulations using CUDA. 


## Configuration files

Simulations follow the instructions provided by a configuration file, which inlcudes the parameters indicated in the following table. The configuration file is on .TOML format. 


| Section     | Parameter             | Description                                         | Type         |
| ------------|-----------------------|-----------------------------------------------------|--------------|
| general     | neuron_type           | neuron type for simulations                         | int (n > 0)  |
|             | n_process             | number of parallel threads                          | int (n > 0)  |
|             | learn                 | inference or learning                               | 0 / 1        | int             |
|             | cuda                  | cuda not used or used                               | 0 / 1        | int             |
|             | multi_gpu             | number of GPUs: all / 1 / n                         | n >= 0       | int             |
|             | batch_size            | batch size                                          | n > 0        | int             |
|             | thrN                  | TODO                                                | n > 0        | int             |
|             | load_network          | TODO                                                | TODO         | TODO            |
|             | load_dataset          | TODO                                                | TODO         | TODO            |
|-------------|-----------------------|-----------------------------------------------------|--------------|-----------------|
| simulation  | time_steps            | time steps for simulating each sample               | n > 0        | int             |
|             | max_input_spikes      | maximum number spikes on each input spike train     | n > 0        | int             |
|-------------|-----------------------|-----------------------------------------------------|--------------|-----------------|
| dataset     | train_provided        | train set not provided or provided                  | 0 / 1        | int             |
|             | train_set             | path to the file with the train set                 | String       | char[]          |
|             | train_labels          | number of classes in the dataset                    | n >= 0       | int             |
|             | n_train               | number of epochs to train the network               | n >= 0       | int             |
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