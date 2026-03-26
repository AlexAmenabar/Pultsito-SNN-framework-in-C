#include "networks/snn_generator.h"

#include <stddef.h>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>

int main(int argc, char *argv[]) {
    
    // randomize exection to create different networks each time
    srand(time(NULL));

    // load configuration file
    generator_conf_t *conf = read_configuration_file(argv[1]);

    // initialize topology
    topology_t topology = generate_topology(conf);

    // initialize neurons and synapses
    topology.neurons = initialize_neurons(conf);
    topology.synapses = initialize_synapses(conf);

    // store generated network
    store_network(&topology, conf, 0);
}