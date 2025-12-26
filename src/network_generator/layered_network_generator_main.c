//#include "network_generator/network_generator.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>

int main(int argc, char *argv[]) {
    
    // randomize exection to create different networks each time
    srand(time(NULL));


    // read how much neurons are on each layer
    char *file_name;
    char *original_file_name, *file_name_neurons, *file_name_synapses;
    FILE *f, *f_neurons, *f_synapses;
    int i;




    file_name = argv[1];
    int n_layers = strtoul(argv[2], NULL, 10);
    printf(" > Number of layers = %u\n", n_layers);

    int *n_neuron_per_layer = (int*)malloc(n_layers * sizeof(int));
    int *first_neuron_in_layer = (int*)malloc(n_layers * sizeof(int));
    int n_neurons = 0;
    int n_synapses = 0;
    for(size_t i=0; i<n_layers; i++){

        first_neuron_in_layer[i] = n_neurons;
        n_neuron_per_layer[i] = strtoul(argv[i+3], NULL, 10);
        n_neurons += n_neuron_per_layer[i];

        printf(" n neurons in layer %d: %d\n", i, n_neuron_per_layer[i]);

        if(i > 0){
            n_synapses += (n_neuron_per_layer[i] * n_neuron_per_layer[i-1]);
        }
    }

    // add input and output synapses
    // TMP
    n_synapses += n_neuron_per_layer[0];
    n_synapses += n_neuron_per_layer[n_layers-1];

    printf(" n_neurons = %d, n_synapses = %d\n", n_neurons, n_synapses);


    f = fopen(file_name, "w"); // write
    if (f == NULL){
        perror("Error opening the file\n");
        exit(1);
    }    
    
    // write general information
    printf("Storing general information...\n");

    fprintf(f, "[general]\n");
    fprintf(f, "    neurons = %d\n", n_neurons);
    fprintf(f, "    input_neurons = %d\n", n_neuron_per_layer[0]);
    fprintf(f, "    output_neurons = %d\n", n_neuron_per_layer[n_layers - 1]);
    fprintf(f, "    synapsis = %d\n", n_synapses);
    fprintf(f, "    input_synapsis = %d\n", n_neuron_per_layer[0]);
    fprintf(f, "    output_synapsis = %d\n", n_neuron_per_layer[n_layers - 1]);
    fprintf(f, "    network_is_separated = %d\n\n", 1);

    printf("General information stored!\n");

    /*
    If output is separaed is 1, then all lists are stored in apart files with the name of the
    original file and an extra indicating what parameter is stored into that file. Else, all
    is stored in only one file. 
    */

    original_file_name = strtok(file_name, "."); // split file name by extension
    int len = strlen(original_file_name);

    file_name_neurons = malloc((len + 20) * sizeof(char));
    file_name_synapses = malloc((len + 25) * sizeof(char));

    strcpy(file_name_neurons, original_file_name);  // copy original file name
    strcpy(file_name_synapses, original_file_name);
    
    strcat(file_name_neurons, "_neurons.toml"); // add extension to original file name for neurons and synapses
    strcat(file_name_synapses, "_synapses.toml");

    printf("%s\n", file_name_neurons);
    printf("%s\n", file_name_synapses);


    // open files if information must be stored separated
    //if(conf->output_is_separated == 1){
    f_neurons = fopen(file_name_neurons, "w"); // write
    if (f_neurons == NULL){
        perror("Error opening the file\n");
        exit(1);
    }  

    f_synapses = fopen(file_name_synapses, "w"); // write
    if (f_synapses == NULL){
        perror("Error opening the file\n");
        exit(1);
    }  

    // write neurons information
    printf("Storing neurons information...\n");

    fprintf(f, "[neurons]\n");
    fprintf(f, "    behaviour = 1\n"); // TODO: THIS IS TEMPORAL



    // set input and output synapse indexes



    fprintf(f, "\n");

    // if is separated == 1, store neurons information in plain text
    //if(output_is_separated == 1){
    for(i=0; i<n_neurons; i++){
        fprintf(f_neurons, "%d ", 1);
    }
    fprintf(f_neurons, "\n");

    for(i=0; i<n_neurons; i++){
        fprintf(f_neurons, "%lf ", 1.0);
    }
    fprintf(f_neurons, "\n");

    for(i=0; i<n_neurons; i++){
        fprintf(f_neurons, "%lf ", 0.0);
    }
    fprintf(f_neurons, "\n");

    for(i=0; i<n_neurons; i++){
        fprintf(f_neurons, "%d ", 1);
    }
    fprintf(f_neurons, "\n");

    /*for(i=0; i<n_neurons; i++){
        fprintf(f_neurons, "%d ", n_input_synapses_per_neuron[i]);
    }
    fprintf(f_neurons, "\n");

    for(i=0; i<network_data->n_neurons; i++){
        fprintf(f_neurons, "%d ", network_data->n_output_synapses_per_neuron[i]);
    }*/

    for(int i = 0; i<n_layers; i++){

        // input layer
        if(i == 0){
        
            for(int j = 0; j<n_neuron_per_layer[i]; j++){
            
                fprintf(f_neurons, "%d ", 1);
            }
        }
        else{

            for(int j = 0; j<n_neuron_per_layer[i]; j++){
                
                fprintf(f_neurons, "%d ", n_neuron_per_layer[i-1]);
            }
        }
    }
    fprintf(f_neurons, "\n");

    // set input and output synapse indexes
    for(int i = 0; i<n_layers; i++){

        // input layer
        if(i < n_layers - 1){
            
            for(int j = 0; j<n_neuron_per_layer[i]; j++){
                
                fprintf(f_neurons, "%d ", n_neuron_per_layer[i + 1]);
            }
        }
        else{

            // output layer
            for(int j = 0; j<n_neuron_per_layer[i]; j++){
                
                fprintf(f_neurons, "%d ", 1);
            }
        }
    }
    fprintf(f_neurons, "\n");

    printf("Neurons information stored!\n");



    // write synapses information
    printf("Storing synapses information...\n");
    
    /*fprintf(f, "[synapsis]\n");

    fprintf(f, "    latency = 1\n"); // TODO: THIS IS TEMPORAL
    
    if(conf->output_is_separated != 1){
        fprintf(f, "    latency_list = [");
        for(i=0; i<network_data->n_synapses-1; i++){
            fprintf(f, "%d, ", network_data->latency_list[i]);
        }
        fprintf(f, "%d]\n", network_data->latency_list[network_data->n_synapses-1]);
    }

    if(conf->output_is_separated != 1){
        fprintf(f, "    weights = [");
        for(i=0; i<network_data->n_synapses-1; i++){
            fprintf(f, "%lf, ", network_data->weights[i]);
        }
        fprintf(f, "%lf]\n", network_data->weights[network_data->n_synapses-1]);
    }

    fprintf(f, "    training_zones = 0\n"); // TODO: THIS IS TEMPORAL
    if(conf->output_is_separated != 1){
        fprintf(f, "    training_zones_list = [");
        for(i=0; i<network_data->n_synapses-1; i++){
            fprintf(f, "%d, ", network_data->training_zones_list[i]);
        }
        fprintf(f, "%d]\n", network_data->training_zones_list[network_data->n_synapses-1]);
    }

    if(conf->output_is_separated != 1){
        fprintf(f, "    connections = [");
        for(int i=0; i<(network_data->n_neurons+1); i++){
            fprintf(f, "[");
            for(int j=0; j<(network_data->connections[i][0] * 2); j++)
                fprintf(f, "%d, ", network_data->connections[i][j]);

            if(i < network_data->n_neurons)
                fprintf(f, "%d], ", network_data->connections[i][network_data->connections[i][0] * 2]);
            else
                fprintf(f, "%d]", network_data->connections[i][network_data->connections[i][0]* 2]);
        }
        fprintf(f, "]\n");
    }*/
    printf("Synapses information stored!\n");

    // if is separated == 1, store neurons information in plain text
    //if(output_is_separated == 1){
        
    for(i=0; i<n_synapses; i++){
        fprintf(f_synapses, "%d ", 1);
    }
    fprintf(f_synapses, "\n");

    for(i=0; i<n_synapses; i++){
        fprintf(f_synapses, "%lf ", 0.5);
    }
    fprintf(f_synapses, "\n");

    for(i=0; i<n_synapses; i++){
        fprintf(f_synapses, "%d ", 0);
    }
    fprintf(f_synapses, "\n\n");


    // COnnectivity
    fprintf(f_synapses, "%d ", n_neuron_per_layer[0]);
    for(int j = 0; j<n_neuron_per_layer[0]; j++){
        fprintf(f_synapses, "%d %d ", j, 1);
    }
    fprintf(f_synapses, "\n");

    for(int i = 0; i<n_layers; i++){

        // hidden layer neurons
        for(int j = 0; j<n_neuron_per_layer[i]; j++){

            // number of neurons to which is connected
            if(i == n_layers - 1){
                fprintf(f_synapses, "%d ", n_neuron_per_layer[i+1] + 1);
            }
            else{
                fprintf(f_synapses, "%d ", n_neuron_per_layer[i+1]);
            }

            for(int l = 0; l<n_neuron_per_layer[i+1]; l++){
                
                fprintf(f_synapses, "%d %d ", first_neuron_in_layer[i+1] + l, 1);
            }
            if(i == n_layers - 1){

                fprintf(f_synapses, "%d %d ", -1, 1);
            }
            fprintf(f_synapses,"\n");
        }

        //fprintf(f_synapses, "\n");
    }

    fclose(f);
    fclose(f_neurons);
    fclose(f_synapses);
}