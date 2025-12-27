//#include "network_generator/network_generator.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>


size_t rand_lim(int limit) {
/* return a random number between 0 and limit inclusive.
 */

    int divisor = RAND_MAX/(limit+1);
    int retval;

    do { 
        retval = rand() / divisor;
    } while (retval > limit);

    return (size_t)retval;
}

void selectionSort(size_t* arr) {
    
    size_t n = arr[0] * 2;

    for (size_t i = 1; i < n; i+=2) {
        
        size_t min = i; // i index minimum
        for (size_t j = i + 2; j < n; j+=2) {
            if (arr[j] < arr[min])
                min = j;
        }

        if (min != i) {
            size_t temp = arr[min];
            size_t temp2 = arr[min + 1];

            arr[min] = arr[i];
            arr[min + 1] = arr[i + 1];

            arr[i] = temp;
            arr[i + 1] = temp2;
        }
    }
}


int main(int argc, char *argv[]) {
    
    // randomize exection to create different networks each time
    srand(time(NULL));

    // read how much neurons are on each layer

    size_t n_neurons = strtoul(argv[1], NULL, 10);
    size_t n_input = strtoul(argv[2], NULL, 10); // virtual neurons...
    size_t n_output_neurons = strtoul(argv[3], NULL, 10);
    size_t n_synapses = strtoul(argv[4], NULL, 10);

    // compute number of input synapses for each neuron
    size_t connect_per_neuron = n_synapses / (n_neurons + n_input);
    size_t max_connections_per_neuron = n_synapses / 2;

    // the first n_input neurons are the input ones
    size_t **input_neurons_per_neuron = (size_t**)malloc(n_neurons * sizeof(size_t*));

    for(size_t i=0; i<n_neurons; i++){
        
        input_neurons_per_neuron[i] = (size_t*)malloc(max_connections_per_neuron * sizeof(size_t));
    }


    n_synapses = 0;
    // set input neurons for each neuron (i -> i)
    for(size_t i=0; i<n_neurons; i++){
        
        size_t neuron_connections = connect_per_neuron;
        input_neurons_per_neuron[i][0] = neuron_connections;

        for(size_t j=0; j<input_neurons_per_neuron[i][0]; j++){
            
            // neuron i receives at least input spike train i
            if(i <= n_input && j == i){

                input_neurons_per_neuron[i][0] ++;
                input_neurons_per_neuron[i][1] = j;
                input_neurons_per_neuron[i][2] = 1;
                n_synapses ++;
            }
            else{
                size_t in_neuron = rand_lim(n_neurons + n_input);
                if(in_neuron == 0) 
                    in_neuron = 1;

                size_t n = rand_lim(4);
                if(n == 0) 
                    n = 1;

                input_neurons_per_neuron[i][j * 2 + 1] = in_neuron;
                input_neurons_per_neuron[i][j * 2 + 2] = n;
                n_synapses += n;
            }
        }

        // sort array
        selectionSort(input_neurons_per_neuron[i]);
    }

    char *file_name = argv[5];
    FILE *f;
    f = fopen(argv[5], "w"); // write
    if (f == NULL){
        perror("Error opening the file\n");
        exit(1);
    }    
    
    // write general information
    printf("Storing general information...\n");

    fprintf(f, "[general]\n");
    fprintf(f, "    neurons = %d\n", n_neurons);
    fprintf(f, "    input_neurons = %d\n", n_input);
    fprintf(f, "    output_neurons = %d\n", n_output_neurons);
    fprintf(f, "    synapsis = %d\n", n_synapses);
    fprintf(f, "    input_synapsis = %d\n", n_input);
    fprintf(f, "    output_synapsis = %d\n", 0);
    fprintf(f, "    network_is_separated = %d\n\n", 1);

    printf("General information stored!\n");

    /*
    If output is separaed is 1, then all lists are stored in apart files with the name of the
    original file and an extra indicating what parameter is stored into that file. Else, all
    is stored in only one file. 
    */

    char *original_file_name;
    original_file_name = strtok(file_name, "."); // split file name by extension
    int len = strlen(original_file_name);

    char* file_name_neurons = malloc((len + 20) * sizeof(char));
    char* file_name_synapses = malloc((len + 25) * sizeof(char));

    strcpy(file_name_neurons, original_file_name);  // copy original file name
    strcpy(file_name_synapses, original_file_name);
    
    strcat(file_name_neurons, "_neurons.toml"); // add extension to original file name for neurons and synapses
    strcat(file_name_synapses, "_synapses.toml");

    printf("%s\n", file_name_neurons);
    printf("%s\n", file_name_synapses);


    // open files if information must be stored separated
    //if(conf->output_is_separated == 1){
    FILE *f_neurons, *f_synapses;
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
    size_t i;
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

    // store connectivity
    for(size_t i = 0; i<n_neurons; i++){

        fprintf(f_synapses, "%zu ", input_neurons_per_neuron[i][0]);
        
        for(size_t j = 0; j<input_neurons_per_neuron[i][0]; j++){

            fprintf(f_synapses, "%zu ", input_neurons_per_neuron[i][j * 2 + 1]);
            fprintf(f_synapses, "%zu ", input_neurons_per_neuron[i][j * 2 + 2]);
        }

        fprintf(f_synapses,"\n");

    }

    fclose(f);
    fclose(f_neurons);
    fclose(f_synapses);
}