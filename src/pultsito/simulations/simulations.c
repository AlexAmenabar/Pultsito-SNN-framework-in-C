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


void simulate(spiking_nn_t *snn, simulation_configuration_t *conf, simulation_results_t *results, int learning){

    // initialize several control variables
    int time_step = 0, i; 
    int n_process = conf->n_process;
    struct timespec start, end; // to measure simulation complete time
    struct timespec start_neurons, end_neurons; // to measure neurons simulation time
    struct timespec start_neurons_input, end_neurons_input; // to measure input synapses simulation time
    struct timespec start_neurons_output, end_neurons_output; // to measure output synapses simulation time
    struct timespec start_synapses, end_synapses; // to measure synapses simulation time
    struct timespec start_synapses_input, end_synapses_input; // to measure synapses simulation time
    struct timespec start_synapses_output, end_synapses_output; // to measure synapses simulation time
    struct timespec start_learning, end_learning; // to measure learning processing time
    
    // TODO: revise this
    simulation_results_per_sample_t *results_per_sample = &(results->results_per_sample[0]);

    // check simulation schema
    clock_gettime(CLOCK_MONOTONIC, &start);

    if(conf->simulation_type == 0) // clock-based
    { 
        // start measuring time

        // if CUDA is defined, simulate on the GPU
    #ifdef CUDA
         {

         }
    #else
        // simulate over simulation time steps
        while(time_step < conf->time_steps)
        {

            #ifdef DEBUG 
                printf("Processing time step %d\n", time_step);
            #endif       

            // process all neurons (input and output steps)
            #pragma omp parallel num_threads(n_process)
            {

                // if OpenMP is used, simulate first input synapses and then output
                #ifdef OPENMP
                #pragma omp master
                {
                clock_gettime(CLOCK_MONOTONIC, &start_neurons);
                clock_gettime(CLOCK_MONOTONIC, &start_neurons_input);
                }

                // simulate input step
                #pragma omp for schedule(static, 10) private(i) 
                for(i=0; i<snn->n_neurons; i++)
                    snn->input_step(snn, time_step, i, results_per_sample);
                
                #pragma omp master
                {
                clock_gettime(CLOCK_MONOTONIC, &end_neurons_input);
                clock_gettime(CLOCK_MONOTONIC, &start_neurons_output);
                }

                // simulate output step
                #pragma omp for schedule(static, 10) private(i)
                for(i=0; i<snn->n_neurons; i++)
                    snn->output_step(snn, time_step, i, results_per_sample);
     

                // if it is serially simulated, do all on one step
                #else  
                for(i = 0; i<snn->n_neurons; i++)
                    snn->complete_step(snn, time_step, i, results_per_sample);
                #endif

                // time measuring
                #pragma omp master
                {
                    clock_gettime(CLOCK_MONOTONIC, &end_neurons_output);
                    clock_gettime(CLOCK_MONOTONIC, &end_neurons);
                    clock_gettime(CLOCK_MONOTONIC, &start_learning);
                }

                #ifdef DEBUG 
                printf("Processing stdp on time step %d\n", time_step);
                #endif      

                // process learning if it is necessary
                #ifndef NOLEARN 
                if(learning == 1){
                    #pragma omp for schedule(static, 50) private(i)
                    for(i = 0; i<snn->n_synapses; i++)
                            snn->synapses[i].learning_rule(&(snn->synapses[i])); 
                }
                #endif

                #pragma omp master
                {
                    clock_gettime(CLOCK_MONOTONIC, &end_learning);
                }

            }
            

            // time step computed
            time_step++;
        
            // store neuron simulation and training rule simulation times // TODO: revise
            results_per_sample->elapsed_time_neurons += (end_neurons.tv_sec - start_neurons.tv_sec) + (end_neurons.tv_nsec - start_neurons.tv_nsec) / 1e9;
            results_per_sample->elapsed_time_neurons_input += (end_neurons_input.tv_sec - start_neurons_input.tv_sec) + (end_neurons_input.tv_nsec - start_neurons_input.tv_nsec) / 1e9;
            results_per_sample->elapsed_time_neurons_output += (end_neurons_output.tv_sec - start_neurons_output.tv_sec) + (end_neurons_output.tv_nsec - start_neurons_output.tv_nsec) / 1e9;
            //results_per_sample->elapsed_time_synapses += (end_synapses.tv_sec - start_synapses.tv_sec) + (end_synapses.tv_nsec - start_synapses.tv_nsec) / 1e9;
            //results_per_sample->elapsed_time_synapses_input += (end_synapses_input.tv_sec - start_synapses_input.tv_sec) + (end_synapses_input.tv_nsec - start_synapses_input.tv_nsec) / 1e9;
            //results_per_sample->elapsed_time_synapses_output += (end_synapses_output.tv_sec - start_synapses_output.tv_sec) + (end_synapses_output.tv_nsec - start_synapses_output.tv_nsec) / 1e9;
            results_per_sample->elapsed_time_learning += (end_learning.tv_sec - start_learning.tv_sec) + (end_learning.tv_nsec - start_learning.tv_nsec) / 1e9;

            // Print info about the simulation
            if(time_step % 1000 == 0)
                printf("Time step: %d/%d\n", time_step, conf->time_steps);

            #ifdef DEBUG 
                printf(" - Printing synapses weights: \n");
                for(i = 0; i<snn->n_synapses; i++)
                    printf(" -- Synapse %d: %f\n", i, snn->synapses[i].w);

                printf("\n=======================================\n\n");
            #endif      

        }
    #endif
    }
    // event-driven
    else 
    {
        //TODO
        printf(" > NOT IMPLEMENTED YET!\n");
    }
    
    clock_gettime(CLOCK_MONOTONIC, &end);


    results_per_sample->elapsed_time += (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;
}

// TODO
void sample_simulation(spiking_nn_t *snn, simulation_configuration_t *conf, simulation_results_t *results){

}

// TODO
void stream_simulation(spiking_nn_t *snn, simulation_configuration_t *conf, simulation_results_t *results){

}