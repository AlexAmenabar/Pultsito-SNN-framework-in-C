/// FUNCTIONS WITH SIMULATION TYPES
#include <omp.h>

#include "snn_library.h"
#include "load_data.h"
#include "helpers.h"
#include "training_rules/stdp.h"

#include "neuron_models/lif_neuron.h"
#include "metrics.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>



// DEPRECATED???
void simulate(spiking_nn_t *snn, simulation_configuration_t *conf, simulation_results_t *results, input_data_t *dataset){

    // initialize several control variables
    int n = snn->n_neurons / conf->n_process * 0.1;

    int time_step = 0, i; 
    int n_process = conf->n_process;
    struct timespec start, end; // to measure simulation complete time
    struct timespec start_neurons, end_neurons; // to measure neurons simulation time
    struct timespec start_neurons_input, end_neurons_input; // to measure input synapses simulation time
    struct timespec start_neurons_output, end_neurons_output; // to measure output synapses simulation time
    struct timespec start_learning, end_learning; // to measure learning processing time
    
    // TODO: revise this
    simulation_results_per_sample_t *results_per_sample = &(results->results_per_sample[0]);

    printf(" > n_processes = %d\n", conf->n_process);
    fflush(stdout);

    // check simulation schema
    clock_gettime(CLOCK_MONOTONIC, &start);

    if(conf->simulation_type == 0) // clock-based
    { 
        // simulate over simulation time steps
        for(time_step = 0; time_step < conf->time_steps; time_step++)
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
                #pragma omp for schedule(guided, n) private(i) 
                for(i=0; i<snn->n_neurons; i++)
                    snn->input_step(snn, time_step, i, results_per_sample);
                
                #pragma omp master
                {
                clock_gettime(CLOCK_MONOTONIC, &end_neurons_input);
                clock_gettime(CLOCK_MONOTONIC, &start_neurons_output);
                }

                // simulate output step
                #pragma omp for schedule(guided, n) private(i)
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
                //#ifndef NOLEARN 
                if(conf->learn == 1){
                    #pragma omp for schedule(static, 50) private(i)
                    for(i = 0; i<snn->n_synapses; i++)
                            snn->synapses[i].learning_rule(&(snn->synapses[i]), time_step, 5); 
                }
                //#endif

                #pragma omp master
                {
                    clock_gettime(CLOCK_MONOTONIC, &end_learning);
                }

            }
                    
            // store neuron simulation and training rule simulation times // TODO: revise
            results_per_sample->elapsed_time_neurons += (end_neurons.tv_sec - start_neurons.tv_sec) + (end_neurons.tv_nsec - start_neurons.tv_nsec) / 1e9;
            results_per_sample->elapsed_time_neurons_input += (end_neurons_input.tv_sec - start_neurons_input.tv_sec) + (end_neurons_input.tv_nsec - start_neurons_input.tv_nsec) / 1e9;
            results_per_sample->elapsed_time_neurons_output += (end_neurons_output.tv_sec - start_neurons_output.tv_sec) + (end_neurons_output.tv_nsec - start_neurons_output.tv_nsec) / 1e9;
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
    }
    // event-driven
    else 
    {
        //TODO
        printf(" > NOT IMPLEMENTED YET!\n");
    }
    
    clock_gettime(CLOCK_MONOTONIC, &end);


    results->elapsed_time += (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;
}

void simulate_samples(spiking_nn_t *snn, simulation_configuration_t *conf, simulation_results_t *results, input_data_t *dataset){

    // initialize several control variables
    int epochs, n_samples, time_steps, batch_size;
    spiking_nn_t *pr_snns;
    spiking_nn_t *tmp_snn;
    

    if(conf->n_process > conf->batch_size)
        conf->n_process = conf->batch_size;

    // cpy important information
    epochs = conf->epochs;
    n_samples = dataset->n_samples;
    time_steps = conf->time_steps;
    batch_size = conf->batch_size; // batch size


    int i, t, s, e, b; 
    int n_process = conf->n_process;
    struct timespec start, end; // to measure simulation complete time
    struct timespec start_epoch, end_epoch; // to measure simulation time for epochs
    struct timespec start_sample, end_sample; // to measure simulation time for samples
    struct timespec start_neurons_input, end_neurons_input; // to measure input synapses simulation time
    struct timespec start_neurons_output, end_neurons_output; // to measure output synapses simulation time
    struct timespec start_learning, end_learning; // to measure learning processing time
    struct timespec start_load_sample, end_load_sample; // to measure learning processing time
    struct timespec start_re_neurons, end_re_neurons; // to measure learning processing time
    struct timespec start_re_synapses, end_re_synapses; // to measure learning processing time

    // create array for batches
    int n_batches;
    int *batches = create_batches(&n_batches, dataset, n_samples, batch_size);

    // TODO: revise this
    simulation_results_per_sample_t *results_per_sample;

    // start measuring time
    clock_gettime(CLOCK_MONOTONIC, &start);

    // configure the OpenMP parameters depending on parallelization strategy
    int p_outer = 0, p_inner = 0, n_outer = 0, n_inner = 0;
    #ifdef NESTED
        omp_set_dynamic(0); // disable dynamic threads (I do not understand this very well)
        omp_set_nested(1); // allow nested parallelism
        omp_set_max_active_levels(2);

        // number of processes
        p_inner = conf->n_inner_process;
        n_inner = snn->n_neurons / p_inner * 0.1;
        if(n_inner == 0) n_inner = 1;
        
        p_outer = n_process / p_inner;
        n_outer = 1;//snn->n_neurons / p_outer * 0.1;

    #elif PAR_SAMPLES
        p_outer = n_process;
        n_outer = 1;//snn->n_neurons / p_outer * 0.1;
    #else
        p_inner = n_process;
        n_inner = snn->n_neurons / p_inner * 0.1;
        if(n_inner == 0) n_inner = 1;
    #endif


    // create network copies if necessary
    #if defined PAR_SAMPLES || defined NESTED 
        pr_snns = (spiking_nn_t *)malloc(p_outer * sizeof(spiking_nn_t));
        for(i = 0; i<p_outer; i++){
            cp_network(&(pr_snns[i]), snn, conf);
        }
    #endif



    // open files for results storage
    open_results_files(conf);


    // simulate epochs
    for(e = 0; e<epochs; e++){

        printf(" > Epoch %d\n", e);
        fflush(stdout);


        // shuffle samples in batch
        if(conf->shuffle_samples == 1)
            shuffle_sample_indexes(batches, n_samples);

        /*for(b = 0; b<n_batches; b++){

            printf("[");
            for(s = 0; s<batch_size; s++){

                printf("%d ", batches[b * batch_size + s]);
            }
            printf("]\n");
        }*/

        // I think all below should be moved to another function

        // loop over batches
        for(b = 0; b<n_batches; b++){
            
            // batch simulation started
            clock_gettime(CLOCK_MONOTONIC, &start_epoch);
            
            //printf(" In batch %d\n", b);

            // loop over batch samples
            // func simulate_batch()

            #if defined PAR_SAMPLES || defined NESTED
            #pragma omp parallel for num_threads(p_outer) schedule(dynamic, n_outer) private(i, s, t, tmp_snn, results_per_sample, start_sample, end_sample, start_re_neurons, end_re_neurons, start_re_synapses, end_re_synapses, start_load_sample, end_load_sample, start_neurons_input, end_neurons_input, start_neurons_output, end_neurons_output, start_learning, end_learning)
            #endif
            for(s = 0; s<batch_size; s++){

                //printf(" > In sample %d\n", s);

                
                // simulate if there is a sample to be simulated
                if(batches[batch_size * b + s] != -1){
                    
                    int sample_index = batches[batch_size * b + s];

                    // get network
                    #if defined PAR_SAMPLES || defined NESTED
                    //printf(" Thread = %d, sample %d\n", omp_get_thread_num(), s);
                    tmp_snn = &(pr_snns[omp_get_thread_num()]);
                    #else
                    tmp_snn = snn;
                    #endif


                    // reinitialize network properties for the sample: move to a function??

                    // sample started
                    clock_gettime(CLOCK_MONOTONIC, &start_sample);

                    // get results struct and reinitialize necessary data // TODO: actually is temporal
                    results_per_sample = &(results->results_per_sample[sample_index]);
                    
                    #if !defined PAR_SAMPLES || defined NESTED
                    #pragma omp parallel for num_threads(p_inner) schedule(guided, n_inner) private(i) 
                    #endif
                    for(i = 0; i<tmp_snn->n_neurons; i++){
                        results_per_sample->n_spikes_per_neuron[i] = 0;
                    }


                    // reinitialize neurons O(n)
                    clock_gettime(CLOCK_MONOTONIC, &start_re_neurons);

                    #if !defined PAR_SAMPLES || defined NESTED
                    #pragma omp parallel for num_threads(p_inner) schedule(guided, n_inner) private(i) 
                    #endif
                    for(i = 0; i<tmp_snn->n_neurons; i++){
                        tmp_snn->neuron_re_initializer(tmp_snn, i);
                    }
                    clock_gettime(CLOCK_MONOTONIC, &end_re_neurons);


                    // reinitialize synapses O(m)
                    clock_gettime(CLOCK_MONOTONIC, &start_re_synapses);

                    #if !defined PAR_SAMPLES || defined NESTED
                    #pragma omp parallel for num_threads(p_inner) schedule(guided, n_inner) private(i) 
                    #endif
                    for(i = 0; i<tmp_snn->n_synapses; i++){
                        re_initialize_synapse(&(tmp_snn->synapses[i]));
                    }
                    clock_gettime(CLOCK_MONOTONIC, &end_re_synapses);

                    // load sample in network 
                    clock_gettime(CLOCK_MONOTONIC, &start_load_sample);
                    tmp_snn->load_sample(tmp_snn, &(dataset->samples[sample_index]));
                    clock_gettime(CLOCK_MONOTONIC, &end_load_sample);

                    // simulate time steps for each sample
                    for(t = 0; t<time_steps; t++){

                        // input step
                        clock_gettime(CLOCK_MONOTONIC, &start_neurons_input);

                        #if !defined PAR_SAMPLES || defined NESTED
                        #pragma omp parallel for num_threads(p_inner) schedule(guided, n_inner) private(i) 
                        #endif
                        for(i = 0; i<tmp_snn->n_neurons; i++) // O(n)
                            tmp_snn->input_step(tmp_snn, t, i, results_per_sample);

                        clock_gettime(CLOCK_MONOTONIC, &end_neurons_input);

                        // output step
                        clock_gettime(CLOCK_MONOTONIC, &start_neurons_output);
                        
                        #if !defined PAR_SAMPLES || defined NESTED
                        #pragma omp parallel for num_threads(p_inner) schedule(guided, n_inner) private(i) 
                        #endif
                        for(i = 0; i<tmp_snn->n_neurons; i++) // O(n)
                            tmp_snn->output_step(tmp_snn, t, i, results_per_sample);

                        clock_gettime(CLOCK_MONOTONIC, &end_neurons_output);

                        // learning rule
                        clock_gettime(CLOCK_MONOTONIC, &start_learning);
                        
                        if(conf->learn == 1){

                            #if !defined PAR_SAMPLES || defined NESTED
                            #pragma omp parallel for num_threads(p_inner) schedule(guided, n_inner) private(i) 
                            #endif
                            for(i = 0; i<tmp_snn->n_synapses; i++) // O(m)
                                tmp_snn->synapses[i].learning_rule(&(tmp_snn->synapses[i]), t, 3); // TODO: Change this!! 
                        }

                        clock_gettime(CLOCK_MONOTONIC, &end_learning);
                    
                        // sum times
                        results_per_sample->elapsed_time_neurons_input += (end_neurons_input.tv_sec - start_neurons_input.tv_sec) + (end_neurons_input.tv_nsec - start_neurons_input.tv_nsec) / 1e9;
                        results_per_sample->elapsed_time_neurons_output += (end_neurons_output.tv_sec - start_neurons_output.tv_sec) + (end_neurons_output.tv_nsec - start_neurons_output.tv_nsec) / 1e9;
                        results_per_sample->elapsed_time_learning += (end_learning.tv_sec - start_learning.tv_sec) + (end_learning.tv_nsec - start_learning.tv_nsec) / 1e9;
                    }

                    // sample finished
                    clock_gettime(CLOCK_MONOTONIC, &end_sample);

                    results_per_sample->elapsed_time_re_neurons += (end_re_neurons.tv_sec - start_re_neurons.tv_sec) + (end_re_neurons.tv_nsec - start_re_neurons.tv_nsec) / 1e9;
                    results_per_sample->elapsed_time_re_synapses += (end_re_synapses.tv_sec - start_re_synapses.tv_sec) + (end_re_synapses.tv_nsec - start_re_synapses.tv_nsec) / 1e9;
                    results_per_sample->elapsed_time_load_sample += (end_load_sample.tv_sec - start_load_sample.tv_sec) + (end_load_sample.tv_nsec - start_load_sample.tv_nsec) / 1e9;
                    results_per_sample->elapsed_time_sample += (end_sample.tv_sec - start_sample.tv_sec) + (end_sample.tv_nsec - start_sample.tv_nsec) / 1e9;
                }
            }


            // func store batch results()
            // reinitialize 
            results->elapsed_time_neurons_input = 0; // elapsed time processing neurons input 
            results->elapsed_time_neurons_output = 0; // elapsed time processing neurons output
            results->elapsed_time_learning = 0; // elapsed time processing learning rules
            results->elapsed_time_sample = 0;
            results->elapsed_time_load_sample = 0;
            results->elapsed_time_re_neurons = 0;
            results->elapsed_time_re_synapses = 0;


            // update weights (in case there are copies, first sum all copies)
            // copies // FUNCTION???
            
            // func update_weighs()
            if(conf->learn == 1){
                clock_gettime(CLOCK_MONOTONIC, &start_learning);
                #if defined PAR_SAMPLES || defined NESTED

                    // sum weight differences of all copies
                    for(s = 0; s<p_outer; s++){
                        #pragma omp parallel for num_threads(p_outer) schedule(dynamic, n_outer) private(i, t, tmp_snn, results_per_sample, start_sample, end_sample, start_re_neurons, end_re_neurons, start_re_synapses, end_re_synapses, start_load_sample, end_load_sample, start_neurons_input, end_neurons_input, start_neurons_output, end_neurons_output, start_learning, end_learning)
                        for(i = 0; i < snn->n_synapses; i++){
                            snn->synapses[i].dw += pr_snns[s].synapses[i].dw;
                        }
                    }
                    
                    // update weights in original snn
                    #pragma omp parallel for num_threads(p_outer) schedule(dynamic, n_outer) private(i, s, t, tmp_snn, results_per_sample, start_sample, end_sample, start_re_neurons, end_re_neurons, start_re_synapses, end_re_synapses, start_load_sample, end_load_sample, start_neurons_input, end_neurons_input, start_neurons_output, end_neurons_output, start_learning, end_learning)
                    for(i = 0; i<snn->n_synapses; i++){    
                        snn->synapses[i].w += snn->synapses[i].dw / batch_size;
                        snn->synapses[i].init_w = snn->synapses[i].w;
                        snn->synapses[i].dw = 0;
                    }    

                    // update weights in copies
                    #pragma omp parallel for num_threads(p_outer) schedule(dynamic, n_outer) private(i, s, t, tmp_snn, results_per_sample, start_sample, end_sample, start_re_neurons, end_re_neurons, start_re_synapses, end_re_synapses, start_load_sample, end_load_sample, start_neurons_input, end_neurons_input, start_neurons_output, end_neurons_output, start_learning, end_learning)
                    for(s = 0; s<p_outer; s++){
                        for(i = 0; i < snn->n_synapses; i++){
                            pr_snns[s].synapses[i].w = snn->synapses[i].w;
                            pr_snns[s].synapses[i].init_w = snn->synapses[i].init_w;
                            pr_snns[s].synapses[i].dw = 0;
                        }
                    }

                #else
                    #pragma omp parallel for num_threads(p_inner) schedule(guided, n_inner) private(i) 
                    for(i = 0; i<snn->n_synapses; i++){
                        
                        // update w
                        snn->synapses[i].w += snn->synapses[i].dw / batch_size;
                        snn->synapses[i].init_w = snn->synapses[i].w;
                        snn->synapses[i].dw = 0;
                    }    
                #endif
                clock_gettime(CLOCK_MONOTONIC, &end_learning);
                results->elapsed_time_learning += (end_learning.tv_sec - start_learning.tv_sec) + (end_learning.tv_nsec - start_learning.tv_nsec) / 1e9;
            }

            // sum all execution times
            for(i = 0; i<batch_size; i++){
                
                if(batches[batch_size * b + i] != -1){
                    int sample = batches[batch_size * b + i];

                    results->elapsed_time_neurons_input += results->results_per_sample[sample].elapsed_time_neurons_input; // elapsed time processing neurons input 
                    results->elapsed_time_neurons_output += results->results_per_sample[sample].elapsed_time_neurons_output; // elapsed time processing neurons output
                    results->elapsed_time_learning += results->results_per_sample[sample].elapsed_time_learning; // elapsed time processing learning rules
                    results->elapsed_time_sample += results->results_per_sample[sample].elapsed_time_sample;
                    results->elapsed_time_load_sample += results->results_per_sample[sample].elapsed_time_load_sample;
                    results->elapsed_time_re_neurons += results->results_per_sample[sample].elapsed_time_re_neurons;
                    results->elapsed_time_re_synapses += results->results_per_sample[sample].elapsed_time_re_synapses;
                }
            }
            // batch finished
            clock_gettime(CLOCK_MONOTONIC, &end_epoch);
            results->elapsed_time_epoch = (end_epoch.tv_sec - start_epoch.tv_sec) + (end_epoch.tv_nsec - start_epoch.tv_nsec) / 1e9;
    
            // store results of the batch
            store_results(results, conf, snn, dataset, batches, b);
        }
        
    }

    clock_gettime(CLOCK_MONOTONIC, &end);

    results->elapsed_time = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;
    
    printf(" > Total elapsed time %lf\n", results->elapsed_time);

    // close results files
    close_results_files(conf);
}

void train_network(spiking_nn_t *snn, simulation_configuration_t *conf, simulation_results_t *results, input_data_t *dataset){

    // store the configuration
    int learn  = conf->learn;
    
    // it's training. so learn = 1
    conf->learn = 1;

    simulate_samples(snn, conf, results, dataset);

    // recover the previous value
    conf->learn = learn;
}

void test_network(spiking_nn_t *snn, simulation_configuration_t *conf, simulation_results_t *results, input_data_t *dataset){

    int learn, epochs; 

    // store information
    learn = conf->learn;
    epochs = conf->epochs;
    
    // update information
    conf->learn = 0; // no learn
    conf->epochs = 1; // only one epoch in test, we are not training

    simulate_samples(snn, conf, results, dataset);

    // recover previous values
    conf->learn = learn;
    conf->epochs = epochs;
}