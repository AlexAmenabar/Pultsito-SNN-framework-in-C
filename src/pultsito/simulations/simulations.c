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


#include <immintrin.h>  // AVX intrinsics



// define masks
#ifdef AVX512
    // batch size = 2
    const __mmask16 mask2  = 0x0003;  // 0000 0000 0000 0011

    // batch size = 4
    const __mmask16 mask4  = 0x000F;  // 0000 0000 0000 1111

    // batch size = 8
    const __mmask16 mask8  = 0x00FF;  // 0000 0000 1111 1111

    // batch size = 12
    const __mmask16 mask12 = 0x0FFF;  // 0000 1111 1111 1111

    // batch size = 16 (full)
    const __mmask16 mask16 = 0xFFFF;  // 1111 1111 1111 1111

    const __mmask16 get_mask(size_t batch_size){

        __mmask16 m;

        switch (batch_size) {

            case 2:  
                m = mask2;  
                break;
            case 4:  
                m = mask4;  
                break;
            case 8:  
                m = mask8;  
                break;
            case 12: 
                m = mask12; 
                break;
            default: 
                m = mask16; 
                break;
        }

        return m;
    }
#endif



void simulate_samples(spiking_nn_t *snn, simulation_configuration_t *conf, simulation_results_t *results, input_data_t *dataset){

    // initialize several control variables
    int epochs, n_samples, time_steps, batch_size;
    spiking_nn_t *pr_snns;
    spiking_nn_t *tmp_snn;
    

    /*#if defined PAR_SAMPLES || defined NESTED
    if(conf->n_process > conf->batch_size)
        conf->n_process = conf->batch_size;
    #endif*/

    // cpy important information
    epochs = conf->epochs;
    n_samples = dataset->n_samples;
    time_steps = conf->time_steps;
    batch_size = conf->batch_size; // batch size
    results->elapsed_time_total_inference = 0.0;
    results->elapsed_time_total_learning = 0.0;

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
    struct timespec start_total_learning, end_total_learning; // to measure learning processing time
    struct timespec start_total_inference, end_total_inference; // to measure learning processing time

    // create array for batches
    int n_batches;
    int extra_samples;
    int *batches = create_batches(&n_batches, &extra_samples, dataset, n_samples, batch_size);

    // TODO: revise this
    simulation_results_per_sample_t *results_per_sample;



    // configure the OpenMP parameters depending on parallelization strategy
    int p_outer = 0, p_inner = 0, n_outer = 0, n_inner = 0, s_inner = 0, s_rem = 0;
    #ifdef NESTED
    
        // allow nedted parallelism
        omp_set_dynamic(0); // disable dynamic threads (I do not understand this very well)
        omp_set_nested(1); // allow nested parallelism
        omp_set_max_active_levels(2);

        // number of processes
        p_inner = conf->n_inner_process;
        n_inner = snn->n_neurons / p_inner * 0.1;
        if(n_inner == 0) 
            n_inner = 1;
        
        p_outer = n_process / p_inner;
        n_outer = 1;//snn->n_neurons / p_outer * 0.1;

        // it makes no sense to have more outer processes than batch_size
        if(p_outer > batch_size)
            p_outer = batch_size;

    #elif PAR_SAMPLES

        // it makes no sense to have more outer processes than batch size
        int tmp_n_process = n_process;
        if(tmp_n_process > batch_size)
            tmp_n_process = batch_size;

        p_outer = tmp_n_process;
        n_outer = 1;//snn->n_neurons / p_outer * 0.1;

        // initialize whether it is necessary to use them
        p_inner = n_process;
        n_inner = 1;

    #else

        p_inner = n_process;
        n_inner = snn->n_neurons / p_inner * 0.1;
        if(n_inner == 0) 
            n_inner = 1;
        
        s_inner = snn->n_synapses / p_inner * 0.1;

        // initialize whether it is necessary to use them
        p_outer = n_process;
        n_outer = 1;
    
    #endif


    // if it is necessary, create network copies
    #if defined PAR_SAMPLES || defined NESTED 
        pr_snns = (spiking_nn_t *)malloc(p_outer * sizeof(spiking_nn_t));
        for(i = 0; i<p_outer; i++){
            
            cp_network(&(pr_snns[i]), snn, conf);
            
            #ifdef REORDER

                printf(" > Reordering array of synapses...\n");

                #ifdef CUDA
                // if cuda is used, initialized delays to 0
                if(conf.cuda != 0)
                    for(int i = 0; i<snn->n_input_synapses; i++){
                        pr_snns[i].synapses[i].delay = 0;
                    }
                #endif

                // reorder array
                reorder_synapse_list(&(pr_snns[i]));

                printf(" > Array of synapses reordered!\n");
                fflush(stdout);

            #endif
        }
    #endif

    #ifdef REORDER

        printf(" > Reordering array of synapses...\n");

        #ifdef CUDA
        // if cuda is used, initialized delays to 0
        if(conf->cuda != 0)
            for(int i = 0; i<snn->n_input_synapses; i++){
                snn->synapses[i].delay = 0;
            }
        #endif

        // reorder array
        reorder_synapse_list(snn);

        printf(" > Array of synapses reordered!\n");
        fflush(stdout);

    #endif

    // start measuring time
    clock_gettime(CLOCK_MONOTONIC, &start);


    // open files for results storage
    open_results_files(conf);


    // simulate epochs
    for(e = 0; e<epochs; e++){

        printf(" > Epoch %d\n", e);
        fflush(stdout);

        // shuffle samples in batch if it is necessary
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


            // loop over batch samples and accumulate weight changes in dw

            clock_gettime(CLOCK_MONOTONIC, &start_total_inference);


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
                    
                    if(conf->store_n_spikes == 1){
                        #if !defined PAR_SAMPLES || defined NESTED
                        #pragma omp parallel for num_threads(p_inner) schedule(guided, n_inner) private(i) 
                        #endif
                        for(i = 0; i<tmp_snn->n_neurons; i++){
                            results_per_sample->n_spikes_per_neuron[i] = 0;
                        }
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
                        re_initialize_synapse(&(tmp_snn->synapses[i])); // reinitialize w to init_w
                    }
                    clock_gettime(CLOCK_MONOTONIC, &end_re_synapses);

                    // load sample in network 
                    clock_gettime(CLOCK_MONOTONIC, &start_load_sample);
                    tmp_snn->load_sample(tmp_snn, &(dataset->samples[sample_index]));
                    clock_gettime(CLOCK_MONOTONIC, &end_load_sample);

                    // simulate time steps for each sample // t starts from 1 since in t=0 no neuron can spike
                    for(t = 0; t<time_steps; t++){

                        // input step
                        clock_gettime(CLOCK_MONOTONIC, &start_neurons_input);

                        #if !defined PAR_SAMPLES || defined NESTED
                        #pragma omp parallel for num_threads(p_inner) schedule(guided, n_inner) private(i) 
                        #endif
                        for(i = 0; i<tmp_snn->n_neurons; i++) // O(n)
                            tmp_snn->input_step(tmp_snn, t, i, results_per_sample, conf);

                        clock_gettime(CLOCK_MONOTONIC, &end_neurons_input);


                        // output step
                        clock_gettime(CLOCK_MONOTONIC, &start_neurons_output);
                        
                        #if !defined PAR_SAMPLES || defined NESTED
                        #pragma omp parallel for num_threads(p_inner) schedule(guided, n_inner) private(i) 
                        #endif
                        for(i = 0; i<tmp_snn->n_neurons; i++) // O(n)
                            tmp_snn->output_step(tmp_snn, t, i, results_per_sample, conf);

                        clock_gettime(CLOCK_MONOTONIC, &end_neurons_output);

                        // learning rule
                        clock_gettime(CLOCK_MONOTONIC, &start_learning);
                        
                        if(conf->learn == 1){

                            #if !defined PAR_SAMPLES || defined NESTED
                            #pragma omp parallel for num_threads(p_inner) schedule(guided, 50) private(i) 
                            #endif
                            for(i = 0; i<tmp_snn->n_synapses; i++){ // O(m)
                                //printf(" > Synapse %d, lr %d, w = %lf\n", i, tmp_snn->synapses[i].lr, tmp_snn->synapses[i].w);
                                tmp_snn->synapses[i].learning_rule(&(tmp_snn->synapses[i]), t, 3); // TODO: Change this!! 
                                //fflush(stdout);
                            }
                            //printf("\n");
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

            clock_gettime(CLOCK_MONOTONIC, &end_total_inference);
                            //printf("\n");


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
            clock_gettime(CLOCK_MONOTONIC, &start_total_learning);


            // get tmp_batch_size
            int tmp_batch_size = batch_size;

            // if it is the last batch, rest extra samples
            if(b == n_batches - 1)
                tmp_batch_size = batch_size - extra_samples;

            if(conf->learn == 1){
                
                clock_gettime(CLOCK_MONOTONIC, &start_learning);

                // if there are several network copies, compute the average dw and update weights
                #if defined PAR_SAMPLES || defined NESTED

                    // sum all the dw of the copies in the original SNN structure
                    /*#pragma omp parallel for num_threads(p_inner) schedule(guided, n_inner) private(i)
                    for(i = 0; i < snn->n_synapses; i++){
                        snn->synapses[i].dw = pr_snns[0].synapses[i].dw;
                    }
                    for(s = 1; s<p_outer; s++){
                        #pragma omp parallel for num_threads(p_inner) schedule(guided, n_inner) private(i)
                        for(i = 0; i < snn->n_synapses; i++){
                            snn->synapses[i].dw += pr_snns[s].synapses[i].dw;
                        }
                    }
                    
                    // update weights in original snn
                    #pragma omp parallel for num_threads(p_inner) schedule(guided, n_inner) private(i, s)
                    for(i = 0; i<snn->n_synapses; i++){    
                        //printf("Synapse %d, init_w = %lf, dw = %lf, b = %d, fw = %lf\n", i, snn->synapses[i].init_w, (double)(snn->synapses[i].dw / (double)tmp_batch_size), tmp_batch_size, snn->synapses[i].init_w + (double)(snn->synapses[i].dw / (double)tmp_batch_size)); 
                        snn->synapses[i].init_w += (double)(snn->synapses[i].dw / (double)tmp_batch_size); // update init_w
                        snn->synapses[i].w = snn->synapses[i].init_w; // set w to init_w
                        snn->synapses[i].dw = 0; // reinit dw

                        for(s = 0; s<p_outer; s++){
                            // update w
                            pr_snns[s].synapses[i].w = snn->synapses[i].w;
                            pr_snns[s].synapses[i].init_w = snn->synapses[i].init_w;
                            pr_snns[s].synapses[i].dw = 0;
                        }
                    }    */

                    #pragma omp parallel for num_threads(p_inner) schedule(guided, n_inner) private(i, s)
                    for(i = 0; i< snn->n_synapses; i++){

                        for(s = 0; s<p_outer; s++){

                            snn->synapses[i].dw += pr_snns[s].synapses[i].dw;
                        }

                        snn->synapses[i].init_w += (double)(snn->synapses[i].dw / (double)tmp_batch_size); // update init_w
                        snn->synapses[i].w = snn->synapses[i].init_w; // set w to init_w
                        snn->synapses[i].dw = 0; // reinit dw


                        for(s = 0; s<p_outer; s++){
                            // update w
                            pr_snns[s].synapses[i].w = snn->synapses[i].w;
                            pr_snns[s].synapses[i].init_w = snn->synapses[i].init_w;
                            pr_snns[s].synapses[i].dw = 0;
                        }
                    }
                
                #else
                
                    // no copies
                    #pragma omp parallel for num_threads(p_inner) private(i) 
                    for(i = 0; i<snn->n_synapses; i++){
                        // update w
                        //printf("Synapse %d, init_w = %lf, dw = %lf, b = %d, fw = %lf\n", i, snn->synapses[i].init_w, (double)(snn->synapses[i].dw / (double)tmp_batch_size), tmp_batch_size, snn->synapses[i].init_w + (double)(snn->synapses[i].dw / (double)tmp_batch_size)); 
                        snn->synapses[i].init_w += (double)(snn->synapses[i].dw / (double)tmp_batch_size); // update init_w
                        //printf(" dw[%d] = %f\n", i, snn->synapses[i].dw / (double)tmp_batch_size);
                        snn->synapses[i].w = snn->synapses[i].init_w; // set w to init_w
                        snn->synapses[i].dw = 0; // reinit dw
                    }
                    //printf("\n\n");
                

                #endif
                clock_gettime(CLOCK_MONOTONIC, &end_total_learning);

                clock_gettime(CLOCK_MONOTONIC, &end_learning);

                results->elapsed_time_total_inference += (end_total_inference.tv_sec - start_total_inference.tv_sec) + (end_total_inference.tv_nsec - start_total_inference.tv_nsec) / 1e9;
                results->elapsed_time_total_learning += (end_total_learning.tv_sec - start_total_learning.tv_sec) + (end_total_learning.tv_nsec - start_total_learning.tv_nsec) / 1e9;

                results->elapsed_time_learning += (end_learning.tv_sec - start_learning.tv_sec) + (end_learning.tv_nsec - start_learning.tv_nsec) / 1e9;
            }


            // sum all execution times of the batch
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
            results->elapsed_time_neurons_input = (double)(results->elapsed_time_neurons_input / (double)(tmp_batch_size)); // elapsed time processing neurons input 
            results->elapsed_time_neurons_output = (double)(results->elapsed_time_neurons_output / (double)(tmp_batch_size)); // elapsed time processing neurons output
            results->elapsed_time_learning = (double)(results->elapsed_time_learning / (double)(tmp_batch_size)); // elapsed time processing learning rules
            results->elapsed_time_sample = (double)(results->elapsed_time_sample / (double)(tmp_batch_size));
            results->elapsed_time_load_sample = (double)(results->elapsed_time_load_sample / (double)(tmp_batch_size));
            results->elapsed_time_re_neurons = (double)(results->elapsed_time_re_neurons / (double)(tmp_batch_size));
            results->elapsed_time_re_synapses = (double)(results->elapsed_time_re_synapses / (double)(tmp_batch_size));


            // batch finished
            clock_gettime(CLOCK_MONOTONIC, &end_epoch);
            results->elapsed_time_epoch = (end_epoch.tv_sec - start_epoch.tv_sec) + (end_epoch.tv_nsec - start_epoch.tv_nsec) / 1e9;
    
            // store results of the samples simulated in the batch
            store_results(results, conf, snn, dataset, batches, b);
        }
        
    }

    clock_gettime(CLOCK_MONOTONIC, &end);

    results->elapsed_time = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    printf(" > Total elapsed time %lf (%lf, %lf)\n", results->elapsed_time, results->elapsed_time_total_inference, results->elapsed_time_total_learning);

    // store total execution time
    if(conf->store_execution_times == 1){
        fprintf(conf->f_et, "%lf %lf %lf\n", results->elapsed_time, results->elapsed_time_total_inference, results->elapsed_time_total_learning);
    }

    // close results files
    close_results_files(conf);
}


void simulate_samples_synapses(spiking_nn_t *snn, simulation_configuration_t *conf, simulation_results_t *results, input_data_t *dataset){


    // initialize several control variables
    int epochs, n_samples, time_steps, batch_size;
    spiking_nn_t *pr_snns;
    spiking_nn_t *tmp_snn;
    

    /*#if defined PAR_SAMPLES || defined NESTED
    if(conf->n_process > conf->batch_size)
        conf->n_process = conf->batch_size;
    #endif*/

    // cpy important information
    epochs = conf->epochs;
    n_samples = dataset->n_samples;
    time_steps = conf->time_steps;
    batch_size = conf->batch_size; // batch size
    results->elapsed_time_total_inference = 0.0;
    results->elapsed_time_total_learning = 0.0;

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
    struct timespec start_total_learning, end_total_learning; // to measure learning processing time
    struct timespec start_total_inference, end_total_inference; // to measure learning processing time

    // create array for batches
    int n_batches;
    int extra_samples;
    int *batches = create_batches(&n_batches, &extra_samples, dataset, n_samples, batch_size);

    // TODO: revise this
    simulation_results_per_sample_t *results_per_sample;



    // configure the OpenMP parameters depending on parallelization strategy
    int p_outer = 0, p_inner = 0, n_outer = 0, n_inner = 0, s_inner = 0, s_rem = 0;
    #ifdef NESTED
    
        // allow nedted parallelism
        omp_set_dynamic(0); // disable dynamic threads (I do not understand this very well)
        omp_set_nested(1); // allow nested parallelism
        omp_set_max_active_levels(2);

        // number of processes
        p_inner = conf->n_inner_process;
        n_inner = snn->n_neurons / p_inner * 0.1;
        if(n_inner == 0) 
            n_inner = 1;
        
        p_outer = n_process / p_inner;
        n_outer = 1;//snn->n_neurons / p_outer * 0.1;

        // it makes no sense to have more outer processes than batch_size
        if(p_outer > batch_size)
            p_outer = batch_size;

    #elif PAR_SAMPLES

        // it makes no sense to have more outer processes than batch size
        int tmp_n_process = n_process;
        if(tmp_n_process > batch_size)
            tmp_n_process = batch_size;

        p_outer = tmp_n_process;
        n_outer = 1;//snn->n_neurons / p_outer * 0.1;

        // initialize whether it is necessary to use them
        p_inner = n_process;
        n_inner = 1;

    #else

        p_inner = n_process;
        n_inner = snn->n_neurons / p_inner * 0.1;
        if(n_inner == 0) 
            n_inner = 1;
        
        s_inner = snn->n_synapses / p_inner * 0.1;

        // initialize whether it is necessary to use them
        p_outer = n_process;
        n_outer = 1;
    
    #endif


    // if it is necessary, create network copies
    #if defined PAR_SAMPLES || defined NESTED 
        pr_snns = (spiking_nn_t *)malloc(p_outer * sizeof(spiking_nn_t));
        for(i = 0; i<p_outer; i++){
            
            cp_network(&(pr_snns[i]), snn, conf);
            
            #ifdef REORDER

                printf(" > Reordering array of synapses...\n");

                #ifdef CUDA
                // if cuda is used, initialized delays to 0
                if(conf.cuda != 0)
                    for(int i = 0; i<snn->n_input_synapses; i++){
                        pr_snns[i].synapses[i].delay = 0;
                    }
                #endif

                // reorder array
                reorder_synapse_list(&(pr_snns[i]));

                printf(" > Array of synapses reordered!\n");
                fflush(stdout);

            #endif
        }
    #endif

    #ifdef REORDER

        printf(" > Reordering array of synapses...\n");

        #ifdef CUDA
        // if cuda is used, initialized delays to 0
        if(conf->cuda != 0)
            for(int i = 0; i<snn->n_input_synapses; i++){
                snn->synapses[i].delay = 0;
            }
        #endif

        // reorder array
        reorder_synapse_list(snn);

        printf(" > Array of synapses reordered!\n");
        fflush(stdout);

    #endif

    // start measuring time
    clock_gettime(CLOCK_MONOTONIC, &start);


    // open files for results storage
    open_results_files(conf);


    // simulate epochs
    for(e = 0; e<epochs; e++){

        printf(" > Epoch %d\n", e);
        fflush(stdout);

        // shuffle samples in batch if it is necessary
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


            // loop over batch samples and accumulate weight changes in dw

            clock_gettime(CLOCK_MONOTONIC, &start_total_inference);


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
                    
                    if(conf->store_n_spikes == 1){
                        #if !defined PAR_SAMPLES || defined NESTED
                        #pragma omp parallel for num_threads(p_inner) schedule(guided, n_inner) private(i) 
                        #endif
                        for(i = 0; i<tmp_snn->n_neurons; i++){
                            results_per_sample->n_spikes_per_neuron[i] = 0;
                        }
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
                        re_initialize_synapse(&(tmp_snn->synapses[i])); // reinitialize w to init_w
                    }
                    clock_gettime(CLOCK_MONOTONIC, &end_re_synapses);

                    // load sample in network 
                    clock_gettime(CLOCK_MONOTONIC, &start_load_sample);
                    tmp_snn->load_sample(tmp_snn, &(dataset->samples[sample_index]));
                    clock_gettime(CLOCK_MONOTONIC, &end_load_sample);

                    // simulate time steps for each sample // t starts from 1 since in t=0 no neuron can spike
                    for(t = 0; t<time_steps; t++){

                        // input step
                        clock_gettime(CLOCK_MONOTONIC, &start_neurons_input);

                        #if !defined PAR_SAMPLES || defined NESTED
                        #pragma omp parallel for num_threads(p_inner) schedule(guided, n_inner) private(i) 
                        #endif
                        for(i = 0; i<tmp_snn->n_neurons; i++) // O(n)
                            tmp_snn->input_step(tmp_snn, t, i, results_per_sample, conf);

                        clock_gettime(CLOCK_MONOTONIC, &end_neurons_input);


                        // output step
                        clock_gettime(CLOCK_MONOTONIC, &start_neurons_output);
                        
                        #if !defined PAR_SAMPLES || defined NESTED
                        #pragma omp parallel for num_threads(p_inner) schedule(guided, n_inner) private(i) 
                        #endif
                        for(i = 0; i<tmp_snn->n_neurons; i++) // O(n)
                            tmp_snn->output_step(tmp_snn, t, i, results_per_sample, conf);

                        clock_gettime(CLOCK_MONOTONIC, &end_neurons_output);

                        // learning rule
                        clock_gettime(CLOCK_MONOTONIC, &start_learning);
                        
                        if(conf->learn == 1){

                            #if !defined PAR_SAMPLES || defined NESTED
                            #pragma omp parallel for num_threads(p_inner) schedule(guided, 50) private(i) 
                            #endif
                            for(i = 0; i<tmp_snn->n_synapses; i++){ // O(m)
                                //printf(" > Synapse %d, lr %d, w = %lf\n", i, tmp_snn->synapses[i].lr, tmp_snn->synapses[i].w);
                                tmp_snn->synapses[i].learning_rule(&(tmp_snn->synapses[i]), t, 3); // TODO: Change this!! 
                                //fflush(stdout);
                            }
                            //printf("\n");
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

            clock_gettime(CLOCK_MONOTONIC, &end_total_inference);
                            //printf("\n");


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
            clock_gettime(CLOCK_MONOTONIC, &start_total_learning);


            // get tmp_batch_size
            int tmp_batch_size = batch_size;

            // if it is the last batch, rest extra samples
            if(b == n_batches - 1)
                tmp_batch_size = batch_size - extra_samples;

            if(conf->learn == 1){
                
                clock_gettime(CLOCK_MONOTONIC, &start_learning);

                // if there are several network copies, compute the average dw and update weights
                #if defined PAR_SAMPLES || defined NESTED

                    // sum all the dw of the copies in the original SNN structure
                    /*#pragma omp parallel for num_threads(p_inner) schedule(guided, n_inner) private(i)
                    for(i = 0; i < snn->n_synapses; i++){
                        snn->synapses[i].dw = pr_snns[0].synapses[i].dw;
                    }
                    for(s = 1; s<p_outer; s++){
                        #pragma omp parallel for num_threads(p_inner) schedule(guided, n_inner) private(i)
                        for(i = 0; i < snn->n_synapses; i++){
                            snn->synapses[i].dw += pr_snns[s].synapses[i].dw;
                        }
                    }
                    
                    // update weights in original snn
                    #pragma omp parallel for num_threads(p_inner) schedule(guided, n_inner) private(i, s)
                    for(i = 0; i<snn->n_synapses; i++){    
                        //printf("Synapse %d, init_w = %lf, dw = %lf, b = %d, fw = %lf\n", i, snn->synapses[i].init_w, (double)(snn->synapses[i].dw / (double)tmp_batch_size), tmp_batch_size, snn->synapses[i].init_w + (double)(snn->synapses[i].dw / (double)tmp_batch_size)); 
                        snn->synapses[i].init_w += (double)(snn->synapses[i].dw / (double)tmp_batch_size); // update init_w
                        snn->synapses[i].w = snn->synapses[i].init_w; // set w to init_w
                        snn->synapses[i].dw = 0; // reinit dw

                        for(s = 0; s<p_outer; s++){
                            // update w
                            pr_snns[s].synapses[i].w = snn->synapses[i].w;
                            pr_snns[s].synapses[i].init_w = snn->synapses[i].init_w;
                            pr_snns[s].synapses[i].dw = 0;
                        }
                    }    */

                    #pragma omp parallel for num_threads(p_inner) schedule(guided, n_inner) private(i, s)
                    for(i = 0; i< snn->n_synapses; i++){

                        for(s = 0; s<p_outer; s++){

                            snn->synapses[i].dw += pr_snns[s].synapses[i].dw;
                        }

                        snn->synapses[i].init_w += (double)(snn->synapses[i].dw / (double)tmp_batch_size); // update init_w
                        snn->synapses[i].w = snn->synapses[i].init_w; // set w to init_w
                        snn->synapses[i].dw = 0; // reinit dw


                        for(s = 0; s<p_outer; s++){
                            // update w
                            pr_snns[s].synapses[i].w = snn->synapses[i].w;
                            pr_snns[s].synapses[i].init_w = snn->synapses[i].init_w;
                            pr_snns[s].synapses[i].dw = 0;
                        }
                    }
                
                #else
                
                    // no copies
                    #pragma omp parallel for num_threads(p_inner) private(i) 
                    for(i = 0; i<snn->n_synapses; i++){
                        // update w
                        //printf("Synapse %d, init_w = %lf, dw = %lf, b = %d, fw = %lf\n", i, snn->synapses[i].init_w, (double)(snn->synapses[i].dw / (double)tmp_batch_size), tmp_batch_size, snn->synapses[i].init_w + (double)(snn->synapses[i].dw / (double)tmp_batch_size)); 
                        snn->synapses[i].init_w += (double)(snn->synapses[i].dw / (double)tmp_batch_size); // update init_w
                        //printf(" dw[%d] = %f\n", i, snn->synapses[i].dw / (double)tmp_batch_size);
                        snn->synapses[i].w = snn->synapses[i].init_w; // set w to init_w
                        snn->synapses[i].dw = 0; // reinit dw
                    }
                    //printf("\n\n");
                

                #endif
                clock_gettime(CLOCK_MONOTONIC, &end_total_learning);

                clock_gettime(CLOCK_MONOTONIC, &end_learning);

                results->elapsed_time_total_inference += (end_total_inference.tv_sec - start_total_inference.tv_sec) + (end_total_inference.tv_nsec - start_total_inference.tv_nsec) / 1e9;
                results->elapsed_time_total_learning += (end_total_learning.tv_sec - start_total_learning.tv_sec) + (end_total_learning.tv_nsec - start_total_learning.tv_nsec) / 1e9;

                results->elapsed_time_learning += (end_learning.tv_sec - start_learning.tv_sec) + (end_learning.tv_nsec - start_learning.tv_nsec) / 1e9;
            }


            // sum all execution times of the batch
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
            results->elapsed_time_neurons_input = (double)(results->elapsed_time_neurons_input / (double)(tmp_batch_size)); // elapsed time processing neurons input 
            results->elapsed_time_neurons_output = (double)(results->elapsed_time_neurons_output / (double)(tmp_batch_size)); // elapsed time processing neurons output
            results->elapsed_time_learning = (double)(results->elapsed_time_learning / (double)(tmp_batch_size)); // elapsed time processing learning rules
            results->elapsed_time_sample = (double)(results->elapsed_time_sample / (double)(tmp_batch_size));
            results->elapsed_time_load_sample = (double)(results->elapsed_time_load_sample / (double)(tmp_batch_size));
            results->elapsed_time_re_neurons = (double)(results->elapsed_time_re_neurons / (double)(tmp_batch_size));
            results->elapsed_time_re_synapses = (double)(results->elapsed_time_re_synapses / (double)(tmp_batch_size));


            // batch finished
            clock_gettime(CLOCK_MONOTONIC, &end_epoch);
            results->elapsed_time_epoch = (end_epoch.tv_sec - start_epoch.tv_sec) + (end_epoch.tv_nsec - start_epoch.tv_nsec) / 1e9;
    
            // store results of the samples simulated in the batch
            store_results(results, conf, snn, dataset, batches, b);
        }
        
    }

    clock_gettime(CLOCK_MONOTONIC, &end);

    results->elapsed_time = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    printf(" > Total elapsed time %lf (%lf, %lf)\n", results->elapsed_time, results->elapsed_time_total_inference, results->elapsed_time_total_learning);

    // store total execution time
    if(conf->store_execution_times == 1){
        fprintf(conf->f_et, "%lf %lf %lf\n", results->elapsed_time, results->elapsed_time_total_inference, results->elapsed_time_total_learning);
    }

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



void large_simulation(GPU_SNN_t *cpu_snn, GPU_dataset_t *cpu_dataset, simulation_configuration_t *conf){
    
    
    struct timespec start, end; // to measure simulation complete time
    struct timespec start_step1, end_step1; // to measure simulation complete time
    struct timespec start_step2, end_step2; // to measure simulation complete time
    struct timespec start_step3, end_step3; // to measure simulation complete time
    struct timespec start_step4, end_step4; // to measure simulation complete time
    struct timespec start_step5, end_step5; // to measure simulation complete time

    double elapsed_time, et1 =0.0, et2=0.0, et3=0.0, et4=0.0, et5=0.0;

    const float alpha = 0.95f;
    const float beta = 0.05f;

    // start timing
    clock_gettime(CLOCK_MONOTONIC, &start);

    size_t s, sidx;
    int *n_spikes = (int*)calloc(cpu_snn->n_neurons, sizeof(int));
    float *arrI = (float*)calloc(cpu_snn->n_neurons, sizeof(float));
    int *inR = (int*)calloc(cpu_snn->n_neurons, sizeof(int));

    /* Restricts for vectorization */
    float * restrict rV      = cpu_snn->v;
    float * restrict rRest = cpu_snn->v_rest;
    float * restrict rRef = inR;
    float * restrict rI      = arrI;


    for(s = 0; s<(size_t)(cpu_dataset->n_samples); s = s + (size_t)(conf->batch_size)){

        printf(" Simulating sample %lu\n", s);
        fflush(stdout);
        for(sidx = 0; sidx < (size_t)(cpu_dataset->n_samples); sidx++){

            printf(" Simulating sample %lu\n", sidx);
            fflush(stdout);

            clock_gettime(CLOCK_MONOTONIC, &start_step1);

            // reinitialize spike matrix
            #pragma omp parallel for num_threads(conf->n_process)
            for(size_t i=0; i<(size_t)(cpu_snn->n_neurons + cpu_snn->n_input_neurons); i++){
                for(size_t j = 0; j<(size_t)(cpu_snn->max_spikes); j++){
                    
                    cpu_snn->spk_matrix[(size_t)(cpu_snn->max_spikes) * i + j] = 0;
                }
            }

            // reinitialize neurons
            #pragma omp parallel for num_threads(conf->n_process)
            for(size_t i = 0; i<(size_t)(cpu_snn->n_neurons); i++){

                cpu_snn->v[i] = cpu_snn->v_rest[i];
                cpu_snn->r_period_remain[i] = 0;
            }


            // reinitialize synapses 
            #pragma omp parallel for num_threads(conf->n_process)
            for(size_t i=0; i<(size_t)(cpu_snn->n_synapses); i++){

                cpu_snn->w[i] = cpu_snn->init_w[i];
                //cpu_snn->next_pre_spike[i] = 0;
                //cpu_snn->next_post_spike[i] = 0;
            }

            // load the sample
            #pragma omp parallel for num_threads(conf->n_process)
            for(size_t i=0; i<(size_t)(cpu_snn->n_input_neurons); i++){

                if(sidx < cpu_dataset->n_samples){

                    size_t fidx = i;

                    size_t start_sample = cpu_dataset->sample_offset[sidx];
                    size_t start_feature = cpu_dataset->feature_offset[sidx * (size_t)(cpu_dataset->n_features) + fidx];
                    size_t n_spikes_per_feature = (size_t)(cpu_dataset->n_spikes_per_feature[sidx * (size_t)(cpu_dataset->n_features) + fidx]);
                
                    for(size_t j = 0; j<n_spikes_per_feature; j++){

                        // matrix [N * T]
                        //cpu_snn->spk_matrix[(size_t)(cpu_snn->max_spikes) * i + (size_t)(cpu_dataset->spikes[start_feature + j])] = 1;//cpu_dataset->spikes[start_feature + j];
                        
                        // matrix [T * N]
                        cpu_snn->spk_matrix[(size_t)(cpu_dataset->spikes[start_feature + j]) * (size_t)(cpu_snn->n_neurons + cpu_snn->n_input_neurons) + i] = 1;
                    }
                }
            }

            clock_gettime(CLOCK_MONOTONIC, &end_step1);
            et1+=(end_step1.tv_sec - start_step1.tv_sec) + (end_step1.tv_nsec - start_step1.tv_nsec) / 1e9;


            // simulate the sample 
            for(size_t t=0; t<(size_t)(conf->time_steps); t++){
                
                clock_gettime(CLOCK_MONOTONIC, &start_step2);


                // process neurons
                #pragma omp parallel for num_threads(conf->n_process) 
                for(size_t i=0; i<(size_t)(cpu_snn->n_neurons); i++){

                    //printf(" > > Neuron %lu: \n", i);
                    //fflush(stdout);
                    float I = 0.0;
                    
                    // convert conditional to multiplication

                    if(cpu_snn->r_period_remain[i] <= 0){

                        size_t synapse_index, in_neuron_index, spk_time_index;
                        int delay, spk_time; 
                        float spk;    

                        size_t base_synapse = (size_t)(cpu_snn->neuron_input_synapses_offset[i]);

                        // correction here: I am probably suming more input synapses than only those of the neuron
                        size_t j;

                        #if defined AVX512
                        for(j = 0; j + 15 < (size_t)(cpu_snn->n_neuron_input_synapses[i]); j+=16){

                            __m512i sInd_vec    = _mm512_loadu_si512(&(cpu_snn->neuron_input_synapses_offset[base_synapse + j]));
                            __m512i d_vec = _mm512_loadu_si512(&(cpu_snn->delay[base_synapse + j]));                           
                            __m512 w_vec = _mm512_loadu_ps(&(cpu_snn->w[base_synapse + j]));    
                            
                            // get spike times vector
                            __m512i spk_time_vec = _mm512_sub_epi32(_mm512_set1_epi32(t), d_vec);
                            __m512i spk_time_clamped =_mm512_max_epi32(spk_time_vec, _mm512_setzero_si512());


                            // valid mask: lower or equal to 0
                            __mmask16 valid_mask = _mm512_cmp_epi32_mask(spk_time_vec, _mm512_setzero_si512(), _MM_CMPINT_GE);

                            // gather in_neuron_index (still indirect)
                            __m512i in_vec = _mm512_loadu_si512(&(cpu_snn->pre_neuron_index[base_synapse + j]));  

                            // compute offsets into spk_matrix
                            // matrix [N * T]
                            //__m512i idx = _mm512_add_epi32(_mm512_mullo_epi32(in_vec, _mm512_set1_epi32(cpu_snn->max_spikes)), spk_time_clamped);

                            // matrix [T * N]
                            __m512i idx = _mm512_add_epi32(_mm512_mullo_epi32(spk_time_clamped, _mm512_add_epi32(_mm512_set1_epi32(cpu_snn->n_neurons), _mm512_set1_epi32(cpu_snn->n_input_neurons))), in_vec);


                            // gather spikes with mask
                            __m512i spikes_i = _mm512_mask_i32gather_epi32(_mm512_setzero_si512(), valid_mask, idx, cpu_snn->spk_matrix, 4);
                            __m512 spikes = _mm512_cvtepi32_ps(spikes_i);
                            
                            // multiply by weights
                            __m512 contrib = _mm512_mul_ps(w_vec, spikes);

                            // horizontal sum into scalar I
                            float I_vec = _mm512_reduce_add_ps(contrib);

                            I += I_vec;
                        }
                        #endif
                        
                        for(; j<(size_t)(cpu_snn->n_neuron_input_synapses[i]); j++){

                            synapse_index = base_synapse + j; // no copies
                            in_neuron_index = (size_t)(cpu_snn->pre_neuron_index[synapse_index]); // absurd, equal to 0
                            delay = cpu_snn->delay[synapse_index]; // no copies
                            spk_time = t - delay; // actual position

                            //int valid = (spk_time >=0); // think how to eliminate this if (larger matrix?)
                            if(spk_time >=0){
                                // [N * T]
                                //spk = (float)(cpu_snn->spk_matrix[in_neuron_index * (size_t)(cpu_snn->max_spikes) + (size_t)(spk_time)]);    
                                
                                // T * N
                                spk = (float)(cpu_snn->spk_matrix[in_neuron_index + (size_t)(spk_time)  * (size_t)(cpu_snn->n_neurons + cpu_snn->n_input_neurons)]);    
                                
                                I += cpu_snn->w[synapse_index] * spk; // spk = 0 / 1
                            }
                            //cpu_snn->v[i] = (1.0-0.05) * cpu_snn->v[i] + cpu_snn->v_rest[i] * 0.05 + I;
                        }
                    }
                    arrI[i] = I;
                    inR[i] = (cpu_snn->r_period_remain[i] <= 0);
                }
                clock_gettime(CLOCK_MONOTONIC, &end_step2);
                et2+=(end_step2.tv_sec - start_step2.tv_sec) + (end_step2.tv_nsec - start_step2.tv_nsec) / 1e9;
                
                clock_gettime(CLOCK_MONOTONIC, &start_step3);

                                size_t i=0;
                // constants in registers
                #if defined AVX512
                __m512 alpha_v = _mm512_set1_ps(alpha);
                __m512 beta_v  = _mm512_set1_ps(beta);

                for (i = 0; i + 15 < cpu_snn->n_neurons; i += 16) {
                    // load 16 neurons
                    __m512 v_vec    = _mm512_loadu_ps(&rV[i]);
                    __m512 vrest_vec= _mm512_loadu_ps(&rRest[i]);
                    __m512 I_vec    = _mm512_loadu_ps(&rI[i]);
                    __m512 ref_vec  = _mm512_loadu_ps(&rRef[i]);

                    // compare ref <= 0 → mask
                    __mmask16 mask = _mm512_cmp_ps_mask(ref_vec, _mm512_setzero_ps(), _CMP_LE_OS);

                    // compute newv
                    //__m512 newv = _mm512_add_ps(_mm512_mul_ps(alpha_v, v_vec),
                    //                            _mm512_add_ps(_mm512_mul_ps(beta_v, vrest_vec), I_vec));
                    // compute new values (vectorized)
                    __m512 newv = _mm512_add_ps(
                            _mm512_mul_ps(alpha_v, v_vec),
                            _mm512_add_ps(_mm512_mul_ps(beta_v, vrest_vec), I_vec)
                        );


                    // combine old & new: v = mask ? newv : v_old
                    __m512 v_combined = _mm512_mask_mov_ps(v_vec, mask, newv);

                    // masked store: only write to v[i] where mask is true
                    //_mm512_mask_storeu_ps(&rV[i], mask, newv);
                    // store the full vector
                    _mm512_storeu_ps(&rV[i], v_combined);
                
                }
                #endif
                // handle remaining neurons
                for (; i < (size_t)(cpu_snn->n_neurons); i++) {
                    if (rRef[i] == 1)
                        rV[i] = alpha * rV[i] + beta * rRest[i] + rI[i];
                }

                clock_gettime(CLOCK_MONOTONIC, &end_step3);
                et3+=(end_step3.tv_sec - start_step3.tv_sec) + (end_step3.tv_nsec - start_step3.tv_nsec) / 1e9;

                clock_gettime(CLOCK_MONOTONIC, &start_step4);
      
                // process neuron output // vectorize?
                
                
                // N * T



                // T * N without vectorization
                /*size_t row = (size_t)(cpu_snn->n_input_neurons + cpu_snn->n_neurons) * (size_t)(t) + (size_t)(cpu_snn->n_input_neurons);

                #pragma omp parallel for num_threads(conf.n_process)
                for(size_t i = 0; i<(size_t)(cpu_snn->n_neurons); i++){
                    
                    cpu_snn->r_period_remain[i] --;
                    if(cpu_snn->v[i] >= cpu_snn->v_thresh[i]){
                        
                        // N * T
                        //size_t idx = ((size_t)(cpu_snn->n_input_neurons) + i) * (size_t)(cpu_snn->max_spikes) + (size_t)(t); 
                        // T * N
                        size_t idx = row + i; // row start at


                        //global_index = matrix_index + row + col;
                        // add spike to matrix // global memory write, bad in GPU :(
                        cpu_snn->spk_matrix[idx] = 1;

                        // reinit neuron values
                        cpu_snn->r_period_remain[i] = cpu_snn->r_period[i];
                        
                        cpu_snn->v[i] = cpu_snn->v_rest[i]; // reinit v_rest

                        n_spikes[i] += 1;
                    }
                }*/


                // T* N vectorized
                size_t idx = (size_t)(cpu_snn->n_input_neurons + cpu_snn->n_neurons) * (size_t)(t) + (size_t)(cpu_snn->n_input_neurons);
                i=0;

                #if defined AVX512
                __m512i ones_vec = _mm512_set1_epi32(1);
                for(i = 0; i + 15 < cpu_snn->n_neurons; i+=16){

                    // update refractory period
                    __m512i ref_vec = _mm512_loadu_si512(&cpu_snn->r_period_remain[i]);
                    ref_vec = _mm512_sub_epi32(ref_vec, ones_vec);
                    _mm512_storeu_si512(&cpu_snn->r_period_remain[i], ref_vec);

                    // load voltages and thresholds
                    __m512 v_vec   = _mm512_loadu_ps(&cpu_snn->v[i]);
                    __m512 vth_vec = _mm512_loadu_ps(&cpu_snn->v_thresh[i]);

                    // firing mask
                    __mmask16 fire_mask = _mm512_cmp_ps_mask(v_vec, vth_vec, _CMP_GE_OS);


                    __m512i spike_vec =
                        _mm512_mask_mov_epi32(_mm512_setzero_si512(),
                                            fire_mask,
                                            _mm512_set1_epi32(1));  
                    _mm512_storeu_si512(&cpu_snn->spk_matrix[idx + i], spike_vec);
                    
                    // reset refractory where fired
                    __m512i ref_reset_vec = _mm512_loadu_si512(&cpu_snn->r_period[i]);
                    ref_vec = _mm512_mask_mov_epi32(ref_vec, fire_mask, ref_reset_vec);
                    _mm512_storeu_si512(&cpu_snn->r_period_remain[i], ref_vec);

                    // reset voltage where fired 
                    __m512 vrest_vec = _mm512_loadu_ps(&cpu_snn->v_rest[i]);
                    v_vec = _mm512_mask_mov_ps(v_vec, fire_mask, vrest_vec);
                    _mm512_storeu_ps(&cpu_snn->v[i], v_vec);

                    // increment spike count
                    __m512i nspk_vec = _mm512_loadu_si512(&n_spikes[i]);
                    nspk_vec = _mm512_mask_add_epi32(nspk_vec, fire_mask, nspk_vec, ones_vec);
                    _mm512_storeu_si512(&n_spikes[i], nspk_vec);

                }
                #endif
                //#pragma omp parallel for num_threads(conf.n_process)
                for(; i<(size_t)(cpu_snn->n_neurons); i++){
                    
                    cpu_snn->r_period_remain[i] --;
                    if(cpu_snn->v[i] >= cpu_snn->v_thresh[i]){
                        
                        // N * T
                        //size_t idx = ((size_t)(cpu_snn->n_input_neurons) + i) * (size_t)(cpu_snn->max_spikes) + (size_t)(t); 
                        // T * N
                        size_t idx_neuron = idx + i; // row start at


                        //global_index = matrix_index + row + col;
                        // add spike to matrix // global memory write, bad in GPU :(
                        cpu_snn->spk_matrix[idx_neuron] = 1;

                        // reinit neuron values
                        cpu_snn->r_period_remain[i] = cpu_snn->r_period[i];
                        
                        cpu_snn->v[i] = cpu_snn->v_rest[i]; // reinit v_rest

                        n_spikes[i] += 1;
                    }
                }

                /*#pragma omp parallel for num_threads(conf.n_process)
                for(size_t i = 0; i<(size_t)(cpu_snn->n_neurons); i++){
                    
                    cpu_snn->r_period_remain[i] --;
                    if(cpu_snn->v[i] >= cpu_snn->v_thresh[i]){
                        
                        // N * T
                        //size_t idx = ((size_t)(cpu_snn->n_input_neurons) + i) * (size_t)(cpu_snn->max_spikes) + (size_t)(t); 
                        // T * N
                        size_t idx = row + i; // row start at


                        //global_index = matrix_index + row + col;
                        // add spike to matrix // global memory write, bad in GPU :(
                        cpu_snn->spk_matrix[idx] = 1;

                        // reinit neuron values
                        cpu_snn->r_period_remain[i] = cpu_snn->r_period[i];
                        
                        cpu_snn->v[i] = cpu_snn->v_rest[i]; // reinit v_rest

                        n_spikes[i] += 1;
                    }
                }*/

                clock_gettime(CLOCK_MONOTONIC, &end_step4);
                et4+=(end_step4.tv_sec - start_step4.tv_sec) + (end_step4.tv_nsec - start_step4.tv_nsec) / 1e9;

            }
        }

    }
}








/* NEW */

void print_spk_matrix_2D(GPU_SNN_t *snn){
  
    size_t iN, N, LT;

    iN = snn->n_input_neurons;
    N = snn->n_neurons;

    LT = snn->LT;

    for(size_t t = 0; t<LT; t++){

        for(size_t i = 0; i<iN + N; i++){

            printf("%d ", snn->spk_matrix[(iN + N) * t + i]);
        }
        printf("\n");
    }
    printf("\n");
}

void print_spk_matrix_3D(GPU_SNN_t *snn, size_t B){
  
    size_t iN, N, LT;
    size_t t, i, b;

    iN = snn->n_input_neurons;
    N = snn->n_neurons;
    LT = snn->LT;
    
    for(t = 0; t<LT; t++){

        for(i = 0; i<iN + N - 1; i++){

            printf("[");
            for(b = 0; b < B-1; b++){
                
                printf("%d ", snn->spk_matrix[(iN + N) * B * t + i * B + b]);
            }
            printf("%d], ", snn->spk_matrix[(iN + N) * B * t + i * B + b]);
        }

        printf("[");
        for(b = 0; b < B-1; b++){
            
            printf("%d ", snn->spk_matrix[(iN + N) * B * t + i * B + b]);
        }
        printf("%d]", snn->spk_matrix[(iN + N) * B * t + i * B + b]);
        
        printf("\n");
    }
    printf("\n");
}


void print_weights_batch(GPU_SNN_t *snn, size_t B, size_t b){
  
    size_t S, LT;
    size_t t, i;

    S = snn->n_synapses;
    LT = snn->LT;
    
    for(i = 0; i < S; i++){

        //printf("[");
        //for(b = 0; b < B-1; b++){
            
            printf("%f ", snn->w[i * B + b]);
        //}
        //printf("%f], ", snn->w[i  * B + B - 1]);
    }
    printf("\n");
}

void print_weights_3D(GPU_SNN_t *snn, size_t B){
  
    size_t S, LT;
    size_t t, i, b;

    S = snn->n_synapses;
    LT = snn->LT;
    
    printf(" W = ");
    for(i = 0; i < S; i++){

        printf("[");
        for(b = 0; b < B-1; b++){
            
            printf("%f ", snn->w[i * B + b]);
        }
        printf("%f], ", snn->w[i  * B + B - 1]);
    }
    printf("\n");
}

void print_dw_3D(GPU_SNN_t *snn, size_t B){
  
    size_t S, LT;
    size_t t, i, b;

    S = snn->n_synapses;
    LT = snn->LT;
    
    printf(" dW = ");
    for(i = 0; i < S; i++){

        printf("[");
        for(b = 0; b < B-1; b++){
            
            printf("%f ", snn->dw[i * B + b]);
        }
        printf("%f], ", snn->dw[i  * B + B - 1]);
    }
    printf("\n");
}


void update_weights_cpu(GPU_SNN_t *snn, simulation_configuration_t *conf){

    size_t S, P, B;
    size_t i, b;

    S = snn->n_synapses;
    B = conf->batch_size;
    P = conf->n_process;

    #if defined AVX512
    {

        __m512 vZero = _mm512_set1_ps(0.0);

        if(B == 8){
            
            __mmask16 m = 0xFF; // only 256 bits used

            #pragma omp parallel for num_threads(P) private(i)
            for(i = 0; i<S; i++){

                float dw = 0.0;
                
                // sum Dw of each batch element (256 bits)
                __m512 vDw = _mm512_maskz_loadu_ps(m, &(snn->dw[i * B + 0])); // need 256 bits, 512 are loaded (if it is the final, padding is loaded)
                for(b = 8; b+7 < B; b+=8){

                    vDw = _mm512_mask_add_ps(vDw, m, vDw, _mm512_maskz_loadu_ps(m, &(snn->dw[i * B + b])));
                }
                
                // reduce to a value
                dw = _mm512_mask_reduce_add_ps(m, vDw);

                // compute mean dw
                dw = dw / (float)B;

                // update initial w and w for all copies
                vDw = _mm512_set1_ps(dw);
                for(b = 0; b+7<B; b+=8){

                    // reinitialize dw
                    _mm512_mask_storeu_ps(&(snn->dw[i * B + b]), m, vZero);

                    // init_w = init_w + dw
                    __m512 vInitW = _mm512_maskz_loadu_ps(m, &(snn->init_w[i * B + b]));
                    vInitW = _mm512_mask_add_ps(vInitW, m, vInitW, vDw);
                    _mm512_mask_storeu_ps(&(snn->init_w[i * B + b]), m, vInitW);

                    // w = init_w
                    _mm512_mask_storeu_ps(&(snn->w[i * B + b]), m, vInitW);
                }    
            }
        }
        else if(B % 16 == 0)
        {
            #pragma omp parallel for num_threads(P) private(i)
            for(i = 0; i<S; i++){

                float dw = 0.0;
                
                // sum Dw of each batch element
                __m512 vDw = _mm512_loadu_ps(&(snn->dw[i * B + 0]));
                for(b = 16; b+15 < B; b+=16){

                    vDw = _mm512_add_ps(vDw, _mm512_loadu_ps(&(snn->dw[i * B + b])));
                }
                
                // reduce to a value
                dw = _mm512_reduce_add_ps(vDw);

                // compute mean dw
                dw = dw / (float)B;

                // update initial w and w for all copies
                vDw = _mm512_set1_ps(dw);
                for(b = 0; b+15<B; b+=16){

                    // dw = 0.0
                    _mm512_storeu_ps(&(snn->dw[i * B + b]), vZero); // dw = 0.0, reinitialize

                    // init_w = init_w + dw
                    __m512 vInitW = _mm512_loadu_ps(&(snn->init_w[i * B + b]));
                    vInitW = _mm512_add_ps(vInitW, vDw);
                    _mm512_storeu_ps(&(snn->init_w[i * B + b]), vInitW);

                    // w = init_w
                    _mm512_storeu_ps(&(snn->w[i * B + b]), vInitW);
                }    
            }        
        }
        else{
            // TODO

        }
    }
    #else
    {
        
        #pragma omp parallel for num_threads(P) private(i)
        for(i = 0; i<S; i++){

            float dw = 0.0;

            // sum dw of the batch for each synapse
            for(size_t b = 0; b<B; b++){

                dw += snn->dw[i * B + b];
            }

            // compute mean dw
            dw = dw / (float)B;

            // update initial w and w for all copies
            for(size_t b = 0; b<B; b++){

                snn->dw[i * B + b] = 0.0; // reinitialize dw
                snn->init_w[i * B + b] += dw; // update initial w for incoming batches
                snn->w[i * B + b] = snn->init_w[i * B + b]; // update w
            }    
        }
    }
    #endif
}


/// Function to loop over all input synapses of a neuron and calculate the input current
float compute_neuron_input_current(GPU_SNN_t *snn, size_t t, size_t i){

    size_t max_spikes = (size_t)snn->max_spikes;
    size_t N = (size_t)snn->n_neurons;
    size_t iN = (size_t)snn->n_input_neurons;

    float I = 0.0;
         
    
    // convert conditional to multiplication
    // inlcude STDP here looping over all synapses
    if(snn->r_period_remain[i] <= 0){

        // variables declaration
        size_t synapse_index, in_neuron_index, spk_time_index;
        int delay, spk_time; 
        float spk;    
        size_t base_synapse = (size_t)(snn->neuron_input_synapses_offset[i]);
        size_t iS = (size_t)snn->n_neuron_input_synapses[i];

        // var to iterate
        size_t j = 0; 

        // vectorize each neuron processing using AVX512
        #if defined AVX512

        for(j = 0; j + 15 < iS; j+=16){

            // load delays and weights
            __m512i d_vec = _mm512_loadu_si512(&(snn->delay[base_synapse + j]));                           
            __m512 w_vec = _mm512_loadu_ps(&(snn->w[base_synapse + j]));    
            
            // compute spike times for each synapse (t-delay)
            __m512i spk_time_vec = _mm512_sub_epi32(_mm512_set1_epi32(t), d_vec);

            // mask: 1 if the spk_time is spk_time >= 0
            __mmask16 valid_mask = _mm512_cmp_epi32_mask(spk_time_vec, _mm512_setzero_si512(), _MM_CMPINT_GE);

            // clamp spike times to avoid spk_time < 0 values
            spk_time_vec =_mm512_max_epi32(spk_time_vec, _mm512_setzero_si512());


            // load input neuron indeces for each synapse
            __m512i in_vec = _mm512_loadu_si512(&(snn->pre_neuron_index[base_synapse + j]));  
            
            // matrix [T * N] // This could be fixed to avoid the computation
            // compute indices to gather spikes [in_index + spk_time * (n_neurons + n_input_neurons)]
            __m512i idx = 
                _mm512_add_epi32(
                    _mm512_mullo_epi32(spk_time_vec, 
                        _mm512_add_epi32(_mm512_set1_epi32(snn->n_neurons), 
                        _mm512_set1_epi32(snn->n_input_neurons))), 
                    in_vec);


            // gather spikes from spike matrix using the indices (int)
            __m512i spikes_i = _mm512_i32gather_epi32(idx, snn->spk_matrix, 4);

            // convert to float
            __m512 spikes = _mm512_cvtepi32_ps(spikes_i);
            
            // multiply by weights to compute I(t)
            __m512 contrib = _mm512_mask_mul_ps(_mm512_setzero_ps(), valid_mask, w_vec, spikes); // multiply
            float I_vec = _mm512_reduce_add_ps(contrib); // reduce
            I += I_vec; 


            // store wether the presynaptic neuron fired (spikes reached now)
            _mm512_storeu_si512(&(snn->pre_fired[base_synapse + j]), spikes_i); // spikes_i stores which neurons fired
        }
        #endif

        // loop over remaining neuron input synapses
        for(; j<iS; j++){

            synapse_index = base_synapse + j; // no copies
            in_neuron_index = (size_t)(snn->pre_neuron_index[synapse_index]); // absurd, equal to 0
            delay = snn->delay[synapse_index]; // no copies
            spk_time = t - delay; // actual position

            // check if spike is bigger than 0
            if(spk_time >=0){ 
                
                // T * N
                // get wether the neuron fired in the matrix and compute I(t)
                spk = (float)(snn->spk_matrix[(size_t)(spk_time)  * (N + iN) + in_neuron_index]);    
                I += snn->w[synapse_index] * spk; // spk = 0 / 1
                snn->pre_fired[synapse_index] = (int)spk;
            }
        }
    }

    // store refractory period and I
    return I;


}

void compute_input_current(GPU_SNN_t *snn, simulation_configuration_t *conf, size_t t){
    
    size_t N = (size_t)snn->n_neurons;
    size_t P = (size_t)conf->n_process;

    size_t i = 0;

    #pragma omp parallel for num_threads(P) 
    for(i=0; i<N; i++){

        snn->arrI[i] = compute_neuron_input_current(snn, t, i);
        snn->inR[i] = (snn->r_period_remain[i] <= 0); // store wether the neuron is in refractary period // necessary?
    }
}


void compute_input_current_batch(GPU_SNN_t *snn, simulation_configuration_t *conf, size_t t, size_t gt){

    size_t i, j, b;
    size_t N, iN, P, B, LT;

    // store general variables
    N = snn->n_neurons;
    iN = snn->n_input_neurons;
    P = conf->n_process;
    B = conf->batch_size;
    LT = snn->LT;
    

    /*for(i = 0; i<N; i++){

        printf("[");
        for(b = 0; b<B; b++){
        
            printf(" %f, ", snn->arrI[i * B + b]);
        }
        printf("%f],", snn->arrI[i * B + B - 1]);
    }
    printf("\n");*/

    // loop over neurons
    #pragma omp parallel for num_threads(P) private (i, j, b) 
    for(i=0; i<N; i++){
        
        // variables declaration
        size_t synapse_index, g_synapse_index, in_neuron_index;
        int delay, spk_time; 
        float spk;  

        size_t neuron_index = i; // for not copied variables
        size_t g_neuron_index = i * B; // for copied variables

        // get neuron information
        size_t iS = snn->n_neuron_input_synapses[neuron_index]; 
        size_t base_synapse = snn->neuron_input_synapses_offset[neuron_index]; // there are several copies of each synapse, so the offset depends

        #if defined AVX512
        {

            int preIndex = (iN + N) * B;

            // vectorize batch
            if(B == 8){ // masked operations

                size_t n_blks = B / 8;

                // initialize vector of I (for all batch elements)
                __m256 vI_arr[n_blks];
                __m256i vRefPer_arr[n_blks];
                __mmask8 mRefPer_arr[n_blks];
                __mmask8 mRefPer = 0;
                __m256i vOnes = _mm256_set1_epi8(1);
                __m256i vZeros = _mm256_set1_epi8(0);
                __m128i vZeros128 = _mm512_castsi512_si128(_mm512_set1_epi8(0));

                // the same as above but masked
                for(size_t blk = 0; blk<n_blks; blk++){

                    vI_arr[blk] = _mm256_set1_ps(0.0); // I = 0.0;
                    vRefPer_arr[blk] = _mm256_loadu_epi32(&(snn->r_period_remain[g_neuron_index + blk * 8])); 
                    mRefPer_arr[blk] = _mm256_cmp_epi32_mask(vRefPer_arr[blk], vZeros, _MM_CMPINT_LE); 

                    mRefPer |= mRefPer_arr[blk];
                }

                if(mRefPer || conf->learn){ // if at least a bit is 1

                    // loop over input synapses of the neuron
                    for(j = 0; j<iS; j++){

                        // get synapse index. In not copied variables use synapse_index, in copied g_synapse_index 
                        synapse_index = base_synapse + j;
                        g_synapse_index = synapse_index * B; // scale base synapse, since now there are B copies of each one
                        
                        // get synapse delay (not copied)
                        delay = snn->delay[synapse_index];

                        // get spike time
                        spk_time = (int)t - delay;
                    
                        // correct spk_time
                        if(spk_time < 0 && gt >= LT){ // CHECK
                            spk_time = (int)LT + spk_time;
                        }

                        // process spike 
                        if(spk_time >= 0){
                                                            
                            // get input neuron index (not copied)
                            in_neuron_index = snn->pre_neuron_index[synapse_index]; 

                            // index to load the spike from
                            size_t index = (preIndex * (size_t)spk_time) + ((in_neuron_index * B));
                            size_t b_index;

                            for(size_t blk = 0; blk < n_blks; blk++){

                                // check if this part of the batch is in refractory period
                                if(mRefPer_arr[blk] || conf->learn){
                                    
                                    b_index = index + blk * 8;
                                    
                                    // load spikes
                                    //__m128i vCharSpikes = _mm_loadl_epi64((const __m128i*)&snn->spk_matrix[b_index]); // load 8 bytes                                
                                    //__mmask16 spk_mask16 = _mm_cmpgt_epi8_mask(vCharSpikes, _mm_setzero_si128());
                                    //__mmask8 spk_mask = (__mmask8)(spk_mask16 & 0xFF); // store the neurons that received an spike
                                    
                                    // AVX512
                                    __m128i vCharSpikes = _mm_loadu_epi8((const __m128i*)&snn->spk_matrix[b_index]); // 128 bits loaded, but we only need 64 (char(8) * batch_size(8))                                
                                    __mmask16 spk_mask16 = _mm_cmpgt_epi8_mask(vCharSpikes, vZeros128);
                                    __mmask8 spk_mask = (__mmask8)(spk_mask16 & 0xFF); // store the neurons that received an spike


                                    // update I for neurons that are not in refractory period
                                    if(mRefPer_arr[blk] && spk_mask){
                                        
                                        // load weights
                                        __m256 vW = _mm256_loadu_ps(&(snn->w[g_synapse_index + blk * 8]));

                                        // update I for those neurons where r <= 0 and received a spike
                                        vI_arr[blk] = _mm256_mask_add_ps(vI_arr[blk], spk_mask & mRefPer_arr[blk], vI_arr[blk], vW);
                                    }

                                    // if TR-STDP is used, store whether the neuron fired
                                    if(conf->learn && spk_mask){
                                        
                                        // 1 if a spike is received and neuron is not in refractary period
                                        //_mm512_mask_storeu_epi8(&(snn->pre_fired[g_synapse_index + blk * 8]), spk_mask ,_mm512_set1_epi8(1));
                                        _mm_mask_storeu_epi8(&(snn->pre_fired[g_synapse_index + blk * 8]), spk_mask, _mm256_castsi256_si128(vOnes));
                                    }
                                }
                            }
                        }
                    }
                }    
                for(size_t blk = 0; blk < n_blks; blk++){
                    
                    _mm256_storeu_ps(&(snn->arrI[g_neuron_index + 8 * blk]), vI_arr[blk]);
                }                    
            }

            // [B % 16 == 0]
            else if(B % 16 == 0){

                size_t n_blks = B / 16, blk;

                // initialize vector of I (for all batch elements) and refractary periods
                __m512 vI_arr[n_blks];
                __m512i vRefPer_arr[n_blks];
                __mmask16 mRefPer_arr[n_blks];
                __mmask16 mRefPer = 0;


                for(blk = 0; blk<n_blks; blk++){

                    vI_arr[blk] = _mm512_set1_ps(0.0); // I = 0.0;

                    // refractary periods can be used to avoid some computations
                    vRefPer_arr[blk] = _mm512_loadu_si512(&(snn->r_period_remain[g_neuron_index + blk * 16])); 

                    // create mask (epi32 -- int -- 32 bits)
                    mRefPer_arr[blk] = _mm512_cmp_epi32_mask(vRefPer_arr[blk], _mm512_setzero_si512(), _MM_CMPINT_LE); 

                    mRefPer |= mRefPer_arr[blk];
                }

                // loop over input synapses only if there is a neuron in the block of 16 which is not in refractory period
                if(mRefPer || conf->learn){

                    // loop over input synapses of the neuron
                    for(j = 0; j<iS; j++){

                        // get synapse index. In not copied variables use synapse_index, in copied g_synapse_index 
                        synapse_index = base_synapse + j;
                        g_synapse_index = synapse_index * B; // scale base synapse, since now there are B copies of each one
                        
                        // get synapse delay (not copied)
                        delay = snn->delay[synapse_index];

                        // get spike time
                        spk_time = (int)t - delay;
                    
                        // correct spk_time (circular array, negative indexes the last)
                        if(spk_time < 0 && gt >= LT){ // CHECK
                            spk_time = (int)LT + spk_time;
                        }

                        // process spike if necessary
                        if(spk_time >= 0){
                                                            
                            // get input neuron index (not copied)
                            in_neuron_index = snn->pre_neuron_index[synapse_index]; 

                            // index to load the spike from in the spike matrix
                            size_t index = (preIndex * (size_t)spk_time) + ((in_neuron_index * B));
                            size_t b_index; // neuron index

                            // loop over blocks: set of 16 neurons (batched)
                            for(blk = 0; blk < n_blks; blk++){

                                // if all neurons in this mini-batch are in refractary period, avoid computation
                                if(mRefPer_arr[blk] || conf->learn){
                                    
                                    // compute first neuron index
                                    b_index = index + blk * 16;
                                    
                                    // load spikes and get mask
                                    // SSE - AVX512
                                    //__m128i vCharSpikes = _mm_loadu_si128((const __m128i*)&snn->spk_matrix[b_index]);                                
                                    //__mmask16 spk_mask = _mm_cmpgt_epi8_mask(vCharSpikes, _mm_setzero_si128()); // mask the received spikes

                                    // AVX512
                                    __m128i vCharSpikes = _mm_loadu_epi8((const __m128i*)&snn->spk_matrix[b_index]);                                
                                    __mmask16 spk_mask = _mm_cmpgt_epi8_mask(vCharSpikes, _mm_setzero_si128()); // mask the received spikes

                                    // update I only if there are neurons out of refractory period
                                    if(mRefPer_arr[blk] & spk_mask){

                                        // load weights   
                                        __m512 vW = _mm512_loadu_ps(&(snn->w[g_synapse_index + blk * 16]));

                                        // compute I: store in spk_mask & mRefPer: received an spike and not in refractory period
                                        vI_arr[blk] = _mm512_mask_add_ps(vI_arr[blk], spk_mask & mRefPer_arr[blk] , vI_arr[blk], vW);
                                    }

                                    // if TR-STDP is used, store whether the neuron fired
                                    if(conf->learn && spk_mask){
                                        
                                        // 1 if a spike is received and neuron is not in refractary period
                                        _mm512_mask_storeu_epi8(&(snn->pre_fired[g_synapse_index + blk * 16]), spk_mask ,_mm512_set1_epi8(1));
                                        //_mm_mask_storeu_epi8(&(snn->pre_fired[g_synapse_index + blk * 16]), spk_mask,_mm512_castsi512_si128(vOnes));
                                    }
                                } 
                            }
                        }
                    }
                }  

                // store I for all neurons
                for(blk = 0; blk < n_blks; blk++){
                    
                    _mm512_storeu_ps(&(snn->arrI[g_neuron_index + 16 * blk]), vI_arr[blk]);
                }                    
            }
        }   
        #else
        {
            // loop over batches
            for(b = 0; b<B; b++){

                float I = 0.0;

                // check wether the neuron is in refractary period: pre_fired flag not activated although it should be
                if(snn->r_period_remain[g_neuron_index + b] <= 0 || conf->learn){
                
                    // loop over input synapses
                    for(j=0; j<iS; j++){

                        // get synapse index. In not copied variables use synapse_index, in copied g_synapse_index 
                        synapse_index = base_synapse + j;
                        g_synapse_index = synapse_index * B; // scale base synapse, since now there are B copies of each one
                        
                        // get input neuron index (not copied)
                        in_neuron_index = snn->pre_neuron_index[synapse_index]; 

                        // get synapse delay (not copied)
                        delay = snn->delay[synapse_index];

                        // get spike time
                        spk_time = (int)t - delay;
                    
                        // correct spk_time
                        if(spk_time < 0 && gt >= LT){ // CHECK
                            spk_time = (int)LT + spk_time;
                        }

                        // process spike 
                        if(spk_time >= 0){
                            
                            // index to load the spike from
                            size_t index = ((iN + N) * B * (size_t)spk_time) + ((in_neuron_index * B));

                            // get spike
                            spk = (float)(snn->spk_matrix[index + b]); 

                            // update I if r period remain <= 0 and spike received
                            if(snn->r_period_remain[g_neuron_index + b] <= 0 && (int)spk == 1){
                                
                                I += snn->w[g_synapse_index + b] * spk; // spk = 0 / 1
                            }

                            // store if the presynaptic neuron fired for TR-STDP
                            if(conf->learn && (int)spk == 1){

                                snn->pre_fired[g_synapse_index + b] = (char)spk;
                            }
                        }
                    }
                }

                snn->arrI[g_neuron_index + b] = I;
            }
        }
        #endif
    }

    /*for(i = 0; i<N; i++){

        printf("[");
        for(b = 0; b<B; b++){
        
            printf(" %f, ", snn->arrI[i * B + b]);
        }
        printf("%f],", snn->arrI[i * B + B - 1]);
    }
    printf("\n");*/
}


// compute LIF neuron model dynamics for all neurons
void compute_LIF_V(GPU_SNN_t *snn, simulation_configuration_t *conf, size_t t){

    // get general information
    size_t max_spikes = (size_t)snn->max_spikes;
    size_t N = (size_t)snn->n_neurons;
    size_t iN = (size_t)snn->n_input_neurons;
    size_t P = (size_t)conf->n_process;

    size_t n_tasks = N / 16, tsk, i = 0;

    // helpers
    const float alpha = 0.95f;
    const float beta = 0.05f;

    // parallelize in tasks and vectorize if openMP and AVX512 are defined
    /*#if defined AVX512

    // load alpha and beta
    __m512 alpha_v = _mm512_set1_ps(alpha);
    __m512 beta_v  = _mm512_set1_ps(beta);

    #pragma omp parallel for num_threads(P) private(tsk, i)
    for(tsk = 0; tsk<n_tasks; tsk++){
        
        // get first element to process (each task processes 16 elements)
        i = tsk * 16;

        // load neuron information
        __m512 v_vec    = _mm512_loadu_ps(&snn->v[i]);
        __m512 vrest_vec= _mm512_loadu_ps(&snn->v_rest[i]);
        __m512 I_vec    = _mm512_loadu_ps(&snn->arrI[i]);
        __m512i ref_vec  = _mm512_loadu_si512(&snn->inR[i]);

        // mask which will be computed: compute if inR > 0
        __mmask16 mask = _mm512_cmp_epi32_mask(ref_vec, _mm512_setzero_si512(), _MM_CMPINT_GT); // lower or equal

        // compute new v-s for all neurons
        __m512 newv = _mm512_add_ps(
                _mm512_mul_ps(alpha_v, v_vec),
                _mm512_add_ps(_mm512_mul_ps(beta_v, vrest_vec), I_vec)
            );


        // combine old & new: v = mask ? newv : v_old
        __m512 v_combined = _mm512_mask_mov_ps(v_vec, mask, newv);

        // store V[t]
        _mm512_storeu_ps(&snn->v[i], v_combined);
    }

    // handle remaining neurons (or all)
    #pragma omp parallel for num_threads(P)
    for (i = n_tasks * 16; i < N; i++) {

        // check if the neuron is in refractary period
        if (snn->inR[i] == 1){
        
            // compute v[t]
            snn->v[i] = alpha * snn->v[i] + beta * snn->v_rest[i] + snn->arrI[i];
        }
    }

    #else*/

    #pragma omp parallel for num_threads(P)
    for (i = 0; i < N; i++) {

        // check if the neuron is in refractary period
        if (snn->r_period[i] == 1){
        
            // compute v[t]
            snn->v[i] = alpha * snn->v[i] + beta * snn->v_rest[i] + snn->arrI[i];
        }
    }

    //#endif
}


void compute_LIF_V_batch(GPU_SNN_t *snn, simulation_configuration_t *conf, size_t t, size_t gt){

    // get general information
    size_t N, P, B;
    size_t i, b;
    size_t neuron_index, g_neuron_index;

    N = snn->n_neurons;
    P = conf->n_process;
    B = conf->batch_size;

    // helpers
    const float alpha = 0.95f;
    const float beta = 0.05f;

    #if defined AVX512
    {
        if(B == 8){
            
            // initialize constants
            __m256 avec = _mm256_set1_ps(alpha);
            __m256 bvec = _mm256_set1_ps(beta);

            #pragma omp parallel for num_threads(P) private(i, b, neuron_index, g_neuron_index)
            for(i = 0; i<N; i++){   
                
                neuron_index = i;
                __m256 v_rest_vec = _mm256_set1_ps(snn->v_rest[neuron_index]);

                for(b = 0; b+7<B; b+=8){

                    g_neuron_index = i * B + b;

                    // load refractary period and create mask
                    __m256i refract_vec = _mm256_loadu_epi32(&(snn->r_period_remain[g_neuron_index]));
                    __mmask8 rft_mask = _mm256_cmp_epi32_mask(refract_vec, _mm256_setzero_si256(), _MM_CMPINT_LE); // if lower or eual, process

                    // load I, rest and v, if any neuron is not in refractory period
                    if(rft_mask){            

                        __m256 v_vec = _mm256_loadu_ps(&(snn->v[g_neuron_index]));
                        __m256 I_vec = _mm256_loadu_ps(&(snn->arrI[g_neuron_index]));

                        // compute dynamics
                        __m256 tmp = _mm256_add_ps(I_vec, _mm256_mul_ps(v_rest_vec, bvec));
                        __m256 rhs = _mm256_fmadd_ps(v_vec, avec, tmp);  // alpha*v + beta*v_rest + I
                        v_vec = _mm256_mask_mov_ps(v_vec, rft_mask, rhs); // only update active lanes

                        // store v
                        _mm256_storeu_ps(&(snn->v[g_neuron_index]), v_vec);
                    }
                }
            }
        }
        else if(B % 16 == 0){

            // initialize constants
            __m512 avec = _mm512_set1_ps(alpha);
            __m512 bvec = _mm512_set1_ps(beta);

            #pragma omp parallel for num_threads(P) private(i, b, neuron_index, g_neuron_index)
            for(i = 0; i<N; i++){   
                
                neuron_index = i;

                // load v_rest
                __m512 v_rest_vec = _mm512_set1_ps(snn->v_rest[neuron_index]);

                for(b = 0; b+15<B; b+=16){

                    // compute neuron index
                    g_neuron_index = i * B + b;

                    // load refractary period and create mask
                    __m512i refract_vec = _mm512_loadu_si512(&(snn->r_period_remain[g_neuron_index]));
                    __mmask16 rft_mask = _mm512_cmp_epi32_mask(refract_vec, _mm512_setzero_si512(), _MM_CMPINT_LE); // if lower or eual, process

                    // avoid computation if all neurons are in refractory period
                    if(rft_mask){

                        // load v and I
                        __m512 v_vec = _mm512_loadu_ps(&(snn->v[g_neuron_index]));
                        __m512 I_vec = _mm512_loadu_ps(&(snn->arrI[g_neuron_index]));

                        // compute dynamics v = rft ? v : compute_lif()  <-- this is done by the next lines
                        __m512 tmp = _mm512_add_ps(I_vec, _mm512_mul_ps(v_rest_vec, bvec));
                        __m512 rhs = _mm512_fmadd_ps(v_vec, avec, tmp);  // alpha*v + beta*v_rest + I
                        v_vec = _mm512_mask_mov_ps(v_vec, rft_mask, rhs); // only update active lanes

                        // store v (all values are stored, mask is not necessary since the previous operation stores the final values)
                        _mm512_storeu_ps(&(snn->v[g_neuron_index]), v_vec); 
                    }
                }
            }
        }
    }
    #else
    {
        // loop over neurons
        #pragma omp parallel for num_threads(P) private(i, b, neuron_index, g_neuron_index)
        for (i = 0; i < N; i++) {
            
            neuron_index = i;
            g_neuron_index = i * B;

            // loop over batches
            for(b = 0; b<B; b++){
                
                if (snn->r_period_remain[g_neuron_index + b] <= 0){
                    
                    snn->v[g_neuron_index + b] = alpha * snn->v[g_neuron_index + b] + beta * snn->v_rest[neuron_index] + snn->arrI[g_neuron_index + b];
                }
            }
        }
    }
    #endif
}


void process_neuron_firing(GPU_SNN_t *snn, simulation_configuration_t *conf, GPU_results_t *results, size_t t){

    // get general information
    size_t max_spikes = (size_t)snn->max_spikes;
    size_t N = (size_t)snn->n_neurons;
    size_t iN = (size_t)snn->n_input_neurons;
    size_t P = (size_t)conf->n_process;

    // T * N vectorized
    size_t idx = (iN + N) * t + iN; // index of the first neuron in the spike matrix in time step t

    // compute the number of tasks
    size_t n_tasks = N / 16, tsk, i = 0;


    #if defined AVX512
    
    __m512i ones_vec = _mm512_set1_epi32(1);

    #pragma omp parallel for num_threads(P) private(tsk, i)
    for(tsk = 0; tsk<n_tasks; tsk++){
        
        // get first element to process (each task processes 16 elements)
        i = tsk * 16;
        
        // update refractory period
        __m512i ref_vec = _mm512_loadu_si512(&(snn->r_period_remain[i])); // load
        ref_vec = _mm512_sub_epi32(ref_vec, ones_vec); // r_time = r_time - 1, stored later


        // load v and v_thresh
        __m512 v_vec   = _mm512_loadu_ps(&snn->v[i]);
        __m512 vth_vec = _mm512_loadu_ps(&snn->v_thresh[i]);

        // firing mask (if v >= v_thresh by mask)
        __mmask16 fire_mask = _mm512_cmp_ps_mask(v_vec, vth_vec, _CMP_GE_OS);

        // create the array of spikes (1 if v >= v_thresh, 0 else)
        __m512i spike_vec =
            _mm512_mask_mov_epi32(_mm512_setzero_si512(),
                                fire_mask,
                                _mm512_set1_epi32(1));

        // store which neurons fired in the spk matrix [T * N], neurons are in sequential positions
        _mm512_storeu_si512(&snn->spk_matrix[idx + i], spike_vec);
        
        // reset refractory for neurons that fired
        __m512i ref_reset_vec = _mm512_loadu_si512(&snn->r_period[i]); // load r_period
        ref_vec = _mm512_mask_mov_epi32(ref_vec, fire_mask, ref_reset_vec); // r_period_remain = r_period (if fired)
        _mm512_storeu_si512(&snn->r_period_remain[i], ref_vec); // store

        // reset 
        __m512 vrest_vec = _mm512_loadu_ps(&snn->v_rest[i]); // load reset
        v_vec = _mm512_mask_mov_ps(v_vec, fire_mask, vrest_vec); // v = reset (if fired)
        _mm512_storeu_ps(&snn->v[i], v_vec); // store

        // increment spike count
        __m512i nspk_vec = _mm512_loadu_si512(&results->n_spks[i]); // load n_spikes
        nspk_vec = _mm512_mask_add_epi32(nspk_vec, fire_mask, nspk_vec, ones_vec); // n_spikes += 1 (if fired)
        _mm512_storeu_si512(&results->n_spks[i], nspk_vec); // store

        // store which neurons fired
        _mm512_storeu_si512(&(snn->post_fired[i]), spike_vec);
    }

    #pragma omp parallel for num_threads(P)
    for(i = n_tasks * 16; i<N; i++){
            
        snn->r_period_remain[i] --;
        if(snn->v[i] >= snn->v_thresh[i]){
            
            // T * N
            size_t idx_neuron = idx + i; // row start at

            snn->spk_matrix[idx_neuron] = 1;

            // reinit neuron values
            snn->r_period_remain[i] = snn->r_period[i];
            snn->v[i] = snn->v_rest[i]; // reinit v_rest

            results->n_spks[i] += 1;
        }

        snn->post_fired[i] = snn->v[i] >= snn->v_thresh[i];
    }

    #else

    #pragma omp parallel for num_threads(P)
    for(i = 0; i<N; i++){
            
        snn->r_period_remain[i] --;
        if(snn->v[i] >= snn->v_thresh[i]){
            
            // T * N
            size_t idx_neuron = idx + i; // row start at

            snn->spk_matrix[idx_neuron] = 1;

            // reinit neuron values
            snn->r_period_remain[i] = snn->r_period[i];
            snn->v[i] = snn->v_rest[i]; // reinit v_rest

            results->n_spks[i] += 1;
        }

        snn->post_fired[i] = snn->v[i] >= snn->v_thresh[i];
    }

    #endif
}


void process_neuron_firing_batch(GPU_SNN_t *snn, simulation_configuration_t *conf, GPU_results_t *results, size_t t, size_t gt){

    // get general information
    size_t N, iN, P, B;
    size_t neuron_index, g_neuron_index;
    size_t i, b;

    N = (size_t)snn->n_neurons;
    iN = (size_t)snn->n_input_neurons;
    P = (size_t)conf->n_process;
    B = conf->batch_size;

    #if defined AVX512
    {
        if(B == 8){

            // constants
            const __m256i vOnes = _mm256_set1_epi8(1);

            #pragma omp parallel for num_threads(P) private(i, b, neuron_index, g_neuron_index)
            for(i = 0; i<N; i++){
            
                // get global and local neuron indexes
                neuron_index = i;
                g_neuron_index = i * B;
                size_t b_neuron;

                // index of the first neuron in the spike matrix in time step t
                size_t idx = (iN + N) * B * t + B * (iN + neuron_index); 
                size_t batch_idx;

                // load threshold and resting potential
                __m256 vThresh = _mm256_set1_ps(snn->v_thresh[neuron_index]);
                __m256 vRest = _mm256_set1_ps(snn->v_rest[neuron_index]);
                __m256i vRef = _mm256_set1_epi32(snn->r_period[neuron_index]);

                for(b = 0; b+7<B; b+=8){

                    batch_idx = idx + b;
                    b_neuron = g_neuron_index + b; // first neuron index to be loaded in register

                    // reduce refractary period
                    __m256i vRefPer = _mm256_loadu_epi32(&(snn->r_period_remain[b_neuron]));
                    vRefPer = _mm256_sub_epi32(vRefPer, _mm256_set1_epi32(1));

                    // check if neuron fired, and in that case, add spike to matrix
                    __m256 vV = _mm256_loadu_ps(&(snn->v[b_neuron]));

                    // mask for v >= v_thresh
                    __mmask8 fire_mask = _mm256_cmp_ps_mask(vV, vThresh, _MM_CMPINT_GE);

                    if(fire_mask){

                        // v = v_rest for those that fired
                        vV = _mm256_mask_mov_ps(vV, fire_mask, vRest);         
                        _mm256_storeu_ps(&(snn->v[b_neuron]), vV);

                        // r_period_remain = r_period for those that fired
                        vRefPer = _mm256_mask_mov_epi32(vRefPer, fire_mask, vRef); // store r_period where the neuron fired
                        
                        // load, update and store generated spikes
                        __m256i vNSpikes = _mm256_loadu_epi32(&(results->n_spks[b_neuron]));
                        vNSpikes = _mm256_mask_add_epi32(vNSpikes, fire_mask, vNSpikes, _mm256_set1_epi32(1));
                        _mm256_storeu_epi32(&(results->n_spks[b_neuron]), vNSpikes);

                        // store generated spikes in matrix (64 bits)
                        _mm256_mask_storeu_epi8(&(snn->spk_matrix[batch_idx]), fire_mask, _mm256_set1_epi8(1));


                        // if TB-STDP is used, store which neurons fired
                        if(conf->learn){
                            
                            _mm256_mask_storeu_epi8(&(snn->post_fired[b_neuron]), fire_mask, _mm256_set1_epi8(1));
                        }
                    }
                    _mm256_storeu_epi32(&(snn->r_period_remain[b_neuron]), vRefPer);
                }
            }
        }   
        else if(B % 16 == 0){     

            // constants
            const __m512i vOnes = _mm512_set1_epi8(1);

            #pragma omp parallel for num_threads(P) private(i, b, neuron_index, g_neuron_index)
            for(i = 0; i<N; i++){
            
                // get global and local neuron indexes
                neuron_index = i;
                g_neuron_index = i * B;
                size_t b_neuron;

                // index of the first neuron in the spike matrix in time step t
                size_t idx = (iN + N) * B * t + B * (iN + neuron_index); 
                size_t batch_idx;

                // load threshold, resting potential and refractory period
                __m512 vThresh = _mm512_set1_ps(snn->v_thresh[neuron_index]);
                __m512 vRest = _mm512_set1_ps(snn->v_rest[neuron_index]);
                __m512i vRef = _mm512_set1_epi32(snn->r_period[neuron_index]);

                // loop over batch
                for(b = 0; b+15<B; b+=16){

                    // get indexes
                    batch_idx = idx + b; // matrix index
                    b_neuron = g_neuron_index + b; // first neuron index to be loaded in register

                    // reduce refractary period
                    __m512i vRefPer = _mm512_loadu_si512(&(snn->r_period_remain[b_neuron]));
                    vRefPer = _mm512_sub_epi32(vRefPer, _mm512_set1_epi32(1)); // r_remain = r_remain - 1

                    // load V
                    __m512 vV = _mm512_loadu_ps(&(snn->v[b_neuron]));

                    // mask for v >= v_thresh (fired)
                    __mmask16 fire_mask = _mm512_cmp_ps_mask(vV, vThresh, _MM_CMPINT_GE);

                    // avoid computation if no one neuron fired
                    if(fire_mask){
                        
                        // v = v_rest for those that fired
                        vV = _mm512_mask_mov_ps(vV, fire_mask, vRest);         
                        _mm512_storeu_ps(&(snn->v[b_neuron]), vV); // v = v_rest

                        // r_period_remain = r_period for those that fired
                        vRefPer = _mm512_mask_mov_epi32(vRefPer, fire_mask, vRef); // store r_period where the neuron fired, r_remain = r_period

                        // load, update and store generated spikes
                        __m512i vNSpikes = _mm512_loadu_si512(&(results->n_spks[b_neuron])); 
                        vNSpikes = _mm512_mask_add_epi32(vNSpikes, fire_mask, vNSpikes, _mm512_set1_epi32(1));
                        _mm512_storeu_si512(&(results->n_spks[b_neuron]), vNSpikes); 

                        // store generated spikes in matrix
                        _mm512_mask_storeu_epi8(&(snn->spk_matrix[batch_idx]), fire_mask, _mm512_set1_epi8(1));

                        // if TB-STDP is used, store which neurons fired
                        if(conf->learn){
                            
                            _mm512_mask_storeu_epi8(&(snn->post_fired[b_neuron]), fire_mask, _mm512_set1_epi8(1));
                        }
                    }

                    _mm512_storeu_si512(&(snn->r_period_remain[b_neuron]), vRefPer); // store r_remain
                }
            }
        }
    }
    #else
    {
        
        #pragma omp parallel for num_threads(P) private(i, b, neuron_index, g_neuron_index)
        for(i = 0; i<N; i++){
                
            neuron_index = i;
            g_neuron_index = i * B;
            
            // compute index to write in the spk matrix
            // [T * B * (iN + N)]
            size_t idx = (iN + N) * B * t + B * (iN + neuron_index); // index of the first neuron in the spike matrix in time step t

            for(b = 0; b<B; b++){
            
                // reduce refractary period
                snn->r_period_remain[g_neuron_index + b] --;

                // check whether fired
                if(snn->v[g_neuron_index + b] >= snn->v_thresh[neuron_index]){
            
                    snn->spk_matrix[idx + b] = 1;

                    // reinit neuron values
                    snn->r_period_remain[g_neuron_index + b] = snn->r_period[neuron_index];
                    snn->v[g_neuron_index + b] = snn->v_rest[neuron_index]; // reinit v_rest

                    results->n_spks[g_neuron_index + b] += 1;

                    // store that the neuron fired for TB-STDP and update the trace
                    if(conf->learn){

                        snn->post_fired[g_neuron_index + b] = 1;
                    }
                }
            }
        }
    }
    #endif
}


void process_trace_based_STDP_batch(GPU_SNN_t *snn, simulation_configuration_t *conf){

    size_t N, S, P, B;
    size_t i, b;

    N = snn->n_neurons;
    S = snn->n_synapses;
    B = conf->batch_size;
    P = conf->n_process;

    float pA = 0.1f, mA = 0.1f, mu = 1.0f, decay = 0.6f;
    
    /*printf(" Post_fired = ");
    for(i = 0; i<N; i++){

        printf("[");
        for(b = 0; b<B-1; b++){
            printf("%d, ", snn->post_fired[i * B + b]);
        }
        
        printf("%d], ", snn->post_fired[i * B + B-1]);
    }
    printf("\n");

    printf(" Pre_fired = ");
    for(i = 0; i<S; i++){

        printf("[");
        for(b = 0; b<B-1; b++){
            printf("%d, ", snn->pre_fired[i * B + b]);
        }
        
        printf("%d], ", snn->pre_fired[i * B + B-1]);
    }
    printf("\n");*/

    #if defined AVX512
    {

        if(B == 8){

            // define constants
            const __m256 vOne   = _mm256_set1_ps(1.0f);
            const __m256 vDecay = _mm256_set1_ps(decay);
            const __m256 vpA = _mm256_set1_ps(pA);
            const __m256 vmA = _mm256_set1_ps(mA);
            const __m128i vZero = _mm_set1_epi8(0); // 16 bytes of zeros

            size_t n_blks = B / 8, blk;


            #pragma omp parallel for num_threads(P) private(i, b, blk)
            for(i = 0; i<N; i++){

                // loop over neuron input synapses
                size_t base_synapse = snn->neuron_input_synapses_offset[i];
                size_t n_iS = snn->n_neuron_input_synapses[i];

                // global neuron index
                size_t g_post_neuron = i * B;

                __mmask8 mPostFired[n_blks];
                __m256 vPostTrace[n_blks];

                // loop over batch to update post synaptic trace (N * B)
                for(blk = 0; blk<n_blks; blk++){

                    // compute neuron index
                    size_t post_idx = g_post_neuron + blk * 16;
                    
                    // load if fired and trace
                    //__m128i vPostFired = _mm_loadu_si128((__m128i*)(&(snn->post_fired[post_idx]))); // char
                    __m128i vPostFired = _mm_loadu_epi8((__m128i*)(&(snn->post_fired[post_idx])));
                    mPostFired[blk] = _mm_cmpgt_epi8_mask(vPostFired, vZero); // mask
                    vPostTrace[blk] = _mm256_loadu_ps(&(snn->post_trace[post_idx])); // load trace
                    
                    // update trace and store trace
                    __m256 vPostDecay = _mm256_mul_ps(vPostTrace[blk], vDecay);
                    vPostTrace[blk] = _mm256_mask_blend_ps(mPostFired[blk], vPostDecay, vOne); // select conditionally
                    _mm256_storeu_ps(&(snn->post_trace[post_idx]), vPostTrace[blk]); // store trace
                }

                // loop over input synapses
                for(size_t j = 0; j<n_iS; j++){

                    // get synapse index
                    size_t synapse_index = base_synapse + j;
                    size_t g_synapse_index = synapse_index * B;

                    // loop over batch
                    for(blk = 0; blk<n_blks; blk++){

                        // compute post and pre indexes
                        size_t post_idx = g_post_neuron + blk * 16;
                        size_t pre_idx  = g_synapse_index + blk * 16;

                        // load wether pre fired and trace
                        //__m128i vPreFired = _mm_loadu_si128((__m128i*)(&(snn->pre_fired[pre_idx]))); // SSE
                        __m128i vPreFired = _mm_loadu_epi8(&(snn->pre_fired[pre_idx])); // AVX512
                        __mmask8 mPreFired = _mm_cmpgt_epi8_mask(vPreFired, vZero); // mask

                        // update pre trace
                        __m256 vPreTrace = _mm256_loadu_ps(&(snn->pre_trace[pre_idx])); // load trace
                        __m256 vPreDecay = _mm256_mul_ps(vPreTrace, vDecay);
                        vPreTrace = _mm256_mask_blend_ps(mPreFired, vPreDecay, vOne);
                        _mm256_storeu_ps(&(snn->pre_trace[pre_idx]), vPreTrace); // store trace

                        // avoid updating v and dw if no one neuron fired
                        if(mPreFired || mPostFired[blk]){ // time reduced to 1/2 (more or less)
                           
                            // update dw
                            __m256 vW = _mm256_loadu_ps(&(snn->w[pre_idx])); // load weight
                            __m256 vPot = _mm256_maskz_mul_ps(mPostFired[blk], vpA, _mm256_mul_ps(_mm256_sub_ps(vOne, vW), vPreTrace)); // LTP
                            __m256 vDep = _mm256_maskz_mul_ps(mPreFired, vmA, _mm256_mul_ps(vW, vPostTrace[blk])); // LTD
                            __m256 vdW = _mm256_sub_ps(vPot, vDep); // dw


                            // update and store w and dw
                            _mm256_storeu_ps(&(snn->w[pre_idx]), _mm256_add_ps(vW, vdW));
                            _mm256_storeu_ps(&(snn->dw[pre_idx]), _mm256_add_ps(_mm256_loadu_ps(&(snn->dw[pre_idx])), vdW));
                        }

                        // reinit pre fired to 0
                        if(mPreFired){
                            //_mm_storeu_si128((__m128i*)&(snn->pre_fired[pre_idx]), vZero); // SSE
                            _mm_mask_storeu_epi8(&(snn->pre_fired[pre_idx]), mPreFired, vZero); // AVX 512           
                        }         
                    }
                }

                // reinit post fired
                for(blk = 0; blk<n_blks; blk++){

                    if(mPostFired[blk]){
                        //_mm_storeu_si128((__m128i*)&(snn->post_fired[g_post_neuron + blk * 16]), vZero); // SSE  
                        _mm_mask_storeu_epi8(&(snn->post_fired[g_post_neuron + blk * 16]), mPostFired[blk], vZero); // AVX 512           
                    }
                }
            }        
        }
        else if(B % 16 == 0){

            // define constants
            const __m512 vOne   = _mm512_set1_ps(1.0f);
            const __m512 vDecay = _mm512_set1_ps(decay);
            const __m512 vpA = _mm512_set1_ps(pA);
            const __m512 vmA = _mm512_set1_ps(mA);
            const __m128i vZero = _mm_set1_epi8(0); // 16 bytes of zeros

            size_t n_blks = B / 16, blk;


            #pragma omp parallel for num_threads(P) private(i, b, blk)
            for(i = 0; i<N; i++){

                // loop over neuron input synapses
                size_t base_synapse = snn->neuron_input_synapses_offset[i];
                size_t n_iS = snn->n_neuron_input_synapses[i];

                // global neuron index
                size_t g_post_neuron = i * B;

                __mmask16 mPostFired[n_blks];
                __m512 vPostTrace[n_blks];

                // loop over batch to update post synaptic trace (N * B)
                for(blk = 0; blk<n_blks; blk++){

                    // compute neuron index
                    size_t post_idx = g_post_neuron + blk * 16;
                    
                    // load if fired and trace
                    //__m128i vPostFired = _mm_loadu_si128((__m128i*)(&(snn->post_fired[post_idx]))); // char
                    __m128i vPostFired = _mm_loadu_epi8((__m128i*)(&(snn->post_fired[post_idx])));
                    mPostFired[blk] = _mm_cmpgt_epi8_mask(vPostFired, vZero); // mask
                    vPostTrace[blk] = _mm512_loadu_ps(&(snn->post_trace[post_idx])); // load trace
                    
                    // update trace and store trace
                    __m512 vPostDecay = _mm512_mul_ps(vPostTrace[blk], vDecay);
                    vPostTrace[blk] = _mm512_mask_blend_ps(mPostFired[blk], vPostDecay, vOne); // select conditionally
                    _mm512_storeu_ps(&(snn->post_trace[post_idx]), vPostTrace[blk]); // store trace
                }

                // loop over input synapses
                for(size_t j = 0; j<n_iS; j++){

                    // get synapse index
                    size_t synapse_index = base_synapse + j;
                    size_t g_synapse_index = synapse_index * B;

                    // loop over batch
                    for(blk = 0; blk<n_blks; blk++){

                        // compute post and pre indexes
                        size_t post_idx = g_post_neuron + blk * 16;
                        size_t pre_idx  = g_synapse_index + blk * 16;

                        // load wether pre fired and trace
                        //__m128i vPreFired = _mm_loadu_si128((__m128i*)(&(snn->pre_fired[pre_idx]))); // SSE
                        __m128i vPreFired = _mm_loadu_epi8(&(snn->pre_fired[pre_idx])); // AVX512
                        __mmask16 mPreFired = _mm_cmpgt_epi8_mask(vPreFired, vZero); // mask

                        // update pre trace
                        __m512 vPreTrace = _mm512_loadu_ps(&(snn->pre_trace[pre_idx])); // load trace
                        __m512 vPreDecay = _mm512_mul_ps(vPreTrace, vDecay);
                        vPreTrace = _mm512_mask_blend_ps(mPreFired, vPreDecay, vOne);
                        _mm512_storeu_ps(&(snn->pre_trace[pre_idx]), vPreTrace); // store trace

                        // avoid updating v and dw if no one neuron fired
                        if(mPreFired || mPostFired[blk]){ // time reduced to 1/2 (more or less)
                           
                            // update dw
                            __m512 vW = _mm512_loadu_ps(&(snn->w[pre_idx])); // load weight
                            __m512 vPot = _mm512_maskz_mul_ps(mPostFired[blk], vpA, _mm512_mul_ps(_mm512_sub_ps(vOne, vW), vPreTrace)); // LTP
                            __m512 vDep = _mm512_maskz_mul_ps(mPreFired, vmA, _mm512_mul_ps(vW, vPostTrace[blk])); // LTD
                            __m512 vdW = _mm512_sub_ps(vPot, vDep); // dw


                            // update and store w and dw
                            _mm512_storeu_ps(&(snn->w[pre_idx]), _mm512_add_ps(vW, vdW));
                            _mm512_storeu_ps(&(snn->dw[pre_idx]), _mm512_add_ps(_mm512_loadu_ps(&(snn->dw[pre_idx])), vdW));
                        }

                        // reinit pre fired to 0
                        if(mPreFired){
                            //_mm_storeu_si128((__m128i*)&(snn->pre_fired[pre_idx]), vZero); // SSE
                            _mm_storeu_epi8(&(snn->pre_fired[pre_idx]), vZero); // AVX 512           
                        }         
                    }
                }

                // reinit post fired
                for(blk = 0; blk<n_blks; blk++){

                    if(mPostFired[blk]){
                        //_mm_storeu_si128((__m128i*)&(snn->post_fired[g_post_neuron + blk * 16]), vZero); // SSE  
                        _mm_storeu_epi8(&(snn->post_fired[g_post_neuron + blk * 16]), vZero); // AVX 512           
                    }
                }
            }
        }
    }
    #else
    {
        // loop over synapses, update traces and weights
        #pragma omp parallel for num_threads(P) private(i, b)
        for(i = 0; i<N; i++){

            // loop over neuron input synapses
            size_t base_synapse = snn->neuron_input_synapses_offset[i];
            size_t n_iS = snn->n_neuron_input_synapses[i];
            size_t post_idx, pre_idx;
            size_t synapse_index;
            size_t g_synapse_index;
            
            // compute neuron global index
            size_t g_post_neuron = i * B;

            // update post trace (depends only on the neuron)
            for(b = 0; b<B; b++){
                
                post_idx = g_post_neuron + b;
                snn->post_trace[post_idx] = snn->post_fired[post_idx] ? 1.0f : snn->post_trace[post_idx] * decay;
            }

            // loop over input synapses, update presynaptic traces (depends on synapses) and update weights
            for(size_t j = 0; j<n_iS; j++){

                synapse_index = base_synapse + j;
                g_synapse_index = synapse_index * B;

                for(b = 0; b<B; b++){

                    post_idx = g_post_neuron + b;
                    pre_idx  = g_synapse_index + b;

                    // update pre_trace
                    snn->pre_trace[pre_idx] = snn->pre_fired[pre_idx] ? 1.0f : snn->pre_trace[pre_idx] * decay;


                    // update w: avoid computation if no neuron fired
                    float dPot = 0.0, dDep = 0.0;
                    if(snn->post_fired[post_idx]){

                        dPot = pA * (1.0f - snn->w[pre_idx]) * snn->pre_trace[pre_idx];
                    }
                    if(snn->pre_fired[pre_idx]){

                        dDep = mA * (snn->w[pre_idx]) * snn->post_trace[post_idx];
                        
                        // reinit pre_fired
                        snn->pre_fired[pre_idx] = 0;
                    }
                    
                    if(dPot != 0 || dDep != 0){

                        float dw = dPot - dDep;
                            
                        snn->w[pre_idx] += dw;
                        snn->dw[pre_idx] += dw;
                    }
                }
            }

            // reiniti neuron fired
            for(b = 0; b<B; b++){
                
                snn->post_fired[g_post_neuron + b] = 0;
            }        
        }
    }
    #endif

    /*printf(" Post_fired = ");
    for(i = 0; i<N; i++){

        printf("[");
        for(b = 0; b<B-1; b++){
            printf("%d, ", snn->post_fired[i * B + b]);
        }
        
        printf("%d], ", snn->post_fired[i * B + B-1]);
    }
    printf("\n");

    printf(" Pre_fired = ");
    for(i = 0; i<S; i++){

        printf("[");
        for(b = 0; b<B-1; b++){
            printf("%d, ", snn->pre_fired[i * B + b]);
        }
        
        printf("%d], ", snn->pre_fired[i * B + B-1]);
    }
    printf("\n");*/

}


void reinitialize_LIF_neurons(GPU_SNN_t *snn, simulation_configuration_t *conf){
    
    size_t i, N, P;

    N = snn->n_neurons;
    P = conf->n_process;
    
    // reinitialize neurons
    #pragma omp parallel for num_threads(P) private(i)
    for(i = 0; i < N; i++){
        
        snn->v[i] = snn->v_rest[i];
        snn->r_period_remain[i] = 0;
        snn->post_fired[i] = 0;
    }
}

void reinitialize_LIF_neurons_batch(GPU_SNN_t *snn, simulation_configuration_t *conf){
    
    size_t i, b;
    size_t B, N, P;
    size_t g_index;

    N = snn->n_neurons;
    P = conf->n_process;
    B = conf->batch_size;
    

    // TODO: vectorize
    
    // loop over neurons
    #pragma omp parallel for num_threads(P) private(i, b, g_index)
    for(i = 0; i < N; i++){

        g_index = i * B;
        
        // loop over copies in batch (serial)
        for(b = 0; b<B; b++){

            snn->v[g_index + b] = snn->v_rest[i];
            snn->r_period_remain[g_index + b] = 0;
            snn->post_fired[g_index + b] = 0;
            snn->post_trace[g_index + b] = 0.0;
        }
    }
}


void reinitialize_synapses(GPU_SNN_t *snn, simulation_configuration_t *conf){
    
    size_t i, S, P;

    S = snn->n_synapses;
    P = conf->n_process;

    // reinitialize synapses 
    #pragma omp parallel for num_threads(P) private(i)
    for(i = 0; i<S; i++){
        
        snn->w[i] = snn->init_w[i];
        snn->pre_fired[i] = 0;
    }
}


void reinitialize_synapses_batch(GPU_SNN_t *snn, simulation_configuration_t *conf){
    
    size_t i, j, b, S, P, B, gindex, lindex;

    S = snn->n_synapses;
    P = conf->n_process;
    B = conf->batch_size;


    // reinitialize synapses 
    #pragma omp parallel for num_threads(P) private(i, b)
    for(i = 0; i<S; i++){

        size_t g_index = i * B;

        // loop over the batch
        for(b = 0; b<B; b++){

            snn->w[g_index + b] = snn->init_w[i];
            snn->pre_fired[g_index + b] = 0;
            snn->pre_trace[g_index + b] = 0.0;
        }
    }
}

void reinitialize_spk_matrix(GPU_SNN_t *snn, simulation_configuration_t *conf){
    
    size_t i, j, N, iN, P, max_spikes;
    
    N = snn->n_neurons;
    iN = snn->n_input_neurons;
    P = conf->n_process;
    max_spikes = snn->max_spikes;

    // reinitialize spike matrix
    #pragma omp parallel for num_threads(P) private(i, j)
    for(i=0; i < N + iN; i++){
        
        for(j = 0; j<max_spikes; j++){
                
            snn->spk_matrix[max_spikes * i + j] = 0;        
        }
    }
}

void reinitialize_spk_matrix_batch(GPU_SNN_t *snn, simulation_configuration_t *conf){
    
    size_t i, j, b, lindex, gindex;
    size_t N, iN, P, LT, B;

    N = snn->n_neurons;
    iN = snn->n_input_neurons;
    P = conf->n_process;
    B = conf->batch_size;
    LT = snn->LT;

    // reinitialize spike matrix
    #pragma omp parallel for num_threads(P) private(i, j, b, gindex, lindex)
    for(i=0; i < N + iN; i++){
        
        for(j = 0; j<LT; j++){
        
            gindex = i * LT * B + j * B;

            // loop over batch
            for(b = 0; b<B; b++){
                
                lindex = gindex + b;
                snn->spk_matrix[lindex] = 0;
            }    
        }
    }
}


void load_sample_in_SNN(GPU_SNN_t *snn, GPU_dataset_t *dataset, simulation_configuration_t *conf, size_t sidx){
    
    size_t i, j, N, iN, n_features, P;    
    
    N = snn->n_neurons;
    iN = snn->n_input_neurons;
    n_features = dataset->n_features;
    P = conf->n_process;

    // load the sample
    #pragma omp parallel for num_threads(P) private(i, j)
    for(i=0; i < iN; i++){

        if(sidx < dataset->n_samples){

            size_t fidx = i; // the feature index in the sample
            size_t start_feature = (size_t)(dataset->feature_offset[sidx * n_features + fidx]);
            size_t n_spikes_per_feature = (size_t)(dataset->n_spikes_per_feature[sidx * n_features + fidx]);
        
            for(j = 0; j<n_spikes_per_feature; j++){

                // matrix [T * N]
                snn->spk_matrix[((size_t)(dataset->spikes[start_feature + j]) % snn->max_spikes) * (N + iN) + i] = 1;
            }
        }
    }
}

void load_sample_in_SNN_batch(GPU_SNN_t *snn, GPU_dataset_t *dataset, simulation_configuration_t *conf, size_t bidx){
    
    size_t i, j, N, iN, n_features, P, B, b;    
    
    N = snn->n_neurons;
    iN = snn->n_input_neurons;
    n_features = dataset->n_features;
    P = conf->n_process;
    B = conf->batch_size;

    size_t fsidx = bidx * B; // first sample in batch
    size_t sidx;

    // load the sample
    #pragma omp parallel for num_threads(P) private(i, j, b, sidx)
    for(i=0; i < iN; i++){

        // b controls the sample to be loaded
        for(b = 0; b<B; b++){

            // get sample index in batch
            sidx = fsidx + b;

            // check if we reached dataset final
            if(sidx < dataset->n_samples){

                size_t fidx = i; // the feature index in the sample
                size_t start_feature = dataset->feature_offset[sidx * n_features + fidx];
                size_t n_spikes_per_feature = dataset->n_spikes_per_feature[sidx * n_features + fidx];
            
                for(j = 0; j<n_spikes_per_feature; j++){

                    size_t t = dataset->spikes[start_feature + j];

                    // [(iN + N) * B * T]
                    snn->spk_matrix[((iN + N) * B * t) + (B * i) + b] = 1;
                }   
            }
        }
    }
}

// initialize temporal struct to store the spike sproduced by the samples in a batch in a matrix of [iN * T * B] 
// Not vectorized. Only executed once per each batch.
tmp_batch_cpu_t* initialize_batch_matrix(GPU_SNN_t *snn, GPU_dataset_t *dataset, simulation_configuration_t *conf, size_t bidx){

    size_t i, j, n_features, b;    
    size_t N, iN, P, B, T;

    N = snn->n_neurons;
    iN = snn->n_input_neurons;
    P = conf->n_process;
    B = conf->batch_size;
    T = conf->max_input_spikes;

    // get number of features in samples
    n_features = dataset->n_features;

    // index of the first sample in the batch
    size_t fsidx = bidx * B; 
    size_t sidx; // var to store the global index of the sample

    // allocate memory for temporal batch samples matrix // TODO: refactorize to function
    tmp_batch_cpu_t *tmp_batch = (tmp_batch_cpu_t*)malloc(sizeof(tmp_batch_cpu_t));
    tmp_batch->spikes = (char*)calloc(iN * B * T, sizeof(char));


    // load the samples in the batch in the matrix: loop over sample features
    #pragma omp parallel for num_threads(P) private(i, j, b, sidx)
    for(i=0; i < iN; i++){

        // loop over batch (sample = fsidx + b)
        for(b = 0; b<B; b++){

            // get sample index in batch
            sidx = fsidx + b;

            // check if the sample is in the dataset
            if(sidx < dataset->n_samples){

                // get the feature index in the sample
                size_t fidx = i; 
                
                // get the feature offset in the array of spikes
                size_t feature_offset = dataset->feature_offset[sidx * n_features + fidx];

                // get the number of spikes of the feature
                size_t n_spikes_feature = dataset->n_spikes_per_feature[sidx * n_features + fidx];
        
                // loop over spikes and store in the matrix
                for(j = 0; j<n_spikes_feature; j++){

                    // get time step of the spike
                    size_t t = dataset->spikes[feature_offset + j];

                    // store the spike in the matrix of [iN * B * T]
                    tmp_batch->spikes[(iN * B * t) + (B * i) + b] = 1; 
                }
            }
        }
    }

    return tmp_batch;
}

// load the next spikes of the input neurons from the temporal batch dataset
void load_batch_time_step_in_SNN_batch(GPU_SNN_t *snn, GPU_dataset_t *dataset, simulation_configuration_t *conf, tmp_batch_cpu_t *batch_data, size_t bidx, size_t t){
    
    size_t i, b;    
    size_t N, iN, P, B, L;
    size_t fsidx, sidx, lt, spk_mtr_idx, batch_mtr_idx;
    
    // get general information
    N = snn->n_neurons;
    iN = snn->n_input_neurons;
    P = conf->n_process;
    B = conf->batch_size;
    L = snn->LT;
    
    // get the index of the first sample in the batch
    fsidx = bidx * B;
    lt = t % L; // compute local time step for the matrix


    #if defined AVX512
    {
        size_t n_tasks, tsk, idx;

        // spike matrix has chars (8 bits), so in 512 bits, 512 / 8 == 64 chars

        // compute matrix indexes first time step index
        spk_mtr_idx = (iN + N) * B * lt;
        batch_mtr_idx = (iN * B * t);
        
        if(B == 8){

            #pragma omp parallel for num_threads(P) private(i, b, spk_mtr_idx, batch_mtr_idx)
            for(i = 0; i < iN; i++){
                
                // compute matrix indexes
                spk_mtr_idx = (iN + N) * B * lt + (B * i);
                batch_mtr_idx = (iN * B * t) + (B * i);

                for(b = 0; b+7<B; b += 8){
                    
                    // load spikes from batch_data
                    __m128i vSpikes = _mm_loadl_epi64((__m128i*)(&(batch_data->spikes[batch_mtr_idx + b])));

                    // store spikes in spk matrix
                    _mm_storel_epi64((__m128i*)(&(snn->spk_matrix[spk_mtr_idx + b])), vSpikes);
                }
            }

            // reinit the rest of the row
            #pragma omp parallel for num_threads(P) private(i, b, spk_mtr_idx)
            for(i = 0; i<N; i++){

                // compute matrix general index
                spk_mtr_idx = (iN + N) * B * lt + (B * (iN + i));

                for(b = 0; b+7<B; b += 8){
                
                    _mm_storel_epi64((__m128i*)(&(snn->spk_matrix[spk_mtr_idx + b])), _mm_setzero_si128());
                }
            }
        }
        else{

            #pragma omp parallel for num_threads(P) private(i, b, spk_mtr_idx, batch_mtr_idx)
            for(i = 0; i < iN; i++){
                
                // compute matrix indexes
                spk_mtr_idx = (iN + N) * B * lt + (B * i);
                batch_mtr_idx = (iN * B * t) + (B * i);

                for(b = 0; b+15<B; b+= 16){
                    
                    // load spikes from batch_data
                    __m128i vSpikes = _mm_loadu_si128((__m128i*)(&(batch_data->spikes[batch_mtr_idx + b])));

                    // store spikes in spk matrix
                    _mm_storeu_si128((__m128i*)(&(snn->spk_matrix[spk_mtr_idx + b])), vSpikes);
                }
            }

            // reinit the rest of the row
            #pragma omp parallel for num_threads(P) private(i, b, spk_mtr_idx)
            for(i = 0; i<N; i++){

                // compute matrix general index
                spk_mtr_idx = (iN + N) * B * lt + (B * (iN + i));

                for(b = 0; b+15<B; b+= 16){
                
                    _mm_storeu_si128((__m128i*)(&(snn->spk_matrix[spk_mtr_idx + b])), _mm_setzero_si128());
                }
            }
        }
    }
    #else
    {
        // load the batch samples from time step [t] to [t+L]
        #pragma omp parallel for num_threads(P) private(i, b, sidx, spk_mtr_idx, batch_mtr_idx)
        for(i=0; i < iN; i++){ // load features

            // compute matrix indexes
            spk_mtr_idx = (iN + N) * B * lt + (B * i);
            batch_mtr_idx = (iN * B * t) + (B * i);

            // loop over batch
            for(b = 0; b<B; b++){

                // get sample index 
                sidx = fsidx + b;

                // check if the sample is in the dataset
                if(sidx < dataset->n_samples){ // not really necessary, but it's okay
                    
                    //printf(" > Loading spike in neuron %zu in %zu (from %zu): %d\n", i, (iN + N) * B * (t % L) + (B * i) + b, iN * B * t + (B * i) + b, batch_data->spikes[iN * B * t + (B * i) + b]);
                    // store if the neuron fired in the spk matrix [(iN + N) * B * L] from the batch data matrix [iN * B * T]
                    snn->spk_matrix[spk_mtr_idx + b] = batch_data->spikes[batch_mtr_idx + b]; 
                }
            }
        }

        // reinit the rest of the row
        #pragma omp parallel for num_threads(P) private(i, b, spk_mtr_idx)
        for(i = 0; i<N; i++){

            // compute matrix general index
            spk_mtr_idx = (iN + N) * B * lt + (B * (iN + i));

            for(b = 0; b<B; b++){
                snn->spk_matrix[spk_mtr_idx + b] = 0; 
            }
        }
    }
    #endif
}



// TODO: revise
void STDP_learning(GPU_SNN_t *snn, simulation_configuration_t *conf, size_t t){

    size_t n_tasks, tsk, i, S, P;

    S = snn->n_synapses;
    P = conf->n_process;

    // helpers
    const float decay = 0.95f;//;expf(-1/); // -dt / tau

    // parallelized and vectorized
    #if defined AVX512

        n_tasks = S / 16; 

        const __m512 decay_vec = _mm512_set1_ps(decay);
        const __m512 A_vec = _mm512_set1_ps(0.01);

        #pragma omp parallel for num_threads(P) private(tsk, i)
        for(tsk = 0; tsk<n_tasks; tsk++){

            // load traces
            __m512 preT_vec = _mm512_loadu_ps(&(snn->pre_trace[i]));
            __m512 postT_vec = _mm512_loadu_ps(&(snn->post_trace[i]));

            // update the postsynaptic trace

            // [IF (post_fired) THEN update_post_trace() ENDIF]

            // load postsynaptic neurons indexes for synapses
            __m512i out_vec = _mm512_loadu_si512(&(snn->post_neuron_index[i])); // out neurons

            // index = index - n_input_neurons
            out_vec = _mm512_sub_epi32(out_vec, _mm512_set1_epi32(snn->n_input_neurons));
            
            // get if the neuron [index] fired using gather
            __m512i post_fired_vec = _mm512_i32gather_epi32(out_vec, snn->post_fired, 4); // get which output neurons fired

            // mask depending on which neurons fired [post_fired[i] == 1]
            __mmask16 post_valid_mask = _mm512_cmp_epi32_mask(post_fired_vec, _mm512_setzero_epi32(), _MM_CMPINT_GT);

            // update postsynaptic trace if neuron fired
            postT_vec = _mm512_mask_add_ps(postT_vec, post_valid_mask, postT_vec, _mm512_set1_ps(1.0));
            

            // pre_fired = pre_fired - post_fired, if pre_fired == 1, compute LTD, if post_fired == 1, compute LTP
            // load pre_fired and mask
            __m512i pre_fired_vec = _mm512_loadu_si512(&(snn->pre_fired[i]));
            pre_fired_vec = _mm512_sub_epi32(pre_fired_vec, post_fired_vec);
            __mmask16 pre_valid_mask = _mm512_cmp_epi32_mask(pre_fired_vec, _mm512_setzero_epi32(), _MM_CMPINT_GT);


            // update weights 

            // load dw and w
            __m512 dw_vec = _mm512_loadu_ps(&(snn->dw[i]));
            __m512 w_vec = _mm512_loadu_ps(&(snn->w[i]));


            //
            //    USING MASKS
            //    if post_fired == 1:
            //        LTP
            //    else if pre_fired == 1:
            //        LTD

            // update w and dw (depending on the trace) [w += lr * pre_trace] --> LTP
            dw_vec = 
                _mm512_mask_add_ps(dw_vec, post_valid_mask, dw_vec, _mm512_mul_ps(preT_vec, _mm512_set1_ps(0.01)));
            w_vec = 
                _mm512_mask_add_ps(w_vec, post_valid_mask, w_vec, _mm512_mul_ps(preT_vec, _mm512_set1_ps(0.01)));

            // update w and dw (depending on the trace) [w -= lr * post_trace] --> LTD
            dw_vec = 
                _mm512_mask_sub_ps(dw_vec, pre_valid_mask, dw_vec, _mm512_mul_ps(postT_vec, _mm512_set1_ps(0.01)));
            w_vec = 
                _mm512_mask_sub_ps(w_vec, pre_valid_mask, w_vec, _mm512_mul_ps(postT_vec, _mm512_set1_ps(0.01)));

            // store w and dw
            _mm512_storeu_ps(&(snn->w[i]), w_vec);
            _mm512_storeu_ps(&(snn->dw[i]), dw_vec);


            // compute the traces decays (for all, not only masked)
            preT_vec =_mm512_mul_ps(preT_vec, _mm512_set1_ps(decay));
            postT_vec = _mm512_mul_ps(postT_vec, _mm512_set1_ps(decay));

            // store traces
            _mm512_storeu_ps(&(snn->pre_trace[i]), preT_vec);
            _mm512_storeu_ps(&(snn->post_trace[i]), postT_vec);
        }

        // process remaining
        #pragma omp parallel for num_threads(P) private(i)
        for(i=n_tasks * 16; i<S; i++){

            if(snn->post_fired[snn->post_neuron_index[i]] == 1){
                
                snn->pre_trace[i] += 1.0;

                // traces updated, compute 

                // compute STDP
                snn->dw[i] += 0.01 * snn->pre_trace[i];
                snn->w[i] += 0.01 * snn->pre_trace[i];
            }
            else if(snn->pre_fired[i] == 1){
                
                snn->post_trace[i] += 1.0;

                // traces updated, compute 

                // compute STDP
                snn->dw[i] -= 0.01 * snn->post_trace[i];
                snn->w[i] -= 0.01 * snn->post_trace[i];
            }

            // update traces using the decay
            snn->pre_trace[i] *= decay;
            snn->post_trace[i] *= decay;
        }

        #else
        
        // process synapses
        #pragma omp parallel for num_threads(P)
        for(i=0; i<S; i++){

            // this is slow: non contiguous memory accesses: integrate STDP previously
            if(snn->post_fired[snn->post_neuron_index[i]] == 1){
                
                snn->pre_trace[i] += 1.0;

                // traces updated, compute 

                // compute STDP
                snn->dw[i] += 0.01 * snn->pre_trace[i];
                snn->w[i] += 0.01 * snn->pre_trace[i];
            }
            else if(snn->pre_fired[i] == 1){
                
                snn->post_trace[i] += 1.0;

                // traces updated, compute 

                // compute STDP
                snn->dw[i] -= 0.01 * snn->post_trace[i];
                snn->w[i] -= 0.01 * snn->post_trace[i];
            }

            // update traces using the decay
            snn->pre_trace[i] *= decay;
            snn->post_trace[i] *= decay;
        }

    #endif    
}

void simulate_sample_CPU(GPU_SNN_t *snn, GPU_dataset_t *dataset, simulation_configuration_t *conf, GPU_results_t *results, size_t sidx){

    struct timespec start, end; 
    struct timespec start_step1, end_step1; 
    struct timespec start_step2, end_step2; 
    struct timespec start_step3, end_step3; 
    struct timespec start_step4, end_step4; 
    struct timespec start_step5, end_step5; 
    double elapsed_time, et1 =0.0, et2=0.0, et3=0.0, et4=0.0, et5=0.0;

    // store information in local variables
    size_t N = snn->n_neurons;
    size_t T = conf->time_steps;

    printf(" > Starting sample simulation\n");

    clock_gettime(CLOCK_MONOTONIC, &start);
    clock_gettime(CLOCK_MONOTONIC, &start_step1);

    /* Simulation step 1: initialize network state and load sample */
    reinitialize_LIF_neurons(snn, conf);
    reinitialize_synapses(snn, conf);
    reinitialize_spk_matrix(snn, conf);
    
    load_sample_in_SNN(snn, dataset, conf, sidx);

    clock_gettime(CLOCK_MONOTONIC, &end_step1);
    et1+=(end_step1.tv_sec - start_step1.tv_sec) + (end_step1.tv_nsec - start_step1.tv_nsec) / 1e9;


    /* Simulate time steps */
    for(size_t t=0; t<T; t++){

        
        /* Step 2: process neuron's incoming spikes (input step, independent of the neuron model) */
        clock_gettime(CLOCK_MONOTONIC, &start_step2);
        
        compute_input_current(snn, conf, t);
        
        clock_gettime(CLOCK_MONOTONIC, &end_step2);
        et2+=(end_step2.tv_sec - start_step2.tv_sec) + (end_step2.tv_nsec - start_step2.tv_nsec) / 1e9;
        

        /* Step3: compute V[t] for all neurons */
        clock_gettime(CLOCK_MONOTONIC, &start_step3);
        
        compute_LIF_V(snn, conf, t);
        
        clock_gettime(CLOCK_MONOTONIC, &end_step3);
        et3+=(end_step3.tv_sec - start_step3.tv_sec) + (end_step3.tv_nsec - start_step3.tv_nsec) / 1e9;


        /* Step 4: Fire spikes */
        clock_gettime(CLOCK_MONOTONIC, &start_step4);
        
        process_neuron_firing(snn, conf, results, t);
        
        clock_gettime(CLOCK_MONOTONIC, &end_step4);
        et4+=(end_step4.tv_sec - start_step4.tv_sec) + (end_step4.tv_nsec - start_step4.tv_nsec) / 1e9;


 
        /* Step 5: learning / training */
        clock_gettime(CLOCK_MONOTONIC, &start_step5);
        // run if it is training (REVISE)
        if(conf->learn == 1){
         
            STDP_learning(snn, conf, t);
        }  

        clock_gettime(CLOCK_MONOTONIC, &end_step5);
        et5+=(end_step5.tv_sec - start_step5.tv_sec) + (end_step5.tv_nsec - start_step5.tv_nsec) / 1e9;
    }

    clock_gettime(CLOCK_MONOTONIC, &end);
    elapsed_time = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // print execution time    
    printf(" > Finished in %lf! (s1 = %lf s2 = %lf s3 = %lf s4 = %lf, s5 = %lf)\n", 
            elapsed_time, et1, et2, et3, et4, et5);

    printf(" > Generated number of spikes per neuron: ");
    for(size_t i = 0; i<N; i++){
        printf("%d ", results->n_spks[i]);
    }
    printf("\n");
    fflush(stdout);
}

void simulate_samples_CPU(GPU_SNN_t *snn, GPU_dataset_t *dataset, simulation_configuration_t *conf, GPU_results_t *results){

    size_t s;

    for(s = 0; s<dataset->n_samples; s++){

        simulate_sample_CPU(snn, dataset, conf, results, s);
    }
}


/* Same as above but for batches */

void simulate_batch_CPU(GPU_SNN_t *snn, GPU_dataset_t *dataset, simulation_configuration_t *conf, GPU_results_t *results, size_t bidx, int print_data){

    // variables to store execution times
    struct timespec start, end; 
    struct timespec start_step1, end_step1; 
    struct timespec start_step2, end_step2; 
    struct timespec start_step3, end_step3; 
    struct timespec start_step4, end_step4; 
    struct timespec start_step5, end_step5; 
    struct timespec start_load_batch, end_load_batch;
    double elapsed_time, et1 =0.0, et2=0.0, et3=0.0, et4=0.0, et5=0.0, et_load = 0.0;

    // store information in local variables
    size_t T = (size_t)conf->time_steps;
    size_t B = conf->batch_size;
    size_t LT = snn->LT;

    // indices for looping
    size_t i, t, lt;

    clock_gettime(CLOCK_MONOTONIC, &start);
    clock_gettime(CLOCK_MONOTONIC, &start_step1);

    /* Simulation step 1: initialize network state and load sample */
    reinitialize_LIF_neurons_batch(snn, conf);
    reinitialize_synapses_batch(snn, conf);
    reinitialize_spk_matrix_batch(snn, conf);
    

    clock_gettime(CLOCK_MONOTONIC, &end_step1);
    et1+=(end_step1.tv_sec - start_step1.tv_sec) + (end_step1.tv_nsec - start_step1.tv_nsec) / 1e9;

    clock_gettime(CLOCK_MONOTONIC, &start_load_batch);

    // initialize batch matrix
    tmp_batch_cpu_t *batch_data = initialize_batch_matrix(snn, dataset, conf, bidx);
    clock_gettime(CLOCK_MONOTONIC, &end_load_batch);
    et_load+=(end_load_batch.tv_sec - start_load_batch.tv_sec) + (end_load_batch.tv_nsec - start_load_batch.tv_nsec) / 1e9;

    /* Simulate time steps */
    for(t=0; t<T; t++){

        // compute local time step
        lt = t % LT;

        //print_weights_3D(snn, B);

        // load batch partially
        clock_gettime(CLOCK_MONOTONIC, &start_load_batch);

        load_batch_time_step_in_SNN_batch(snn, dataset, conf, batch_data, bidx, t);
            
        clock_gettime(CLOCK_MONOTONIC, &end_load_batch);
        et_load+=(end_load_batch.tv_sec - start_load_batch.tv_sec) + (end_load_batch.tv_nsec - start_load_batch.tv_nsec) / 1e9;

        /* Step 2: process neuron's incoming spikes (input step, independent of the neuron model) */
        clock_gettime(CLOCK_MONOTONIC, &start_step2);
        
        compute_input_current_batch(snn, conf, lt, t);
        
        clock_gettime(CLOCK_MONOTONIC, &end_step2);
        et2+=(end_step2.tv_sec - start_step2.tv_sec) + (end_step2.tv_nsec - start_step2.tv_nsec) / 1e9;
        

        /* Step3: compute V[t] for all neurons */
        clock_gettime(CLOCK_MONOTONIC, &start_step3);
        
        compute_LIF_V_batch(snn, conf, lt, t);
        
        clock_gettime(CLOCK_MONOTONIC, &end_step3);
        et3+=(end_step3.tv_sec - start_step3.tv_sec) + (end_step3.tv_nsec - start_step3.tv_nsec) / 1e9;


        /* Step 4: Fire spikes */
        clock_gettime(CLOCK_MONOTONIC, &start_step4);
        
        process_neuron_firing_batch(snn, conf, results, lt, t);
        
        clock_gettime(CLOCK_MONOTONIC, &end_step4);
        et4+=(end_step4.tv_sec - start_step4.tv_sec) + (end_step4.tv_nsec - start_step4.tv_nsec) / 1e9;


        /* Step 5: learning / training */
        clock_gettime(CLOCK_MONOTONIC, &start_step5);
        // run if it is training (REVISE)
        if(conf->learn){
         
            process_trace_based_STDP_batch(snn, conf);
        } 

        clock_gettime(CLOCK_MONOTONIC, &end_step5);
        et5+=(end_step5.tv_sec - start_step5.tv_sec) + (end_step5.tv_nsec - start_step5.tv_nsec) / 1e9;
        
    }

    // update weights if learning
    clock_gettime(CLOCK_MONOTONIC, &start_step5);
    if(conf->learn)
        update_weights_cpu(snn, conf);
    clock_gettime(CLOCK_MONOTONIC, &end_step5);
    et5+=(end_step5.tv_sec - start_step5.tv_sec) + (end_step5.tv_nsec - start_step5.tv_nsec) / 1e9;


    clock_gettime(CLOCK_MONOTONIC, &end);
    elapsed_time = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;


    // store execution times and number of spikes
    results->t = elapsed_time;
    results->t_in = et2;
    results->t_v = et3;
    results->t_out = et4;
    results->t_reinit = et1;
    results->t_load = et_load;
    results->t_learn = et5;

    // print execution time
    //if((bidx + 1) % 10 == 0){    
    if(print_data == 1){
        // print execution times
        printf(" > Finished in %lf! (s_rnt = %lf s_load = %lf s_in = %lf s_v = %lf s_out = %lf, s_tr = %lf)\n", 
            elapsed_time, et1, et_load, et2, et3, et4, et5);
        fflush(stdout);

        // print generated spikes
        printf(" > Generated number of spikes per neuron: ");
        for(size_t b = 0; b<conf->batch_size; b++){

            for(i = 0; i<snn->n_neurons; i++){
                printf("%d ", results->n_spks[i * B + b]);
                results->n_spks[i * B + b] = 0;
            }
            printf("\n");
        }
        printf("\n");
    }
}

void simulate_batches(GPU_SNN_t *snn, GPU_dataset_t *dataset, simulation_configuration_t *conf, GPU_results_t *results){

    struct timespec start, end; 
    struct timespec start_uw, end_uw; 
    double elapsed_time;

    // TODO: refactorize number of batches computation
    // compute number of batches
    size_t n_batches = dataset->n_samples / conf->batch_size;
    size_t r_samples = dataset->n_samples % conf->batch_size;
    // there are remaining sampels
    if(r_samples > 0){

        n_batches += 1; // last batch contains less samples
    }

    // copy necessary data of the network to compute by batches
    cpy_snn(snn, conf);

    // initialize results struct to store batch results
    GPU_results_t *batch_results = initialize_batch_results_cpu(conf, snn->n_neurons, conf->batch_size, 1);
    simulate_batch_CPU(snn, dataset, conf, batch_results, 0, 1);


    // simulate and print first batch results

    // start timer
    clock_gettime(CLOCK_MONOTONIC, &start);

    // loop over batches
    for(size_t b = 0; b<n_batches; b++){
        
        if((b+1) % 100 == 0){
            printf(" In Batch %zu\n", b+1);
            fflush(stdout);
        }

        // reinitialize results struct
        reinitialize_batch_results_cpu(batch_results, conf, snn->n_neurons, conf->batch_size, 1);

        // simulate batch
        simulate_batch_CPU(snn, dataset, conf, batch_results, b, 0);

        // accumulate batch execution times in general results struct
        acc_batch_execution_times(results, batch_results);
    }

    clock_gettime(CLOCK_MONOTONIC, &end);
    elapsed_time = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;
    printf(" > Finished in %lf!\n", elapsed_time); 


    // print information in results struct
    printf("   > Total execution time: %.3lf | Mean per batch %.3lf | Mean per time step %.3lf\n", 
                results->t, results->t / (float)n_batches, results->t / (float)n_batches / (float)conf->time_steps);
    printf("   > Total in time:        %.3lf | Mean per batch %.3lf | Mean per time step %.3lf\n", 
                results->t_in, results->t_in / (float)n_batches, results->t_in / (float)n_batches / (float)conf->time_steps);
    printf("   > Total v time:         %.3lf | Mean per batch %.3lf | Mean per time step %.3lf\n", 
                results->t_v, results->t_v / (float)n_batches, results->t_v / (float)n_batches / (float)conf->time_steps);
    printf("   > Total out time:       %.3lf | Mean per batch %.3lf | Mean per time step %.3lf\n", 
                results->t_out, results->t_out / (float)n_batches, results->t_out / (float)n_batches / (float)conf->time_steps);
    printf("   > Total reinit time:    %.3lf | Mean per batch %.3lf | Mean per time step %.3lf\n", 
                results->t_reinit, results->t_reinit / (float)n_batches, results->t_reinit / (float)n_batches / (float)conf->time_steps);
    printf("   > Total load time:      %.3lf | Mean per batch %.3lf | Mean per time step %.3lf\n", 
                results->t_load, results->t_load / (float)n_batches, results->t_load / (float)n_batches / (float)conf->time_steps);
    printf("   > Total learn time:     %.3lf | Mean per batch %.3lf | Mean per time step %.3lf\n", 
                results->t_learn, results->t_learn / (float)n_batches, results->t_learn / (float)n_batches / (float)conf->time_steps);
}