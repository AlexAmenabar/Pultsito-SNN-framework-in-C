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

                // constants in registers
                __m512 alpha_v = _mm512_set1_ps(alpha);
                __m512 beta_v  = _mm512_set1_ps(beta);

                size_t i;
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
                __m512i ones_vec = _mm512_set1_epi32(1);
                for(size_t i = 0; i + 15 < cpu_snn->n_neurons; i+=16){

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


static inline __mmask16 get_batch_mask(int batch_size)
{
    // batch_size is 8 or 16
    return (batch_size == 8) ? 0x00FF : 0xFFFF;
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




void compute_input_current_batch(GPU_SNN_t *snn, simulation_configuration_t *conf, size_t t){
    
    size_t N = snn->n_neurons;
    size_t iN = snn->n_input_neurons;
    size_t P = conf->n_process;
    size_t B = conf->batch_size;
    size_t max_spikes = snn->max_spikes;
    size_t i = 0;
    size_t b = 0;
    size_t gindex, lindex;
    

    #if defined AVX512
    {

        size_t n_tasks = B / 16;
        if(n_tasks == 0) // if B < 16
            n_tasks = 1;
        
        // create mask if batch size is smaller than 16
        const __mmask16 m = get_mask(B);

        #pragma omp parallel for num_threads(P) private(gindex)
        for(i = 0; i<N; i++){

            // first neuron index
            gindex = i * B;

            // initialize array I with zeros
            __m512 I_vec = _mm512_set1_ps(0.0);

            // check if neuron is in refractary period and create mask
            __m512i refract_vec = _mm512_maskz_loadu_epi32(m, &(snn->r_period_remain[gindex]));
            __mmask16 rft_mask = _mm512_cmp_epi32_mask(refract_vec, _mm512_setzero_si512(), _MM_CMPINT_LE) & m; // if lower or eual, process

        
            // load number of input synapses and first input synapse offset 
            size_t iS = snn->n_neuron_input_synapses[gindex]; 
            size_t base_synapse = snn->neuron_input_synapses_offset[gindex] * B; // there are several copies of each synapse, so the offset depends

            //printf(" > IN neuron %zu, starting processing synapses\n", i);
            //fflush(stdout);

            for(size_t j = 0; j<iS; j++){

                size_t synapse_index = base_synapse + j * B; // next synapse index
                size_t in_neuron_index = snn->pre_neuron_index[synapse_index]; // this is the in neuron index scaled (the index for the first in the batch)

                // compute spike time: the same for all in the batch
                int spk_time = (int)t - snn->delay[synapse_index];
                int position = spk_time % (int)snn->max_spikes; // get position in matrix

                //printf(" > Spike time loaded\n");
                //fflush(stdout);
                if(spk_time >= 0){

                    // load whether neuron received a spike
                    
                    // [T * B * (iN + N)]
                    //size_t index = (max_spikes * in_neuron_index * B) + (B * (size_t)spk_time); // spk_time -> position
                    
                    // [(iN + N) * B * T]
                    size_t index = ((iN + N) * B * (size_t)spk_time) + (in_neuron_index * B);

                    __m512i spikes_i = _mm512_maskz_loadu_epi32(m, &(snn->spk_matrix[index]));
                    
                    // convert to float
                    __m512 spikes = _mm512_maskz_cvtepi32_ps(m, spikes_i);

                    // load weights
                    __m512 w_vec = _mm512_maskz_loadu_ps(m, &(snn->w[synapse_index]));
                    
                    // update I
                    I_vec = _mm512_mask_add_ps(I_vec, m, I_vec, _mm512_mask_mul_ps(w_vec, m, w_vec, spikes));  

                    // update pre_foired
                }
            }

            _mm512_mask_storeu_ps(&(snn->arrI[gindex]), m, I_vec);
            //_mm512_mask_storeu_epi32(&(snn->inR[gindex]), m, refract_vec);
        }
    }
    #else
    {

        #pragma omp parallel for num_threads(P) private (gindex, lindex) 
        for(i=0; i<N; i++){
            gindex = i * B; // there are B copies of each neuron, so multiply to get the first neuron copy index

            //if(b==1)
            //    printf(" > > Processing neuron %zu (gindex = %zu)\n", i, gindex);

            for(b = 0; b<B; b++){

                float I = 0.0;
                lindex = gindex + b; // index of the neuron copy

                if(snn->r_period_remain[lindex] <= 0){

                    // variables declaration
                    size_t synapse_index, in_neuron_index, spk_time_index;
                    int delay, spk_time; 
                    float spk;    
                    size_t j = 0; 

                    // get if it is in refractary period
                    snn->inR[lindex] = (snn->r_period_remain[lindex] <= 0); // store wether the neuron is in refractary period // necessary?

                    // get number of input synapses and offset
                    size_t iS = snn->n_neuron_input_synapses[lindex]; 
                    size_t base_synapse = snn->neuron_input_synapses_offset[lindex] * B + b; // there are several copies of each synapse, so the offset depends

                    //if(b==1)
                    //    printf(" > >>> lindex %zu, iS %zu, base_synapse %zu\n", lindex, iS, base_synapse);

                    // loop over input synapses
                    for(j=0; j<iS; j++){

                        synapse_index = base_synapse + j * B; // scale base synapse, since now there are B copies of each one
                        
                        in_neuron_index = snn->pre_neuron_index[synapse_index]; 
                        //if(b==1)
                        //    printf(" > >>> >>> synapse index = %zu / In neuron index = %zu\n", synapse_index, in_neuron_index);

                        delay = snn->delay[synapse_index]; // no copies
                        spk_time = (int)t - delay; // actual position

                        // check if spike is bigger than 0
                        if(spk_time >=0){ 
                            
                            // T * N
                            // get wether the neuron fired in the 3D matrix (used as 2D matrix)
                            //if(b==1)
                            //    printf(" > >>> >>> Accessing matrix position %zu\n", (max_spikes * in_neuron_index) + (B * (size_t)spk_time) + b);
                            
                            // [T * B * (iN + N)]
                            spk = (float)(snn->spk_matrix[(max_spikes * B * in_neuron_index) + (B * (size_t)spk_time) + b]);    
                            
                            // [(iN + N) * B * T]
                            //spk = (float)(snn->spk_matrix[((iN + N) * B * (size_t)spk_time) + (B * in_neuron_index) + b]);    
                            

                            I += snn->w[synapse_index] * spk; // spk = 0 / 1
                            snn->pre_fired[synapse_index] = (int)spk;
                        }
                    }
                }                
                snn->arrI[lindex] = I;
                snn->inR[lindex] = snn->r_period_remain[lindex] <= 0;
            }
        }
    }
    #endif
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
    #if defined AVX512

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

    #else

    #pragma omp parallel for num_threads(P)
    for (i = 0; i < N; i++) {

        // check if the neuron is in refractary period
        if (snn->r_period[i] == 1){
        
            // compute v[t]
            snn->v[i] = alpha * snn->v[i] + beta * snn->v_rest[i] + snn->arrI[i];
        }
    }

    #endif
}


void compute_LIF_V_batch(GPU_SNN_t *snn, simulation_configuration_t *conf, size_t t){

    // get general information
    size_t max_spikes = snn->max_spikes;
    size_t N = snn->n_neurons;
    size_t iN = snn->n_input_neurons;
    size_t P = conf->n_process;
    size_t B = conf->batch_size;

    size_t n_tasks = N / 16, tsk, i = 0, b = 0;

    // helpers
    const float alpha = 0.95f;
    const float beta = 0.05f;


    size_t gindex, lindex;


    #if defined AVX512
    {

        size_t n_tasks = N * B / 16;
        if(n_tasks == 0) // if B < 16
            n_tasks = 1;


        // initialize constants
        __m512 avec = _mm512_set1_ps(alpha);
        __m512 bvec = _mm512_set1_ps(beta);

        #pragma omp parallel for num_threads(P) private(i)
        for(size_t tsk = 0; tsk<n_tasks; tsk++){   
            
            i = tsk * 16;

            // load refractary period and create mask
            __m512i refract_vec = _mm512_loadu_si512(&(snn->r_period_remain[i]));
            __mmask16 rft_mask = _mm512_cmp_epi32_mask(refract_vec, _mm512_setzero_si512(), _MM_CMPINT_LE); // if lower or eual, process

            // load I, rest and v
            __m512 v_vec = _mm512_loadu_ps(&(snn->v[i]));
            __m512 v_rest_vec = _mm512_loadu_ps(&(snn->v_rest[i]));
            __m512 I_vec = _mm512_loadu_ps(&(snn->arrI[i]));

            // compute dynamics
            __m512 tmp = _mm512_add_ps(I_vec, _mm512_mul_ps(v_rest_vec, bvec));
            __m512 rhs = _mm512_fmadd_ps(v_vec, avec, tmp);  // alpha*v + beta*v_rest + I
            v_vec = _mm512_mask_mov_ps(v_vec, rft_mask, rhs); // only update active lanes

            // store v
            _mm512_storeu_ps(&(snn->v[i]), v_vec);
        }

        // process remaining
        #pragma omp parallel for num_threads(P)
        for(i = n_tasks * 16; i<N*B; i++){
            
            if (snn->r_period_remain[i] <= 0){

                snn->v[i] = alpha * snn->v[i] + beta * snn->v_rest[i] + snn->arrI[i];
            }
        }
    }
    #else
    {
        #pragma omp parallel for num_threads(P) private(lindex, gindex)
        for (i = 0; i < N; i++) {

            gindex = i * B;
            for(b = 0; b<B; b++){

                lindex = gindex + b;
                if (snn->r_period_remain[lindex] <= 0){
                    
                    //printf(" > Neuron %zu I = %f\n", i, snn->arrI[lindex]);
                    snn->v[lindex] = alpha * snn->v[lindex] + beta * snn->v_rest[lindex] + snn->arrI[lindex];
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


void process_neuron_firing_batch(GPU_SNN_t *snn, simulation_configuration_t *conf, GPU_results_t *results, size_t t){

    // get general information
    size_t max_spikes = (size_t)snn->max_spikes;
    size_t N = (size_t)snn->n_neurons;
    size_t iN = (size_t)snn->n_input_neurons;
    size_t P = (size_t)conf->n_process;
    size_t B = conf->batch_size;


    size_t gindex, lindex, i, b;


    #if defined AVX512
    {

        size_t n_tasks = B / 16;
        if(n_tasks == 0) // if B < 16
            n_tasks = 1;
        
        // create mask if batch size is smaller than 16
        const __mmask16 m = get_mask(B);

        // process neurons in parallel and vectorize batch inside each neuron
        #pragma omp parallel for num_threads(P) private(i, gindex)
        for(i = 0; i<N; i++){

            gindex = i * B;
            //size_t idx = (((i + iN) * B) * max_spikes) + (t * B); // index for the first neuron to write
            size_t idx = ((iN + N) * B * t) + (iN * B) + (i * B);

            // load and update refractary period, update and store
            __m512i r_period_vec = _mm512_maskz_loadu_epi32(m, &(snn->r_period_remain[gindex]));
            r_period_vec = _mm512_sub_epi32(r_period_vec, _mm512_set1_epi32(1));            
            _mm512_mask_storeu_epi32(&(snn->r_period_remain[gindex]), m, r_period_vec);

            // load v and v_thresh
            __m512 v_vec = _mm512_maskz_loadu_ps(m, &(snn->v[gindex]));
            __m512 v_thresh_vec = _mm512_maskz_loadu_ps(m, &(snn->v_thresh[gindex]));
            __m512 v_rest_vec = _mm512_maskz_loadu_ps(m, &(snn->v_rest[gindex]));

            // compare v and v_thresh
            __mmask16 fire_mask = _mm512_cmp_ps_mask(v_vec, v_thresh_vec, _MM_CMPINT_GE) & m;

            // store 1 for neurons that fired
            _mm512_mask_storeu_epi32(&(snn->spk_matrix[idx]), fire_mask, _mm512_set1_epi32(1));

            // v = v_rest
            _mm512_mask_storeu_ps(&(snn->v[gindex]), fire_mask, v_rest_vec);

            // r_period_remain = r_period
            _mm512_mask_storeu_epi32(&(snn->r_period_remain[gindex]), fire_mask, _mm512_maskz_loadu_epi32(fire_mask, &(snn->r_period[gindex])));

            // store results too
            __m512i n_spk_vec = _mm512_maskz_loadu_epi32(fire_mask, &(results->n_spks[gindex]));
            n_spk_vec = _mm512_add_epi32(n_spk_vec, _mm512_set1_epi32(1));
            _mm512_mask_storeu_epi32(&(results->n_spks[gindex]), fire_mask, n_spk_vec);
        }
    }
    #else
    {
        
        #pragma omp parallel for num_threads(P) private(gindex, lindex)
        for(i = 0; i<N; i++){
                
            gindex = i * B; // batch first neuron index

            // fst.elem.id.: pass input neurons +  
            
            // [T * B * (iN + N)]
            size_t idx = ((gindex + iN * B) * max_spikes) + (t * B); // index of the first neuron in the spike matrix in time step t

            // [(iN + N) * B * T]
            //size_t idx = ((iN + N) * B * t) + (gindex);

            for(b = 0; b<B; b++){
            
                lindex = gindex + b; // neuron in batch index

                // reduce refractary period
                snn->r_period_remain[lindex] --;

                // check whether fired
                if(snn->v[lindex] >= snn->v_thresh[lindex]){
            
                    // T * N
                    size_t idx_neuron = idx + b; // row start at

                    snn->spk_matrix[idx_neuron] = 1;

                    //printf(" > >>> Neuron %zu (gindex %zu / lindex %zu) fired in T = %zu at position = %zu\n", i, gindex, lindex, t, idx_neuron);

                    // reinit neuron values
                    snn->r_period_remain[lindex] = snn->r_period[lindex];
                    snn->v[lindex] = snn->v_rest[lindex]; // reinit v_rest

                    results->n_spks[lindex] += 1;

                    //if(b == 1)
                    //    printf(" > Neuron %zu fired!\n", i);
                }

                snn->post_fired[lindex] = snn->v[lindex] >= snn->v_thresh[lindex];
            }
        }
    }
    #endif
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
    
    size_t i, j, b, N, P;
    size_t B, gindex, lindex;

    N = snn->n_neurons;
    P = conf->n_process;
    B = conf->batch_size;
    

    #ifdef AVX512
    {
        // compute number of tasks (N * B, since all neurons of all copies are reinitialized)
        size_t n_tasks = N * B / 16;
        size_t task;

        if(n_tasks == 0) // if n_tasks < 16
            n_tasks = 1;         

        // process tasks in parallel if openMP is defined
        #pragma omp parallel for num_threads(P) private(i, j, b, task)
        for(task = 0; task<n_tasks; task++){

            // get the index of the batch element
            b = task * 16;

            // load v_rest
            __m512 v_rest_vec = _mm512_loadu_ps(&(snn->v_rest[b])); 

            // store v_rest in v
            _mm512_storeu_ps(&(snn->v[b]), v_rest_vec);

            // reinit r_period_remain and post_fired
            _mm512_storeu_si512(&(snn->r_period_remain[b]), _mm512_set1_epi32(0));
            _mm512_storeu_si512(&(snn->post_fired[b]), _mm512_set1_epi32(0));
        }

        // process remaining
        for(j = n_tasks * 16; j < N * B; j++){
            
            // loop over copies in batch (serial)
            snn->v[j] = snn->v_rest[j];
            snn->r_period_remain[j] = 0;
            snn->post_fired[j] = 0;
        }
    }
    #else
    {
        // loop over neurons
        #pragma omp parallel for num_threads(P) private(i, b, gindex, lindex)
        for(i = 0; i < N; i++){

            gindex = i * B;
            
            // loop over copies in batch (serial)
            for(b = 0; b<B; b++){

                lindex = gindex + b;

                snn->v[lindex] = snn->v_rest[lindex];
                snn->r_period_remain[lindex] = 0;
                snn->post_fired[lindex] = 0;
            }
        }
    }
    #endif
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


    #ifdef AVX512
    {
        // compute number of tasks (N * B, since all neurons of all copies are reinitialized)
        size_t n_tasks = S * B / 16;
        size_t task;

        if(n_tasks == 0) // if n_tasks < 16
            n_tasks = 1;         

        // process tasks in parallel if openMP is defined
        #pragma omp parallel for num_threads(P) private(i, j, b, task)
        for(task = 0; task<n_tasks; task++){

            // get the index of the batch element
            b = task * 16;

            // load init_w
            __m512 init_w_vec = _mm512_loadu_ps(&(snn->init_w[b])); 

            // store init_w in w
            _mm512_storeu_ps(&(snn->w[b]), init_w_vec);

            // reinit pre_fired 
            _mm512_storeu_si512(&(snn->pre_fired[b]), _mm512_set1_epi32(0));
        }

        // process remaining
        for(j = n_tasks * 16; j < S * B; j++){
            
            // loop over copies in batch (serial)
            snn->w[j] = snn->init_w[j];
            snn->pre_fired[j] = 0;
        }
    }
    #else
    {
        // reinitialize synapses 
        #pragma omp parallel for num_threads(P) private(i, b, gindex, lindex)
        for(i = 0; i<S; i++){
            
            gindex = i * B;

            // loop over the batch
            for(b = 0; b<B; b++){

                lindex = gindex + b;

                snn->w[lindex] = snn->init_w[lindex];
                snn->pre_fired[lindex] = 0;
            }
        }
    }
    #endif
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
    
    size_t i, j, N, iN, P, max_spikes, b, B, lindex, gindex;
    
    N = snn->n_neurons;
    iN = snn->n_input_neurons;
    P = conf->n_process;
    max_spikes = snn->max_spikes;
    B = conf->batch_size;

    // reinitialize spike matrix
    #pragma omp parallel for num_threads(P) private(i, j, b, gindex, lindex)
    for(i=0; i < N + iN; i++){
        
        for(j = 0; j<max_spikes; j++){
        
            gindex = i * max_spikes * B + j * B;

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

    // NO VECTORIZABLE

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

                    // matrix [T * B * (iN + N)]
                    //snn->spk_matrix[(snn->max_spikes * B * i) + (B * t) + b] = 1;

                    // [(iN + N) * B * T]
                    snn->spk_matrix[((iN + N) * B * t) + (B * i) + b] = 1;

                    //if(b == 1){
                    //    printf(" > > Loading spike in %zu\n", (snn->max_spikes * B * i) + (B * dataset->spikes[start_feature + j]) + b);
                    //}
                }   
            }
        }
    }
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

void simulate_batch_CPU(GPU_SNN_t *snn, GPU_dataset_t *dataset, simulation_configuration_t *conf, GPU_results_t *results, size_t bidx){

    struct timespec start, end; 
    struct timespec start_step1, end_step1; 
    struct timespec start_step2, end_step2; 
    struct timespec start_step3, end_step3; 
    struct timespec start_step4, end_step4; 
    struct timespec start_step5, end_step5; 
    double elapsed_time, et1 =0.0, et2=0.0, et3=0.0, et4=0.0, et5=0.0;

    // store information in local variables
    size_t T = (size_t)conf->time_steps;
    size_t B = conf->batch_size;

    // indices for looping
    size_t i, t;

    //printf(" > Starting sample simulation\n");

    clock_gettime(CLOCK_MONOTONIC, &start);
    clock_gettime(CLOCK_MONOTONIC, &start_step1);

    /* Simulation step 1: initialize network state and load sample */
    reinitialize_LIF_neurons_batch(snn, conf);
    reinitialize_synapses_batch(snn, conf);
    reinitialize_spk_matrix_batch(snn, conf);
    
    load_sample_in_SNN_batch(snn, dataset, conf, bidx);

    clock_gettime(CLOCK_MONOTONIC, &end_step1);
    et1+=(end_step1.tv_sec - start_step1.tv_sec) + (end_step1.tv_nsec - start_step1.tv_nsec) / 1e9;


    /* Simulate time steps */
    for(t=0; t<T; t++){

        /* Step 2: process neuron's incoming spikes (input step, independent of the neuron model) */
        clock_gettime(CLOCK_MONOTONIC, &start_step2);
        
        compute_input_current_batch(snn, conf, t);
        
        clock_gettime(CLOCK_MONOTONIC, &end_step2);
        et2+=(end_step2.tv_sec - start_step2.tv_sec) + (end_step2.tv_nsec - start_step2.tv_nsec) / 1e9;
        

        /* Step3: compute V[t] for all neurons */
        clock_gettime(CLOCK_MONOTONIC, &start_step3);
        
        compute_LIF_V_batch(snn, conf, t);
        
        clock_gettime(CLOCK_MONOTONIC, &end_step3);
        et3+=(end_step3.tv_sec - start_step3.tv_sec) + (end_step3.tv_nsec - start_step3.tv_nsec) / 1e9;


        /* Step 4: Fire spikes */
        clock_gettime(CLOCK_MONOTONIC, &start_step4);
        
        process_neuron_firing_batch(snn, conf, results, t);
        
        clock_gettime(CLOCK_MONOTONIC, &end_step4);
        et4+=(end_step4.tv_sec - start_step4.tv_sec) + (end_step4.tv_nsec - start_step4.tv_nsec) / 1e9;


 
        /* Step 5: learning / training */
        clock_gettime(CLOCK_MONOTONIC, &start_step5);
        // run if it is training (REVISE)
        /*if(conf->learn == 1){
         
            STDP_learning(snn, conf, t);
        }  */

        clock_gettime(CLOCK_MONOTONIC, &end_step5);
        et5+=(end_step5.tv_sec - start_step5.tv_sec) + (end_step5.tv_nsec - start_step5.tv_nsec) / 1e9;
    }

    clock_gettime(CLOCK_MONOTONIC, &end);
    elapsed_time = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // print execution time
    if(bidx % 10 == 0){    
        printf(" > Finished in %lf! (s1 = %lf s2 = %lf s3 = %lf s4 = %lf, s5 = %lf)\n", 
            elapsed_time, et1, et2, et3, et4, et5);
        fflush(stdout);
    }
    //printf(" > Generated number of spikes per neuron: ");
    //for(size_t b = 0; b<conf->batch_size; b++){

    //    for(i = 0; i<snn->n_neurons; i++){
    //        printf("%d ", results->n_spks[i * B + b]);
    //        results->n_spks[i * B + b] = 0;
    //    }
    //    printf("\n");
    //}
    //printf("\n");
}


void update_weights_cpu(GPU_SNN_t *snn, size_t batch_size){

    size_t i;

    for(i = 0; i<snn->n_synapses; i++){

        float dw = 0.0;

        // sum dw of all network copies
        for(size_t b = 0; b<batch_size; b++){

            dw += snn->dw[i * batch_size + b];
        }

        // compute mean dw
        dw /= batch_size;

        // update initial w and w for all copies
        for(size_t b = 0; b<batch_size; b++){

            snn->dw[i * batch_size + b] = 0.0; // reinitialize dw
            snn->init_w[i * batch_size + b] += dw; // update initial w for incoming batches
            snn->w[i * batch_size + b] = snn->init_w[i * batch_size + b]; // update w
        }    
    }
}

void simulate_batches(GPU_SNN_t *snn, GPU_dataset_t *dataset, simulation_configuration_t *conf, GPU_results_t *results){

    struct timespec start, end; 

    // TODO: improve in the future

    // compute number of batches
    size_t n_batches = dataset->n_samples / conf->batch_size;
    size_t r_samples = dataset->n_samples % conf->batch_size;

    // there are remaining sampels
    if(r_samples > 0){

        n_batches += 1; // last batch contains less samples
    }

    // create network copies
    cpy_snn(snn, conf);
    //printf(" > >> Network copied %zu times\n", conf->batch_size);
    //fflush(stdout);

    //print_networks(snn, conf);

    clock_gettime(CLOCK_MONOTONIC, &start);

    // loop over batches
    for(size_t b = 0; b<n_batches; b++){

        // simulate batch
        simulate_batch_CPU(snn, dataset, conf, results, b);

        if(b % 100 == 0){
            printf(" In Batch %zu\n", b);
            fflush(stdout);
        }
        // update weights
        // update_weights();
    }

    clock_gettime(CLOCK_MONOTONIC, &end);
    double elapsed_time = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;
    printf(" > Finished in %lf!\n", elapsed_time); 
}