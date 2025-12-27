#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

#include "snn_library.h"
#include "load_data.h"

#include "simulations/simulations.h"

#include "neuron_models/lif_neuron.h"

#include "training_rules/stdp.h"

#include "helpers.h"

#ifdef CUDA
    #include "cuda/GPU_simulations.cuh"

    #include "neuron_models/GPU_lif_neuron.cuh"
    #include "cuda/cuda_helpers.cuh"
#endif

#include <immintrin.h>  // AVX intrinsics
#include <omp.h>



/*
    Functions to be moved

*/




void cpy_snn(GPU_SNN_t *snn, int n_cps, simulation_configuration_t *conf){

    size_t i, j, t;

    // temporal variables to store values
    float *tmp_v, *tmp_thresh, *tmp_rest, *tmp_w, *tmp_init_w, *tmp_dw, *tmp_pre_trace, *tmp_post_trace;
    int *tmp_rperiod, *tmp_rperiod_remain, *tmp_res, *tmp_pstfired, *tmp_prefired, *tmp_delay, *tmp_lr, *tmp_spk_mtrx;
    size_t *tmp_in_off, *tmp_pre_neuron_index, *tmp_post_neuron_index, *tmp_n_insyn;
    // allocate neural properties
    tmp_v                         = snn->v;
    tmp_thresh                    = snn->v_thresh;
    tmp_rest                      = snn->v_rest;
    tmp_rperiod                   = snn->r_period;
    tmp_rperiod_remain            = snn->r_period_remain;
    tmp_res                       = snn->res;
    tmp_pstfired                  = snn->post_fired;
    tmp_n_insyn                   = snn->n_neuron_input_synapses;
    tmp_in_off                    = snn->neuron_input_synapses_offset;

    tmp_w                         = snn->w;
    tmp_init_w                    = snn->init_w;
    tmp_dw                        = snn->dw;
    tmp_pre_trace                 = snn->pre_trace;
    tmp_post_trace                = snn->post_trace;
    tmp_prefired                  = snn->pre_fired;
    tmp_delay                     = snn->delay;
    tmp_lr                        = snn->lr;
    tmp_pre_neuron_index          = snn->pre_neuron_index;
    tmp_post_neuron_index         = snn->post_neuron_index;

    tmp_spk_mtrx                  = snn->spk_matrix;


    /* Memory allocation */
    snn->v                            = (float*)malloc(snn->n_neurons * n_cps * sizeof(float));
    snn->v_thresh                     = (float*)malloc(snn->n_neurons * n_cps * sizeof(float));
    snn->v_rest                       = (float*)malloc(snn->n_neurons * n_cps * sizeof(float));
    snn->r_period                     = (int*)malloc(snn->n_neurons * n_cps * sizeof(float));
    snn->r_period_remain              = (int*)malloc(snn->n_neurons * n_cps * sizeof(int));
    snn->res                          = (int*)malloc(snn->n_neurons * n_cps * sizeof(int));
    snn->post_fired                   = (int*)malloc(snn->n_neurons * n_cps * sizeof(int));
    snn->n_neuron_input_synapses      = (size_t*)malloc(snn->n_neurons * n_cps * sizeof(size_t));
    snn->neuron_input_synapses_offset = (size_t*)malloc(snn->n_neurons * n_cps * sizeof(size_t));

    snn->w                            = (float*)malloc(snn->n_synapses * n_cps * sizeof(float));
    snn->init_w                       = (float*)malloc(snn->n_synapses * n_cps * sizeof(float));
    snn->dw                           = (float*)malloc(snn->n_synapses * n_cps * sizeof(float));
    snn->pre_trace                    = (float*)malloc(snn->n_synapses * n_cps * sizeof(float));
    snn->post_trace                   = (float*)malloc(snn->n_synapses * n_cps * sizeof(float));
    snn->pre_fired                    = (int*)malloc(snn->n_synapses * n_cps * sizeof(int));
    snn->delay                        = (int*)malloc(snn->n_synapses * n_cps * sizeof(int));
    snn->lr                           = (int*)malloc(snn->n_synapses * n_cps * sizeof(int));
    snn->pre_neuron_index             = (size_t*)malloc(snn->n_synapses * n_cps * sizeof(size_t));
    snn->post_neuron_index            = (size_t*)malloc(snn->n_synapses * n_cps * sizeof(size_t));

    snn->spk_matrix                   = (int*)malloc((snn->n_neurons + snn->n_input_neurons) * conf->time_steps * n_cps * sizeof(int));


    /* Initialization */
    for(i = 0; i<snn->n_neurons; i++){

        for(j = 0; j<n_cps; j++){

            snn->v_thresh[n_cps * i + j] = tmp_thresh[i];
            snn->v_rest[n_cps * i + j] = tmp_rest[i];
            snn->r_period[n_cps * i + j] = tmp_rperiod[i];
            snn->r_period_remain[n_cps * i + j] = tmp_rperiod_remain[i];
            snn->res[n_cps * i + j] = tmp_res[i];
            snn->post_fired[n_cps * i + j] = tmp_pstfired[i];
            snn->n_neuron_input_synapses[n_cps * i + j] = tmp_n_insyn[i];
            snn->neuron_input_synapses_offset[n_cps * i + j] = tmp_in_off[i];
        }
    }

    for(i = 0; i<snn->n_synapses; i++){

        for(j = 0; j<n_cps; j++){

            snn->w[n_cps * i + j] = tmp_w[i];
            snn->init_w[n_cps * i + j] = tmp_init_w[i];
            snn->dw[n_cps * i + j] = tmp_dw[i];
            snn->pre_trace[n_cps * i + j] = tmp_pre_trace[i];
            snn->post_trace[n_cps * i + j] = tmp_post_trace[i];
            snn->pre_fired[n_cps * i + j] = tmp_prefired[i];
            snn->delay[n_cps * i + j] = tmp_delay[i];
            snn->lr[n_cps * i + j] = tmp_lr[i];
            snn->pre_neuron_index[n_cps * i + j] = tmp_pre_neuron_index[i];
            snn->post_neuron_index[n_cps * i + j] = tmp_post_neuron_index[i];
        }
    }

    for(t = 0; t<conf->time_steps; t++){
        
        for(i = 0; i<snn->n_neurons + snn->n_input_neurons; i++){
            
            for(j = 0; j<n_cps; j++){

                snn->spk_matrix[
                    (snn->n_neurons + snn->n_input_neurons) * n_cps * t + n_cps * i + j
                ] = tmp_spk_mtrx[t * snn->n_neurons + snn->n_input_neurons + j];
            }
        }
    }
}

static inline __mmask16 get_batch_mask(int batch_size)
{
    // batch_size is 8 or 16
    return (batch_size == 8) ? 0x00FF : 0xFFFF;
}

/* main.c */
int main(int argc, char *argv[]) {

    // I think that too much structures are used, probably this should be refactorized
    //spiking_nn_t snn; // base SNN structure and copies to process samples in parallel
    simulation_results_t train_results, test_results; // simulation results
    //network_construction_lists_t lists; // structures to 
    input_data_t train_dataset, test_dataset; // train and test datasets



    // randomize execution
    srand(time(NULL));

    /*
    Load and initialize
    */

    printf(" ============================= \n Loading and initializing data \n ============================= \n");

    // load configuration parameters from input file
    printf(" > Loading configuration file...\n");
    simulation_configuration_t *conf = load_configuration_params_from_toml(argv[1]);
    printf(" > Configuration file loaded!\n\n");
    fflush(stdout);

    // load network initialization arrays
    printf(" > Loading network information in lists...\n");
    network_construction_lists_t *data = load_network_information_in_lists(conf);
    printf(" > Network information loaded!\n\n");
    fflush(stdout);

    // initialize network
    printf(" > Initializing network...\n");
    GPU_SNN_t *snn = initialize_network_cpu(conf, data);
    printf(" > Network initialized!\n");
    fflush(stdout);

    // load dataset
    printf(" > Loading dataset... \n");
    GPU_dataset_t *dataset = load_dataset_from_file_cpu(conf->train_set, conf->train_labels, conf->n_train, conf);
    printf(" > Dataset loaded!\n");
    fflush(stdout);

    // TODO
    GPU_results_t *cpu_results = (GPU_results_t*)calloc(1, sizeof(GPU_results_t));
    cpu_results->n_spks = (int*)calloc(snn->n_neurons, sizeof(int));


    print_network(snn);


    printf("\n ============================= \n ==== Starting simulation ==== \n ============================= \n");

    simulate_sample_CPU(snn, dataset, conf, cpu_results, 0);


    /*
    // load information about the snn from the network file //
    printf(" > Loading network data from file...\n");
    load_network_information(conf->network_file, &snn, &lists, conf); // I don't like that this function loads some data into the SNN structure directly, it's confusing
    printf(" > Network data loaded!\n\n");
    fflush(stdout);

    // initialize the network
    printf(" > Initializing network...\n"); 
    initialize_network(&snn, conf, &lists);
    printf(" > Network initialized!\n\n");
    fflush(stdout);

    // load input spike train from file (different depending on execution type) // ESTO DEBERÍA CAMBIARLO; NO ME TERMINA DE GUSTAR COMO ESTÁ PASANDO DIRECTAMENTE EL PARÁMETRO DE ENTRADA
    printf(" > Loading dataset...\n");
    load_dataset_from_file(&train_dataset, conf->train_set, conf->train_labels, conf->n_train, conf);
    if(conf->test_provided == 1)
        load_dataset_from_file(&test_dataset, conf->test_set, conf->test_labels, conf->n_test, conf);
    printf(" > Dataset loaded!\n");
    fflush(stdout);


    // initialize struct to store results // REVISE THIS
    printf(" > Initializing results struct...\n");
    initialize_results_struct(&train_results, conf, train_dataset.n_samples, snn.n_neurons);
    if(conf->test_provided == 1)
        initialize_results_struct(&test_results, conf, test_dataset.n_samples, snn.n_neurons);

    printf(" > Results strcut initialized!\n");
    fflush(stdout);*/


 

/*
    struct timespec start, end; // to measure simulation complete time
    struct timespec start_step1, end_step1; // to measure simulation complete time
    struct timespec start_step2, end_step2; // to measure simulation complete time
    struct timespec start_step3, end_step3; // to measure simulation complete time
    struct timespec start_step4, end_step4; // to measure simulation complete time
    struct timespec start_step5, end_step5; // to measure simulation complete time

    double elapsed_time, et1 =0.0, et2=0.0, et3=0.0, et4=0.0, et5=0.0;

    #ifdef REORDER

        printf(" > Reordering array of synapses...\n");
        fflush(stdout);
        
        #ifdef CUDA
        // if cuda is used, initialized delays to 0
        if(conf.cuda != 0)
            for(int i = 0; i<snn.n_input_synapses; i++){
                snn.synapses[i].delay = 0;
            }
        #endif

        // reorder array
        reorder_synapse_list(&snn);

        printf(" > Array of synapses reordered!\n");
        fflush(stdout);

    #endif


    printf(" > Mapping structures to new struct type\n");
    GPU_SNN_t *cpu_snn = SNN_CPU2GPU_mapping(&snn, conf);
    GPU_dataset_t *cpu_dataset = dataset_CPU2GPU_mapping(&train_dataset, conf);
    GPU_results_t *cpu_results = (GPU_results_t*)calloc(1, sizeof(GPU_results_t));

    // initialize results struct
    cpu_results->n_spks = (int*)calloc(cpu_snn->n_neurons, sizeof(int));
    //cpu_results->gnt_spks = (char*)calloc(cpu_snn->n_neurons * conf->time_steps, sizeof(char));

    printf(" > Structs mapped!\n");
    fflush(stdout);



    printf(" > Printing network information\n");
    
    for(size_t i = 0; i<cpu_snn->n_synapses; i++){
        printf(" Synapse %ld: In(%ld) Out(%ld)\n", i, cpu_snn->pre_neuron_index[i], cpu_snn->post_neuron_index[i]);
    }



    simulate_sample_CPU(cpu_snn, cpu_dataset, conf, cpu_results, 0);
*/
    // start timing
/*    clock_gettime(CLOCK_MONOTONIC, &start);

    size_t s, sidx;

    float decay = 0.95;//;expf(-1/); // -dt / tau
    int *n_spikes = (int*)calloc(cpu_snn->n_neurons, sizeof(int));
    float *arrI = (float*)calloc(cpu_snn->n_neurons, sizeof(float));
    int *inR = (int*)calloc(cpu_snn->n_neurons, sizeof(int));

    // Restricts for vectorization
    float * restrict rV      = cpu_snn->v;
    float * restrict rRest = cpu_snn->v_rest;
    int * restrict rRef = inR;
    float * restrict rI      = arrI;

    const float alpha = 0.95f;
    const float beta = 0.05f;

    // preload vectors
    // reinitialize post and pre fired


    // Copy network 
    int batch_size=16, n_cps = 16;
    printf(" > Copying network...\n");
    fflush(stdout);
    cpy_snn(cpu_snn, n_cps, &conf);
    printf(" > Network copied!\n");
    fflush(stdout);


    // create batches
    int n_b;
    int e_samples;
    int *batches = create_batches(&n_batches, &extra_samples, dataset, n_samples, batch_size);

    size_t n_batches, extra_samples;
    n_batches = (size_t)n_b;
    extra_samples = (size_t)e_samples;

    // matrix dimensions
    int n_rows = (cpu_snn->n_neurons + cpu_snn->n_input_neurons);
    int n_cols = n_cps;

    

    for(s = 0; s<(size_t)(cpu_dataset->n_samples); s = s + (size_t)(conf.batch_size)){

        printf(" Simulating sample %lu\n", s);
        fflush(stdout);
        for(sidx = 0; sidx < (size_t)(cpu_dataset->n_samples); sidx++){

            printf(" Simulating sample %lu\n", sidx);
            fflush(stdout);

            clock_gettime(CLOCK_MONOTONIC, &start_step1);

            // reinitialize spike matrix
            __mmask16 batch_mask = get_batch_mask();
            __512 zero_row = _mm512_setzero_ps();
            #pragma omp parallel for num_threads(conf.n_process)
            for(size_t t = 0; t<(size_t)(cpu_snn->max_spikes); t++){

                for(size_t i=0; i<(size_t)(cpu_snn->n_neurons + cpu_snn->n_input_neurons); i++){
                    
                    // vectorized
                    for(size_t c = 0; c+batch_size-1<(size_t)n_cps; c+=batch_size){

                        //__m512 spk_row = _mm512_maskz_load_ps(batch_mask, &(cpu_snn->spk_matrix[t * (cpu_snn->n_neurons + cpu_snn->n_input_neurons) * n_cps + i * n_cps + c]));  
                        
                        // store spike row
                        _mm512_mask_store_ps(&(cpu_snn->spk_matrix[t * (cpu_snn->n_neurons + cpu_snn->n_input_neurons) * n_cps + i * n_cps + c]), batch_mask, zero_row);
                    }
                    //cpu_snn->spk_matrix[(size_t)(cpu_snn->max_spikes) * i + j] = 0; // vectorize
                }
            }

            // reinitialize neurons
            #pragma omp parallel for num_threads(conf.n_process)
            for(size_t i = 0; i<(size_t)(cpu_snn->n_neurons); i++){

                // vectorized
                for(size_t c = 0; c+batch_size-1<(size_t)n_cps; c+=batch_size){
                    
                    // load data
                    __m512 v_rest_vec = _mm512_maskz_load_ps(batch_mask, &(cpu_snn->v_rest[i * n_cps + c]));

                    // store
                    _mm512_mask_store_ps(&(cpu_snn->v[i * n_cps + c]), batch_mask, v_rest_vec); // v = v_rest
                    _mm512_mask_store_ps(&(cpu_snn->r_period_remain[i * n_cps + c]), batch_mask, zero_row); // r remain = 0
                }
                //cpu_snn->v[i] = cpu_snn->v_rest[i];
                //cpu_snn->r_period_remain[i] = 0;
            }


            // reinitialize synapses 
            #pragma omp parallel for num_threads(conf.n_process)
            for(size_t i=0; i<(size_t)(cpu_snn->n_synapses); i++){

                // vectorized
                for(size_t c = 0; c+batch_size-1<(size_t)n_cps; c+=batch_size){
                    
                    // load data
                    __m512 init_w_vec = _mm512_maskz_load_ps(batch_mask, &(cpu_snn->v_rest[i * n_cps + c]));

                    // store
                    _mm512_mask_store_ps(&(cpu_snn->w[i * n_cps + c]), batch_mask, init_w_vec); // v = v_rest
                    _mm512_mask_store_ps(&(cpu_snn->next_pre_spike[i * n_cps + c]), batch_mask, zero_row); // r remain = 0
                    _mm512_mask_store_ps(&(cpu_snn->next_post_spike[i * n_cps + c]), batch_mask, zero_row); // r remain = 0

                }

                //cpu_snn->w[i] = cpu_snn->init_w[i];
                //cpu_snn->next_pre_spike[i] = 0;
                //cpu_snn->next_post_spike[i] = 0;
            }

            // load the sample
            #pragma omp parallel for num_threads(conf.n_process)
            for(size_t i=0; i+batch_size-1<(size_t)(cpu_snn->n_input_neurons); i+=batch_size){

                for(size_t c = 0; c<batch_size; c++){
                    
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
            }

            clock_gettime(CLOCK_MONOTONIC, &end_step1);
            et1+=(end_step1.tv_sec - start_step1.tv_sec) + (end_step1.tv_nsec - start_step1.tv_nsec) / 1e9;


            // simulate the sample 
            for(int t=1; t<conf.time_steps; t++){
                
                clock_gettime(CLOCK_MONOTONIC, &start_step2);


                // process neurons
                #pragma omp parallel for num_threads(conf.n_process) 
                for(size_t i=0; i<(size_t)(cpu_snn->n_neurons); i++){

                    //printf(" > > Neuron %lu: \n", i);
                    //fflush(stdout);
                    float I = 0.0;
                    
                    // convert conditional to multiplication
                    // inlcude STDP here looping over all synapses
                    if(cpu_snn->r_period_remain[i] <= 0){

                        size_t synapse_index, in_neuron_index, spk_time_index;
                        int delay, spk_time; 
                        float spk;    

                        size_t base_synapse = (size_t)(cpu_snn->neuron_input_synapses_offset[i]);


                        // correction here: I am probably suming more input synapses than only those of the neuron --> this will happen with low synapse values
                        size_t j;
                        for(j = 0; j + 15 < (size_t)(cpu_snn->n_neuron_input_synapses[i]); j+=16){

                            // load delays and weights
                            __m512i d_vec = _mm512_loadu_si512(&(cpu_snn->delay[base_synapse + j]));                           
                            __m512 w_vec = _mm512_loadu_ps(&(cpu_snn->w[base_synapse + j]));    
                            
                            // compute spike times for each synapse (t-delay)
                            __m512i spk_time_vec = _mm512_sub_epi32(_mm512_set1_epi32(t), d_vec);

                            // mask: 1 if the spk_time is spk_time >= 0
                            __mmask16 valid_mask = _mm512_cmp_epi32_mask(spk_time_vec, _mm512_setzero_si512(), _MM_CMPINT_GE);

                            // clamp spike times to avoid spk_time < 0 values
                            spk_time_vec =_mm512_max_epi32(spk_time_vec, _mm512_setzero_si512());


                            // load input neuron indeces for each synapse
                            __m512i in_vec = _mm512_loadu_si512(&(cpu_snn->pre_neuron_index[base_synapse + j]));  

                            // compute offsets into spk_matrix
                            // matrix [N * T]
                            //__m512i idx = _mm512_add_epi32(_mm512_mullo_epi32(in_vec, _mm512_set1_epi32(cpu_snn->max_spikes)), spk_time_clamped);

                            
                            // matrix [T * N] // This could be fixed to avoid the computation
                            // compute indices to gather spikes [in_index + spk_time * (n_neurons + n_input_neurons)]
                            __m512i idx = 
                                _mm512_add_epi32(
                                    _mm512_mullo_epi32(spk_time_vec, 
                                        _mm512_add_epi32(_mm512_set1_epi32(cpu_snn->n_neurons), 
                                        _mm512_set1_epi32(cpu_snn->n_input_neurons))), 
                                    in_vec);


                            // gather spikes from spike matrix using the indices (int)
                            __m512i spikes_i = _mm512_i32gather_epi32(idx, cpu_snn->spk_matrix, 4);

                            // convert to float
                            __m512 spikes = _mm512_cvtepi32_ps(spikes_i);
                            
                            // multiply by weights to compute I(t)
                            __m512 contrib = _mm512_mask_mul_ps(_mm512_setzero_ps(), valid_mask, w_vec, spikes); // multiply
                            float I_vec = _mm512_reduce_add_ps(contrib); // reduce
                            I += I_vec; 


                            // store wether the presynaptic neuron fired (spikes reached now)
                            _mm512_storeu_si512(&(cpu_snn->pre_fired[base_synapse + j]), spikes_i); // spikes_i stores which neurons fired
                        }
                        



                        for(; j<(size_t)(cpu_snn->n_neuron_input_synapses[i]); j++){

                            synapse_index = base_synapse + j; // no copies
                            in_neuron_index = (size_t)(cpu_snn->pre_neuron_index[synapse_index]); // absurd, equal to 0
                            delay = cpu_snn->delay[synapse_index]; // no copies
                            spk_time = t - delay; // actual position

                            // check if spike is bigger than 0
                            if(spk_time >=0){ // this can be removed? spk_time >= 0 --> mlt = spk_time >= 0??
                                // [N * T]
                                //spk = (float)(cpu_snn->spk_matrix[in_neuron_index * (size_t)(cpu_snn->max_spikes) + (size_t)(spk_time)]);    
                                
                                // T * N
                                // get wether the neuron fired in the matrix and compute I(t)
                                spk = (float)(cpu_snn->spk_matrix[in_neuron_index + (size_t)(spk_time)  * (size_t)(cpu_snn->n_neurons + cpu_snn->n_input_neurons)]);    
                                I += cpu_snn->w[synapse_index] * spk; // spk = 0 / 1
                                cpu_snn->pre_fired[synapse_index] = (int)spk;
                            }
                            else{
                                cpu_snn->pre_fired[synapse_index] = 0;
                            }
                        }
                    }
                    arrI[i] = I;
                    inR[i] = (cpu_snn->r_period_remain[i] <= 0); // store wether the neuron is in refractary period
                }
                clock_gettime(CLOCK_MONOTONIC, &end_step2);
                et2+=(end_step2.tv_sec - start_step2.tv_sec) + (end_step2.tv_nsec - start_step2.tv_nsec) / 1e9;
                
                clock_gettime(CLOCK_MONOTONIC, &start_step3);


                // constants in registers
                __m512 alpha_v = _mm512_set1_ps(alpha);
                __m512 beta_v  = _mm512_set1_ps(beta);

                size_t i;
                for (i = 0; i + 15 < cpu_snn->n_neurons; i += 16) {
                    
                    // load neuron data in vectors
                    __m512 v_vec    = _mm512_loadu_ps(&rV[i]);
                    __m512 vrest_vec= _mm512_loadu_ps(&rRest[i]);
                    __m512 I_vec    = _mm512_loadu_ps(&rI[i]);
                    __m512i ref_vec  = _mm512_loadu_si512(&inR[i]);

                    // mask which will be computed: inR > 0, compute
                    __mmask16 mask = _mm512_cmp_epi32_mask(ref_vec, _mm512_setzero_si512(), _MM_CMPINT_GT); // lower or equal

                    // compute new v-s for all neurons
                    __m512 newv = _mm512_add_ps(
                            _mm512_mul_ps(alpha_v, v_vec),
                            _mm512_add_ps(_mm512_mul_ps(beta_v, vrest_vec), I_vec)
                        );


                    // combine old & new: v = mask ? newv : v_old
                    __m512 v_combined = _mm512_mask_mov_ps(v_vec, mask, newv);

                    // store V[t]
                    _mm512_storeu_ps(&rV[i], v_combined);
                
                }

                // handle remaining neurons
                for (; i < (size_t)(cpu_snn->n_neurons); i++) {

                    // check if the neuron is in refractary period
                    if (inR[i] == 1){
                    
                        // compute v[t]
                        rV[i] = alpha * rV[i] + beta * rRest[i] + rI[i];
                    }
                }

                clock_gettime(CLOCK_MONOTONIC, &end_step3);
                et3+=(end_step3.tv_sec - start_step3.tv_sec) + (end_step3.tv_nsec - start_step3.tv_nsec) / 1e9;

                clock_gettime(CLOCK_MONOTONIC, &start_step4);
      

                // T* N vectorized
                size_t idx = (size_t)(cpu_snn->n_input_neurons + cpu_snn->n_neurons) * (size_t)(t) + (size_t)(cpu_snn->n_input_neurons);
                __m512i ones_vec = _mm512_set1_epi32(1);
                
                i = 0;
                for(; i + 15 < (size_t)(cpu_snn->n_neurons); i+=16){

                    // update refractory period
                    __m512i ref_vec = _mm512_loadu_si512(&(cpu_snn->r_period_remain[i])); // load
                    ref_vec = _mm512_sub_epi32(ref_vec, ones_vec); // r_time = r_time - 1
                    _mm512_storeu_si512(&cpu_snn->r_period_remain[i], ref_vec); // store


                    // load v and v_thresh
                    __m512 v_vec   = _mm512_loadu_ps(&cpu_snn->v[i]);
                    __m512 vth_vec = _mm512_loadu_ps(&cpu_snn->v_thresh[i]);

                    // firing mask (if v >= v_thresh by mask)
                    __mmask16 fire_mask = _mm512_cmp_ps_mask(v_vec, vth_vec, _CMP_GE_OS);

                    // create the array of spikes (1 if v >= v_thresh, 0 else)
                    __m512i spike_vec =
                        _mm512_mask_mov_epi32(_mm512_setzero_si512(),
                                            fire_mask,
                                            _mm512_set1_epi32(1));

                    // store which neurons fired in the spk matrix [T * N], neurons are in sequential positions
                    _mm512_storeu_si512(&cpu_snn->spk_matrix[idx + i], spike_vec);
                    
                    // reset refractory for neurons that fired
                    __m512i ref_reset_vec = _mm512_loadu_si512(&cpu_snn->r_period[i]); // load r_period
                    ref_vec = _mm512_mask_mov_epi32(ref_vec, fire_mask, ref_reset_vec); // r_period_remain = r_period (if fired)
                    _mm512_storeu_si512(&cpu_snn->r_period_remain[i], ref_vec); // store

                    // reset 
                    __m512 vrest_vec = _mm512_loadu_ps(&cpu_snn->v_rest[i]); // load reset
                    v_vec = _mm512_mask_mov_ps(v_vec, fire_mask, vrest_vec); // v = reset (if fired)
                    _mm512_storeu_ps(&cpu_snn->v[i], v_vec); // store

                    // increment spike count
                    __m512i nspk_vec = _mm512_loadu_si512(&n_spikes[i]); // load n_spikes
                    nspk_vec = _mm512_mask_add_epi32(nspk_vec, fire_mask, nspk_vec, ones_vec); // n_spikes += 1 (if fired)
                    _mm512_storeu_si512(&n_spikes[i], nspk_vec); // store

                    // store which neurons fired
                    _mm512_storeu_si512(&(cpu_snn->post_fired[i]), spike_vec);

                }
                //#pragma omp parallel for num_threads(conf.n_process)
                for(; i<(size_t)(cpu_snn->n_neurons); i++){
                    
                    cpu_snn->r_period_remain[i] --;
                    if(cpu_snn->v[i] >= cpu_snn->v_thresh[i]){
                        
                        // T * N
                        size_t idx_neuron = idx + i; // row start at

                        cpu_snn->spk_matrix[idx_neuron] = 1;

                        // reinit neuron values
                        cpu_snn->r_period_remain[i] = cpu_snn->r_period[i];
                        cpu_snn->v[i] = cpu_snn->v_rest[i]; // reinit v_rest

                        n_spikes[i] += 1;

                    }

                    cpu_snn->post_fired[i] = cpu_snn->v[i] >= cpu_snn->v_thresh[i];
                }

                clock_gettime(CLOCK_MONOTONIC, &end_step4);
                et4+=(end_step4.tv_sec - start_step4.tv_sec) + (end_step4.tv_nsec - start_step4.tv_nsec) / 1e9;

                // learn
                clock_gettime(CLOCK_MONOTONIC, &start_step5);

                if(conf.learn == 1){
                    
                    i = 0;
                    /*
                    
                    // load traces
                    for(; i+15<(size_t)(cpu_snn->n_synapses); i+=16){

                        // load traces
                        __m512 preT_vec = _mm512_loadu_ps(&(cpu_snn->pre_trace[i]));
                        __m512 postT_vec = _mm512_loadu_ps(&(cpu_snn->post_trace[i]));

                        // update the postsynaptic trace

                        // [IF (post_fired) THEN update_post_trace() ENDIF]

                        // load postsynaptic neurons indexes for synapses
                        __m512i out_vec = _mm512_loadu_si512(&(cpu_snn->post_neuron_index[i])); // out neurons

                        // index = index - n_input_neurons
                        out_vec = _mm512_sub_epi32(out_vec, _mm512_set1_epi32(cpu_snn->n_input_neurons));
                        
                        // get if the neuron [index] fired using gather
                        __m512i post_fired_vec = _mm512_i32gather_epi32(out_vec, cpu_snn->post_fired, 4); // get which output neurons fired

                        // mask depending on which neurons fired [post_fired[i] == 1]
                        __mmask16 post_valid_mask = _mm512_cmp_epi32_mask(post_fired_vec, _mm512_setzero_epi32(), _MM_CMPINT_GT);

                        // update postsynaptic trace if neuron fired
                        postT_vec = _mm512_mask_add_ps(postT_vec, post_valid_mask, postT_vec, _mm512_set1_ps(1.0));
                        

                        // pre_fired = pre_fired - post_fired, if pre_fired == 1, compute LTD, if post_fired == 1, compute LTP
                        // load pre_fired and mask
                        __m512i pre_fired_vec = _mm512_loadu_si512(&(cpu_snn->pre_fired[i]));
                        pre_fired_vec = _mm512_sub_epi32(pre_fired_vec, post_fired_vec);
                        __mmask16 pre_valid_mask = _mm512_cmp_epi32_mask(pre_fired_vec, _mm512_setzero_epi32(), _MM_CMPINT_GT);


                        // update weights 

                        // load dw and w
                        __m512 dw_vec = _mm512_loadu_ps(&(cpu_snn->dw[i]));
                        __m512 w_vec = _mm512_loadu_ps(&(cpu_snn->w[i]));


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
                        _mm512_storeu_ps(&(cpu_snn->w[i]), w_vec);
                        _mm512_storeu_ps(&(cpu_snn->dw[i]), dw_vec);


                        // compute the traces decays (for all, not only masked)
                        preT_vec =_mm512_mul_ps(preT_vec, _mm512_set1_ps(decay));
                        postT_vec = _mm512_mul_ps(postT_vec, _mm512_set1_ps(decay));

                        // store traces
                        _mm512_storeu_ps(&(cpu_snn->pre_trace[i]), preT_vec);
                        _mm512_storeu_ps(&(cpu_snn->post_trace[i]), postT_vec);
                    }*/

                    // manual openmp + vectorization
                    
                    //size_t n_tasks = 0;

                    /*
                    n_tasks = cpu_snn->n_synapses / 16;

                    #pragma omp parallel for num_threads(conf.n_process) schedule(static) private(i)
                    for(size_t tsk = 0; tsk<(size_t)n_tasks; tsk++){

                        size_t i = tsk * 16;


                        // load traces
                        __m512 preT_vec = _mm512_loadu_ps(&(cpu_snn->pre_trace[i]));
                        __m512 postT_vec = _mm512_loadu_ps(&(cpu_snn->post_trace[i]));

                        // update the postsynaptic trace

                        // [IF (post_fired) THEN update_post_trace() ENDIF]

                        // load postsynaptic neurons indexes for synapses
                        __m512i out_vec = _mm512_loadu_si512(&(cpu_snn->post_neuron_index[i])); // out neurons

                        // index = index - n_input_neurons
                        out_vec = _mm512_sub_epi32(out_vec, _mm512_set1_epi32(cpu_snn->n_input_neurons));
                        
                        // get if the neuron [index] fired using gather
                        __m512i post_fired_vec = _mm512_i32gather_epi32(out_vec, cpu_snn->post_fired, 4); // get which output neurons fired

                        // mask depending on which neurons fired [post_fired[i] == 1]
                        __mmask16 post_valid_mask = _mm512_cmp_epi32_mask(post_fired_vec, _mm512_setzero_epi32(), _MM_CMPINT_GT);

                        // update postsynaptic trace if neuron fired
                        postT_vec = _mm512_mask_add_ps(postT_vec, post_valid_mask, postT_vec, _mm512_set1_ps(1.0));
                        

                        // pre_fired = pre_fired - post_fired, if pre_fired == 1, compute LTD, if post_fired == 1, compute LTP
                        // load pre_fired and mask
                        __m512i pre_fired_vec = _mm512_loadu_si512(&(cpu_snn->pre_fired[i]));
                        pre_fired_vec = _mm512_sub_epi32(pre_fired_vec, post_fired_vec);
                        __mmask16 pre_valid_mask = _mm512_cmp_epi32_mask(pre_fired_vec, _mm512_setzero_epi32(), _MM_CMPINT_GT);


                        // update weights 

                        // load dw and w
                        __m512 dw_vec = _mm512_loadu_ps(&(cpu_snn->dw[i]));
                        __m512 w_vec = _mm512_loadu_ps(&(cpu_snn->w[i]));


                        //
                        //    USING MASKS
                        //    if post_fired == 1:
                        //        LTP
                        //    else if pre_fired == 1:
                        //        LTD
            
                        // update w and dw (depending on the trace) [w += lr * pre_trace] --> LTP
                        const __m512 A_vec = _mm512_set1_ps(0.01);
                        dw_vec = _mm512_mask_add_ps(dw_vec, post_valid_mask, dw_vec, _mm512_mul_ps(preT_vec, A_vec));
                        w_vec = _mm512_mask_add_ps(w_vec, post_valid_mask, w_vec, _mm512_mul_ps(preT_vec, A_vec));

                        // update w and dw (depending on the trace) [w -= lr * post_trace] --> LTD
                        dw_vec = _mm512_mask_sub_ps(dw_vec, pre_valid_mask, dw_vec, _mm512_mul_ps(postT_vec, A_vec));
                        w_vec = _mm512_mask_sub_ps(w_vec, pre_valid_mask, w_vec, _mm512_mul_ps(postT_vec, A_vec));
                       
                        __m512 delta = _mm512_setzero_ps();    
                        delta = _mm512_mask_add_ps(delta, post_valid_mask, delta, preT_vec);
                        delta = _mm512_mask_sub_ps(delta, pre_valid_mask,  delta, postT_vec);

                        w_vec  = _mm512_fmadd_ps(delta, A_vec, w_vec);
                        dw_vec = _mm512_fmadd_ps(delta, A_vec, dw_vec);


                        // store w and dw
                        _mm512_storeu_ps(&(cpu_snn->w[i]), w_vec);
                        _mm512_storeu_ps(&(cpu_snn->dw[i]), dw_vec);


                        // compute the traces decays (for all, not only masked)
                        const __m512 decay_vec = _mm512_set1_ps(decay);
                        preT_vec =_mm512_mul_ps(preT_vec, decay_vec);
                        postT_vec = _mm512_mul_ps(postT_vec, decay_vec);

                        //// store traces
                        _mm512_storeu_ps(&(cpu_snn->pre_trace[i]), preT_vec);
                        _mm512_storeu_ps(&(cpu_snn->post_trace[i]), postT_vec);

                    }



                    i = n_tasks * 16;                    

                    // manual openmp + vectorization
                    #pragma omp parallel for num_threads(conf.n_process)
                    for(size_t j=i; j<(size_t)cpu_snn->n_synapses; j++){

                        // loop postsynaptic neurons and update trace

                        int post_fired = cpu_snn->post_fired[cpu_snn->post_neuron_index[j]] == 1;
                        cpu_snn->pre_trace[j] += post_fired * 1.0;

                        // traces updated, compute 

                        // compute STDP
                        float dw = post_fired * 0.01 * cpu_snn->pre_trace[j];
                        cpu_snn->dw[j] += dw;
                        cpu_snn->w[j] += dw;




                        int pre_fired = cpu_snn->pre_fired[j] == 1;
                        cpu_snn->post_trace[j] += pre_fired * 1.0;

                        // traces updated, compute 

                        // compute STDP
                        dw = pre_fired * 0.01 * cpu_snn->post_trace[j];
                        cpu_snn->dw[j] -= dw;
                        cpu_snn->w[j] -= dw;

                        // decay
                        cpu_snn->pre_trace[j] *= decay;
                        cpu_snn->post_trace[j] *= decay;

                        /*if(cpu_snn->post_fired[cpu_snn->post_neuron_index[j]] == 1){
                            
                            cpu_snn->pre_trace[j] += 1.0;

                            // traces updated, compute 

                            // compute STDP
                            cpu_snn->dw[j] += 0.01 * cpu_snn->pre_trace[j];
                            cpu_snn->w[j] += 0.01 * cpu_snn->pre_trace[j];
                        }
                        else if(cpu_snn->pre_fired[j] == 1){
                            
                            cpu_snn->post_trace[j] += 1.0;

                            // traces updated, compute 

                            // compute STDP
                            cpu_snn->dw[j] -= 0.01 * cpu_snn->post_trace[j];
                            cpu_snn->w[j] -= 0.01 * cpu_snn->post_trace[j];
                        }

                        // update traces using the decay
                        cpu_snn->pre_trace[j] *= decay;
                        cpu_snn->post_trace[j] *= decay;
                    }

                    
                }

                clock_gettime(CLOCK_MONOTONIC, &end_step5);
                et5+=(end_step5.tv_sec - start_step5.tv_sec) + (end_step5.tv_nsec - start_step5.tv_nsec) / 1e9;
            }
        }


        // guide the learning after the batch finished
    }   

    clock_gettime(CLOCK_MONOTONIC, &end);

    elapsed_time = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    printf(" FInished in %lf! (s1 = %lf s2 = %lf s3 = %lf s4 = %lf, s5 = %lf)\n", 
            elapsed_time, et1, et2, et3, et4, et5);
    printf(" Generated number of spikes per neuron: ");
    for(int i = 0; i<cpu_snn->n_neurons; i++){
        printf("%d ", n_spikes[i]);
    }
    printf("\n");
    fflush(stdout);
*/
/* second str
    size_t s, sidx;
    int *n_spikes = (int*)calloc(cpu_snn->n_neurons, sizeof(int));

    for(s = 0; s<(size_t)(cpu_dataset->n_samples); s = s + (size_t)(conf.batch_size)){

        printf(" Simulating sample %lu\n", s);
        fflush(stdout);
        for(sidx = 0; sidx < (size_t)(cpu_dataset->n_samples); sidx++){

            printf(" Simulating sample %lu\n", sidx);
            fflush(stdout);
            // reinitialize spike matrix
            #pragma omp parallel for num_threads(conf.n_process)
            for(size_t i=0; i<(size_t)(cpu_snn->n_neurons + cpu_snn->n_input_neurons); i++){
                for(size_t j = 0; j<(size_t)(cpu_snn->max_spikes); j++){
                    
                    cpu_snn->spk_matrix[(size_t)(cpu_snn->max_spikes) * i + j] = 0;
                }
            }

            // reinitialize neurons
            #pragma omp parallel for num_threads(conf.n_process)
            for(size_t i = 0; i<(size_t)(cpu_snn->n_neurons); i++){

                cpu_snn->v[i] = cpu_snn->v_rest[i];
                cpu_snn->r_period_remain[i] = 0;
            }


            // reinitialize synapses 
            #pragma omp parallel for num_threads(conf.n_process)
            for(size_t i=0; i<(size_t)(cpu_snn->n_synapses); i++){

                cpu_snn->w[i] = cpu_snn->init_w[i];
                cpu_snn->next_pre_spike[i] = 0;
                cpu_snn->next_post_spike[i] = 0;
            }

            // load the sample
            #pragma omp parallel for num_threads(conf.n_process)
            for(size_t i=0; i<(size_t)(cpu_snn->n_input_neurons); i++){

                if(sidx < cpu_dataset->n_samples){

                    size_t fidx = i;

                    size_t start_sample = cpu_dataset->sample_offset[sidx];
                    size_t start_feature = cpu_dataset->feature_offset[sidx * (size_t)(cpu_dataset->n_features) + fidx];
                    size_t n_spikes_per_feature = (size_t)(cpu_dataset->n_spikes_per_feature[sidx * (size_t)(cpu_dataset->n_features) + fidx]);
                
                    for(size_t j = 0; j<n_spikes_per_feature; j++){

                        cpu_snn->spk_matrix[(size_t)(cpu_snn->max_spikes) * i + (size_t)(cpu_dataset->spikes[start_feature + j])] = 1;//cpu_dataset->spikes[start_feature + j];
                    }
                }
            }


            // simulate the sample 
            for(int t=0; t<conf.time_steps; t++){
                
                //printf(" > Time step %d\n", t);

                // process neurons
                #pragma omp parallel for num_threads(conf.n_process)
                for(size_t i=0; i<(size_t)(cpu_snn->n_neurons); i++){

                    //printf(" > > Neuron %lu: \n", i);
                    //fflush(stdout);
                    float I = 0.0;

                    if(cpu_snn->r_period_remain[i] <= 0){

                        size_t synapse_index, in_neuron_index, spk_time_index;
                        int delay, spk_time; 
                        float spk;    

                        for(size_t j=0; j<(size_t)(cpu_snn->n_neuron_input_synapses[i]); j++){

                            synapse_index = (size_t)(cpu_snn->neuron_input_synapses_offset[i]) + j; // no copies
                            in_neuron_index = (size_t)(cpu_snn->pre_neuron_index[synapse_index]); // no copies
                            delay = cpu_snn->delay[synapse_index]; // no copies
                            spk_time = t - delay; // actual position

                            if(spk_time >=0){
                            
                                spk = (float)(cpu_snn->spk_matrix[in_neuron_index * (size_t)(cpu_snn->max_spikes) + (size_t)(spk_time)]);    
                                //printf(" > > > spk = %f, w = %f\n", spk, cpu_snn->w[synapse_index]);
                                I += cpu_snn->w[synapse_index] * spk; // spk = 0 / 1
                            }
                        }
                    }
                    //printf("v = %f , I = %f ",cpu_snn->v[i], I);
                    cpu_snn->v[i] = (1-0.05) * cpu_snn->v[i] + cpu_snn->v_rest[i] * 0.05 + I;
                    //printf("v = %f\n", cpu_snn->v[i]);
                }

                
                // process neuron output 
                #pragma omp parallel for num_threads(conf.n_process)
                for(size_t i = 0; i<(size_t)(cpu_snn->n_neurons); i++){
                    
                    cpu_snn->r_period_remain[i] --;
                    if(cpu_snn->v[i] >= cpu_snn->v_thresh[i]){

                        //printf(" SPike generated\n");
                        //fflush(stdout);
                        size_t row, col, matrix_index, global_index;

                        //matrix_index = (size_t)(snn->n_neurons + snn->n_input_neurons) * (size_t)(cpu_snn->max_spikes); // cpy starts at
                        row = ((size_t)(cpu_snn->n_input_neurons) + i) * (size_t)(cpu_snn->max_spikes); // row start at
                        col = (size_t)(t) ;//% (size_t)(cpu_snn->max_spikes); // column to locate

                        //global_index = matrix_index + row + col;
                        global_index = row + col;

                        // add spike to matrix
                        cpu_snn->spk_matrix[global_index] = 1;

                        // reinit neuron values
                        cpu_snn->r_period_remain[i] = cpu_snn->r_period[i];
                        
                        cpu_snn->v[i] = cpu_snn->v_rest[i]; // reinit v_rest

                        n_spikes[i] += 1;
                    }
                }
            }
        }
    }   

    clock_gettime(CLOCK_MONOTONIC, &end);

    double elapsed_time = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    printf(" FInished in %lf!\n", elapsed_time);
    printf(" Generated number of spikes per neuron: ");
    for(int i = 0; i<cpu_snn->n_neurons; i++){
        printf("%d ", n_spikes[i]);
    }
    printf("\n");
    fflush(stdout);
*/

/*

    // train the network
    if(conf.cuda == 0){

        // biological simulation
        if(conf.simulation_obj == 0){

            simulate_samples(&snn, &conf, &train_results, &train_dataset);
        }
        // ML
        else{
            if(conf.train_provided == 1)
                train_network(&snn, &conf, &train_results, &train_dataset);

            // test the network
            if(conf.test_provided == 1)
                test_network(&snn, &conf, &test_results, &test_dataset);

            // another function to integrate train and test?
        }
    }
    else{
        #ifdef CUDA
        
        printf(" > Running for Cuda\n");
        fflush(stdout);

        double elpt;
        struct timespec start, end; // to measure simulation complete time
        struct timespec start_cpy, end_cpy;
        struct timespec start_simulation, end_simulation;
        struct timespec start_mapping, end_mapping;

        #ifdef REORDER

            printf(" > Reordering array of synapses...\n");
            fflush(stdout);
            
            #ifdef CUDA
            // if cuda is used, initialized delays to 0
            if(conf.cuda != 0)
                for(int i = 0; i<snn.n_input_synapses; i++){
                    snn.synapses[i].delay = 0;
                }
            #endif

            // reorder array
            reorder_synapse_list(&snn);

            printf(" > Array of synapses reordered!\n");
            fflush(stdout);

        #endif


        // start timing
        clock_gettime(CLOCK_MONOTONIC, &start);

        // get gpu properties
        cuda_info_t *cuda_info = getProperties();
        // compute_cuda_sim_conf(cuda_info);


        clock_gettime(CLOCK_MONOTONIC, &start_mapping);


        // map data structures to GPU structures in CPU
        GPU_SNN_t *gpu_snn = SNN_CPU2GPU_mapping(&snn, &conf);
        //printf(" > SNN mapped in CPU\n");
        //fflush(stdout);

        GPU_dataset_t *gpu_dataset = dataset_CPU2GPU_mapping(&train_dataset, &conf);
        //printf(" > Dataset mapped in CPU\n");
        
        clock_gettime(CLOCK_MONOTONIC, &end_mapping);
        elpt = (end_mapping.tv_sec - start_mapping.tv_sec) + (end_mapping.tv_nsec - start_mapping.tv_nsec) / 1e9;
        printf(" > > Elapsed time mapping data: %lf\n", elpt);
        fflush(stdout);

        // configure cuda simulation
        configure_cuda_simulation(cuda_info, gpu_snn, gpu_dataset, &conf);
        //printf(" > Cuda simulation configured CPU\n");
        //fflush(stdout);

        // copy data structures to GPU memory
        clock_gettime(CLOCK_MONOTONIC, &start_cpy); 

        GPU_SNN_t **d_gpu_snn = cpy_SNN2GPU(gpu_snn, cuda_info); // structure in GPU
        //printf(" > SNN copied to GPU\n");
        //fflush(stdout);

        GPU_dataset_t **d_gpu_dataset = cpy_dataset2GPU(gpu_dataset, cuda_info);
        //printf(" > Dataset copied to GPU\n");

        clock_gettime(CLOCK_MONOTONIC, &end_cpy);
        elpt = (end_cpy.tv_sec - start_cpy.tv_sec) + (end_cpy.tv_nsec - start_cpy.tv_nsec) / 1e9;
        printf(" > > Elapsed time copying data: %lf\n", elpt);
        fflush(stdout);

        // free memory of GPU structs allocated in CPU memory
        //printf(" > Deallocating memory\n");
        //fflush(stdout);
        free_gpu_snn_in_CPU(gpu_snn);
        free_gpu_dataset_in_CPU(gpu_dataset);
        //printf(" > CPU intermediate structs deallocated!\n");
        //fflush(stdout);

        // perform the simulation
        clock_gettime(CLOCK_MONOTONIC, &start_simulation);
        GPU_results_t *cpu_results = simulate_in_GPU(d_gpu_snn, d_gpu_dataset, &conf, cuda_info, &snn, &train_dataset); // refactorize
        clock_gettime(CLOCK_MONOTONIC, &end_simulation);

        elpt = (end_simulation.tv_sec - start_simulation.tv_sec) + (end_simulation.tv_nsec - start_simulation.tv_nsec) / 1e9;
        printf(" > > Elapsed time simulating in GPU: %lf\n", elpt);
        fflush(stdout);

        // store results // TEMPORAL
        


        // end timing
        clock_gettime(CLOCK_MONOTONIC, &end);
        elpt = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;
    
        printf(" > Elapsed total time: %lf\n", elpt);

        #endif
    }

    */
    // free memory






    return 0;
}
