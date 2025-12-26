for(s = 0; s<(size_t)(cpu_dataset->n_samples); s = s + (size_t)(conf.batch_size)){

        printf(" Simulating sample %lu\n", s);
        fflush(stdout);
        for(sidx = 0; sidx < (size_t)(cpu_dataset->n_samples); sidx++){

            printf(" Simulating sample %lu\n", sidx);
            fflush(stdout);

            clock_gettime(CLOCK_MONOTONIC, &start_step1);

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
            for(int t=1; t<conf.time_steps; t++){
                
                clock_gettime(CLOCK_MONOTONIC, &start_step2);


                /* INPUT STEP */

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


                        /* CONSTANT DELAYS */
                        /*__m512i spk_time_vec = _mm512_sub_epi32(_mm512_set1_epi32(t), _mm512_set1_epi32(1));
                        __m512i n_vec = _mm512_add_epi32(_mm512_set1_epi32(cpu_snn->n_neurons), _mm512_set1_epi32(cpu_snn->n_input_neurons));
                        __m512i row_vec = _mm512_mullo_epi32(spk_time_vec, n_vec);
                        int *row = &cpu_snn->spk_matrix[spk_time * (cpu_snn->n_neurons + cpu_snn->n_input_neurons)];

                        for(j = 0; j + 15 < (size_t)(cpu_snn->n_neuron_input_synapses[i]); j+=16){

                            // load delays and weights
                            __m512 w_vec = _mm512_loadu_ps(&(cpu_snn->w[base_synapse + j]));    
                            
                            // compute spike times for each synapse (t-delay)
                            //__m512i spk_time_vec = _mm512_sub_epi32(_mm512_set1_epi32(t), d_vec);

                            // load input neuron indeces for each synapse
                            __m512i in_vec = _mm512_loadu_si512(&(cpu_snn->pre_neuron_index[base_synapse + j]));  

                            // matrix [T * N] // This could be fixed to avoid the computation
                            // compute indices to gather spikes [in_index + spk_time * (n_neurons + n_input_neurons)]

                            // gather spikes from spike matrix using the indices (int)
                            __m512i spikes_i = _mm512_i32gather_epi32(in_vec, row, 4);
                            __m512i d_vec = _mm512_loadu_si512(&(cpu_snn->delay[base_synapse + j]));                           

                            // convert to float
                            __m512 spikes = _mm512_cvtepi32_ps(spikes_i);
                            
                            // multiply by weights to compute I(t)
                            __m512 contrib = _mm512_mul_ps(w_vec, spikes); // multiply
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
                        }*/


                    }
                    arrI[i] = I;
                    inR[i] = (cpu_snn->r_period_remain[i] <= 0); // store wether the neuron is in refractary period
                }
                clock_gettime(CLOCK_MONOTONIC, &end_step2);
                et2+=(end_step2.tv_sec - start_step2.tv_sec) + (end_step2.tv_nsec - start_step2.tv_nsec) / 1e9;
                
                clock_gettime(CLOCK_MONOTONIC, &start_step3);


                /* V[t] COMPUTATION */

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
      


                /* OUTPUT STEP */

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
                    
                    size_t n_tasks = 0;

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
*/


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
                        cpu_snn->post_trace[j] *= decay;*/
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
