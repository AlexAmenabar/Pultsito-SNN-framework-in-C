#include "snn_library.h"
#include "training_rules/stdp.h"

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

// TODO: This must be revised, since should be introduced as input or in the network, I guess?
#define TAU_PLUS 5.0 // ms
#define TAU_MINUS 5.0 // ms
#define A_PLUS 1.0 // ms
#define A_MINUS 1.0 // ms
#define A 0.25 // modulation magnitude for STDP
#define P_WINDOW 50.0
#define N_WINDOW -75.0
#define EPSILON 0.15

int mod(int a, int m) {
    return ( (a % m) + m ) % m;
}

/*
Functions to do the simple computations of STDP based on the time difference between two spikes
*/

double addSTDP_comp(synapse_t *synapse, double tdiff){

    double dw = 0.0;

    if(tdiff > 0 && tdiff < P_WINDOW){
        
        //synapse->w += A_PLUS * exp(-time_diff / TAU_PLUS);
        dw += A_PLUS * exp(-tdiff / TAU_PLUS);

        //if(initial_weight > 0 && synapse->w < 0)
        //    synapse->w = 0.0001;
    }
    else if(tdiff < 0 && tdiff > N_WINDOW){ // time window to stdp be considered
        
        //synapse->w -= A_MINUS * exp(time_diff / TAU_MINUS);
        dw -= A_MINUS * exp(tdiff / TAU_MINUS);

        //if(initial_weight < 0 && synapse->w > 0)
        //    synapse->w = -0.0001;
    }

    return dw;
}

// TODO: not correctly implemented: max-min bounds?
double mltSTDP_comp(synapse_t *synapse, double tdiff){

    double dw = 0.0;

    if(tdiff > 0 && tdiff < P_WINDOW){
        
        //synapse->w += A_PLUS * exp(-time_diff / TAU_PLUS);
        //dw += A_PLUS * synapse->w * (1 - synapse->w) * exp(-tdiff / TAU_PLUS);
        dw += A_PLUS * exp(-tdiff / TAU_PLUS);


        //if(initial_weight > 0 && synapse->w < 0)
        //    synapse->w = 0.0001;
    }
    else if(tdiff < 0 && tdiff > N_WINDOW){ // time window to stdp be considered
        
        //synapse->w -= A_MINUS * exp(time_diff / TAU_MINUS);
        //dw -= A_MINUS * synapse->w * (1 - synapse->w) * exp(tdiff / TAU_MINUS);
        dw -= A_MINUS * exp(tdiff / TAU_MINUS);


        //if(initial_weight < 0 && synapse->w > 0)
        //    synapse->w = -0.0001;
    }

    return dw;
}

double antiSTDP_comp(synapse_t *synapse, int tdiff){

    double dw = 0.0;

    if(tdiff > 0 && tdiff < P_WINDOW){
        //synapse->w -= A_PLUS * exp(-time_diff / TAU_PLUS);
        dw -= A_PLUS * exp(-tdiff / TAU_PLUS);

        //if(initial_weight > 0 && synapse->w < 0)
        //    synapse->w = 0.0001;
    }
    else if(tdiff < 0 && tdiff > N_WINDOW){ // time window to stdp be considered
        //synapse->w += A_MINUS * exp(time_diff / TAU_MINUS);
        dw += A_MINUS * exp(tdiff / TAU_MINUS);

        //if(initial_weight < 0 && synapse->w > 0)
        //    synapse->w = -0.0001;
    }

    return dw;
}

double tbSTDP_comp(synapse_t *synapse_t, int tdiff){

    double dw = 0.0;
}


/*
STDP general functions
*/
int cond_stdp(synapse_t *synapse, int t){


    // conditions:
    // - presynaptic and postsynaptic neurons exists (in other words, synapse is neither an input or output synapse)
    // - both, postsynaptic and presynaptic neurons have done at least one spike to compute STDP
    // - last spikes timestamps are different
    // - at least one of both spikes is on t timestamp, otherwise it has been already computed

    lif_neuron_t *pre, *post;
    pre = synapse->pre_synaptic_lif_neuron;
    post = synapse->post_synaptic_lif_neuron;

    //return pre != NULL && post != NULL && // is not input or output synapse
    //    post->t_last_spikes[0] != -1 && pre->t_last_spikes[0] != -1 && // both computed
    //    post->t_last_spikes[(post->next_last_spike - 1) % post->n_last_spikes] != pre->t_last_spikes[(pre->next_last_spike - 1) % pre->n_last_spikes] && // if both are equal STDP results is 0
    //    (post->t_last_spikes[(post->next_last_spike - 1) % post->n_last_spikes] == t || pre->t_last_spikes[(pre->next_last_spike - 1) % pre->n_last_spikes] == t); // one of both has to be actual t


    // check that both neurons exist
    if(pre == NULL || post == NULL)
        return 0;

    
    // check that at least one of both is t
    if(pre->spike_times_arr[synapse->next_pre_spike] + synapse->delay != t && post->spike_times_arr[synapse->next_post_spike] != t)
        return 0;


    return 1;
}


void stdp(synapse_t *synapse, int t, int n, double (*stdp_func)(synapse_t *synapse, double tdiff_double)){

    double dw = 0;
    
    // check conditions to compute STDP
    if(cond_stdp(synapse, t) == 1){

        int delay;
        double tdiff_double;
        int tpost, tpre, prev_tpost, prev_tpre, tmp = 0;
        lif_neuron_t *post, *pre;
        int max_spikes, last_spike;

        // store pre and postsynaptic neurons
        post = synapse->post_synaptic_lif_neuron;
        pre = synapse->pre_synaptic_lif_neuron;
        delay = synapse->delay;

        // compute tdiffs
        tpost = post->spike_times_arr[(synapse->next_post_spike) % post->max_spikes];
        tpre = pre->spike_times_arr[(synapse->next_pre_spike) % pre->max_spikes] + delay;

        // get previous spikes for both
        //prev_tpost = post->spike_times_arr[mod(synapse->next_post_spike - 1, post->max_spikes)]; // can be useful in the future?
        //prev_tpre = pre->spike_times_arr[mod(synapse->next_pre_spike - 1, pre->max_spikes)] + delay;

        // if tdiff == 0 --> random: incrase or decrease

        //int lower_bound;
        int tmp_tdiff = tpost - tpre;

        // LTP
        if(tpost == t){

            //lower_bound = synapse->next_pre_spike;
            max_spikes = pre->max_spikes;
            last_spike = pre->last_spike;
            tmp_tdiff = 100000;

            // conditions to keep looping:
            // - tpre - delay > -1, there is an spike
            // - tpost >= tpre, there are spikes to be processed in STDP
            // - lower_bound is smaller or equal to tpre
            while((tpre - delay) > -1 && tpost >= tpre && tpost - tpre < tmp_tdiff){

                // compute stdp
                tmp_tdiff = tpost - tpre;
                tdiff_double = (double)(tpost - tpre) + EPSILON;
                dw += stdp_func(synapse, tdiff_double);
                //printf(" >> tpost = %d, tpre = %d; tdiff %lf; dw = %lf\n", tpost, tpre, tdiff_double, dw);

                // update tmp
                tmp++;
                tpre = pre->spike_times_arr[(synapse->next_pre_spike + tmp) % pre->max_spikes] + delay;
            }

            synapse->next_pre_spike = (synapse->next_pre_spike + tmp) % pre->max_spikes;

        }
        // LTD: if tpost != t, then tpre == t due to the cond_stdp(...) function
        else{
            
            //lower_bound = tpost - 1;
            max_spikes = post->max_spikes;
            last_spike = post->last_spike;
            tmp_tdiff = -1000000;

            // conditions to keep looping:
            // - tpost > -1, there is an spike
            // - tpost < tpre, there are spikes to be processed in STDP
            // - lower_bound is smaller or equal to tpre
            while(tpost > -1 && tpost < tpre && tpost - tpre > tmp_tdiff){

                // compute stdp
                tmp_tdiff = tpost - tpre;
                tdiff_double = (double)(tpost - tpre) - EPSILON; // I can ignore this epsilon?
                dw += stdp_func(synapse, tdiff_double);
                //printf(" >> tpre = %d, tpost = %d; tdiff %lf; dw = %lf\n", tpre, tpost, tdiff_double, dw);

                // update tmp
                tmp++;
                tpost = post->spike_times_arr[(synapse->next_post_spike + tmp) % max_spikes];
            }

            synapse->next_post_spike = (synapse->next_post_spike + tmp) % max_spikes;

        }
    }

    // update weight and store weight change
    synapse->dw += dw;
    synapse->w += dw;

    //printf(" > Sample ?, dw[?] = %f, final dw[?] = %f\n", dw, synapse->dw);

}


void addSTDP(synapse_t *synapse, int t, int n){

    stdp(synapse, t, n, &addSTDP_comp);
}


void mltSTDP(synapse_t *synapse, int t, int n){

    stdp(synapse, t, n, &mltSTDP_comp);
}


void antiSTDP(synapse_t *synapse, int t, int n){

    stdp(synapse, t, n, &addSTDP_comp);
}


