#include "snn_library.h"
#include "training_rules/stdp.h"

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

// TODO: This must be revised, since should be introduced as input or in the network, I guess?
#define TAU_PLUS 5 // ms
#define TAU_MINUS 5 // ms
#define A_PLUS 0.01 // ms
#define A_MINUS 0.01 // ms
#define A 0.25 // modulation magnitude for STDP
#define P_WINDOW 50
#define N_WINDOW -75


/*
Functions to do the simple computations of STDP based on the time difference between two spikes
*/

double addSTDP_comp(synapse_t *synapse, int tdiff){

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

double mltSTDP_comp(synapse_t *synapse, int tdiff){

    double dw = 0.0;

    if(tdiff > 0 && tdiff < P_WINDOW){
        
        //synapse->w += A_PLUS * exp(-time_diff / TAU_PLUS);
        dw += A_PLUS * synapse->w * (1 - synapse->w) * exp(-tdiff / TAU_PLUS);

        //if(initial_weight > 0 && synapse->w < 0)
        //    synapse->w = 0.0001;
    }
    else if(tdiff < 0 && tdiff > N_WINDOW){ // time window to stdp be considered
        
        //synapse->w -= A_MINUS * exp(time_diff / TAU_MINUS);
        dw -= A_MINUS * synapse->w * (1 - synapse->w) * exp(tdiff / TAU_MINUS);

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

    return pre != NULL && post != NULL && // is not input or output synapse
        post->t_last_spikes[0] != -1 && pre->t_last_spikes[0] != -1 && // at least one spike done by both synapses
        (post->t_last_spikes[(post->next_last_spike - 1) % post->n_last_spikes] == t || 
        pre->t_last_spikes[(pre->next_last_spike - 1) % pre->n_last_spikes] == t); // one of both last spike is in timestamp t
}


void stdp(synapse_t *synapse, int t, int n, double (*stdp_func)(synapse_t *synapse, int time_diff)){

    double dw = 0;
    
    // check conditions to compute STDP
    if(cond_stdp(synapse, t) == 1){

        int d, tdiff1, tdiff2;
        lif_neuron_t *post, *pre, *base = 0, *compute = 0;

        // store pre and postsynaptic neurons
        post = synapse->post_synaptic_lif_neuron;
        pre = synapse->pre_synaptic_lif_neuron;

        // compute tdiffs
        tdiff1 = (post->t_last_spikes[(post->next_last_spike - 1) % n]) - (pre->t_last_spikes[(pre->next_last_spike - 1) % n]);
        tdiff2 = (post->t_last_spikes[(post->next_last_spike - 2) % n]) - (pre->t_last_spikes[(pre->next_last_spike - 2) % n]); 

        // decide how to compute STDP
        if(tdiff1 > 0 || (tdiff1 == 0 && tdiff2 > 0)){
           
            base = post;
            compute = pre;
            d = 1;
        }
        else if(tdiff1 < 0 || (tdiff1 == 0 && tdiff2 < 0)){

            base = pre;
            compute = post;
            d = 0; 
        }

        // check if there is something to compute
        if(tdiff1 != 0 || tdiff2 != 0){

            int t1_base, t2_base, t_compute, step, tdiff;
            t1_base = base->t_last_spikes[(base->next_last_spike - 1) % n];
            t2_base = base->t_last_spikes[(base->next_last_spike - 2) % n];

            // loop spikes and compute stdp
            step = 0; 
            while(step < synapse->stdp_steps && compute->t_last_spikes[(compute->next_last_spike - step - 1) % n] > t2_base){
                
                // get spike timestamp
                t_compute = compute->t_last_spikes[(compute->next_last_spike - step - 1) % n];
            
                // compute compute time difference between spikes
                tdiff = t1_base - t_compute;

                if(d == 0) tdiff = -tdiff;
                
                // compute dw
                dw += stdp_func(synapse, tdiff);

                // update step
                step ++;
            }
        }
    }

    // update weight and store weight change
    synapse->dw += dw;
    synapse->w += dw;
}


void addSTDP(synapse_t *synapse, int t, int n){

    stdp(synapse, t, n, &addSTDP_comp);
}


void mltSTDP(synapse_t *synapse, int t, int n){

    stdp(synapse, t, n, &mltSTDP_comp);
}


void antiSTDP(synapse_t *synapse, int t, int n){

    stdp(synapse, t, n, &antiSTDP_comp);
}


