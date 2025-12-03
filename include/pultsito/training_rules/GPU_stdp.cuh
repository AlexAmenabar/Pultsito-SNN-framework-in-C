#ifndef GPU_STDP_CUH
#define GPU_STDP_CUH

#include "snn_library.h"

/*
STDP general functions
*/
int cond_stdp(synapse_t *synapse, int t);
void stdp(synapse_t *synapse, int t, int n, double (*stdp_func)(synapse_t *synapse, int time_diff));
void addSTDP(synapse_t *synapse, int t, int n);
void mltSTDP(synapse_t *synapse, int t, int n);
void antiSTDP(synapse_t *synapse, int t, int n);

#endif