#ifndef STDP_H
#define STDP_H

#include "snn_library.h"

/// @brief Additive STDP computation
/// @param synapse Synapse to be updated
void add_stdp(synapse_t *synapse, int t);


/// @brief Multiplicative STDP computation
/// @param synapse Synapse to be updated
void mult_stdp(synapse_t *synapse, int t);


/// @brief Anti STDP computation
/// @param synapse Synapse to be updated
void anti_stdp(synapse_t *synapse, int t);


/// @brief Tiplet STDP computation // TODO
/// @param synapse Synapse to be updated
void triplet_stdp(synapse_t *synapse, int t);

#endif