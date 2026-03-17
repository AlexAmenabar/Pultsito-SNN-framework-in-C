#ifndef ARCEUS_H
#define ARCEUS_H

// "umbrella" file to include the entire public API

#include "config/config_loader.h"

#include "datasets/datasets.h"
#include "encoders/image_encoders.h"

#include "networks/snn_generator.h"
#include "networks/snn.h"

#include "neuron_models/neuron_models.h"

#include "simulations/results.h"
#include "simulations/simulations.h"

#include "training_rules/stdp.h"

#include "utils.h"

// include CUDA related files if CUDA is included too
#ifdef CUDA

#endif

#endif