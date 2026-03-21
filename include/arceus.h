#ifndef ARCEUS_H
#define ARCEUS_H

// "umbrella" include file to include the entire public API
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <unistd.h>

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
#include "cuda/cuda_simulations.cuh"
#include "cuda/cuda_networks.cuh"
#include "cuda/cuda_datasets.cuh"
#include "cuda/cuda_results.cuh"
#include "cuda/cuda_simulations_conf.cuh"
#include "cuda/cuda_utils.cuh"
#endif

#endif