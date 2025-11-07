#ifndef GPU_LIF_NEURON_CUH
#define GPU_LIF_NEURON_CUH

#include <stdbool.h>

#define cudaCheckError() {                                          \
    cudaError_t e=cudaGetLastError();                                 \
    if(e!=cudaSuccess) {                                              \
      printf("Cuda failure %s:%d: '%s'\n",__FILE__,__LINE__,cudaGetErrorString(e));           \
      exit(0); \
    }                                                                 \
   }
   
   #define gpuErrchk(call)                                 \
     do {                                                        \
       cudaError_t err = call;                                   \
       if (err != cudaSuccess) {                                 \
         printf("CUDA error at %s %d: %s\n", __FILE__, __LINE__, \
                cudaGetErrorString(err));                        \
         exit(EXIT_FAILURE);                                     \
       }                                                         \
     } while (0)



#ifdef __cplusplus
extern "C" {
#endif


/// @brief LIF neuron model structure
typedef struct {

  int n_neurons;
  int n_input_neurons;
  int n_output_neurons;
  int n_synapses;
  int n_input_synapses;
  int n_output_synapses;

  // neccessary parameters for all neurons
  float *v;
  float *v_thresh;
  float *v_rest;
  int *r_period;
  int *r_period_remain;
  int *res;

  // input and output synapses
  int *in_synapses, *in_off;
  int *out_synapses, *out_off;

  // spike arrays
  int *spk_matrix; 
  int *in_spk_matrix;

  // neuron model dependent variables


  // synapses
  float *w; // synapse weight
  float *dw; // weight change
  int *delay; // latency
  int *lr; // learning rule index
  int *pre_neuron_index; // input neurons...
  int *post_neuron_index;

} GPU_SNN_t;


/// @brief LIF neuron model structure
typedef struct {

  int type;
  int n_classes;

  int n_samples; // number of samples in the dataset
  int n_features; // number of features of the samples
  int *n_spikes; // n spikes per each sample element [n_samples * n_features]

  int *sample_offset; // indicates where each sample starts [n_samples] in the 
  int *feature_offset; // indicates the local offset of each element in the sample [n_samples * ~sample_size]. Sample size can be different for each sample
  //int *sample_offset_in_n_spikes; // offset of the sample in the array of the number of spikes. This is necessary for cases in which samples have different number of features

  int *spikes; // the entire dataset is stored in a 1D array

} GPU_input_info_t;


double process_simulation_lif_neuron(spiking_nn_t *snn, int n, int m, int time_steps);
void simulate_in_GPU(spiking_nn_t *snn, simulation_configuration_t *conf, input_data_t *dataset, simulation_results_t *results);

GPU_SNN_t* copy_snn_structure_to_GPU(spiking_nn_t *snn, int n);
GPU_input_info_t* copy_dataset_to_GPU(input_data_t *dataset);
void simulate_snn_in_GPU(GPU_SNN_t *gpu_snn, spiking_nn_t *snn, simulation_configuration_t *conf, GPU_input_info_t *gpu_dataset, input_data_t *dataset, simulation_results_t *results);


#ifdef __cplusplus
}
#endif


#endif