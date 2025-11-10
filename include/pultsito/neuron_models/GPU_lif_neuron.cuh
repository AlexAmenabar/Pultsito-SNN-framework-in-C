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
  int T;

  // neccessary parameters for all neurons
  float *v; // [n_neurons]
  float *v_thresh; // [n_neurons]
  float *v_rest; // [n_neurons]
  int *r_period; // [n_neurons]
  int *r_period_remain; // [n_neurons]
  int *res; // [n_neurons]
  int *t_last_spike; // [n_neurons]
  
  int *n_neuron_input_synapses; // [n_neurons]
  int *neuron_input_synapses; // [n_synapses - n_output_synapses] // if indexes would continuous, it would be smaller
  int *neuron_input_synapses_offset; // [n_neurons]

  int *next_spk; // next spike in input synapses [n_synapses - n_output_synapses] // the same offsets of the previous
  int *last_spk; // last spike [n_neurons]
  

  // synapses
  float *w; // [n_synapses]
  float *dw; // [n_synapses]
  int *delay; // [n_synapses]
  int *lr; // [n_synapses]
  int *pre_neuron_index; // [n_synapses]
  int *post_neuron_index; // [n_synapses]


  int *spk_matrix; // [(n_neurons + n_input_neurons) * T]
  // FLOAT(3N + 2M) + INT(7N + 5M) + NT = 32(10N + 7M) + 32NT bit

  // N = 1.000; M = 100.000; T = 500; 32(10.000 + 700.000 + 500.000) = 32 * 1.210.000 = 38.720.000 = 0,0045 GB
  // N = 100.000; M = 100.000.000; T = 500; COST = 24.032.000.000 = 2,8 GB



  // it should be helpfull to add output synapses to compute STDP easier?
  // neurons are processed? Those variables could help

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