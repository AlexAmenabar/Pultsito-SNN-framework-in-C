#ifndef CUDA_HELPERS_CUH
#define CUDA_HELPERS_CUH

typedef struct cuda_info_t cuda_info_t;

/// Error handling
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


/// @brief Function to get how much space is available on each GPU
/// @param cuda_info Structure to store the avaiable bytes per device
/// @param dev Device index 
/// @return return error code
int get_memory_info(cuda_info_t *cuda_info, int dev);


/// @brief Function to get GPUs properties
/// @return Structure with the information of the device(s)
cuda_info_t* getGPUProperties();


#ifdef __cplusplus
}
#endif


#endif